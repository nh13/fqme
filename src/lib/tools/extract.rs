use std::{
    fs::File,
    io::{self, BufReader, BufWriter, Read, Seek, SeekFrom},
    path::{Path, PathBuf},
};

use anyhow::{bail, ensure, Result};
use byteorder::{LittleEndian, ReadBytesExt};
use clap::Parser;
use env_logger::Env;
use seq_io::BaseRecord;

use crate::utils::{built_info, BUFFERSIZE};

use bytes::BytesMut;
use gzp::{deflate::Bgzf, BlockFormatSpec, FooterValues, FormatSpec, GzpError, BUFSIZE};

const BGZF_BLOCK_SIZE: usize = 65280;

/// Finds the (uncompressed) offset and length to use for bgzip -b <OFFSET> -s <LENGTH>
#[derive(Parser, Debug)]
#[clap(name = "fq2bam", verbatim_doc_comment, version = built_info::VERSION.as_str())]
pub struct Opts {
    /// The input bgzip'ed FASTQ file.
    #[clap(short = 'f', long, display_order = 1)]
    pub input: PathBuf,

    /// The first record to display (1-based).
    #[clap(short = 's', long, display_order = 2)]
    pub start: Option<u64>,

    /// The last record to display (1-based inclusive).
    #[clap(short = 'e', long, display_order = 2)]
    pub end: Option<u64>,
}

// Run extract
#[allow(clippy::too_many_lines)]
pub fn run(opts: &Opts) -> Result<(), anyhow::Error> {
    let (start, end) = match (opts.start, opts.end) {
        (Some(s), Some(e)) => (s, e),
        (Some(s), None) => (s, s),
        (None, Some(e)) => (e, e),
        (None, None) => bail!("Either --start or --end must be given"),
    };
    ensure!(start <= end, "--start must be less than or equal to --end");

    // Create the file names
    let fqi_path = format!("{}.{}", opts.input.to_string_lossy(), "fqi");
    let gzi_path = format!("{}.{}", opts.input.to_string_lossy(), "gzi");

    // Read the FASTQ index
    let fqi_range = FastqIndexRange::from(&fqi_path, start, end);
    // println!(
    //     "The following command will output {} leading and {} trailing records:",
    //     fqi_range.leading_records, fqi_range.trailing_records
    // );
    // println!("    bgzip -b {} -s {} <file.gz>", fqi_range.start_byte, fqi_range.num_bytes);

    // Read the BGZF index and find the compressed offset
    let gzi = BgzfIndex::from(gzi_path);
    let mut start_entry: BgzfIndexOffset = gzi.entries[0];
    let mut num_blocks: usize = 0;
    for entry in gzi.entries {
        if entry.uncompressed_offset < fqi_range.start_byte as u64 {
            start_entry = entry;
        }
        num_blocks += 1;
        if entry.uncompressed_offset >= fqi_range.num_bytes as u64 {
            break;
        }
    }

    // Open the BGZF file and seek to the block offset
    let mut reader = File::open(opts.input.clone())?;
    reader.seek(SeekFrom::Start(start_entry.compressed_offset))?;

    // Loop through all the blocks
    let bgzf = Bgzf::new();
    let mut header_buf = vec![0; Bgzf::HEADER_SIZE];
    let mut compressed_buffer = BytesMut::with_capacity(BGZF_BLOCK_SIZE);
    let mut uncompressed_buffer = BytesMut::with_capacity(BUFSIZE);
    let mut decompressor = libdeflater::Decompressor::new();
    let mut uncompressed_data: Vec<u8> = vec![];
    for _ in 0..num_blocks {
        // Read the block header
        reader.read_exact(&mut header_buf)?;
        bgzf.check_header(&header_buf).unwrap();

        // Read the compressed block data
        let size = bgzf.get_block_size(&header_buf).unwrap();
        compressed_buffer.clear();
        compressed_buffer.resize(size - Bgzf::HEADER_SIZE, 0);
        reader.read_exact(&mut compressed_buffer)?;
        let check = bgzf.get_footer_values(&compressed_buffer);

        // Decompress the block data
        uncompressed_buffer.clear();
        uncompressed_buffer.resize(check.amount as usize, 0);
        decompress(&compressed_buffer, &mut decompressor, &mut uncompressed_buffer, check).unwrap();

        // Append
        uncompressed_data.extend_from_slice(&uncompressed_buffer);
    }

    // Write the FASTQ entries
    let data_iter: Vec<u8> =
        uncompressed_data.iter().skip_while(|x| *x != &b'@').copied().collect();
    let reader = seq_io::fastq::Reader::new(data_iter.as_slice());
    let mut writer = BufWriter::with_capacity(BUFFERSIZE, io::stdout());
    let mut num_to_write: u64 = end - start + 1;
    for (index, result) in reader.into_records().enumerate() {
        let rec = result?;
        if index as u64 >= fqi_range.leading_records {
            rec.write(&mut writer)?;
            num_to_write -= 1;
        }
        if num_to_write == 0 {
            break;
        }
    }

    Ok(())
}

#[inline]
pub fn decompress(
    input: &[u8],
    decoder: &mut libdeflater::Decompressor,
    output: &mut [u8],
    footer_vals: FooterValues,
) -> Result<(), GzpError> {
    if footer_vals.amount != 0 {
        let _bytes_decompressed = decoder.deflate_decompress(&input[..input.len() - 8], output)?;
    }
    let mut new_check = libdeflater::Crc::new();
    new_check.update(output);

    if footer_vals.sum != new_check.sum() {
        return Err(GzpError::InvalidCheck { found: new_check.sum(), expected: footer_vals.sum });
    }
    Ok(())
}

/// Parse args and set up logging / tracing
pub fn setup() -> Opts {
    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info");
    }
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();

    Opts::parse()
}

#[derive(Debug, Clone, PartialEq, Hash, Eq)]
pub struct FastqIndex {
    pub num_records: u64,

    pub total_bytes: u64,
}

#[derive(Debug, Clone, PartialEq, Hash, Eq)]
pub struct FastqIndexRange {
    pub start_byte: u64,
    pub num_bytes: u64,
    pub leading_records: u64,
    pub trailing_records: u64,
    pub num_records: u64,
}

impl FastqIndexRange {
    pub fn from(fastq_index: &str, start_record: u64, end_record: u64) -> FastqIndexRange {
        let indexes = {
            let path = Path::new(fastq_index);
            let mut reader = BufReader::with_capacity(BUFFERSIZE, File::open(path).unwrap());
            let mut entries: Vec<FastqIndex> = vec![];
            while let Ok(num_records) = reader.read_u64::<LittleEndian>() {
                let total_bytes = reader.read_u64::<LittleEndian>().unwrap();
                let entry = FastqIndex { num_records, total_bytes };
                entries.push(entry);
            }
            entries
        };

        let mut leading_records: u64 = start_record - 1;
        let mut trailing_records: u64 = 0;
        let mut start_byte: u64 = 0;
        let mut end_byte: u64 = 0;
        let mut num_records: u64 = 0;
        for entry in indexes {
            if entry.num_records < start_record {
                leading_records = start_record - entry.num_records - 1;
                start_byte = entry.total_bytes;
            }
            num_records += entry.num_records;
            end_byte = entry.total_bytes;
            if entry.num_records >= end_record {
                trailing_records = entry.num_records - end_record;
                break;
            }
        }

        FastqIndexRange {
            start_byte,
            num_bytes: (end_byte - start_byte),
            leading_records,
            trailing_records,
            num_records,
        }
    }
}

pub struct BgzfIndex {
    pub num_entries: u64,
    pub entries: Vec<BgzfIndexOffset>,
}

#[derive(Clone, Copy, Debug)]
pub struct BgzfIndexOffset {
    pub compressed_offset: u64,
    pub uncompressed_offset: u64,
}

impl BgzfIndex {
    pub fn from(gzi_index: String) -> BgzfIndex {
        let mut reader = BufReader::with_capacity(BUFFERSIZE, File::open(gzi_index).unwrap());

        let num_entries = reader.read_u64::<LittleEndian>().unwrap();

        let mut entries = vec![BgzfIndexOffset { compressed_offset: 0, uncompressed_offset: 0 }];
        for _ in 0..num_entries {
            let compressed_offset = reader.read_u64::<LittleEndian>().unwrap();
            let uncompressed_offset = reader.read_u64::<LittleEndian>().unwrap();
            let entry = BgzfIndexOffset { compressed_offset, uncompressed_offset };
            entries.push(entry);
        }

        BgzfIndex { num_entries: num_entries + 1, entries }
    }
}
