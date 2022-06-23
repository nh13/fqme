use std::{
    fs::File,
    io::{self, BufReader, BufWriter, Read, Seek, SeekFrom},
    path::{Path, PathBuf},
};

use anyhow::{bail, ensure, Result};
use byteorder::{LittleEndian, ReadBytesExt};
use clap::Parser;
use env_logger::Env;
use libdeflater::Decompressor;
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

    // Build a BgzfReader starting at the next FASTQ record
    let file = File::open(opts.input.clone()).unwrap();
    let bgzf_reader: BgzfReader = BgzfReader::new(file, start_entry.compressed_offset, num_blocks);

    // Write the FASTQ entries
    let reader = seq_io::fastq::Reader::new(bgzf_reader);
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

pub struct BgzfReader {
    reader: File,
    bgzf: Bgzf,
    header_buf: Vec<u8>,
    compressed_buffer: BytesMut,
    uncompressed_buffer: BytesMut,
    decompressor: Decompressor,
    uncompressed_data: Vec<u8>,
    uncompressed_data_index: usize,
    num_blocks_left: usize,
}

impl BgzfReader {
    fn new(mut reader: File, compressed_offset: u64, num_blocks: usize) -> Self {
        let bgzf = Bgzf::new();
        let header_buf = vec![0; Bgzf::HEADER_SIZE];
        let compressed_buffer = BytesMut::with_capacity(BGZF_BLOCK_SIZE);
        let uncompressed_buffer = BytesMut::with_capacity(BUFSIZE);
        let decompressor = libdeflater::Decompressor::new();
        let uncompressed_data: Vec<u8> = vec![];

        reader.seek(SeekFrom::Start(compressed_offset)).unwrap();

        let mut bgzf_reader = BgzfReader {
            reader,
            bgzf,
            header_buf,
            compressed_buffer,
            uncompressed_buffer,
            decompressor,
            uncompressed_data,
            uncompressed_data_index: 0,
            num_blocks_left: num_blocks,
        };

        // skip over any data until we hit a '@' or end of data
        'outer: loop {
            // fill the data, stop when we have no more data
            if bgzf_reader.fill().unwrap() == 0 {
                break;
            }
            // skip over any data until we hit a '@' or end of data
            while bgzf_reader.uncompressed_data_index < bgzf_reader.uncompressed_data.len() {
                if bgzf_reader.uncompressed_data[bgzf_reader.uncompressed_data_index] == b'@' {
                    break 'outer;
                }
                bgzf_reader.uncompressed_data_index += 1;
            }
        }
        bgzf_reader
    }

    fn bytes_available(&self) -> usize {
        self.uncompressed_data.len() - self.uncompressed_data_index
    }

    fn fill(&mut self) -> io::Result<usize> {
        let available = self.bytes_available();
        if 0 < available {
            return Ok(available);
        }

        // if no more blocks, we fill nothing
        if self.num_blocks_left == 0 {
            return Ok(0);
        }

        // Read in a block!

        // Read the block header
        self.reader.read_exact(&mut self.header_buf)?;
        self.bgzf.check_header(&self.header_buf).unwrap();

        // Read the compressed block data
        let size = self.bgzf.get_block_size(&self.header_buf).unwrap();
        self.compressed_buffer.clear();
        self.compressed_buffer.resize(size - Bgzf::HEADER_SIZE, 0);
        self.reader.read_exact(&mut self.compressed_buffer)?;
        let check = self.bgzf.get_footer_values(&self.compressed_buffer);

        // Decompress the block data
        self.uncompressed_buffer.clear();
        self.uncompressed_buffer.resize(check.amount as usize, 0);
        decompress(
            &self.compressed_buffer,
            &mut self.decompressor,
            &mut self.uncompressed_buffer,
            check,
        )
        .unwrap();

        // Append
        self.uncompressed_data.clear();
        self.uncompressed_data.extend(&self.uncompressed_buffer);
        self.uncompressed_data_index = 0;

        Ok(self.bytes_available())
    }
}

impl Read for BgzfReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        for buf_index in 0..buf.len() {
            // no more data available, try to fill
            if self.bytes_available() == 0 {
                // no more data was filled, return how many bytes we've read
                if let Ok(0) = self.fill() {
                    return Ok(buf_index);
                }
            }
            // add the byte
            buf[buf_index] = self.uncompressed_data[self.uncompressed_data_index];
            self.uncompressed_data_index += 1;
        }
        Ok(buf.len())
    }
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
