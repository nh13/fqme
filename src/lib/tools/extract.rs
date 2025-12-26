use std::{
    fs::File,
    io::{self, BufWriter, ErrorKind, Read, Seek, SeekFrom},
    path::{Path, PathBuf},
};

use anyhow::{bail, ensure, Result};
use clap::Parser;
use env_logger::Env;
use libdeflater::Decompressor;
use seq_io::BaseRecord;

use crate::{
    tools::{
        bgzf_index::{BgzfIndex, BgzfIndexOffset},
        fastq_index::FastqIndex,
    },
    utils::{built_info, BUFFERSIZE},
};

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
    let fastq_index = FastqIndex::read(Path::new(&fqi_path))?;
    let fqi_range = match fastq_index.range(start, end) {
        Some(range) => range,
        None => {
            log::warn!(
                "Requested range [{}, {}] is out of bounds (file has {} records)",
                start,
                end,
                fastq_index.total_records
            );
            return Ok(());
        }
    };
    // println!("{:?}", fqi_range);
    // println!(
    //     "The following command will output {} leading and {} trailing records:",
    //     fqi_range.leading_records, fqi_range.trailing_records
    // );
    // println!("    bgzip -b {} -s {} {:?}", fqi_range.start_byte, fqi_range.num_bytes(), opts.input);

    // Read the BGZF index and find the compressed offset
    let gzi = BgzfIndex::from(&gzi_path)?;
    ensure!(!gzi.entries.is_empty(), "GZI index file is empty or corrupted: {}", gzi_path);
    let mut start_entry: BgzfIndexOffset = gzi.entries[0];
    let mut num_blocks: usize = 0;
    for entry in gzi.entries {
        if entry.uncompressed_offset < fqi_range.start_byte {
            start_entry = entry;
            num_blocks = 0;
        }
        num_blocks += 1;
        if entry.uncompressed_offset >= fqi_range.end_byte {
            break;
        }
    }

    // Build a BgzfReader starting at the next FASTQ record
    let file = File::open(&opts.input)?;
    let bgzf_reader: BgzfReader =
        BgzfReader::new(file, fqi_range.start_byte, start_entry, num_blocks);

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
    fn new(mut reader: File, start_byte: u64, entry: BgzfIndexOffset, num_blocks: usize) -> Self {
        let bgzf = Bgzf::new();
        let header_buf = vec![0; Bgzf::HEADER_SIZE];
        let compressed_buffer = BytesMut::with_capacity(BGZF_BLOCK_SIZE);
        let uncompressed_buffer = BytesMut::with_capacity(BUFSIZE);
        let decompressor = libdeflater::Decompressor::new();
        let uncompressed_data: Vec<u8> = vec![];

        reader.seek(SeekFrom::Start(entry.compressed_offset)).unwrap();

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

        // move to the start uncompressed byte offset
        let mut cur_uncompressed_offset = entry.uncompressed_offset;
        while cur_uncompressed_offset < start_byte {
            // fill the data, stop when we have no more data
            if bgzf_reader.bytes_available() == 0 && bgzf_reader.fill().unwrap() == 0 {
                break;
            }
            cur_uncompressed_offset += 1;
            bgzf_reader.uncompressed_data_index += 1;
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
        match self.reader.read_exact(&mut self.header_buf) {
            Ok(()) => (),
            Err(ref e) if e.kind() == ErrorKind::UnexpectedEof => return Ok(0),
            e => e.unwrap(),
        }
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
        #[allow(clippy::needless_range_loop)]
        for buf_index in 0..buf.len() {
            // no more data available, try to fill
            if self.bytes_available() == 0 {
                // no more data was filled, return how many bytes we've read
                if self.fill().unwrap() == 0 {
                    return Ok(buf_index);
                }
            }
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
