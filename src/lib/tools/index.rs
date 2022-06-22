use fgoxide::io::Io;
use std::{path::PathBuf, io, io::{Write, BufReader}, io:: BufWriter};

use anyhow::Result;
use clap::Parser;
use env_logger::Env;

use seq_io::{BaseRecord, fastq::OwnedRecord};

use crate::utils::{BUFFERSIZE, built_info};

/// Index a FASTQ
#[derive(Parser, Debug)]
#[clap(name = "fq2bam", verbatim_doc_comment, version = built_info::VERSION.as_str())]
pub struct Opts {

    /// The output index file.
    #[clap(short = 'o', long, display_order = 2)]
    pub output: PathBuf,

    /// Index every Nth entry
    #[clap(short = 'n', long, default_value = "100000", display_order = 3)]
    pub nth: usize,

    /// True to omit emitting the FASTQ to stdout
    #[clap(long, display_order = 4)]
    pub no_stdout: bool,
}

// Run index
#[allow(clippy::too_many_lines)]
pub fn run(opts: &Opts) -> Result<(), anyhow::Error> {
    let reader = {
        let reader =  BufReader::with_capacity(BUFFERSIZE, io::stdin());
        seq_io::fastq::Reader::new(reader).into_records()
    };

    let mut index_writer = Io::default().new_writer(&opts.output).unwrap();

    let mut fastq_writer = {
        if opts.no_stdout {
            None
        } else {
            Some(BufWriter::with_capacity(BUFFERSIZE, io::stdout()))
        }
    };

    let mut total_bytes: usize = 0;
    let mut num_records: usize = 0;
    for result in reader {
        let rec: OwnedRecord = result.unwrap();

        // NB: this is incorrect if there exists comment
        let num_bytes: usize = {
            1 // leading '@'
            + rec.head.len() 
            + 1 // newline
            + rec.seq.len()
            + 1 // newline
            + 2 // + and newline
            + rec.qual.len()
            + 1 // newline
        };

        num_records += 1;

        if num_records % opts.nth == 0 {
            index_writer.write_all(format!("{}\t{}\n", num_records, total_bytes).as_bytes())?;
        }

        total_bytes += num_bytes;

        if let Some(ref mut writer) = fastq_writer {
            rec.write(writer)?;
        }
    }
    if num_records % opts.nth != 0 {
        index_writer.write_all(format!("{}\t{}\n", num_records, total_bytes).as_bytes())?;
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
