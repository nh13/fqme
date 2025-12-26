use std::{io, io::BufReader, io::BufWriter, path::PathBuf};

use anyhow::Result;
use clap::Parser;
use env_logger::Env;

use crate::utils::{built_info, BUFFERSIZE};

use super::fastq_index::FastqIndex;

/// Index a FASTQ
#[derive(Parser, Debug)]
#[clap(name = "fq2bam", verbatim_doc_comment, version = built_info::VERSION.as_str())]
pub struct Opts {
    /// The output index file.
    #[clap(short = 'o', long, display_order = 2)]
    pub output: PathBuf,

    /// Index every Nth entry
    #[clap(short = 'n', long, default_value = "100000", display_order = 3)]
    pub nth: u64,

    /// True to omit emitting the FASTQ to stdout
    #[clap(long, display_order = 4)]
    pub no_stdout: bool,
}

// Run index
#[allow(clippy::too_many_lines)]
pub fn run(opts: &Opts) -> Result<(), anyhow::Error> {
    let reader = BufReader::with_capacity(BUFFERSIZE, io::stdin());

    let mut fastq_writer = if opts.no_stdout {
        None
    } else {
        Some(BufWriter::with_capacity(BUFFERSIZE, io::stdout()))
    };

    FastqIndex::from_reader(reader, opts.nth, &mut fastq_writer)?.write(opts.output.as_path())?;

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
