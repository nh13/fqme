#![allow(clippy::missing_errors_doc, clippy::missing_panics_doc)]
use std::process::exit;

use clap::{Parser, Subcommand};
use env_logger::Env;
use fqme_lib::tools::extract::{run as extract, Opts as ExtractOpts};
use fqme_lib::tools::index::{run as index, Opts as IndexOpts};
use log::error;

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
#[clap(propagate_version = true)]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

// This is where sub-commands are added.  The value of each enum should be the corresponding option
// struct
#[derive(Subcommand)]
enum Commands {
    /// Extracts byte offset and length from a FASTQ index
    Extract(ExtractOpts),
    /// Index a FASTQ
    Index(IndexOpts),
}

#[cfg(not(tarpaulin_include))]
#[allow(clippy::too_many_lines)]
fn main() {
    if std::env::var("RUST_LOG").is_err() {
        std::env::set_var("RUST_LOG", "info");
    }
    env_logger::Builder::from_env(Env::default().default_filter_or("info")).init();

    let cli = Cli::parse();

    let result = match &cli.command {
        Commands::Extract(opts) => extract(opts),
        Commands::Index(opts) => index(opts),
    };

    if let Err(err) = result {
        error!("{:#}", err);
        exit(1);
    }
}
