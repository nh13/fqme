#![deny(unsafe_code)]
#![allow(
    clippy::must_use_candidate,
    clippy::missing_panics_doc,
    clippy::missing_errors_doc,
    clippy::module_name_repetitions
)]
pub mod tools {
    pub mod bgzf_index;
    pub mod extract;
    pub mod fastq_index;
    pub mod index;
}
pub mod utils;
