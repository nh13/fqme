use byteorder::{LittleEndian, ReadBytesExt};
use std::{
    fs::File,
    io::{self, BufReader},
    path::Path,
};

use crate::utils::BUFFERSIZE;

pub struct BgzfIndex {
    pub entries: Vec<BgzfIndexOffset>,
}

#[derive(Clone, Copy, Debug)]
pub struct BgzfIndexOffset {
    pub compressed_offset: u64,
    pub uncompressed_offset: u64,
}

impl BgzfIndex {
    pub fn from<P: AsRef<Path>>(gzi_index: P) -> io::Result<BgzfIndex> {
        let mut reader = BufReader::with_capacity(BUFFERSIZE, File::open(gzi_index)?);

        let file_entry_count = reader.read_u64::<LittleEndian>()?;

        // Prepend a synthetic zero entry for the start of the file
        let mut entries = vec![BgzfIndexOffset { compressed_offset: 0, uncompressed_offset: 0 }];
        for _ in 0..file_entry_count {
            let compressed_offset = reader.read_u64::<LittleEndian>()?;
            let uncompressed_offset = reader.read_u64::<LittleEndian>()?;
            let entry = BgzfIndexOffset { compressed_offset, uncompressed_offset };
            entries.push(entry);
        }

        Ok(BgzfIndex { entries })
    }
}
