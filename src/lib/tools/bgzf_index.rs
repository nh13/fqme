use byteorder::{LittleEndian, ReadBytesExt};
use std::{
    fs::File,
    io::{self, BufReader},
    path::Path,
};

use crate::utils::BUFFERSIZE;

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
    pub fn from<P: AsRef<Path>>(gzi_index: P) -> io::Result<BgzfIndex> {
        let mut reader = BufReader::with_capacity(BUFFERSIZE, File::open(gzi_index)?);

        let num_entries = reader.read_u64::<LittleEndian>()?;

        let mut entries = vec![BgzfIndexOffset { compressed_offset: 0, uncompressed_offset: 0 }];
        for _ in 0..num_entries {
            let compressed_offset = reader.read_u64::<LittleEndian>()?;
            let uncompressed_offset = reader.read_u64::<LittleEndian>()?;
            let entry = BgzfIndexOffset { compressed_offset, uncompressed_offset };
            entries.push(entry);
        }

        Ok(BgzfIndex { num_entries: num_entries + 1, entries })
    }
}
