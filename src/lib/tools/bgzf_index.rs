use byteorder::{LittleEndian, ReadBytesExt};
use std::{fs::File, io::BufReader};

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
