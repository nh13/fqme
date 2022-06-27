use std::{
    fs::File,
    io::{BufReader, BufWriter, Stdout},
    path::Path,
};

use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use fgoxide::io::Io;
use seq_io::{
    fastq::{Error, OwnedRecord},
    BaseRecord,
};

use crate::utils::BUFFERSIZE;

#[derive(Debug, Clone, PartialEq, Hash, Eq)]
pub struct FastqIndexEntry {
    // total # of records seen
    pub total_records: u64,
    // total # of bytes seen
    pub total_bytes: u64,
}

#[derive(Debug, Clone, PartialEq, Hash, Eq)]
pub struct FastqIndex {
    pub total_records: u64,
    pub nth: u64,
    pub entries: Vec<FastqIndexEntry>,
}

impl FastqIndex {
    pub fn read(path: &Path) -> FastqIndex {
        let mut reader = BufReader::with_capacity(BUFFERSIZE, File::open(path).unwrap());
        let mut entries: Vec<FastqIndexEntry> = vec![];
        let total_records = reader.read_u64::<LittleEndian>().unwrap();
        let nth = reader.read_u64::<LittleEndian>().unwrap();
        while let Ok(num_records) = reader.read_u64::<LittleEndian>() {
            let total_bytes = reader.read_u64::<LittleEndian>().unwrap();
            let entry = FastqIndexEntry { total_records: num_records, total_bytes };
            entries.push(entry);
        }
        FastqIndex { total_records, nth, entries }
    }

    pub fn from(
        records: impl IntoIterator<Item = Result<OwnedRecord, Error>>,
        nth: u64,
        fastq_writer: &mut Option<BufWriter<Stdout>>,
    ) -> FastqIndex {
        let mut total_bytes: u64 = 0;
        let mut total_records: u64 = 0;
        let mut entries: Vec<FastqIndexEntry> = vec![];
        for result in records {
            let rec: OwnedRecord = result.unwrap();
            let num_bytes = FastqIndex::record_to_num_bytes(&rec);

            if total_records % nth == 0 {
                entries.push(FastqIndexEntry { total_records, total_bytes });
            }

            total_records += 1;
            total_bytes += num_bytes as u64;

            if let Some(ref mut writer) = fastq_writer {
                rec.write(writer).unwrap();
            }
        }
        entries.push(FastqIndexEntry { total_records, total_bytes });
        FastqIndex { total_records, nth, entries }
    }

    pub fn write(self, output: &Path) {
        let mut writer = Io::default().new_writer(&output).unwrap();
        writer.write_u64::<LittleEndian>(self.total_records).unwrap();
        writer.write_u64::<LittleEndian>(self.nth).unwrap();
        for entry in self.entries {
            writer.write_u64::<LittleEndian>(entry.total_records).unwrap();
            writer.write_u64::<LittleEndian>(entry.total_bytes).unwrap();
        }
    }

    pub fn record_to_num_bytes(rec: &OwnedRecord) -> u64 {
        // NB: this is incorrect if there exists comment
        let num_bytes = 1 // leading '@'
            + rec.head.len()
            + 1 // newline
            + rec.seq.len()
            + 1 // newline
            + 2 // + and newline
            + rec.qual.len()
            + 1; // newline
        num_bytes as u64
    }

    // NB: start_record and end_record are 1-based inclusive
    pub fn range(&self, start_record: u64, end_record: u64) -> Option<FastqIndexRange> {
        if end_record < start_record || end_record < 1 || self.total_records < start_record {
            return None;
        }
        if start_record < 1 {
            return self.range(1, end_record);
        }
        if self.total_records < end_record {
            return self.range(start_record, self.total_records);
        }

        let mut leading_records: u64 = start_record - 1;
        let mut trailing_records: u64 = 0;
        let mut start_byte: u64 = 0;
        let mut end_byte: u64 = 0;
        let mut total_records: u64 = 0;
        let mut last_total_records: u64 = 0;

        // TODO: binary search
        let mut start_entry = 0;
        for entry in &self.entries {
            if start_record <= entry.total_records {
                leading_records = (start_record - 1) % self.nth;
                break;
            }
            start_byte = entry.total_bytes;
            last_total_records = entry.total_records;
            start_entry += 1;
        }

        // TODO: binary search
        for i in start_entry..self.entries.len() {
            let entry = &self.entries[i];

            if end_record <= entry.total_records {
                trailing_records = entry.total_records - end_record;
                end_byte = entry.total_bytes;
                total_records = entry.total_records - last_total_records;
                break;
            }
        }

        Some(FastqIndexRange {
            start_byte,
            end_byte,
            leading_records,
            trailing_records,
            total_records,
        })
    }
}

#[derive(Debug, Clone, PartialEq, Hash, Eq)]
pub struct FastqIndexRange {
    // the uncompressed start byte range
    pub start_byte: u64,
    // the uncompressed end byte range
    pub end_byte: u64,
    // the number of leading records in the range
    pub leading_records: u64,
    // the number of trailing records in the range
    pub trailing_records: u64,
    // the total number of records in the range
    pub total_records: u64,
}

impl FastqIndexRange {
    /// the total number of selected records in the ragne
    pub fn selected_records(self) -> u64 {
        self.total_records - self.leading_records - self.trailing_records
    }

    pub fn num_bytes(&self) -> u64 {
        self.end_byte - self.start_byte
    }
}

#[cfg(test)]
mod test {
    use crate::tools::fastq_index::{FastqIndex, FastqIndexEntry};
    use seq_io::fastq::{Error, OwnedRecord};

    use super::FastqIndexRange;

    fn record() -> OwnedRecord {
        OwnedRecord {
            head: b"some-read-name".to_vec(),
            seq: b"GATTACA".to_vec(),
            qual: b"IIIIIII".to_vec(),
        }
    }

    fn index() -> FastqIndex {
        let records: Vec<Result<OwnedRecord, Error>> = {
            vec![
                Ok(record()),
                Ok(record()),
                Ok(record()),
                Ok(record()),
                Ok(record()),
                Ok(record()),
                Ok(record()),
                Ok(record()),
            ]
        };
        FastqIndex::from(records, 3, &mut None)
    }

    #[test]
    fn test_fastq_index_record_to_num_bytes() {
        assert_eq!(FastqIndex::record_to_num_bytes(&record()), 34);
    }

    fn test_fastq_index_from(
        records: Vec<Result<OwnedRecord, Error>>,
        nth: u64,
        index_entries: usize,
    ) {
        let num_input_records = records.len();
        let index = FastqIndex::from(records, nth, &mut None);
        let record_num_bytes = FastqIndex::record_to_num_bytes(&record());
        assert_eq!(index.entries.len(), index_entries);
        for i in 0..nth as usize {
            let num_records: u64 = std::cmp::min(num_input_records as u64, nth * i as u64);
            let expected: FastqIndexEntry = FastqIndexEntry {
                total_records: num_records,
                total_bytes: record_num_bytes * num_records as u64,
            };
            assert_eq!(index.entries[i], expected);
        }
    }

    #[test]
    fn test_fastq_index_from_even_with_even_nth() {
        // even # of inputs, even nth
        let records: Vec<Result<OwnedRecord, Error>> =
            vec![Ok(record()), Ok(record()), Ok(record()), Ok(record())];
        test_fastq_index_from(records, 2, 3);
    }
    #[test]

    fn test_fastq_index_from_odd_with_even_nth() {
        // even # of inputs, even nth
        let records: Vec<Result<OwnedRecord, Error>> =
            vec![Ok(record()), Ok(record()), Ok(record()), Ok(record()), Ok(record())];
        test_fastq_index_from(records, 2, 4);
    }

    #[test]
    fn test_fastq_index_from_even_with_odd_nth() {
        // even # of inputs, even nth
        let records: Vec<Result<OwnedRecord, Error>> =
            vec![Ok(record()), Ok(record()), Ok(record()), Ok(record())];
        test_fastq_index_from(records, 3, 3);
    }

    #[test]
    fn test_fastq_index_from_odd_with_odd_nth() {
        // even # of inputs, even nth
        let records: Vec<Result<OwnedRecord, Error>> =
            vec![Ok(record()), Ok(record()), Ok(record()), Ok(record()), Ok(record())];
        test_fastq_index_from(records, 3, 3);
    }

    #[test]
    fn test_fastq_index_range_num_bytes() {
        let entry = FastqIndexRange {
            start_byte: 15,
            end_byte: 20,
            leading_records: 0,
            trailing_records: 0,
            total_records: 0,
        };
        assert_eq!(entry.num_bytes(), 5);
    }

    #[test]
    fn test_fastq_index_range_out_of_range() {
        let index: FastqIndex = index();

        // [0, 0] should yield zero selected records
        assert_eq!(index.range(0, 0), None);

        // [4, 3] should yield zero selected records
        assert_eq!(index.range(4, 3), None);

        // 10, 10
        assert_eq!(index.range(9, 9), None);
    }

    #[test]
    fn test_fastq_index_range_in_range() {
        let index: FastqIndex = index();

        // [0, 2] should yield two selected records
        let expected = FastqIndexRange {
            start_byte: 0,
            end_byte: 102,
            leading_records: 0,
            trailing_records: 1,
            total_records: 3,
        };
        let range = index.range(0, 2).unwrap();
        assert_eq!(range, expected);
        assert_eq!(range.selected_records(), 2);

        // [1, 3]
        let expected = FastqIndexRange {
            start_byte: 0,
            end_byte: 102,
            leading_records: 0,
            trailing_records: 0,
            total_records: 3,
        };
        let range = index.range(1, 3).unwrap();
        assert_eq!(range, expected);
        assert_eq!(range.selected_records(), 3);

        // [2, 2] should yield a single selected record
        let expected = FastqIndexRange {
            start_byte: 0,
            end_byte: 102,
            leading_records: 1,
            trailing_records: 1,
            total_records: 3,
        };
        let range = index.range(2, 2).unwrap();
        assert_eq!(range, expected);
        assert_eq!(range.selected_records(), 1);

        // [3, 3] should yield a single selected record
        let expected = FastqIndexRange {
            start_byte: 0,
            end_byte: 102,
            leading_records: 2,
            trailing_records: 0,
            total_records: 3,
        };
        let range = index.range(3, 3).unwrap();
        assert_eq!(range, expected);
        assert_eq!(range.selected_records(), 1);

        // [2, 3]
        let expected = FastqIndexRange {
            start_byte: 0,
            end_byte: 102,
            leading_records: 1,
            trailing_records: 0,
            total_records: 3,
        };
        let range = index.range(2, 3).unwrap();
        assert_eq!(range, expected);
        assert_eq!(range.selected_records(), 2);

        // // [2, 4]
        let expected = FastqIndexRange {
            start_byte: 0,
            end_byte: 204,
            leading_records: 1,
            trailing_records: 2,
            total_records: 6,
        };
        let range = index.range(2, 4).unwrap();
        assert_eq!(range, expected);
        assert_eq!(range.selected_records(), 3);

        // [3, 4]
        let expected = FastqIndexRange {
            start_byte: 0,
            end_byte: 204,
            leading_records: 2,
            trailing_records: 2,
            total_records: 6,
        };
        let range = index.range(3, 4).unwrap();
        assert_eq!(range, expected);
        assert_eq!(range.selected_records(), 2);

        // [1, 8]
        let expected = FastqIndexRange {
            start_byte: 0,
            end_byte: 272,
            leading_records: 0,
            trailing_records: 0,
            total_records: 8,
        };
        let range = index.range(1, 8).unwrap();
        assert_eq!(range, expected);
        assert_eq!(range.selected_records(), 8);

        // [8, 9]
        let expected = FastqIndexRange {
            start_byte: 204,
            end_byte: 272,
            leading_records: 1,
            trailing_records: 0,
            total_records: 2,
        };
        let range = index.range(8, 9).unwrap();
        assert_eq!(range, expected);
        assert_eq!(range.selected_records(), 1);
    }
}
