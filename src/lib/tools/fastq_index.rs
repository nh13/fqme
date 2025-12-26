use std::{
    fs::File,
    io::{self, BufRead, BufReader, BufWriter, Stdout, Write},
    path::Path,
};

use byteorder::{LittleEndian, ReadBytesExt, WriteBytesExt};
use fgoxide::io::Io;
use seq_io::{
    fastq::{Error, OwnedRecord, Reader as FastqReader},
    BaseRecord,
};

use crate::utils::BUFFERSIZE;

/// A writer that counts bytes written without storing them.
/// Used to measure the exact byte size of FASTQ records including separator lines.
pub struct ByteCountingWriter {
    count: u64,
}

impl ByteCountingWriter {
    pub fn new() -> Self {
        Self { count: 0 }
    }

    pub fn count(&self) -> u64 {
        self.count
    }

    pub fn reset(&mut self) {
        self.count = 0;
    }
}

impl Default for ByteCountingWriter {
    fn default() -> Self {
        Self::new()
    }
}

impl Write for ByteCountingWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.count += buf.len() as u64;
        Ok(buf.len())
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

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
    pub fn read(path: &Path) -> io::Result<FastqIndex> {
        let mut reader = BufReader::with_capacity(BUFFERSIZE, File::open(path)?);
        let mut entries: Vec<FastqIndexEntry> = vec![];
        let total_records = reader.read_u64::<LittleEndian>()?;
        let nth = reader.read_u64::<LittleEndian>()?;
        while let Ok(num_records) = reader.read_u64::<LittleEndian>() {
            let total_bytes = reader.read_u64::<LittleEndian>()?;
            let entry = FastqIndexEntry { total_records: num_records, total_bytes };
            entries.push(entry);
        }
        Ok(FastqIndex { total_records, nth, entries })
    }

    /// Creates a FastqIndex from an iterator of records.
    ///
    /// **Note**: This method computes byte sizes from record fields, which is incorrect
    /// for FASTQ files with comments on the `+` line. Use `from_reader` for accurate
    /// byte position tracking.
    #[deprecated(since = "0.1.0", note = "Use from_reader for accurate byte position tracking")]
    #[allow(deprecated)]
    pub fn from(
        records: impl IntoIterator<Item = Result<OwnedRecord, Error>>,
        nth: u64,
        fastq_writer: &mut Option<BufWriter<Stdout>>,
    ) -> io::Result<FastqIndex> {
        let mut total_bytes: u64 = 0;
        let mut total_records: u64 = 0;
        let mut entries: Vec<FastqIndexEntry> = vec![];
        for result in records {
            let rec: OwnedRecord =
                result.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
            let num_bytes = FastqIndex::record_to_num_bytes(&rec);

            #[allow(unknown_lints, clippy::manual_is_multiple_of)]
            if total_records % nth == 0 {
                entries.push(FastqIndexEntry { total_records, total_bytes });
            }

            total_records += 1;
            total_bytes += num_bytes;

            if let Some(ref mut writer) = fastq_writer {
                rec.write(writer)?;
            }
        }
        entries.push(FastqIndexEntry { total_records, total_bytes });
        Ok(FastqIndex { total_records, nth, entries })
    }

    /// Creates a FastqIndex from a reader, tracking actual byte positions.
    ///
    /// This method correctly handles FASTQ files with comments on the `+` line
    /// by using `write_unchanged` to measure exact record sizes.
    pub fn from_reader<R: BufRead>(
        reader: R,
        nth: u64,
        fastq_writer: &mut Option<BufWriter<Stdout>>,
    ) -> io::Result<FastqIndex> {
        let mut fastq_reader = FastqReader::new(reader);
        let mut total_records: u64 = 0;
        let mut total_bytes: u64 = 0;
        let mut entries: Vec<FastqIndexEntry> = vec![];
        let mut byte_counter = ByteCountingWriter::new();

        while let Some(result) = fastq_reader.next() {
            let rec = result.map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

            #[allow(unknown_lints, clippy::manual_is_multiple_of)]
            if total_records % nth == 0 {
                entries.push(FastqIndexEntry { total_records, total_bytes });
            }

            // Measure exact byte size using write_unchanged
            byte_counter.reset();
            rec.write_unchanged(&mut byte_counter)?;
            total_bytes += byte_counter.count();

            total_records += 1;

            if let Some(ref mut writer) = fastq_writer {
                rec.write(writer)?;
            }
        }

        // Final entry with total counts
        entries.push(FastqIndexEntry { total_records, total_bytes });

        Ok(FastqIndex { total_records, nth, entries })
    }

    #[allow(unknown_lints, clippy::io_other_error)]
    pub fn write(self, output: &Path) -> io::Result<()> {
        let mut writer = Io::default()
            .new_writer(&output)
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        writer.write_u64::<LittleEndian>(self.total_records)?;
        writer.write_u64::<LittleEndian>(self.nth)?;
        for entry in self.entries {
            writer.write_u64::<LittleEndian>(entry.total_records)?;
            writer.write_u64::<LittleEndian>(entry.total_bytes)?;
        }
        Ok(())
    }

    /// **Note**: This method assumes the '+' line is exactly "+\n", which is incorrect
    /// for FASTQ files with comments on the '+' line. Use `from_reader` with actual
    /// byte counting for accurate results.
    #[deprecated(
        since = "0.1.0",
        note = "Incorrect for FASTQ with '+' line comments; use from_reader"
    )]
    pub fn record_to_num_bytes(rec: &OwnedRecord) -> u64 {
        // NB: this is incorrect if there exists a comment on the + line
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

    #[allow(deprecated)]
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
        FastqIndex::from(records, 3, &mut None).unwrap()
    }

    #[test]
    #[allow(deprecated)]
    fn test_fastq_index_record_to_num_bytes() {
        assert_eq!(FastqIndex::record_to_num_bytes(&record()), 34);
    }

    #[allow(deprecated)]
    fn test_fastq_index_from(
        records: Vec<Result<OwnedRecord, Error>>,
        nth: u64,
        index_entries: usize,
    ) {
        let num_input_records = records.len();
        let index = FastqIndex::from(records, nth, &mut None).unwrap();
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

    #[test]
    fn test_from_reader_tracks_positions() {
        use std::io::Cursor;

        // Create FASTQ data as bytes
        let fastq_data = b"@read1\nACGT\n+\nIIII\n@read2\nGGGG\n+\nJJJJ\n@read3\nTTTT\n+\nKKKK\n";
        let reader = Cursor::new(&fastq_data[..]);

        let index = FastqIndex::from_reader(reader, 2, &mut None).unwrap();

        assert_eq!(index.total_records, 3);
        assert_eq!(index.entries.len(), 3); // entries at 0, 2, and final

        // First entry: at record 0, byte 0
        assert_eq!(index.entries[0].total_records, 0);
        assert_eq!(index.entries[0].total_bytes, 0);

        // Second entry: at record 2, after 2 records (each 19 bytes: @read\nACGT\n+\nIIII\n)
        assert_eq!(index.entries[1].total_records, 2);
        assert_eq!(index.entries[1].total_bytes, 38); // 19 * 2

        // Final entry: all 3 records
        assert_eq!(index.entries[2].total_records, 3);
        assert_eq!(index.entries[2].total_bytes, 57); // 19 * 3
    }

    #[test]
    fn test_from_reader_with_plus_comment() {
        use std::io::Cursor;

        // FASTQ data WITH comment on + line (the bug case)
        let fastq_data = b"@read1\nACGT\n+this is a comment\nIIII\n@read2\nGGGG\n+\nJJJJ\n";
        let reader = Cursor::new(&fastq_data[..]);

        let index = FastqIndex::from_reader(reader, 1, &mut None).unwrap();

        assert_eq!(index.total_records, 2);
        assert_eq!(index.entries.len(), 3); // entries at 0, 1, and final

        // First record is 36 bytes: @read1\nACGT\n+this is a comment\nIIII\n
        // Second record is 19 bytes: @read2\nGGGG\n+\nJJJJ\n
        assert_eq!(index.entries[0].total_bytes, 0);
        assert_eq!(index.entries[1].total_bytes, 36); // First record with comment
        assert_eq!(index.entries[2].total_bytes, 55); // 36 + 19
    }

    // ==================== Edge Case Tests ====================
    // These tests validate the range calculation logic, particularly:
    // - The modulo calculation: (start_record - 1) % self.nth
    // - The loop termination when end_record is at/near total_records

    /// Helper to create an index with a specific number of records
    fn index_with_n_records(n: usize, nth: u64) -> FastqIndex {
        #[allow(deprecated)]
        let index = FastqIndex::from((0..n).map(|_| Ok(record())), nth, &mut None).unwrap();
        index
    }

    #[test]
    fn test_range_single_record_index() {
        // Index with just 1 record, nth=1
        // Entries: [0, 1] (at bytes [0, 34])
        let index = index_with_n_records(1, 1);
        assert_eq!(index.total_records, 1);
        assert_eq!(index.entries.len(), 2);

        // Query the only record
        let range = index.range(1, 1).unwrap();
        assert_eq!(range.start_byte, 0);
        assert_eq!(range.end_byte, 34);
        assert_eq!(range.leading_records, 0);
        assert_eq!(range.trailing_records, 0);
        assert_eq!(range.selected_records(), 1);
    }

    #[test]
    fn test_range_single_record_index_large_nth() {
        // Index with just 1 record, nth=100 (larger than total records)
        // Entries: [0, 1] (at bytes [0, 34])
        let index = index_with_n_records(1, 100);
        assert_eq!(index.total_records, 1);
        assert_eq!(index.entries.len(), 2); // entry at 0 + final entry

        let range = index.range(1, 1).unwrap();
        assert_eq!(range.start_byte, 0);
        assert_eq!(range.end_byte, 34);
        assert_eq!(range.leading_records, 0);
        assert_eq!(range.trailing_records, 0);
        assert_eq!(range.selected_records(), 1);
    }

    #[test]
    fn test_range_empty_index() {
        // Index with 0 records
        let index = index_with_n_records(0, 3);
        assert_eq!(index.total_records, 0);

        // Any query should return None
        assert_eq!(index.range(1, 1), None);
        assert_eq!(index.range(0, 1), None);
    }

    #[test]
    fn test_range_at_exact_chunk_boundaries() {
        // 9 records, nth=3
        // Entries at total_records: [0, 3, 6, 9] with bytes [0, 102, 204, 306]
        let index = index_with_n_records(9, 3);
        assert_eq!(index.total_records, 9);
        assert_eq!(index.entries.len(), 4);

        // Query exactly at first chunk boundary: record 3 (1-indexed)
        // Record 3 is at position 2 in chunk [1,2,3], so leading_records = 2
        let range = index.range(3, 3).unwrap();
        assert_eq!(range.start_byte, 0);
        assert_eq!(range.end_byte, 102);
        assert_eq!(range.leading_records, 2);
        assert_eq!(range.trailing_records, 0);
        assert_eq!(range.selected_records(), 1);

        // Query at second chunk boundary: record 6 (1-indexed)
        // Record 6 is at position 2 in chunk [4,5,6], so leading_records = 2
        let range = index.range(6, 6).unwrap();
        assert_eq!(range.start_byte, 102);
        assert_eq!(range.end_byte, 204);
        assert_eq!(range.leading_records, 2);
        assert_eq!(range.trailing_records, 0);
        assert_eq!(range.selected_records(), 1);

        // Query record right after chunk boundary: record 4 (1-indexed)
        // Record 4 is at position 0 in chunk [4,5,6], so leading_records = 0
        let range = index.range(4, 4).unwrap();
        assert_eq!(range.start_byte, 102);
        assert_eq!(range.end_byte, 204);
        assert_eq!(range.leading_records, 0);
        assert_eq!(range.trailing_records, 2);
        assert_eq!(range.selected_records(), 1);
    }

    #[test]
    fn test_range_spanning_multiple_chunks() {
        // 9 records, nth=3
        let index = index_with_n_records(9, 3);

        // Query spanning all chunks: [1, 9]
        let range = index.range(1, 9).unwrap();
        assert_eq!(range.start_byte, 0);
        assert_eq!(range.end_byte, 306);
        assert_eq!(range.leading_records, 0);
        assert_eq!(range.trailing_records, 0);
        assert_eq!(range.selected_records(), 9);

        // Query spanning middle: [2, 8]
        let range = index.range(2, 8).unwrap();
        assert_eq!(range.start_byte, 0);
        assert_eq!(range.end_byte, 306);
        assert_eq!(range.leading_records, 1);
        assert_eq!(range.trailing_records, 1);
        assert_eq!(range.selected_records(), 7);
    }

    #[test]
    fn test_range_last_record_non_aligned() {
        // 8 records, nth=3 (same as index())
        // Entries at total_records: [0, 3, 6, 8] with bytes [0, 102, 204, 272]
        // The final chunk [7, 8] only has 2 records, not 3
        let index = index();

        // Query the very last record: [8, 8]
        let range = index.range(8, 8).unwrap();
        assert_eq!(range.start_byte, 204);
        assert_eq!(range.end_byte, 272);
        assert_eq!(range.leading_records, 1); // skip record 7
        assert_eq!(range.trailing_records, 0);
        assert_eq!(range.selected_records(), 1);

        // Query second-to-last: [7, 7]
        let range = index.range(7, 7).unwrap();
        assert_eq!(range.start_byte, 204);
        assert_eq!(range.end_byte, 272);
        assert_eq!(range.leading_records, 0);
        assert_eq!(range.trailing_records, 1);
        assert_eq!(range.selected_records(), 1);

        // Query both records in final chunk: [7, 8]
        let range = index.range(7, 8).unwrap();
        assert_eq!(range.start_byte, 204);
        assert_eq!(range.end_byte, 272);
        assert_eq!(range.leading_records, 0);
        assert_eq!(range.trailing_records, 0);
        assert_eq!(range.selected_records(), 2);
    }

    #[test]
    fn test_range_large_nth_few_records() {
        // 5 records, nth=100 (nth > total records)
        // Entries: [0, 5] at bytes [0, 170]
        let index = index_with_n_records(5, 100);
        assert_eq!(index.total_records, 5);
        assert_eq!(index.entries.len(), 2);

        // Query middle record
        let range = index.range(3, 3).unwrap();
        assert_eq!(range.start_byte, 0);
        assert_eq!(range.end_byte, 170);
        assert_eq!(range.leading_records, 2);
        assert_eq!(range.trailing_records, 2);
        assert_eq!(range.selected_records(), 1);

        // Query all records
        let range = index.range(1, 5).unwrap();
        assert_eq!(range.start_byte, 0);
        assert_eq!(range.end_byte, 170);
        assert_eq!(range.leading_records, 0);
        assert_eq!(range.trailing_records, 0);
        assert_eq!(range.selected_records(), 5);
    }

    #[test]
    fn test_range_nth_equals_one() {
        // 5 records, nth=1 (every record indexed)
        // Entries: [0, 1, 2, 3, 4, 5] at bytes [0, 34, 68, 102, 136, 170]
        let index = index_with_n_records(5, 1);
        assert_eq!(index.total_records, 5);
        assert_eq!(index.entries.len(), 6);

        // Each query should have exactly the bytes for those records
        for i in 1..=5 {
            let range = index.range(i, i).unwrap();
            assert_eq!(range.start_byte, (i - 1) * 34);
            assert_eq!(range.end_byte, i * 34);
            assert_eq!(range.leading_records, 0);
            assert_eq!(range.trailing_records, 0);
            assert_eq!(range.selected_records(), 1);
        }

        // Range [2, 4]
        let range = index.range(2, 4).unwrap();
        assert_eq!(range.start_byte, 34);
        assert_eq!(range.end_byte, 136);
        assert_eq!(range.leading_records, 0);
        assert_eq!(range.trailing_records, 0);
        assert_eq!(range.selected_records(), 3);
    }

    #[test]
    fn test_range_clamping_behavior() {
        // 8 records, nth=3
        let index = index();

        // Start before first record (0) should clamp to 1
        let range = index.range(0, 3).unwrap();
        assert_eq!(range.selected_records(), 3);

        // End after last record should clamp to total_records
        let range = index.range(7, 100).unwrap();
        assert_eq!(range.selected_records(), 2); // records 7 and 8

        // Both clamped
        let range = index.range(0, 100).unwrap();
        assert_eq!(range.selected_records(), 8); // all records
    }

    #[test]
    fn test_range_research_claim_validation() {
        // This test specifically validates the research.md claim about bug 1.1:
        // "If nth=3 and entries are indexed at [0, 3, 6, ...]:
        //  Request range(4, 5) should have leading_records=1"
        //
        // Our analysis: The claim appears incorrect. Record 4 (1-indexed) is
        // at position 0 in the chunk [4, 5, 6] (1-indexed), so leading_records=0.

        // 9 records, nth=3 to have clean chunks
        let index = index_with_n_records(9, 3);

        // range(4, 5): records 4 and 5 are in chunk [4, 5, 6]
        // Record 4 is at position 0, Record 5 is at position 1
        // So we need to skip 0 leading records and 1 trailing record
        let range = index.range(4, 5).unwrap();
        assert_eq!(range.start_byte, 102); // starts at byte 102 (after first 3 records)
        assert_eq!(range.end_byte, 204); // ends at byte 204 (after 6 records)
        assert_eq!(range.leading_records, 0); // NO leading records to skip
        assert_eq!(range.trailing_records, 1); // skip record 6
        assert_eq!(range.selected_records(), 2);
    }
}
