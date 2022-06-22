# FASTQ ME

## Quickstart

Compress and index (two-ways) the FASTQ with [`bgzip`](http://www.htslib.org/doc/bgzip.html):

```bash
cat test.fastq | fqme index --output test.fastq.gz.fqi -n 100 | bgzip -c -i --index-name test.fastq.gz.gzi > test.fastq.gz
```

Extract entries:
```bash
fqme extract --input test.fastq.gz -s 100 -e 102
```

## Help

```bash
fqme --help
```

## To Build

```bash
cargo build --release
```

The executable is located in:

```bash
target/release/fqme
```

## To Test

```bash
cargo test
```

