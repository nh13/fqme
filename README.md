# FASTQ ME

[![Check](https://github.com/nh13/fqme/actions/workflows/build_and_test.yml/badge.svg?branch=main)](https://github.com/nh13/fqme/actions/workflows/build_and_test.yml)
[![License](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/nh13/fqme/blob/main/LICENSE)
[![Language](http://img.shields.io/badge/language-rust-brightgreen.svg)](http://www.https://www.rust-lang.org/)

***Please do not use this for any serious work until unit tests have been completed***

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

