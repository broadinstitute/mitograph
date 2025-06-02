# Himito User Guide

## Overview

Himito is a specialized tool for building mitochondrial anchor-based graphical genomes from long-read sequencing data. It provides a comprehensive pipeline for mitochondrial genome analysis, including filtering nuclear mitochondrial sequences (NUMTs), assembling major haplotypes, calling variants, and analyzing methylation signals.

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Commands](#commands)
- [Workflows](#workflows)
- [Input Requirements](#input-requirements)
- [Output Files](#output-files)
- [Best Practices](#best-practices)
- [Troubleshooting](#troubleshooting)

## Installation

### Prerequisites

- Rust programming language
- Long-read sequencing data (BAM format)
- Mitochondrial reference genome (FASTA format)

### Step 1: Install Rust

### install rust
```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```
### Step 2: Install Himito

### install and run Himito
```
git clone https://github.com/broadinstitute/Himito.git
cd Himito
cargo build --release
```

## Quick Start

Here's a basic workflow to get started with Himito:

```
# filter NUMTs-derived reads
./target/release/Himito filter -i <input.bam> -c <chromosome in bam, e.g. "chrM"> -m <mt_output.bam> -n <numts_output.bam>

# construct graph
./target/release/Himito build -i <mt_output.bam> -k <kmer_size> -r <NC_012920.1.fasta> -o <output.gfa> 

# call variants from graph
./target/release/Himito call -g <output.gfa> -r <NC_012920.1.fasta> -k <kmer_size> -s <sampleid> -o <output.vcf>

# extract major haplotype from graph
./target/release/Himito asm -g <output.gfa>  -o <output.majorhaplotpe.fasta> -s <header string, e.g. "HG002 major haplotype">

# call methylation signals
./target/release/Himito methyl -g <output.annotated.gfa> -p <min_prob> -b <mt_test.bam> -o <methyl.bed>
```

## Commands
### filter - Remove NUMTs Reads
Separates genuine mitochondrial reads from nuclear mitochondrial sequences (NUMTs).

```
./target/release/Himito filter [OPTIONS]
```
Required Parameters:
- -i, --input <FILE>: Input BAM file containing long reads
- -c, --chromosome <STRING>: Mitochondrial chromosome name (e.g., "chrM", "MT", "M")
- -m, --mt-output <FILE>: Output BAM file for mitochondrial reads
- -n, --numts-output <FILE>: Output BAM file for NUMTs reads

Example:
```
./target/release/Himito filter -i HG002.bam -c chrM -m HG002_mt.bam -n HG002_numts.bam
```

### build - Construct Mitochondrial Graph
Creates a graph representation of the mitochondrial genome from filtered reads.

```
./target/release/Himito build [OPTIONS]
```
Required Parameters:
- -i, --input <FILE>: Input BAM file (mitochondrial reads only)
- -k, --kmer-size <INT>: K-mer size for graph construction
- -r, --reference <FILE>: Mitochondrial reference genome (FASTA)
- -o, --output <FILE>: Output graph file (GFA format)

Example:
```
./target/release/Himito build -i HG002_mt.bam -k 21 -r rCRS.fasta -o HG002_mt.gfa
```
### call - Variant Calling
Identifies homoplasmic and heteroplasmic variants from the mitochondrial graph.
