# mitograph

Building Mitochondrial anchor-based graphical genome from long reads. Filter Numts reads, assemble major haplotypes, call homoplasmic and heteroplasmic variants, analyze methylation signals

## Usage
### install rust
```
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

### install and run mitograph
```
git clone https://github.com/broadinstitute/mitograph.git
cd mitograph
cargo build --release

# filter NUMTs-derived reads
./target/release/mitograph filter -i <input.bam> -c <chromosome in bam, e.g. "chrM"> -m <mt_output.bam> -n <numts_output.bam>

# construct graph
./target/release/mito_graph build -k <kmer_size> -r <NC_012920.1.fasta> -o <output.gfa> <mt_output.bam>

# call variants from graph
./target/release/mito_graph call -g <output.gfa> -r <NC_012920.1.fasta> -k <kmer_size> --output-file <output.vcf>

# extract major haplotype from graph
./target/release/mitograph asm -g <output.gfa>  -o <output.majorhaplotpe.fasta> -s <header string, e.g. "HG002 major haplotype">

# call methylation signals
./target/release/mitograph methyl -g <output.annotated.gfa> -b <mt_test.bam> -o <methyl.bed>
```
