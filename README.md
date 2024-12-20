# mitograph

Building Mitochondrial anchor-based graphical genome from long reads. Call very low frequent heteroplasmies from graph.

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

# filter reads from numts based on SA tag
./target/release/mitograph filter -i <input.bam> -c <contig_name> -m <mt_test.bam> -n <numts_test.bam>

# construct graph
./target/release/mito_graph build -k <kmer_size> -r <NC_012920.1.fasta> -o <output.gfa> <mt_test.bam>

# call variants from graph
./target/release/mito_graph call -g <output.gfa> -r <header of reference "NC_012920.1"> -k <kmer_size> --output-file <output.vcf>
```
