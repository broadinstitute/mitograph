use std::path::PathBuf;

use clap::{Parser, Subcommand};

mod build;

#[derive(Debug, Parser)]
#[clap(name = "mito_graph")]
#[clap(about = "Analysis of mitochondrial genome using long reads.", long_about = None)]

struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

const DEFAULT_KMER_SIZE: usize = 17;

#[derive(Debug, Subcommand)]

enum Commands {
    /// Build series-parallel graph from long-read data in multi-sample FASTA file.
    #[clap(arg_required_else_help = true)]
    Build {
        /// Output path for series-parallel graph.
        #[clap(short, long, value_parser, default_value = "/dev/stdout")]
        output: PathBuf,

        /// Kmer-size
        #[clap(short, long, value_parser, default_value_t = DEFAULT_KMER_SIZE)]
        kmer_size: usize,

        /// Name of sequence to use as reference.
        #[clap(short, long, value_parser, required = true)]
        reference_path: PathBuf,

        /// bam or fasta file with reads spanning locus of interest.
        #[clap(required = true, value_parser)]
        read_path: PathBuf,
    },
}

fn main() {
    let args = Cli::parse();
    match args.command {
        Commands::Build {
            output,
            kmer_size,
            read_path,
            reference_path,
        } => {
            build::start(&output, kmer_size, &read_path, &reference_path);
        }
    }
}
