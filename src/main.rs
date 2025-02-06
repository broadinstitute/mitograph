use std::path::PathBuf;

use clap::{Parser, Subcommand};

mod build;
mod agg;
mod call;
mod filter;
mod asm;

#[derive(Debug, Parser)]
#[clap(name = "mitograph")]
#[clap(about = "Analysis of mitochondrial genome using long reads.", long_about = None)]

struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

const DEFAULT_KMER_SIZE: usize = 21;


#[derive(Debug, Subcommand)]

enum Commands {
    /// Filter reads derived from Numts
    #[clap(arg_required_else_help = true)]
    Filter {
        /// input path for bam file.
        #[clap(short, long, value_parser, required = true)]
        input_bam: PathBuf,

        /// contig name in the bam file
        #[clap(short, long, value_parser)]
        chromo: String,

        /// output path for mtDNA bam file
        #[clap(short, long, value_parser, required = true)]
        mt_output: PathBuf,

        /// output path for numts bam file
        #[clap(short, long, required = true, value_parser)]
        numts_output: PathBuf,
    },

    /// Build graph from long-read data in FASTA or Bam file.
    #[clap(arg_required_else_help = true)]
    Build {
        /// Output path for anchor graph.
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

    ///Call Variants from Anchor Graph
    #[clap(arg_required_else_help = true)]
    Call {
        /// path for anchor graph.
        #[clap(short, long, value_parser)]
        graphfile: PathBuf,

        /// Name of sequence to use as reference.
        #[clap(short, long, value_parser, required = true)]
        ref_strain: String,

        /// Kmer-size
        #[clap(short, long, value_parser, default_value_t = DEFAULT_KMER_SIZE)]
        k: usize,
        /// max Length to do alignment
        #[clap(short, long, value_parser, default_value_t = 3000)]
        length_max: usize, 

        /// minimal allele count for variants
        #[clap(short, long, value_parser, default_value_t = 1)]
        minimal_ac: usize,
        /// output file name
        #[clap(short, long, value_parser, required = true)]
        output_file: String,

        /// sample name of the bam file
        #[clap(short, long, value_parser, required = true)]
        sample_id: String,

    },

    /// Extract Major Haplotype as Fasta file from Graph
    #[clap(arg_required_else_help = true)]
    Asm {
        /// path for anchor graph.
        #[clap(short, long, value_parser)]
        graphfile: PathBuf,
        /// path for output fasta file
        #[clap(short, long, value_parser)]
        outputfile: PathBuf,

        /// header for the major haplotype, usually the sample name
        #[clap(short, long, value_parser)]
        sample: String,

    }
}

fn main() {
    let args = Cli::parse();
    match args.command {
        Commands::Filter {
            input_bam,
            chromo,
            mt_output,
            numts_output,
        } => {
            filter::start(&input_bam, &chromo, &mt_output, &numts_output);
        }

        Commands::Build {
            output,
            kmer_size,
            read_path,
            reference_path,
        } => {
            build::start(&output, kmer_size, &read_path, &reference_path);
        }

        Commands::Call {
            graphfile,
            ref_strain,
            k,
            length_max,
            minimal_ac,
            output_file,
            sample_id
        } => {
            call::start(&graphfile, &ref_strain, k, length_max, minimal_ac, &output_file, &sample_id);
        }

        Commands::Asm {
            graphfile,
            outputfile,
            sample
        } => {
            asm::start(&graphfile, &outputfile, &sample);
        }
    }
}
