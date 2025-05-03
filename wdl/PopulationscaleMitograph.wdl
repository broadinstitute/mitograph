version 1.0
workflow PopulationMitograph {
    input {
        Array[File] fastas
        File reference_fa
        String prefix
        Int kmer_size

    }

    call MergeFasta {
        input:        
            fastas = fastas,
            prefix = prefix
    }

    call Build {
        input:
            bam = MergeFasta.merged_fasta,
            reference = reference_fa,
            prefix = prefix,
            kmer_size = kmer_size,
            sampleid = prefix
    }

    output {
        File graph = Build.graph
    }
}

task MergeFasta {

    meta {
        description : "cat fastas file into one"
    }

    parameter_meta {
    }

    input {
        Array[File] fastas
        String prefix
    }

    command <<<
        cat ~{sep=" " fastas} > ~{prefix}.fasta
    >>>

    output {
        File merged_fasta = "~{prefix}.fasta"
        
    }
    
    runtime {
        cpu: 1
        memory: 4 + " GiB"
        disks: "local-disk 375 LOCAL"
        bootDiskSizeGb: 10
        preemptible_tries:     1
        max_retries:           0
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.20"
    }

}

task Build {
    input {
        File bam
        File reference
        String prefix
        String sampleid
        Int kmer_size
    }

    command <<<
        set -euxo pipefail

        /mitograph/target/release/mitograph build -i ~{bam} -k ~{kmer_size} -r ~{reference} -o ~{sampleid}.~{prefix}.gfa 

    >>>

    output {
        File graph = "~{sampleid}.~{prefix}.gfa"
    }

    runtime {
        docker: "hangsuunc/mitograph:v3"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}

task Call {

    input {
        File graph_gfa
        File reference_fa
        String prefix
        String sampleid
        Int kmer_size
    }
    
    

    command <<<
        set -euxo pipefail
        /mitograph/target/release/mitograph call -g ~{graph_gfa} -r ~{reference_fa} -k ~{kmer_size} -s ~{sampleid} -o ~{sampleid}.~{prefix}.vcf
        ls

    >>>

    output {
        File graph = "~{sampleid}.~{prefix}.annotated.gfa"
        File matrix = "~{sampleid}.~{prefix}.matrix.csv"
        File vcf = "~{sampleid}.~{prefix}.vcf"
    }

    runtime {
        docker: "hangsuunc/mitograph:v3"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}

task Methyl {

    input {
        File graph_gfa
        File bam
        String sampleid
        String prefix = "methyl"
        Float min_prob
    }

    command <<<
        set -euxo pipefail
        /mitograph/target/release/mitograph methyl -g ~{graph_gfa} -p ~{min_prob} -b ~{bam} -o ~{sampleid}.~{prefix}.bed

    >>>

    output {
        File graph = "~{sampleid}.~{prefix}.methyl.gfa"
        File matrix = "~{sampleid}.~{prefix}.methylation_per_read.csv"
        File bed = "~{sampleid}.~{prefix}.bed"
    }

    runtime {
        docker: "hangsuunc/mitograph:v3"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}

task Asm {

    input {
        File graph_gfa
        String prefix
        String sampleid
    }

    command <<<
        set -euxo pipefail
        /mitograph/target/release/mitograph asm -g ~{graph_gfa} -s ~{sampleid} -o ~{sampleid}.~{prefix}.fasta

    >>>

    output {
        File fasta = "~{sampleid}.~{prefix}.fasta"
    }

    runtime {
        docker: "hangsuunc/mitograph:v3"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}
