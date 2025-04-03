version 1.0
workflow CompareFasta {
    input {
        File whole_read_bam
        File whole_read_bai
        File reference_fa
        File reference_gb
        String prefix
        String sampleid
        Int kmer_size

    }

    # call mitograph assembly
    call Filter {
        input:
            bam = whole_read_bam,
            bai = whole_read_bai,
            prefix = prefix
    }

    call Build {
        input:
            bam = Filter.mt_bam,
            reference = reference_fa,
            prefix = prefix,
            kmer_size = kmer_size,
            sampleid = sampleid
    }

    call Asm {
        input:
            graph_gfa = Build.graph,
            prefix = "mitograph",
            sampleid=sampleid
    }



    output {
        File mitograph_fasta = Asm.fasta
    }
}

task Filter {
    input {
        File bam
        File bai
        String prefix
    }

    command <<<
        set -euxo pipefail
        /mitograph/target/release/mitograph filter -i ~{bam} -c chrM -m ~{prefix}_mt.bam -n ~{prefix}_numts.bam
    >>>

    output {
        File mt_bam = "~{prefix}_mt.bam"
        File numts_bam = "~{prefix}_numts.bam"
    }

    runtime {
        docker: "hangsuunc/mitograph:v3"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 100 SSD"
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

        /mitograph/target/release/mitograph build -k ~{kmer_size} -r ~{reference} -o ~{sampleid}.~{prefix}.gfa ~{bam}

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
        File reference
        String reference_name
        String prefix
        String sampleid
        Int kmer_size
    }

    command <<<
        set -euxo pipefail
        /mitograph/target/release/mitograph call -g ~{graph_gfa} -r ~{reference_name} -k ~{kmer_size} -s ~{sampleid} -o ~{sampleid}.~{prefix}.vcf

    >>>

    output {
        # File graph = "~{prefix}.annotated.gfa"
        File vcf = "~{sampleid}.~{prefix}.vcf"
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
