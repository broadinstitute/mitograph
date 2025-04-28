version 1.0

workflow mitograph_methyl {
    input {
        File whole_genome_bam
        File whole_genome_bai
        File reference_fa
        String prefix
        String reference_header
        String sampleid
        Int kmer_size

    }
    call Filter {
        input:
            bam = whole_genome_bam,
            bai = whole_genome_bai,
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

    call Call {
        input:
            graph_gfa = Build.graph,
            reference_name = reference_header,
            prefix = prefix,
            kmer_size = kmer_size,
            sampleid=sampleid
    }
    
    call Methyl {
        input:
            graph_gfa = Call.graph,
            bam = Filter.mt_bam,
            sampleid = sampleid,
            prefix = "methyl",
            min_prob = 0.5

    }


    output {
        File graph = Methyl.graph
        File vcf_file = Call.vcf
        File methyl_file = Methyl.bed
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
        String reference_name
        String prefix
        String sampleid
        Int kmer_size
    }

    command <<<
        set -euxo pipefail
        /mitograph/target/release/mitograph call -g ~{graph_gfa} -r ~{reference_name} -k ~{kmer_size} -s ~{sampleid} -o ~{sampleid}.~{prefix}.vcf
        ls

    >>>

    output {
        File graph = "~{sampleid}.~{prefix}.annotated.gfa"
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
        /mitograph/target/release/mitograph methyl -g ~{graph_gfa} -b ~{bam} -m ~{min_prob} -o ~{sampleid}.~{prefix}.bed

    >>>

    output {
        File graph = "~{sampleid}.~{prefix}.methyl.gfa"
        File bed = "~{sampleid}.~{prefix}.bed"
    }

    runtime {
        docker: "hangsuunc/mitograph:v3"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}