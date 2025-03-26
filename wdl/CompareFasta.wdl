version 1.0
workflow CompareFasta {
    input {
        # File whole_hap1_bam
        # File whole_hap1_bai
        # File whole_hap2_bam
        # File whole_hap2_bai
        File whole_read_bam
        File whole_read_bai
        File reference_fa
        File reference_gb
        String prefix
        String reference_header
        String sampleid
        Int kmer_size

    }

    # # extract truth assembly
    # call SubsetBam as Hap1 {
    #     input:
    #         bam = whole_hap1_bam,
    #         bai = whole_hap1_bai,
    #         locus = "chrM",
    #         prefix = sampleid + "hap1"
    # }

    # call SubsetBam as Hap2 {
    #     input:
    #         bam = whole_hap2_bam,
    #         bai = whole_hap2_bai,
    #         locus = "chrM",
    #         prefix = sampleid + "hap2"
    # }

    # call MergeFasta {
    #     input:        
    #         fasta1 = Hap1.subset_fasta,
    #         fasta2 = Hap2.subset_fasta,
    #         locus = "chrM",
    #         prefix = sampleid
    # }

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

    # call mitohifi assembly
    call SubsetBam as SubsetReads {
        input:
            bam = whole_read_bam,
            bai = whole_read_bai,
            locus = "chrM",
            prefix = sampleid + "reads"
    }

    call MitoHifiAsm {
        input:
            reads = SubsetReads.subset_fastq,
            reffa = reference_fa,
            refgb = reference_gb

    }


    output {
        # File truth_fasta = MergeFasta.merged_fasta
        File mitograph_fasta = Asm.fasta
        File mitohifi_fasta = MitoHifiAsm.final_fa
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

task SubsetBam {

    meta {
        description : "Subset a BAM file to a specified locus, and extract sequence into fa"
    }

    parameter_meta {
        bam: {
            description: "bam to subset",
            localization_optional: true
        }
        bai:    "index for bam file"
        locus:  "genomic locus to select"
        prefix: "prefix for output bam and bai file names"
    }

    input {
        File bam
        File bai
        String locus
        String prefix
    }

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        samtools view -bhX ~{bam} ~{bai} ~{locus} > ~{prefix}.bam
        samtools index ~{prefix}.bam

        samtools fasta ~{prefix}.bam > ~{prefix}.fasta
        samtools fastq ~{prefix}.bam > ~{prefix}.fastq
    >>>

    output {
        File subset_fasta = "~{prefix}.fasta"
        File subset_fastq = "~{prefix}.fastq"
        File subset_bam = "~{prefix}.bam"
        File subset_bai = "~{prefix}.bam.bai"
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

task MergeFasta {

    meta {
        description : "cat two fasta file into one"
    }

    parameter_meta {
    }

    input {
        File fasta1
        File fasta2
        String locus
        String prefix
    }

    command <<<
        cat ~{fasta1} ~{fasta2} > ~{prefix}~{locus}.fasta
    >>>

    output {
        File merged_fasta = "~{prefix}~{locus}.fasta"
        
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

task MitoHifiAsm{
    input{
        File reads
        File reffa
        File refgb
        Int num_cpus
    }

    #Int disk_size = 50 + ceil(2 * size(reads, "GiB"))

    command <<<
        #set -euxo pipefail
        mitohifi.py -r ~{reads} -f ~{reffa} -g ~{refgb} -t ~{num_cpus} -o 1 
        ls -l
    >>>

    output{
        File final_fa="final_mitogenome.fasta"
        File final_stats="contigs_stats.tsv"
        File mapping_file = "coverage_mapping/HiFi-vs-final_mitogenome.sorted.bam"
        

    }
    runtime {
        cpu: 4
        memory: "16 GiB"
        disks: "local-disk " + "256" + " HDD" #"local-disk 100 HDD"
        bootDiskSizeGb: 64
        preemptible: 2
        maxRetries: 1
        docker: "ghcr.io/marcelauliano/mitohifi:master"
        #docker:"biocontainers/mitohifi:3.0.0_cv1"
    }    
}