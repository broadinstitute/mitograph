version 1.0

workflow coverage_hap_calculation {
    input {
        File whole_genome_bam
        File whole_genome_bai
        String region
        String sampleid

    }
    call SubsetandSplit {
        input:
            bam = whole_genome_bam,
            bai = whole_genome_bai,
            locus = region,
            prefix = sampleid
    }

    output {
        Float coverage_hap1 = SubsetandSplit.coverage_hap1
        Float coverage_hap2 = SubsetandSplit.coverage_hap2
    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

struct DataTypeParameters {
    Int num_shards
    String map_preset
}

task SubsetandSplit {

    meta {
        description : "Subset a BAM file to a specified locus, and split haplotype and calculate coverage"
    }

    parameter_meta {
        bam: {
            description: "bam to subset",
            localization_optional: true
        }
        bai:    "index for bam file"
        locus:  "genomic locus to select"
        prefix: "prefix for output bam and bai file names"
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        File bam
        File bai
        String locus
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail

        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)

        samtools view -bhX ~{bam} ~{bai} ~{locus} > ~{prefix}.bam
        samtools index ~{prefix}.bam

        samtools view -b -d HP:1 ~{prefix}.bam | samtools depth -r ~{locus} | awk '{sum+=$3} END {print sum/NR}' > coverage_1.txt
        samtools view -b -d HP:2 ~{prefix}.bam | samtools depth -r ~{locus} | awk '{sum+=$3} END {print sum/NR}' > coverage_2.txt
    >>>

    output {
        Float coverage_hap1 = read_float("coverage_1.txt")
        Float coverage_hap2 = read_float("coverage_2.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-utils:0.1.9"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
