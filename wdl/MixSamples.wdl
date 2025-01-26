version 1.0

workflow MixSamples {
    input {
        File first_donor_bam
        File first_donor_bai
        File second_donor_bam
        File second_donor_bai
        Float first_proportion
        String sampleid
    }
    call mix_sample {
        input:
            first_donor_bam = first_donor_bam,
            first_donor_bai = first_donor_bai,
            second_donor_bam = second_donor_bam,
            second_donor_bai = second_donor_bai,
            first_proportion = first_proportion,
            prefix = sampleid
    }

    output {
        File bam = mix_sample.merged_bam
        File bai = mix_sample.merged_bai
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

task mix_sample {

    meta {
        description : "sample reads to a specific proportion and mix the two bam"
    }

    parameter_meta {
        first_donor_bam: {
            description: "first donor bam to subset",
            localization_optional: false
        }
        first_donor_bai:    "index for first donor bam file"
        second_donor_bam: {
            description: "second donor bam to subset",
            localization_optional: false
        }
        second_donor_bai:    "index for second donor bam file"
        first_proportion:  "proportion of reads to select in the first bam"
        second_proportion: "proportion of reads to select in the second bam, should add up to 1"
        prefix: "prefix for output bam and bai file names"
        runtime_attr_override: "Override the default runtime attributes."
    }

    input {
        File first_donor_bam
        File first_donor_bai
        File second_donor_bam
        File second_donor_bai
        Float first_proportion
        Float second_proportion = 1 - first_proportion
        String prefix

        RuntimeAttr? runtime_attr_override
    }


    command <<<
        set -euxo pipefail

        samtools view -s ~{first_proportion} -b ~{first_donor_bam} -o ~{prefix}.first.bam

        samtools view -s ~{second_proportion} -b ~{second_donor_bam} -o ~{prefix}.second.bam

        samtools merge -f ~{prefix}.merged.bam ~{prefix}.first.bam ~{prefix}.second.bam

        samtools sort -o ~{prefix}.merged.sorted.bam ~{prefix}.merged.bam
        
        samtools index ~{prefix}.merged.sorted.bam

    >>>

    output {
        File merged_bam = "~{prefix}.merged.sorted.bam"
        File merged_bai = "~{prefix}.merged.sorted.bam.bai"
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
