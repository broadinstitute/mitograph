version 1.0

workflow MixSamples {
    input {
        File wholegenome_bam
        File wholegenome_bai
        File truth_vcf
        File truth_tbi

        File reference_fa
        File reference_fai
        Int desiredCoverage
        Int kmer_size = 21

        String sampleid
        String region = "chrM"
        String reference_header
        String vcf_score_field_mitograph
        String query_field_mitograph


    }

    call CalculateCoverage {
        input:
            bam = wholegenome_bam,
            bai = wholegenome_bai,
            locus = region,
            prefix = sampleid
    }

    call downsampleBam {input:
        input_bam = wholegenome_bam,
        input_bam_bai = wholegenome_bai,
        basename = sampleid,
        desiredCoverage = desiredCoverage,
        currentCoverage = CalculateCoverage.coverage,
        preemptible_tries = 0
    }

    call Filter {
        input:
            bam = downsampleBam.downsampled_bam,
            bai = downsampleBam.downsampled_bai,
            prefix = sampleid
    }

    call Build {
        input:
            bam = Filter.mt_bam,
            reference = reference_fa,
            prefix = sampleid,
            kmer_size = kmer_size,
            sampleid = sampleid
    }

    call Call {
        input:
            graph_gfa = Build.graph,
            reference = reference_fa,
            reference_name = reference_header,
            prefix = sampleid,
            kmer_size = kmer_size,
            sampleid=sampleid
    }

    call VCFEval as Mitograph_Eval {
        input:
            query_vcf = Call.vcf,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            query_output_sample_name = sampleid,
            base_vcf = truth_vcf,
            base_vcf_index = truth_tbi,
            vcf_score_field = vcf_score_field_mitograph,
            query_field = query_field_mitograph,
            threshold = 0.95
    }

    call CalculateVariantNumber{
        input:
            vcf = Call.vcf
    }

    output {
        File mitograph_summary_file = Mitograph_Eval.summary_statistics
        Float mitograph_variantnumber = CalculateVariantNumber.count
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

task CalculateCoverage {

    meta {
        description : "Subset a BAM file to a specified locus, and calculate coverage"
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

        samtools depth -r ~{locus} ~{bam} | awk '{sum+=$3} END {print sum/NR}' > coverage.txt
    >>>

    output {
        Float coverage = read_float("coverage.txt")
        # File subset_bam = "~{prefix}.bam"
        # File subset_bai = "~{prefix}.bam.bai"
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

task downsampleBam {
  input {
    File input_bam
    File input_bam_bai
    String basename
    Int desiredCoverage
    Float currentCoverage
    Float scalingFactor = desiredCoverage / currentCoverage

    Int? preemptible_tries
  }

  meta {
    description: "Uses Picard to downsample to desired coverage based on provided estimate of coverage."
  }
  parameter_meta {
    basename: "Input is a string specifying the sample name which will be used to locate the file on gs."
    downsampled_bam: "Output is a bam file downsampled to the specified mean coverage."
    downsampled_bai: "Output is the index file for a bam file downsampled to the specified mean coverage."
    desiredCoverage: "Input is an integer of the desired approximate coverage in the output bam file."
  }
  command <<<
    set -eo pipefail
    gatk DownsampleSam -I ~{input_bam} -O ~{basename}_~{desiredCoverage}x.bam -R 7 -P ~{scalingFactor} -S ConstantMemory --VALIDATION_STRINGENCY LENIENT --CREATE_INDEX true


  >>>
  runtime {
    preemptible: select_first([preemptible_tries, 5])
    memory: "8 GB"
    cpu: "2"
    disks: "local-disk 500 HDD"
    docker: "us.gcr.io/broad-gatk/gatk"
  }
  output {
    File downsampled_bam = "~{basename}_~{desiredCoverage}x.bam"
    File downsampled_bai = "~{basename}_~{desiredCoverage}x.bai"
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
        docker: "hangsuunc/mitograph:v2"
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
        docker: "hangsuunc/mitograph:v2"
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
        docker: "hangsuunc/mitograph:v2"
        memory: "2 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}

struct RuntimeAttributes {
    Int disk_size
    Int cpu
    Int memory
}

task VCFEval {
    input {
        # Input VCF Files
        File query_vcf
        File reference_fa
        File reference_fai
        String query_output_sample_name
        File base_vcf
        File base_vcf_index
        String vcf_score_field
        String query_field
        Float threshold

        # Runtime params
        Int? preemptible
        RuntimeAttributes runtimeAttributes = {"disk_size": ceil(2 * size(query_vcf, "GB") + 2 * size(base_vcf, "GB") + size(reference_fa, "GB")) + 50,
                                                  "cpu": 8, "memory": 16}
    }
    String query_info = "${query_field}\\>${threshold}"

    command <<<
        set -xeuo pipefail

        # Compress and Index vcf files
        bcftools view ~{query_vcf} -O z -o ~{query_vcf}.vcf.gz
        bcftools index -t ~{query_vcf}.vcf.gz

        # extract AF from query vcf file
        bcftools view -i  ~{query_info} ~{query_vcf}.vcf.gz -O z -o ~{query_output_sample_name}.query.~{threshold}.vcf.gz
        bcftools index -t ~{query_output_sample_name}.query.~{threshold}.vcf.gz
        
        # split multiallelic sites in the base_vcf
        bcftools norm \
                -f ~{reference_fa} \
                -m -both ~{base_vcf} \
                -O z \
                -o ~{base_vcf}.normed.vcf.gz 
        bcftools index -t ~{base_vcf}.normed.vcf.gz 
        
        # rtg vcfeval
        rtg format -o rtg_ref ~{reference_fa}
        rtg vcfeval \
            -b ~{base_vcf}.normed.vcf.gz  \
            -c ~{query_output_sample_name}.query.~{threshold}.vcf.gz \
            -o reg \
            -t rtg_ref \
            --squash-ploidy \
            --sample ALT,ALT \
            --vcf-score-field ~{vcf_score_field}

        mkdir output_dir
        cp reg/summary.txt output_dir/
        # cp reg/weighted_roc.tsv.gz output_dir/
        # cp reg/*.vcf.gz* output_dir/
        # cp output_dir/output.vcf.gz output_dir/~{query_output_sample_name}.vcf.gz
        # cp output_dir/output.vcf.gz.tbi output_dir/~{query_output_sample_name}.vcf.gz.tbi

    
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.1-tmp"
        preemptible: select_first([preemptible, 0])
        disks: "local-disk " + runtimeAttributes.disk_size + " HDD"
        cpu: runtimeAttributes.cpu
        memory: runtimeAttributes.memory + " GB"
    }

    output {
        File summary_statistics = "output_dir/summary.txt"
    }
}

task merge_vcf {

    meta {
        description : "merge two homoplasmic vcf as truth set"
    }

    input {
        File first_donor_vcf
        File first_donor_tbi
        File second_donor_vcf
        File second_donor_tbi
        String prefix
        

        RuntimeAttr? runtime_attr_override
    }


    command <<<
        set -euxo pipefail



        bcftools merge ~{first_donor_vcf} ~{second_donor_vcf} -O z -o ~{prefix}.merged.vcf.gz
        bcftools index -t ~{prefix}.merged.vcf.gz
        

    >>>

    output {
        File truth_vcf = "~{prefix}.merged.vcf.gz"
        File truth_tbi = "~{prefix}.merged.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "broadinstitute/gatk:4.6.1.0"
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


task CalculateVariantNumber {

    meta {
        description : "calculate variant number in a vcf file"
    }

    parameter_meta {
   }

    input {
        File vcf

        RuntimeAttr? runtime_attr_override
    }


    command <<<
        set -euxo pipefail

        bcftools view -H ~{vcf} | wc -l > count.txt
    >>>

    output {
        Float count = read_float("count.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            10,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.20"
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

