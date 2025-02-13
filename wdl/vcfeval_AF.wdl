version 1.0

workflow Evaluation {
    input {
        File query_vcf
        File? query_vcf_index
        File reference_fa
        File reference_fai
        String query_output_sample_name
        File base_vcf
        File base_vcf_index
        String vcf_score_field
        Float threshold
    }
    call AF_Filter {
        input:
            query_vcf = query_vcf,
            base_vcf = base_vcf,
            threshold = threshold,
            prefix = query_output_sample_name

    }

    call GetIndels as Base_indel {
        input:
            query_vcf = base_vcf,
            query_tbi = base_vcf_index,
            prefix = query_output_sample_name
    }

    call GetIndels as Query_indel {
        input:
            query_vcf = AF_Filter.query_output_vcf,
            query_tbi = AF_Filter.query_output_tbi,
            prefix = query_output_sample_name
    }

    call GetSNPS as Base_snp {
        input:
            query_vcf = base_vcf,
            query_tbi = base_vcf_index,
            prefix = query_output_sample_name
    }

    call GetSNPS as Query_snp {
        input:
            query_vcf = AF_Filter.query_output_vcf,
            query_tbi = AF_Filter.query_output_tbi,
            prefix = query_output_sample_name
    }

    call VCFEval as SNPEval {
        input:
            query_vcf = Query_snp.vcf,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            query_output_sample_name = query_output_sample_name,
            base_vcf = Base_snp.vcf,
            base_vcf_index = Base_snp.tbi,
            vcf_score_field = vcf_score_field,
            prefix = "snps"
    }

    call VCFEval as IndelEval {
        input:
            query_vcf = Query_indel.vcf,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            query_output_sample_name = query_output_sample_name,
            base_vcf = Base_indel.vcf,
            base_vcf_index = Base_indel.tbi,
            vcf_score_field = vcf_score_field,
            prefix = "indels"
    }

    output {
        File summary_snp = SNPEval.summary_statistics
        File summary_indel = IndelEval.summary_statistics
    }

}

struct RuntimeAttributes {
    Int disk_size
    Int cpu
    Int memory
}

task AF_Filter {
    input {
        # Input VCF Files
        File query_vcf
        File base_vcf
        Float threshold
        String prefix
        String base_field
        String query_field

        # Runtime params
        Int? preemptible
        RuntimeAttributes runtimeAttributes = {"disk_size": ceil(2 * size(query_vcf, "GB")) + 50,
                                                  "cpu": 1, "memory": 4}
    }

    String query_info = "${query_field}\\>${threshold}"
    String base_info = "${base_field}\\>${threshold}"

    command <<<
        set -xeuo pipefail
        echo ~{query_info}
        echo ~{base_info}

        # Compress and Index vcf files
        bcftools view ~{query_vcf} -O z -o ~{query_vcf}.vcf.gz
        bcftools index -t ~{query_vcf}.vcf.gz
        bcftools view ~{base_vcf} -O z -o ~{base_vcf}.vcf.gz
        bcftools index -t ~{base_vcf}.vcf.gz

        # extract AF from query vcf file
        bcftools view -i  ~{query_info} ~{query_vcf}.vcf.gz -O z -o ~{prefix}.query.~{threshold}.vcf.gz
        bcftools index -t ~{prefix}.query.~{threshold}.vcf.gz

        # extract AF from base vcf file
        bcftools view -i  ~{base_info} ~{base_vcf}.vcf.gz -O z -o ~{prefix}.base.~{threshold}.vcf.gz
        bcftools index -t ~{prefix}.base.~{threshold}.vcf.gz
        
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.1-tmp"
        preemptible: select_first([preemptible, 0])
        disks: "local-disk " + runtimeAttributes.disk_size + " HDD"
        cpu: runtimeAttributes.cpu
        memory: runtimeAttributes.memory + " GB"
    }

    output {
        File base_output_vcf = "~{prefix}.base.~{threshold}.vcf.gz"
        File base_output_tbi = "~{prefix}.base.~{threshold}.vcf.gz.tbi"
        File query_output_vcf = "~{prefix}.query.~{threshold}.vcf.gz"
        File query_output_tbi = "~{prefix}.query.~{threshold}.vcf.gz.tbi"
    }
}

task GetSNPS {
    input {
        # Input VCF Files
        File query_vcf
        File query_tbi
        String prefix

        # Runtime params
        Int? preemptible
        RuntimeAttributes runtimeAttributes = {"disk_size": ceil(2 * size(query_vcf, "GB")) + 50,
                                                  "cpu": 1, "memory": 4}
    }

    command <<<
        set -xeuo pipefail

        # # Compress and Index vcf files
        # bcftools view ~{query_vcf} -O z -o ~{query_vcf}.vcf.gz
        # bcftools index -t ~{query_vcf}.vcf.gz
        # extract SNPs from vcf file
        bcftools view -v snps ~{query_vcf} -O z -o ~{prefix}.snps.vcf.gz
        bcftools index -t ~{prefix}.snps.vcf.gz
        
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.1-tmp"
        preemptible: select_first([preemptible, 0])
        disks: "local-disk " + runtimeAttributes.disk_size + " HDD"
        cpu: runtimeAttributes.cpu
        memory: runtimeAttributes.memory + " GB"
    }

    output {
        File vcf = "~{prefix}.snps.vcf.gz"
        File tbi = "~{prefix}.snps.vcf.gz.tbi"
    }
}

task GetIndels {
    input {
        # Input VCF Files
        File query_vcf
        File query_tbi
        String prefix

        # Runtime params
        Int? preemptible
        RuntimeAttributes runtimeAttributes = {"disk_size": ceil(2 * size(query_vcf, "GB")) + 50,
                                                  "cpu": 1, "memory": 4}
    }

    command <<<
        set -xeuo pipefail

        # Compress and Index vcf files
        # bcftools view ~{query_vcf} -O z -o ~{query_vcf}.vcf.gz
        # bcftools index -t ~{query_vcf}.vcf.gz
        # extract Indels from vcf file
        bcftools view -v indels ~{query_vcf} -O z -o ~{prefix}.indels.vcf.gz
        bcftools index -t ~{prefix}.indels.vcf.gz
        
    >>>

    runtime {
        docker: "us.gcr.io/broad-dsde-methods/vcfeval_docker:v1.1-tmp"
        preemptible: select_first([preemptible, 0])
        disks: "local-disk " + runtimeAttributes.disk_size + " HDD"
        cpu: runtimeAttributes.cpu
        memory: runtimeAttributes.memory + " GB"
    }

    output {
        File vcf = "~{prefix}.indels.vcf.gz"
        File tbi = "~{prefix}.indels.vcf.gz.tbi"
    }
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
        String prefix

        # Runtime params
        Int? preemptible
        RuntimeAttributes runtimeAttributes = {"disk_size": ceil(2 * size(query_vcf, "GB") + 2 * size(base_vcf, "GB") + size(reference_fa, "GB")) + 50,
                                                  "cpu": 8, "memory": 16}
    }

    command <<<
        set -xeuo pipefail

        # Compress and Index vcf files
        bcftools view ~{query_vcf} -O z -o ~{query_vcf}.vcf.gz
        bcftools index -t ~{query_vcf}.vcf.gz
        
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
            -c ~{query_vcf}.vcf.gz \
            -o reg \
            -t rtg_ref \
            --squash-ploidy \
            --sample ALT,ALT \
            --vcf-score-field ~{vcf_score_field}

        mkdir output_dir
        cp reg/summary.txt output_dir/~{query_output_sample_name}.~{prefix}.summary.txt
        cp reg/weighted_roc.tsv.gz output_dir/
        cp reg/*.vcf.gz* output_dir/
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
        File summary_statistics = "output_dir/~{query_output_sample_name}.~{prefix}.summary.txt"
        File weighted_roc = "output_dir/weighted_roc.tsv.gz"
        Array[File] combined_output = glob("output_dir/*.vcf.gz")
        Array[File] combined_output_index = glob("output_dir/*.vcf.gz.tbi")
    }
}
