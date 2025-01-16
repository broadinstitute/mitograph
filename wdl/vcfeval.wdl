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
    }
    call VCFEval {
        input:
            query_vcf = query_vcf,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            query_output_sample_name = query_output_sample_name,
            base_vcf = base_vcf,
            base_vcf_index = base_vcf_index,
            vcf_score_field = vcf_score_field
    }

    output {
        File summary = VCFEval.summary_statistics
        File roc = VCFEval.weighted_roc
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
            --output-mode combine \
            --decompose \
            --squash-ploidy \
            --sample ALT,ALT \
            --vcf-score-field ~{vcf_score_field}

        mkdir output_dir
        cp reg/summary.txt output_dir/
        cp reg/weighted_roc.tsv.gz output_dir/
        cp reg/*.vcf.gz* output_dir/
        cp output_dir/output.vcf.gz output_dir/~{query_output_sample_name}.vcf.gz
        cp output_dir/output.vcf.gz.tbi output_dir/~{query_output_sample_name}.vcf.gz.tbi

    
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
        File weighted_roc = "output_dir/weighted_roc.tsv.gz"
        File combined_output = "output_dir/~{query_output_sample_name}.vcf.gz"
        File combined_output_index = "output_dir/~{query_output_sample_name}.vcf.gz.tbi"
    }
}
