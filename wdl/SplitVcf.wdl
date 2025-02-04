version 1.0

workflow SplitCohortVcf {

    input {
        File joint_short_vcf
        File joint_short_vcf_tbi
        String sampleid
        String region
        Int memory
    }

    call SplitVcf { input:
        joint_vcf = joint_short_vcf,
        joint_vcf_tbi = joint_short_vcf_tbi,
        sampleid = sampleid,
        region = region, 
        memory = memory
    }

    output {
        File splitted_vcf = SplitVcf.vcf
        File splitted_tbi = SplitVcf.tbi
    }
}


task SplitVcf {

    input {
        File joint_vcf
        File joint_vcf_tbi
        String sampleid
        String region
        Int memory
    }

    command <<<
        set -euxo pipefail
        bcftools view -r ~{region} -s ~{sampleid} ~{joint_vcf} --write-index -O z -o ~{sampleid}.vcf.gz
        bcftools norm -m-any --do-not-normalize ~{sampleid}.vcf.gz \
            --write-index -Oz -o ~{sampleid}.normed.vcf.gz
        bcftools view -e 'GT="0/0" || GT="0|0" || GT~"0/." || GT~"./0" || GT~".|0" || GT~"0|."' ~{sampleid}.normed.vcf.gz -O z -o ~{sampleid}.splitted.vcf.gz 
        bcftools index -t ~{sampleid}.splitted.vcf.gz

    >>>

    output {
        File vcf = "~{sampleid}.splitted.vcf.gz"
        File tbi = "~{sampleid}.splitted.vcf.gz.tbi"
        # Array[File] vcf_by_sample_tbi = glob("output/*vcf.gz.tbi")

    }
    ###################
    runtime {
        cpu: 4
        memory: memory + " GiB"
        disks: "local-disk 375 LOCAL"
        bootDiskSizeGb: 10
        preemptible_tries:     1
        max_retries:           0
        docker:"us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.20"
    }
}