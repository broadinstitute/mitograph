version 1.0

workflow VEPAnnotation {

    input {
        File vcf
        File tbi
        String prefix

    }


    call MergePerChrVcfWithBcftools as Merge { input:
        vcf_input = BcftoolsIndex.output_vcf,
        tbi_input = BcftoolsIndex.output_tbi,
        pref = prefix,
        threads_num = merge_num_threads,
        batch_size = batchsize
    }

    output {
        File merged_vcf = Merge.merged_vcf
        File merged_tbi = Merge.merged_tbi
     
    }
}


task VEPAnnotation {

    input{
        Array[File] vcf_input
        Array[File] tbi_input
        String pref
        Int threads_num
        Int batch_size
    }

    command <<<
        set -eux

        # we do single-sample phased VCFs localization ourselves
        mkdir -p ssp_vcfs
        time \
        gcloud storage cp ~{sep=" " vcf_input} /cromwell_root/ssp_vcfs/ &

        time \
        gcloud storage cp ~{sep=" " tbi_input} /cromwell_root/ssp_vcfs/ &
        wait

        # then merge, and safely assume all ssp-VCFs are sorted in the same order, on one chr
        cd ssp_vcfs
        ls *.vcf.gz | split -l ~{batch_size} - subset_vcfs

        cnt=0
        for i in subset_vcfs*;
        do
            bcftools merge \
                --threads 6 \
                --merge none \
                --force-single \
                -l $i \
                -O z \
                -o ~{pref}.merge.$i.vcf.gz &
            cnt=$((cnt+1))
            if [[ $cnt -eq 18 ]]; then cnt=0; wait; fi
        done
        wait
        for i in ~{pref}.merge.*.vcf.gz;
        do 
            bcftools index --threads 6 -t $i &
            cnt=$((cnt+1))
            if [[ $cnt -eq 18 ]]; then cnt=0; wait; fi
        done
        wait
        ls ~{pref}.merge.*.vcf.gz > merge.txt

        time \
        bcftools merge \
            --threads ~{threads_num} \
            --force-single \
            --merge none \
            -l merge.txt \
            -O z \
            -o ~{pref}.AllSamples.vcf.gz
        bcftools index --threads ~{threads_num} -t ~{pref}.AllSamples.vcf.gz
        # move result files to the correct location for cromwell to de-localize
        mv ~{pref}.AllSamples.vcf.gz ~{pref}.AllSamples.vcf.gz.tbi /cromwell_root/
    >>>

    output{
        File merged_vcf = "~{pref}.AllSamples.vcf.gz"
        File merged_tbi = "~{pref}.AllSamples.vcf.gz.tbi"
    }

    runtime {
        cpu: 8
        memory: "32 GiB"
        disks: "local-disk 3000 LOCAL"
        preemptible: 1
        maxRetries: 0
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.20"
    }
}
