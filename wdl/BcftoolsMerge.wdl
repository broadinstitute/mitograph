version 1.0

workflow HierarchicalMergeVCFs {

    input {
        Array[File] vcfs
        Array[String] sample_ids

        String prefix
        Int merge_num_threads = 1
        Int batchsize
        Float threshold

    }
    Int data_length = length(vcfs)
    Array[Int] indexes= range(data_length)

    scatter (idx in indexes)  {
        File vcf = vcfs[idx]
        String sampleid = sample_ids[idx]
        call BcftoolsIndex { input:
            vcf = vcf,
            sampleid = sampleid,
            base_field = "FORMAT/HF",
            threshold = threshold
        }   
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

task BcftoolsIndex {
    input{
        File vcf
        String sampleid
        String base_field
        Float threshold
    }

    String base_info = "${base_field}\\>${threshold}"

    command <<<
        bcftools view ~{vcf} -O z -o ~{sampleid}.vcf.gz
        bcftools index -t ~{sampleid}.vcf.gz
        bcftools view -i  ~{base_info} ~{sampleid}.vcf.gz -O z -o ~{sampleid}.~{threshold}.vcf.gz
        bcftools index -t ~{sampleid}.~{threshold}.vcf.gz

    >>>
    output{
        File output_vcf = "~{sampleid}.~{threshold}.vcf.gz"
        File output_tbi = "~{sampleid}.~{threshold}.vcf.gz.tbi"
    }
    runtime {
        cpu: 1
        memory: "4 GiB"
        disks: "local-disk 375 LOCAL"
        preemptible: 0
        maxRetries: 0
        docker: "us.gcr.io/broad-dsp-lrma/lr-gcloud-samtools:0.1.20"
    }
}

task MergePerChrVcfWithBcftools {
    parameter_meta {
        vcf_input: {localization_optional: true}
        tbi_input: {localization_optional: true}
    }
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
