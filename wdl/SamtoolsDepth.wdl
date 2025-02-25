version 1.0

workflow CalculateDepth {

    input {
        Array[File] bams
        Array[File] bais
        Array[String] sample_ids

    }
    Int data_length = length(bams)
    Array[Int] indexes= range(data_length)

    scatter (idx in indexes)  {
        File bam = bams[idx]
        File bai = bais[idx]
        String sampleid = sample_ids[idx]

        call SamtoolsDepth { input:
            bam = bam,
            bai = bai,
            sampleid = sampleid,
            region = "chrM"
        }   
    }

    output {
        Array[File] depthfile = SamtoolsDepth.depth_file
     
    }
}

task SamtoolsDepth {
    input{
        File bam
        File bai
        String sampleid
        String region
       
    }


    command <<<
        samtools depth \
            -a \
            -r ~{region} \
            ~{bam} > ~{sampleid}.depth.txt
        
    >>>
    output{
        File depth_file = "~{sampleid}.depth.txt"
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

