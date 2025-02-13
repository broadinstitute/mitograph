version 1.0

workflow HaplogroupDetection {
    input {
        File single_sample_vcf
        String sampleid
    }

    call Haplogrep {
        input:
            vcf = single_sample_vcf,
            prefix = sampleid
    }

    output {
        File haplogroup = Haplogrep.haplogroup
    }
}

task Haplogrep {
    input {
        File vcf
        String prefix
    }

    command <<<
        set -euxo pipefail
        curl -sL haplogrep.now.sh | bash
        ./haplogrep classify --in ~{vcf} --format vcf --out ~{prefix}.haplogroups.txt

    >>>

    output {
        File haplogroup = "~{prefix}.haplogroups.txt"
        
    }

    runtime {
        docker: "hangsuunc/mitograph:v2"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 100 SSD"
    }
}

