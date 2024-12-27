# best practices of
path_python = "/data/home/pengjia/miniconda3/envs/plot/bin/python/"
script = "/data/home/pengjia/HapTR/src/main.py"

bams = {
    "HiFi": {
        "LCL5": "/data/DATA/ChineseQuartet/ref_based_analysis/aligned_reads/ChineseQuartet/LCL5/ChineseQuartet.LCL5.GRCh38.HiFi.minimap2.bam",
        "LCL6": "/data/DATA/ChineseQuartet/ref_based_analysis/aligned_reads/ChineseQuartet/LCL6/ChineseQuartet.LCL6.GRCh38.HiFi.minimap2.bam",
        "LCL7": "/data/DATA/ChineseQuartet/ref_based_analysis/aligned_reads/ChineseQuartet/LCL7/ChineseQuartet.LCL7.GRCh38.HiFi.minimap2.bam",
        "LCL8": "/data/DATA/ChineseQuartet/ref_based_analysis/aligned_reads/ChineseQuartet/LCL8/ChineseQuartet.LCL8.GRCh38.HiFi.minimap2.bam",
    },
    "ONT": {
        "LCL5": "/data/DATA/ChineseQuartet/ref_based_analysis/aligned_reads/ChineseQuartet/LCL5/ChineseQuartet.LCL5.GRCh38.ONT.minimap2.bam",
        "LCL6": "/data/DATA/ChineseQuartet/ref_based_analysis/aligned_reads/ChineseQuartet/LCL6/ChineseQuartet.LCL6.GRCh38.ONT.minimap2.bam",
        "LCL7": "/data/DATA/ChineseQuartet/ref_based_analysis/aligned_reads/ChineseQuartet/LCL7/ChineseQuartet.LCL7.GRCh38.ONT.minimap2.bam",
        "LCL8": "/data/DATA/ChineseQuartet/ref_based_analysis/aligned_reads/ChineseQuartet/LCL8/ChineseQuartet.LCL8.GRCh38.ONT.minimap2.bam",
    }
}

path_regions = {"v0.1.3": "/data/home/pengjia/Project_Data/TR_analysis/regions/GRCh38/v0.1.3/work_dir/release/GRCh38.HapTR.bed"}
dir_work = "/data/home/pengjia/HapTR/test/"
path_ref = "/data/DATA/Reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome/GRCh38_full_analysis_set_plus_decoy_hla.fa"

wildcard_constraints:
    region="v0.1.3"
rule all:
    input:
        expand(dir_work + "{sample}.{tech}.{region}.{item}.tsv",
            sample=["LCL5"],tech="HiFi",region="v0.1.3",item=["train"])


rule run_test:
    input:
        bam=lambda wildcards: bams[wildcards.tech][wildcards.sample],
        bai=lambda wildcards: bams[wildcards.tech][wildcards.sample] + ".bai",
        bed=lambda wildcards: path_regions[wildcards.region],
        ref=path_ref
    output:
        dir_work + "{sample}.{tech}.{region}.{item}.tsv"

    threads: 40
    run:
        if item in ["train"]:
            shell("{path_python} {script} train -i {input.bam} -r {input.bed} -t {wildcards.tech} -o {output} -ref {input.ref}")
