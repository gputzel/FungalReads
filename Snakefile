configfile: "config.json"

rule bbduk:
    input:
        "FASTQ/{sample}_L002_R1_001.fastq.gz"
    output:
        "output/FASTQ_unmapped/{sample}.fq"
    benchmark:
        "benchmarks/bbduk/{sample}.txt"
    threads: 4
    shell:
        "bbduk.sh -Xmx56g in={input} outu={output} ref=" + config["mouse_transcriptome"] + " k=27 "

rule fungal_star_index:
    input:
        "fungal_genomes/MUC.fasta"
    output:
        "fungal_genomes/MUC_star/Genome"
    threads: 4
    shell:
        "STAR --runThreadN {threads} --runMode genomeGenerate --genomeDir fungal_genomes/MUC_star --genomeFastaFiles {input}"

rule fungal_align:
    input:
        fastq="output/FASTQ_unmapped/{sample}.fq",
        index="fungal_genomes/MUC_star"
    output:
        "output/STAR_align_MUC/{sample}_Aligned.out.sam"
    threads: 4
    shell:
        "star --runThreadN {threads} --genomeDir {input.index} --readFilesIn {input.fastq} --outFileNamePrefix output/STAR_align_MUC/{wildcards.sample}_"
