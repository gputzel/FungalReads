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

rule mask_fungal_genomes:
    input:
        "fungal_genomes/MUC.fasta"
    output:
        "fungal_genomes/MUC_masked_regions.txt"
    shell:
        "dustmasker -in {input} -out {output} -outfmt acclist"

rule masked_region_bed:
    input:
        "fungal_genomes/MUC_masked_regions.txt"
    output:
        "fungal_genomes/MUC_masked_regions.bed"
    shell:
        "cat {input} | sed 's/>//g' > {output}"

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

rule filter_masked_regions:
    input:
        bed="fungal_genomes/MUC_masked_regions.bed",
        sam="output/STAR_align_MUC/{sample}_Aligned.out.sam"
    output:
        "output/STAR_align_MUC_filter_masked/{sample}.sam"
    shell:
        "samtools view -h -L {input.bed} {input.sam} -U {output} > /dev/null"
