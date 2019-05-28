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

rule download_SILVA_LSU:
    output:
        temp("resources/SILVA/SILVA_132_LSURef_tax_silva.fasta.gz")
    shell:
        "wget https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_LSURef_tax_silva.fasta.gz -O {output}"

rule SILVA_LSU_DNA:
    input:
        "resources/SILVA/SILVA_132_LSURef_tax_silva.fasta.gz"
    output:
        "resources/SILVA/SILVA_132_LSURef_tax_silva_DNA.fasta"
    shell:
        'cat {input} | seqkit replace -s -p "U" -r "T" - > {output}'

rule download_SILVA_SSU:
    output:
        temp("resources/SILVA/SILVA_132_SSURef_tax_silva.fasta.gz")
    shell:
        "wget https://www.arb-silva.de/fileadmin/silva_databases/release_132/Exports/SILVA_132_SSURef_tax_silva.fasta.gz -O {output}"

rule SILVA_SSU_DNA:
    input:
        "resources/SILVA/SILVA_132_SSURef_tax_silva.fasta.gz"
    output:
        "resources/SILVA/SILVA_132_SSURef_tax_silva_DNA.fasta"
    shell:
        'cat {input} | seqkit replace -s -p "U" -r "T" - > {output}'

rule filter_SSU:
    input:
        fastq="output/FASTQ_unmapped/{sample}.fq",
        db="resources/SILVA/SILVA_132_SSURef_tax_silva_DNA.fasta"
    output:
        "output/FASTQ_filter_SSU/{sample}.fq"
    threads: 4
    shell:
        "cat {input.fastq} | bbduk.sh -Xmx56g in=stdin.fq outu={output} ref={input.db} k=27 int=f"

rule filter_LSU:
    input:
        fastq="output/FASTQ_filter_SSU/{sample}.fq",
        db="resources/SILVA/SILVA_132_LSURef_tax_silva_DNA.fasta"
    output:
        "output/FASTQ_filter_LSU/{sample}.fq"
    threads: 4
    shell:
        "cat {input.fastq} | bbduk.sh -Xmx56g in=stdin.fq outu={output} ref={input.db} k=27 int=f"

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
        fastq="output/FASTQ_filter_LSU/{sample}.fq",
        index="fungal_genomes/MUC_star"
    output:
        "output/STAR_align_MUC/{sample}_Aligned.out.sam"
    threads: 4
    shell:
        "star --runThreadN {threads} --genomeDir {input.index} --readFilesIn {input.fastq} --outFilterMismatchNmax 2 --outFilterMultimapNmax 3 --outFileNamePrefix output/STAR_align_MUC/{wildcards.sample}_"

rule filter_masked_regions:
    input:
        bed="fungal_genomes/MUC_masked_regions.bed",
        sam="output/STAR_align_MUC/{sample}_Aligned.out.sam"
    output:
        "output/STAR_align_MUC_filter_masked/{sample}.sam"
    shell:
        "samtools view -h -L {input.bed} {input.sam} -U {output} > /dev/null"
