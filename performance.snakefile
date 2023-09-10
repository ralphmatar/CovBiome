SAMPLES = ['599', '48951', '128', '95', '50229', '34710', '168', '46183', '2', '41586', SRR14215324, SRR14215325, SRR14215327, SRR14215328, SRR14215331 , SRR14215332, SRR1421533]


rule all:
    input:
        "analysis/benchmarking/genomes/hybrid.fna.amb",
        "analysis/benchmarking/genomes/GRCh38.fna.amb",
        "analysis/benchmarking/genomes/cov_genome.fna.amb",
        "analysis/benchmarking/dummy/load_genome_dummy.txt",
        "analysis/benchmarking/dummy/unload_genome_dummy.txt",
        "analysis/benchmarking/dummy/load_genome_dummy1.txt",
        "analysis/benchmarking/dummy/unload_genome_dummy2.txt",
        "analysis/benchmarking/dummy/load_genome_dummy3.txt",
        "analysis/benchmarking/dummy/unload_genome_dummy4.txt",
        expand("analysis/benchmarking/covmap/{sample}.sam", sample=SAMPLES),
        expand("analysis/benchmarking/unmapped_cov/covum_{sample}_R1.fastq", sample=SAMPLES),
        expand("analysis/benchmarking/unmapped_cov/covum_{sample}_R2.fastq", sample=SAMPLES),
        expand("analysis/benchmarking/hmap/{sample}.sam", sample=SAMPLES),
        expand("analysis/benchmarking/unmapped_hum/humum_{sample}_R1.fastq", sample=SAMPLES),
        expand("analysis/benchmarking/unmapped_hum/humum_{sample}_R2.fastq", sample=SAMPLES),
        expand("analysis/benchmarking/hbmap/{sample}.sam", sample=SAMPLES),
        expand("analysis/benchmarking/unmapped_hyb/hyum_{sample}_R1.fastq", sample=SAMPLES),
        expand("analysis/benchmarking/unmapped_hyb/hyum_{sample}_R2.fastq", sample=SAMPLES)

        
rule index_concat:
    input:
        "analysis/benchmarking/genomes/hybrid.fna"
    output:
        "analysis/benchmarking/genomes/hybrid.fna.amb"
    benchmark:
        "analysis/benchmarking/index_concat.tsv"
    shell:
        """
        bwa index {input}
        """
 
 
rule index_human:
    input:
        "analysis/benchmarking/genomes/GRCh38.fna"
    output:
        "analysis/benchmarking/genomes/GRCh38.fna.amb"
    benchmark:
        "analysis/benchmarking/index_human.tsv"
    shell:
        """
        bwa index {input}
        """
        
        
rule index_cov:
    input:
        "analysis/benchmarking/genomes/cov_genome.fna"
    output:
        "analysis/benchmarking/genomes/cov_genome.fna.amb"
    benchmark:
        "analysis/benchmarking/index_cov.tsv"
    shell:
        """
        bwa index {input}
        """
        
        
rule load_cov_genome:
    input:
        "analysis/benchmarking/genomes/cov_genome.fna"
    output:
        "analysis/benchmarking/dummy/load_genome_dummy.txt"
    priority: 20
    shell:
        """
        bwa shm {input} && \
        touch {output}
        """
    
    
rule map2covid:
    input:
        genome = "analysis/benchmarking/genomes/cov_genome.fna",
        read1 = "data/trimmed_data/{sample}_read_1.fastq.gz",
        read2 = "data/trimmed_data/{sample}_read_2.fastq.gz"
    output:
        "analysis/benchmarking/covmap/{sample}.sam"
    priority: 18
    benchmark:
        "analysis/benchmarking/map2cov_{sample}.tsv"
    shell:
        """
        bwa mem \
            -t 16 \
            {input.genome} \
            {input.read1} \
            {input.read2} \
            -o {output}
        """
        
        
rule unload_cov_genome:
    input: 
        rules.load_cov_genome.output
    output:
        "analysis/benchmarking/dummy/unload_genome_dummy.txt"
    priority: 16
    shell:
        """
        bwa shm -d && \
        touch {input} && \
        touch {output}
        """
        

rule unmapped_cov:
    input:
        rules.map2covid.output
    output:
        read1 = "analysis/benchmarking/unmapped_cov/covum_{sample}_R1.fastq",
        read2 = "analysis/benchmarking/unmapped_cov/covum_{sample}_R2.fastq"
    priority: 14
    benchmark:
        "analysis/benchmarking/umcov_{sample}.tsv"
    shell:
        """
        samtools fastq \
            -@ 16 \
            -f 4 \
            -1 {output.read1} \
            -2 {output.read2} \
            -s /dev/null \
            -0 /dev/null \
            -n {input}
        """        

rule load_hum_genome:
    input:
        "analysis/benchmarking/genomes/GRCh38.fna"
    output:
        "analysis/benchmarking/dummy/load_genome_dummy1.txt"
    priority: 12
    shell:
        """
        bwa shm {input} && \
        touch {output}
        """
    
    
rule map2human:
    input:
        genome = "analysis/benchmarking/genomes/GRCh38.fna",
        read1 = "analysis/benchmarking/unmapped_cov/covum_{sample}_R1.fastq",
        read2 = "analysis/benchmarking/unmapped_cov/covum_{sample}_R2.fastq"
    output:
        "analysis/benchmarking/hmap/{sample}.sam"
    priority: 10
    benchmark:
        "analysis/benchmarking/map2hu_{sample}.tsv"
    shell:
        """
        bwa mem \
            -t 16 \
            {input.genome} \
            {input.read1} \
            {input.read2} \
            -o {output}
        """
        
        
rule unload_hum_genome:
    input: 
        rules.load_hum_genome.output
    output:
        "analysis/benchmarking/dummy/unload_genome_dummy2.txt"
    priority: 8
    shell:
        """
        bwa shm -d && \
        touch {input} && \
        touch {output}
        """      
        
        
rule unmapped_hum:
    input:
        rules.map2human.output
    output:
        read1 = "analysis/benchmarking/unmapped_hum/humum_{sample}_R1.fastq",
        read2 = "analysis/benchmarking/unmapped_hum/humum_{sample}_R2.fastq"
    priority: 6
    benchmark:
        "analysis/benchmarking/umhum_{sample}.tsv"
    shell:
        """
        samtools fastq \
            -@ 16 \
            -f 4 \
            -1 {output.read1} \
            -2 {output.read2} \
            -s /dev/null \
            -0 /dev/null \
            -n {input}
        """
        
        
rule load_hyb_genome:
    input:
        "analysis/benchmarking/genomes/hybrid.fna"
    output:
        "analysis/benchmarking/dummy/load_genome_dummy3.txt"
    priority: 4
    shell:
        """
        bwa shm {input} && \
        touch {output}
        """
    
    
rule map2hybrid:
    input:
        genome = "analysis/benchmarking/genomes/hybrid.fna",
        read1 = "data/trimmed_data/{sample}_read_1.fastq.gz",
        read2 = "data/trimmed_data/{sample}_read_2.fastq.gz"
    output:
        "analysis/benchmarking/hbmap/{sample}.sam"
    priority: 3
    benchmark:
        "analysis/benchmarking/map2hb_{sample}.tsv"
    shell:
        """
        bwa mem \
            -t 16 \
            {input.genome} \
            {input.read1} \
            {input.read2} \
            -o {output}
        """
        
        
rule unload_hyb_genome:
    input: 
        rules.load_hyb_genome.output
    output:
        "analysis/benchmarking/dummy/unload_genome_dummy4.txt"
    priority: 2
    shell:
        """
        bwa shm -d && \
        touch {input} && \
        touch {output}
        """      
        
        
rule unmapped_hyb:
    input:
        rules.map2hybrid.output
    output:
        read1 = "analysis/benchmarking/unmapped_hyb/hyum_{sample}_R1.fastq",
        read2 = "analysis/benchmarking/unmapped_hyb/hyum_{sample}_R2.fastq"
    priority: 1
    benchmark:
        "analysis/benchmarking/umhyb_{sample}.tsv"
    shell:
        """
        samtools fastq \
            -@ 16 \
            -f 4 \
            -1 {output.read1} \
            -2 {output.read2} \
            -s /dev/null \
            -0 /dev/null \
            -n {input}
        """                                