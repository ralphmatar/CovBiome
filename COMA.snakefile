import os
SAMPLES = list(set([i.replace("sample_", "") for i in os.listdir("data/raw_data")]))


rule all:
    input:
        "data/genomes/hybrid/hybrid.fna.amb",
        "analysis/dummy/unload_k2db_dummy.txt",
        "analysis/dummy/load_k2db_dummy.txt",
        "analysis/dummy/load_genome_dummy.txt",
        "analysis/dummy/unload_genome_dummy.txt",
        "analysis/QC/multiqc/merged_multiqc_report.html",
        "analysis/plots/bracken_count_matrix.csv",
        "analysis/fcounts.txt",
        expand("analysis/QC/fastqc/trimmed/{sample}_read_1_fastqc.html", sample=SAMPLES),
        expand("analysis/QC/fastqc/trimmed/{sample}_read_2_fastqc.html", sample=SAMPLES),
        expand("analysis/bracken/{sample}.bracken", sample=SAMPLES),
        expand("analysis/QC/qualimap/{sample}/{sample}.pdf", sample=SAMPLES),
        expand("data/trimmed_data/{sample}_read_1.fastq.gz", sample=SAMPLES),
        expand("data/trimmed_data/{sample}_read_2.fastq.gz", sample=SAMPLES),
        
         
rule concatenate_genomes:
    input:
        human = "data/genomes/human/GRCh38.fna",
        covid = "data/genomes/cov/cov_genome.fna"
    output: 
        "data/genomes/hybrid/hybrid.fna"
    priority: 20
    shell:
        """
        cat {input.human} {input.covid} > {output}
        """
        
        
rule genome_idexing:
    input:
        rules.concatenate_genomes.output
    output:
        "data/genomes/hybrid/hybrid.fna.amb"
    priority: 20
    shell:
        """
        bwa index {input}
        """


rule adapter_barcode_trimming:
    input:
        barcodes = "data/final_bc.fasta",
        read1 = "data/raw_data/sample_{sample}/read1_sample_{sample}.fastq.gz",
        read2 = "data/raw_data/sample_{sample}/read1_sample_{sample}.fastq.gz"
    output:
        "data/trimmed_data/{sample}_read_1.fastq.gz",
        "data/trimmed_data/{sample}_read_2.fastq.gz"
    priority: 19
    shell:
        """
        flexbar \
            -threads 4 \
            -r {input.read1}  \
            -p {input.read2}  \
            -aa Nextera \
            -br {input.barcodes} \
            -t data/trimmed_data/{wildcards.sample}_read  \
            -z GZ
        """
        
        
rule fastqc:
    input:
        read1 = "data/trimmed_data/{sample}_read_1.fastq.gz",
        read2 = "data/trimmed_data/{sample}_read_2.fastq.gz", 
    output:
        "analysis/QC/fastqc/trimmed/{sample}_read_1_fastqc.html",
        "analysis/QC/fastqc/trimmed/{sample}_read_2_fastqc.html"
    priority: 18
    shell:
        """
        fastqc \
            -t 10 \
            {input.read1} {input.read2} \
            -o analysis/QC/fastqc/trimmed
        """
        
        
rule load_genome:
    input:
        "data/genomes/hybrid/hybrid.fna"
    output:
        "analysis/dummy/load_genome_dummy.txt"
    priority: 17
    shell:
        """
        bwa shm {input} && \
        touch {output}
        """
        
        
rule map_to_hybrid_genome:
    input:
        genome = "data/genomes/hybrid/hybrid.fna",
        read1 = "data/trimmed_data/{sample}_read_1.fastq.gz",
        read2 = "data/trimmed_data/{sample}_read_2.fastq.gz"
    output:
        "analysis/mapped/{sample}.sorted.bam"
    priority: 16
    shell:
        """
        bwa mem \
            -t 16 \
            {input.genome} \
            {input.read1} \
            {input.read2} | \
        samtools sort \
            -@ 10 \
            -o {output} && \
        samtools index \
            -@ 10 \
            {output}
        """
        
        
rule unload_genome:
    input: 
        rules.load_genome.output
    output:
        "analysis/dummy/unload_genome_dummy.txt"
    priority: 15
    shell:
        """
        bwa shm -d && \
        touch {input} && \
        touch {output}
        """
        
        
rule get_unmapped_reads:
    input:
        rules.map_to_hybrid_genome.output
    output:
        read1 = "analysis/non_mapped_seqs/{sample}_R1.fastq",
        read2 = "analysis/non_mapped_seqs/{sample}_R2.fastq"
    priority: 14
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
   
       
rule load_database:
    input:
        "k2db16"
    output:
        "analysis/dummy/load_k2db_dummy.txt"
    priority: 13
    shell:
        """
        cp k2db16/*k2d /dev/shm/ && \
        touch {output}
        """
        
             
rule metagenomic_classification_with_kraken2:
    input:
        read1 = "analysis/non_mapped_seqs/{sample}_R1.fastq",
        read2 = "analysis/non_mapped_seqs/{sample}_R2.fastq",
        db = "k2db16"
    output:
        "analysis/classification/{sample}.kraken2"
    priority: 11
    shell:
        """
        kraken2 \
            --db {input.db} \
            --memory-mapping \
            --threads 16 \
            --classified-out analysis/classified_out/{wildcards.sample}_c#.fq \
            --use-names \
            --report {output} \
            --paired {input.read1} {input.read2} \
            --confidence 0.1
        """
        
        
rule unload_databse:
    input: 
        rules.load_database.output
    output:
        "analysis/dummy/unload_k2db_dummy.txt"
    priority: 9
    shell:
        """
        rm /dev/shm/*k2d && \
        touch {input} && \
        touch {output}
        """
 
 
rule abundance_reestimation_with_bracken:
    input:
        file = rules.metagenomic_classification_with_kraken2.output,
        database = "k2db16"
    output:
        "analysis/bracken/{sample}.bracken"
    priority: 7
    shell:
        """
        set +e
        bracken \
            -d {input.database} \
            -i {input.file} \
            -o {output} \
            -r 50 \
            -t 10 
        if [ ! -s "{output}" ]; then
            echo "Warning: Bracken output file {output} is empty or does not exist." >&2
            touch "{output}"
        fi
        set -e
        """
              
                
rule count_matrix:
    input:
        "analysis/bracken"
    output:
        "analysis/plots/bracken_count_matrix.csv"
    priority: 6
    shell:
        """
        python3 scripts/generate_count_matrix.py {input}    
        """
        
        
rule multiqc:
    input: 
        "analysis"
    output:
        "analysis/QC/multiqc/merged_multiqc_report.html"
    priority: 0
    shell:
        """
        multiqc \
            -o analysis/QC/multiqc \
            -i merged \
            -x analysis/benchmarking \
            analysis
        """
        
        
rule gene_count:
    input:
        annotation = "data/genomes/human/annotation.gff",
        bam = expand("analysis/mapped/{sample}.sorted.bam", sample=SAMPLES)
    output:
        "analysis/fcounts.txt"
    priority: 5
    shell:
        """
        featureCounts \
            -a {input.annotation} \
            -p \
            -g gene \
            -T 14 \
            -o {output} \
            {input.bam}
        """
        
        
rule qualimap:
    input:
        rules.map_to_hybrid_genome.output
    output:
        "analysis/QC/qualimap/{sample}/{sample}.pdf"
    priority: 3
    shell:
        """
        qualimap bamqc \
            -bam {input} \
            -c \
            -gd HUMAN \
            -nt 16 \
            -outdir analysis/QC/qualimap/{wildcards.sample}/ \
            -outfile {wildcards.sample}.pdf \
            -outformat PDF
        """