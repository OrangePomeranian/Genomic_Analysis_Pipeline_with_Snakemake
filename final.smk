# All rules contained together for slurm run
sample_names, = glob_wildcards("working_folder/gz/{sample}.fq")
# Create rule all
rule all:
    input:
        "working_folder",
        "working_folder/final.smk",
        "working_folder/reference/chr21.fa",
        "working_folder/gz",
        "working_folder/results_fastqc",
        expand("working_folder/results_bwa/{sample}.bam", sample=sample_names),
        expand("working_folder/results_bwa/{sample}.bam.bai", sample=sample_names),
        "working_folder/results_bcf/output.vcf",
        "working_folder/results_bcf_cleaned/output_cleaned.vcf",
        "working_folder/snpeff/snpEff_summary.html",
        "working_folder/snpeff/snpEff_genes.txt",
        "working_folder/snpeff_final/snps.annotated2.vcf",
        "working_folder/genes.vcf",
        "working_folder/test/fastqc_unzipped",
        "working_folder/pictures/summary.png"




 




# Creation of new working folder 
rule create_folder_copy_smk:
    input:
        "final.smk"
    output:
        folder = directory("working_folder"),
        file = "working_folder/final.smk"
    shell:
        """
        mkdir -p {output.folder} &&
        cp {input} {output.file} &&
        if test -f "{output.file}"; then
            echo "File {output.file} copied successfully."
        else
            echo "Failed to copy {input} to {output.file}."
            exit 1
        fi
        """


# Copy of 21th BWA and indexing
rule copy_21th_BWA:
    input:
        '/lustre1/project/stg_00079/teaching/hg38_21/chr21.fa'
    output:
        reference = "working_folder/reference/chr21.fa"
    shell:
        """
        mkdir -p working_folder/reference
        cp {input} {output.reference}
        bwa index {output.reference} 
        if grep FAIL {output.reference}; then
            echo "FAILED!"
            false
        fi        
        """


# Files ending with 9 are being copied and unpacked to the working_folder/gz 
rule copy_fq_unpack_files:
    input:
        "/staging/leuven/stg_00079/teaching/1000genomes/"
    output:
        directory("working_folder/gz")
    shell:
        """
        mkdir -p working_folder/gz
        cp {input}/*9.GRCh38DH.exome.chr21.fq.gz working_folder/gz &&
        gunzip {output}/*9.GRCh38DH.exome.chr21.fq.gz 

        if grep FAIL {output}/*9.GRCh38DH.exome.chr21.fq; then
            echo "FAILED!"
            false
        else
            echo "Process completed successfully."
        fi
        """

# Running fastqc on files
rule fastqc:
    input:
        "working_folder/gz"
    output:
        directory("working_folder/results_fastqc")
    shell:
        """
        mkdir -p working_folder/results_fastqc &&
        fastqc {input}/*9.GRCh38DH.exome.chr21.fq -o {output} 
        if grep FAIL {output}; then
            echo "FAILED!"
            false
        fi
        """


# unzip
# creation of plot

# Aligning fasta files to reference
rule bwa:
    input:
         reference ='working_folder/reference/chr21.fa',
         files = "working_folder/gz/{sample}.fq"
         #bams=expand("working_folder/results_bwa/{sample}.bam", sample=sample_names)
    output:
        bam = "working_folder/results_bwa/{sample}.bam",
        bai = "working_folder/results_bwa/{sample}.bam.bai"
    shell:
        """
        mkdir -p working_folder/results_bwa &&
        bwa mem {input.reference} {input.files} \
            | samtools sort - \
            > {output.bam} 
        samtools index {output.bam} {output.bai} 
        if grep FAIL {output.bam} || grep FAIL {output.bai}; then
            echo "FAILED!"
            false
        fi
        """

# Calling variants
rule call_variants:
    input:
        reference="working_folder/reference/chr21.fa",
        bams=expand("working_folder/results_bwa/{sample}.bam", sample=sample_names),

    output:
        vcf="working_folder/results_bcf/output.vcf"
    params:
        threads=8
    shell:
        """
        mkdir -p working_folder/results_bcf &&
        bcftools mpileup -Ou -f {input.reference} {input.bams} | bcftools call -mv -o {output.vcf}
        if grep FAIL {output.vcf}; then
            echo "FAILED!"
            false
        fi
        """

# Cleaning of variants
rule variant_cleanup:
    input:
        reference="working_folder/reference/chr21.fa",
        vcf="working_folder/results_bcf/output.vcf"
    output:
        vcf="working_folder/results_bcf_cleaned/output_cleaned.vcf"
    shell:
        """
        mkdir -p working_folder/results_bcf_cleaned &&
        ( cat {input.vcf} \
           | vt decompose - \
           | vt normalize -n -r {input.reference} - \
           | vt uniq - \
           | vt view -f "QUAL>20" -h - \
           > {output.vcf} )
        if grep FAIL {output.vcf}; then
            echo "FAILED!"
            false
        fi
        """

# setting the paths for references
snpeff_jar = "/lustre1/project/stg_00079/teaching/I0U19a_conda_2024/share/snpeff-5.2-0/snpEff.jar"
snpeff_genome = 'hg38'
snpeff_db_folder = '/staging/leuven/stg_00079/teaching/snpeff_db'

# Running snpeff
rule snpeff:
    input:
        vcf = "working_folder/results_bcf_cleaned/output_cleaned.vcf",
    params:
        snpeff_db_folder = snpeff_db_folder,
        snpeff_jar = snpeff_jar,
        snpeff_genome = snpeff_genome,
    log:
        err = "working_folder/snpeff/snakemake.err",
    output:
        vcf = "working_folder/snpeff/snps.annotated.vcf",
        html = "working_folder/snpeff/snpEff_summary.html",
        genetxt = "working_folder/snpeff/snpEff_genes.txt",
    shell:
        """
        mkdir -p working_folder/snpeff &&
        java -Xmx4096m -jar \
            {params.snpeff_jar} eff {params.snpeff_genome} \
            -dataDir {params.snpeff_db_folder} \
            {input.vcf} > {output.vcf}

        # move output files to the snpeff output folder
        mv snpEff_genes.txt snpEff_summary.html working_folder/snpeff/
        if grep FAIL {output.vcf} || grep FAIL {output.html} || grep FAIL {output.genetxt}; then
            echo "FAILED!"
            false
        fi
        """
# extracting SNPs
rule extract_SNP:
    input:
        vcf = "working_folder/snpeff/snps.annotated.vcf",
    output:
        output = "working_folder/snpeff_final/snps.annotated2.vcf"
    shell:
        """
        mkdir -p working_folder/snpeff_final &&
        cp {input.vcf} {output.output}
        if grep FAIL {output.output}; then
            echo "FAILED!"
        false
        fi
        """

# Specification of genes of interests
genes_of_interest = ["APP", "SOD1", "DYRK1A"]

# Rule to filter SNPs associated with specified genes
rule filter_snps:
    input:
        "working_folder/snpeff_final/snps.annotated2.vcf"
    output:
        "working_folder/genes.vcf"
    shell:
        """
        grep -e "{genes_of_interest[0]}" -e "{genes_of_interest[1]}" -e "{genes_of_interest[2]}" {input} > {output}
        if grep FAIL {output}; then
            echo "FAILED!"
        false
        fi
        """
# unziping the results from fastqc to get the information that would be used for creation of png 
rule unzip_fastqc:
    input:

        "working_folder/results_fastqc/HG00099.GRCh38DH.exome.chr21_fastqc.zip"
    output:
        directory("working_folder/test/fastqc_unzipped")
    shell:
        """
        mkdir -p {output} &&
        unzip {input} -d {output}
        """

# Creation of png 
rule fastqc_report_image:
    input:
        summarytxt = "working_folder/fastqc_unzipped/HG00099.GRCh38DH.exome.chr21_fastqc/summary.txt"
    output:
        statuspng = "working_folder/pictures/summary.png"

    run:
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt

        data = pd.read_csv(input.summarytxt, sep="\t", header=None)
        data.columns = ['status', 'test', 'sample']

        data['x'] = 1

        sns.set(style="whitegrid")
        plt.figure(figsize=(6, 6))
        sns.scatterplot(data=data, x='x', y='test', hue='status', s=200)
        plt.xlabel('')
        plt.ylabel('Test')
        plt.title('FastQC Summary')
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.tight_layout()

        plt.savefig(output.statuspng)



# snakemake -s final.smk --report final.html
# snakemake --delete-all-output -c1
# rm -rf .snakemake
# snakemake -s final.smk -c2
