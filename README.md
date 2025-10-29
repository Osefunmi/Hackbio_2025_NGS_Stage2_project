# Hackbio_2025_NGS_Stage2_project
## Transcriptomic Profiling of Staphylococcus aureus During Acute vs Chronic Phases of Periprosthetic Joint Infection (PJI)
Background and Rationale
Periprosthetic joint infections (PJIs) are among the most devastating complications of orthopedic implants. They increase morbidity, prolong hospital stays, and often require costly revision surgeries. Staphylococcus aureus—particularly methicillin-resistant strains (MRSA)—is a leading cause of PJIs.
One critical feature of S. aureus is its ability to switch phenotypes between acute and chronic infection phases.

Acute phase: Bacteria adopt an aggressive, planktonic growth mode, expressing virulence factors such as toxins, adhesins, and immune evasion genes.
Chronic phase: Bacteria adapt to a biofilm-like state, downregulating overt virulence and upregulating persistence pathways (stress response, metabolic rewiring, antibiotic tolerance).
This adaptive flexibility makes chronic PJIs notoriously difficult to eradicate. Antibiotic regimens often fail, and host immune responses are blunted by biofilm shielding.

Why RNA-seq?
RNA sequencing provides a window into the global transcriptional programs that underpin this acute-to-chronic transition. By capturing gene expression profiles directly from S. aureus isolates in different clinical phases of PJI, RNA-seq can:
## The goal
Identify virulence genes uniquely expressed in acute infection.
Detect metabolic and stress-response pathways that dominate during chronic infection.
Reveal regulatory RNAs and transcriptional signatures linked to biofilm persistence.

## Step 1: Preprocessing and Quality Control
### Perform read trimming, alignment to the S. aureus reference genome, and assessment of sequencing quality.

`mkdir data && cd data`

`nano data.sh`

```
#!/usr/bin/env bash
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/079/SRR20959679/SRR20959679_1.fastq.gz -o SRR20959679_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/079/SRR20959679/SRR20959679_2.fastq.gz -o SRR20959679_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/080/SRR20959680/SRR20959680_1.fastq.gz -o SRR20959680_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/080/SRR20959680/SRR20959680_2.fastq.gz -o SRR20959680_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/081/SRR20959681/SRR20959681_1.fastq.gz -o SRR20959681_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/081/SRR20959681/SRR20959681_2.fastq.gz -o SRR20959681_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/077/SRR20959677/SRR20959677_1.fastq.gz -o SRR20959677_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/077/SRR20959677/SRR20959677_2.fastq.gz -o SRR20959677_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/082/SRR20959682/SRR20959682_1.fastq.gz -o SRR20959682_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/082/SRR20959682/SRR20959682_2.fastq.gz -o SRR20959682_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/076/SRR20959676/SRR20959676_1.fastq.gz -o SRR20959676_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/076/SRR20959676/SRR20959676_2.fastq.gz -o SRR20959676_2.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/078/SRR20959678/SRR20959678_1.fastq.gz -o SRR20959678_1.fastq.gz
curl -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/078/SRR20959678/SRR20959678_2.fastq.gz -o SRR20959678_2.fastq.gz

```
save and close the bash script

`bash data.sh`

`cd ../`

`nano preprocess.sh`

```
#!/bin/bash
#Quality control
mkdir -p fastqc
fastqc data/*.fastq.gz -o fastqc

#performing multiqc
multiqc fastqc

#Trimming
#Create output directories if they don't exist
mkdir -p trimmed_reads fastp_reports

for file in data/*_1.fastq.gz
do
    #Extract sample name (everything before "_1.fastq.gz")
    sample=$(basename "$file" _1.fastq.gz)
    
    echo "Processing $sample ..."

    fastp \
        -i data/${sample}_1.fastq.gz \
        -I data/${sample}_2.fastq.gz \
        -o trimmed_reads/${sample}_1.trim.fastq.gz \
        -O trimmed_reads/${sample}_2.trim.fastq.gz \
        -h fastp_reports/${sample}_fastp.html \
        -j fastp_reports/${sample}_fastp.json \
        -q 20 -u 30 -l 50 -w 4
        
    echo "Completed $sample"
done
gunzip data/*.fastq.gz
```

`bash preprocess.sh`

Getting the reference genome: go to the site https://www.ncbi.nlm.nih.gov/datasets/genome/. Search for the organism of interest (Staphylococcus aureus) and select the appropriate assembly. A good choice is a complete, "Reference" or "Representative" genome from a common strain. Use the download button or click FTP. Next, copy the link with ".fna.gz" as your reference genome and ".gff.gz or .gtf.gz" as annotation

`mkdir Genome && cd Genome`

`wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.fna.gz`

`mv GCF_000013425.1_ASM1342v1_genomic.fna.gz s_aureus.fna.gz`

`gunzip s_aureus.fna.gz`

`cd ./`

`nano runstar.sh`

```
#!/bin/bash
#performing star
#creating genome index
cd Genome
STAR --runMode genomeGenerate --genomeDir genomeIndex --genomeFastaFiles S_aureus.fna
cd ../
# Loop through all R1 (forward) reads
mkdir mapped_reads
gunzip trimmed_reads/*.fastq.gz
for infile in trimmed_reads/*_1.trim.fastq; do
    # Extract the sample name (everything before _1.fastq.gz)
    sample=$(basename "$infile" _1.trim.fastq)

    echo "Processing $sample ..."

    # Define paired-end files
    R1=$infile
    R2=trimmed_reads/${sample}_2.trim.fastq

    # Define paired-end files
    R1=$infile
    R2=trimmed_reads/${sample}_2.trim.fastq

    # Define output prefix and directory
    outfile=mapped_reads/${sample}_

    # Run STAR alignment
    STAR --runThreadN 8 \
         --genomeDir Genome/genomeIndex \
         --readFilesIn "$R1" "$R2" \
         --outFileNamePrefix "$outfile" \
         --outSAMtype BAM SortedByCoordinate \
         --outSAMattributes All

    echo "$sample done!"
done
mkdir IGV
cp mapped_reads/*.bam IGV/
cd IGV
#creating indices for all bam files
for infile in *.bam; do
    samtools index "$infile"
done 
echo "All samples processed successfully!"

```

`pwd`

download the IGV folder to your local computer. The next code stores it in the folder Stage3_project in HackbioNGS in Downloads.

`scp -r a_adegite@135.181.163.242:/home/a_adegite/Oluwasefunmi/Project/RNA_Seq/IGV Downloads/HackbioNGS/Stage3_project`

visit (https://igv.org/app/) specifiy your genome and upload the bam and bai files to view.

### Generate count matrices of gene expression for acute and chronic isolates.

`mkdir Counts`

`cd Counts`

`wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/013/425/GCF_000013425.1_ASM1342v1/GCF_000013425.1_ASM1342v1_genomic.gff.gz`

`mv GCF_000013425.1_ASM1342v1_genomic.gff.gz S_aureus.gff.gz`

`gunzip S_aureus.gff.gz`
`ls`
FEATURE COUNTS
If data was unpaired the right way to do feature counts is to use the code "featureCounts -O -t gene -g ID -a S_aureus.gff -o counts.txt ../IGV/*.bam" our data is paired so we have to go through another route.

`featureCounts -p -B -C --primary -T 8 -t gene -g ID -a S_aureus.gff -o g_counts.txt ../IGV/*.bam` 
Then you copy your counts.txt file to your computer.
`scp -r a_adegite@135.181.163.242:/home/a_adegite/Oluwasefunmi/Project/RNA_Seq/Counts/g_counts.txt Downloads/HackbioNGS/Stage3_project`
enter the password and the file would be saved.

## Step 2: Differential Gene Expression Analysis
First create your metadata in Google sheets and store it as a csv file.

```
Sample	State
SRR20959676	Chronic
SRR20959677	Chronic
SRR20959678	Chronic
SRR20959679	Chronic
SRR20959680	Acute
SRR20959681	Acute
SRR20959682	Acute
```
The next steps would be done in your R studio. Make a 


