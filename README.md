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
# The goal
Identify virulence genes uniquely expressed in acute infection.
Detect metabolic and stress-response pathways that dominate during chronic infection.
Reveal regulatory RNAs and transcriptional signatures linked to biofilm persistence.

## Step 1: Preprocessing and Quality Control
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
