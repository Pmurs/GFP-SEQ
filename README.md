<h1>GFP-SEQ</h1>
This project is for the alignment and analysis of our Green Fluorescent Protein (GFP) samples.

<h2>Wetlab Context</h2>

<h3>🧬 GFP Reporter Design</h3>

N-GFP (G) and C-GFP (FP) are split GFP fragments

Combined, they form full GFP only after successful transposition

Artificial intron connects them in silico

Acts as a reporter: only fluoresces if SB transposase integrates C-GFP into the genome with N-GFP

---

<h3>🧪 Integration Capture</h3>

gDNA digested with either:

AluI or

FatI (restriction enzymes targeting LTR-adjacent sites)

Self-ligation performed

Then, LTR-sequence–primed PCR amplifies the junction
![image](https://github.com/user-attachments/assets/50faf5d7-945a-48d0-9b18-58ca5a36bb77)

---


<h3>🔬 Purpose</h3>

Detect integration sites using LTR-based mapping

GFP sequence shouldn't appear in sequencing reads — we're sequencing genome-junction regions, not full transgene
    
---

<h3>🧪 Sample types</h3>
    
| File Label | Meaning                                                 |
|------------|---------------------------------------------------------|
| `Alu`      | Digested with **AluI** restriction enzyme               |
| `Fat`      | Digested with **FatI** restriction enzyme               |
| `SB`       | **Sleeping Beauty** transposase–mediated integration    |
| `cout`     | **Cas9-mediated cut** outside the LTR region            |

---

<h3>🧫 Cell Line</h3>

293T (HEK293T human cells)

Align to HG38

---

<h3>🧬 Sequencing Context Summary</h3>

| Aspect                  | Details                                                                 |
|-------------------------|-------------------------------------------------------------------------|
| **Sequencing Provider** | [Plasmidsaurus](https://plasmidsaurus.com/premium_PCR_sequencing)                         |
| **Service Type**        | Premium PCR sequencing                                                  |
| **Platform**            | Oxford Nanopore Technologies (ONT)                                      |
| **Flow Cell**           | **R10.4.1**                                                             |
| **Read Type**           | Single-end, long-read (Nanopore)                                        |
| **Target**              | PCR amplicons from genomic junction capture                             |

---

<h2>Infrastructure</h2>

<h3>🖥️ Alignment Pipeline Summary</h3>

| Step | Description                                                                 |
|------|-----------------------------------------------------------------------------|
| **1. FastQC** | Runs quality control on raw FASTQ files using `fastqc` |
| **2. Filtlong** | Filters reads by length (100–1000 bp) |
| **3. Minimap2** | Aligns filtered reads to the **human genome (hg38)** |
| **4. Samtools** | Converts SAM to BAM, sorts, and indexes the alignments |
| **5. Medaka (optional)** | Performs consensus polishing and variant calling using filtered FASTQ reads and the reference genome |

> ⚠️ Note: Medaka requires the raw FASTQ input for consensus generation and internally realigns reads, so it does not use the BAM output from Minimap2.

---

## ⚙️ How to Create an Environment and Run the Pipeline

1. **Launch EC2 instance**
   - Use the `SMM-SEQ` launch template on AWS

2. **Download reference genomes**
    ```bash
    sudo bash SMM-SEQ/data/reference_setup.bash
    ```

3. **Configure SLURM job scheduler**
    ```bash
    sudo bash SMM-SEQ/data/slurm_setup.bash
    ```

4. **Clone the pipeline repository**
    ```bash
    git clone https://github.com/Pmurs/GFP-SEQ.git
    ```

5. **Install dependencies and create Python virtual environment**
    ```bash
    sudo bash GFP-SEQ/setup.bash
    ```

6. **Run the alignment and variant calling pipeline**
    ```bash
    sudo bash GFP-SEQ/data/process_files.sh
    ```

