<h1>GFP-SEQ</h1>

This project is for the alignment and analysis of our **Green Fluorescent Protein (GFP)** samples.

This project uses a split-GFP reporter system to validate integration and targeting specificity of DNA delivered via transposase or plasmid-based methods. 

We aim to identify **genomic integration sites**, characterize their **on- vs. off-target nature**, and distinguish between **transposase-mediated** and **random plasmid insertions**.

---

<h3>🔬 Purpose</h3>

- **Determine where GFP integrates in the genome**
    - Record chromosomal location and frequency
    - Focus on **targeted integration sites** near known guide RNA regions (~1 kb window)
- **Quantify on-target vs off-target integration**
    - Create a ratio of integrations near expected loci vs elsewhere
    - Higher precision expected when using guide RNAs
    - Detect integration sites using **Long Term Repeat (LTR)**-based mapping
    
---

### ⚠️ Interpretation Caveats

- **GFP Sequence**
    - Should not appear in reads
    - We're sequencing genome-junction regions, not full transgene
- **Read count ≠ Cell count**
    - PCR amplification and unknown cell input prevent using read counts to estimate number of cells
- **Multiple insertions**
    - Multiple insertions in the same genome can't be distinguished from separate cells
- **Off-target insertions**
    - Any insertions far from the expected target region (e.g. >1 kb) are considered off-target
- **Plasmid contamination**
    - Insertions containing plasmid backbone sequence are likely random/plasmid-driven events
    - High rates of plasmid integration may indicate protocol issues

---

### 🧭 Downstream Analysis Considerations

- Align reads to:
    - Human reference genome (hg38)
    - GFP reporter construct
    - Plasmid backbone sequence
- Classify insertions as:
    - **On-target** (near guide RNA)
    - **Off-target**
    - **Plasmid-based (non-transposase)**

---

<h2>Wetlab Context</h2>

<h3>🧬 GFP Reporter Design</h3>

- N-GFP (G) and C-GFP (FP) are split GFP fragments
- Combined, they form full GFP only after successful transposition
- Artificial intron connects them in silico
- Acts as a reporter: only fluoresces if SB transposase integrates C-GFP into the genome with N-GFP

---

<h3>🧪 Integration Capture</h3>

1. gDNA digested with either:
   - AluI
   - FatI 
   > Restriction enzymes targeting LTR-adjacent sites
2. Self-ligation performed
3. LTR-sequence–primed PCR amplifies the junction

![image](https://github.com/user-attachments/assets/50faf5d7-945a-48d0-9b18-58ca5a36bb77)

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

- **Cells used**: 293T human cells
    - **Align to Hg38**
- **Integration mechanisms**:
    - **Transposase-mediated (canonical)**: inserts only the sequence **between LTRs** into the genome
    - **Plasmid integration (off-target/aberrant)**: occurs randomly, often includes plasmid backbone sequence
- **Electroporation delivery method**:
    - Plasmid and transposase mRNA introduced via electric pulse
    - DNA enters cells due to negative charge and membrane permeabilization
    - Cell uptakes transposase and manufactures it as a protein
    - Buffers and conditions optimized for ~80% cell survival
    - Electroporation preferred over lipid nanoparticles (LNPs) for high-efficiency DNA delivery

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

> ⚠️ Note: Medaka requires the raw FASTQ input for consensus generation and internally realigns reads, so it does not use the BAM output from Minimap2. It can be GPU-accelerated for faster alignment

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

4. **Clone the GFP-SEQ repository**
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
