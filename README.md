# PRAISE: A quantitative pseudoridine sequencing method
This is the bioinformatics guide and scripts for PRAISE, a quantitative pseudouridine sequencing method. Our paper is now available at [BioRxiv](https://www.biorxiv.org/content/10.1101/2022.10.25.513650v1). \
For all the codes shown in this document, the parameters are the same as our published data. You can change suitable parameters for yourself.

## 0 Requirments

- Softwares
  - cutadapt
  - umi_tools
  - seqkit
  - hisat2
  - samtools
- python3
- packages of python3
	- pysam

## 1 Overview of the bioinformatics pipeline of PRAISE
This part is a detailed bioinformatics guide for our quantitative pseudoridine sequencing method: PRAISE.

### 1.1 Prepare your reads
Here, we provide optimized reads pre-processing for three different library: eCLIP, KAPA, and Takara (A modified Takara library method, see experimental details in our [Paper](https://www.biorxiv.org/content/10.1101/2022.10.25.513650v1)).

#### 1.1.1 Takara Library
Takara Library is suitable for quantitative purposes of RNA with regular length (e.g. mRNA sequencing). 
- Cut adaptor, __cutadapt__ is recommended.
```bash
cutadapt -j {CORES} --times 1 -e 0.1 -O 3 --quality-cutoff 25 -m 55 \
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-o {output.fix_R1} \
-p {output.fix_R2} \
 {input.raw_R1} \
 {input.raw_R2}
 # {CORES}: core number used
 # {output.fix_R1}: Output read 1 file (.fq.gz)
 # {output.fix_R2}: Output read 2 file (.fq.gz)
 # {input.raw_R1}: Input read 1 file (.fq.gz)
 # {input.raw_R2}: Input read 2 file (.fq.gz)
```

- Remove PCR duplication, __seqkit__ is recommended.
  - We remove PCR duplication before mapping instead of removing after mapping to decrease the computation of realignment step.
  - We only need Read 2.
```bash
seqkit rmdup {input.fix_R2} -s -j {CORES} -o {output.dedup_R2}
# {input.fix_2}: read 2 file after cutadapt (.fq.gz)
# {CORES}: core number used
# {output.dedup_R2}: read 2 file after deduplication (.fq.gz)
```

- Cut 14 mer of the 5' end of read, __umi_tools__ is recommended.
  - the length of UMI is 8 mer
  - cut another 6 mer of template switch to increase the confidence of the mapping result
  - 14 mer total
```bash
umi_tools extract --extract-method=string --bc-pattern=NNNNNNNNNNNNNN \
-I {input.dedup_R2} -S {output.cut5prime_R2}
# {input.dedup_R2}: read 2 file after deduplication (.fq.gz)
# {output.cut5prime_R2}: read file after cutting 14 mer of 5' end of read (.fq.gz)
```

- Cut 6 mer of the 3' end of read, __umi_tools__ is recommended.
	- Combination of random primer of RT may introduce INDELs and/or mismatches, cut 6 mer to increase the confidence of the mapping result.
	- This is the final step of reads pre-processing of Takara library.
```bash
umi_tools extract --extract-method=string --3prime --bc-pattern=NNNNNN \
-I {input.cut5prime_R2} -S {output.cut3prime_R2}
# {input.cut5prime_R2}: read 2 file after cutting 14 mer of 5' end of read (.fq.gz)
# {output.cut3prime_R2}: read 2 file after cutting 6 mer of 3' end of read (.fq.gz)
```

#### 1.1.2 KAPA Library
Takara Library is suitable for semi-quantitative purposes of RNA with regular length (e.g. mRNA sequencing). From our [data](https://www.biorxiv.org/content/10.1101/2022.10.25.513650v1), the KAPA library result is pretty consistant with the Takara library result in whole transcriptome pseudouridine sequencing using PRAISE. However, if you want to get precise quantitative result, we still recomend Takara library since it has UMI for deduplication. 

- Cut adaptor, __cutadapt__ is recommended.
	- This is the only step for pre-processing KAPA library reads.
```bash
cutadapt -j {CORES} --times 1 -e 0.1 -O 3 --quality-cutoff 25 -m 55 \
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-o {output.fix_R1} \
-p {output.fix_R2} \
 {input.raw_R1} \
 {input.raw_R2}
 # {CORES}: core number used
 # {output.fix_R1}: Output read 1 file (.fq.gz)
 # {output.fix_R2}: Output read 2 file (.fq.gz)
 # {input.raw_R1}: Input read 1 file (.fq.gz)
 # {input.raw_R2}: Input read 2 file (.fq.gz)
```

#### 1.1.3 eCLIP Library
eCLIP Library is suitable for quantitative purposes of RNA with short length (e.g. tRNA sequencing).

- Cut adaptor, __cutadapt__ is recommended.
```bash
cutadapt -j {CORES} --times 1 -e 0.1 -O 3 --quality-cutoff 25 -m 30 \
-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
-o {output.fix_R1} \
-p {output.fix_R2} \
 {input.raw_R1} \
 {input.raw_R2}
 # {CORES}: core number used
 # {output.fix_R1}: Output read 1 file (.fq.gz)
 # {output.fix_R2}: Output read 2 file (.fq.gz)
 # {input.raw_R1}: Input read 1 file (.fq.gz)
 # {input.raw_R2}: Input read 2 file (.fq.gz)
```

- Remove PCR duplication, __seqkit__ is recommended.
  - We remove PCR duplication before mapping instead of removing after mapping to decrease the computation of realignment step.
  - We only need Read 2.
```bash
seqkit rmdup {input.fix_R2} -s -j {CORES} -o {output.dedup_R2}
# {input.fix_2}: read 2 file after cutadapt (.fq.gz)
# {CORES}: core number used
# {output.dedup_R2}: read 2 file after deduplication (.fq.gz)
```

- Cut 10 mer of the 5' end of read, __umi_tools__ is recommended.
  - the length of UMI is 8 mer.
  - cut additional 2 mer of possible template switch to increase the confidence of the mapping result.
	- This is the final step of reads pre-processing of eCLIP library.
```bash
umi_tools extract --extract-method=string --bc-pattern=NNNNNNNNNNNNNN \
-I {input.dedup_R2} -S {output.cut5prime_R2}
# {input.dedup_R2}: read 2 file after deduplication (.fq.gz)
# {output.cut5prime_R2}: read file after cutting 14 mer of 5' end of read (.fq.gz)
```

#### 1.1.4 Other libraries
- Please follow the pre-processing protocols of the specific libraries

### 1.2 Mapping your reads

---
__IMPORTANT MESSAGES:__
- If you want to do realignment, be sure to use references without any intron (i.e. using transcriptome instead of genome). Because our realigment scripts do not support mapping results with intron.
- If you do not want to do realignment, we still recomend using references without any intron and the mapping method without spliced alignment. Because our signal is deletion, it may cause errors with spliced alignment.
- This part is a general tutorial, see details about realignment at Part 3.
---

#### 1.2.1 Mapping with realignment
Realignment will give you a more precise result but require a amount of computational resource. If you want to get a precise quantitative result, especially on genome wide sequencing, it is better to do realignment.

- First mapping with hisat2, __hisat2__ is required
	- Be sure to use parameter '--no-spliced-alignment'
	- Be sure to use parameter '--very-sensitive' because this step is to get all possible mapping results
	- Be sure to use the reference without any intron

```bash
hisat2 -p {CORES} \
-x {INDEX_HISAT2} \
-q --repeat --no-spliced-alignment --very-sensitive \
-U {input.processed_R2} \
-S {output.SAM}
# {CORES}: core number used
# {INDEX_HISAT2}: index needed for hisat2
# {input.processed_R2}: pre-processed read 2 file (.fq.gz)
# {output.SAM}: output SAM file (.SAM)
```

- Sort and change into BAM, __samtools__ is recommended

```bash
samtools sort -O BAM -o {output.bam} -@ {CORES} -m {MEM} -T {output.bam_temp} {input.sam}
# {CORES}: core number used
# {MEM}: Memory used per core (e.g. 2G)
# {output.bam}: output sorted BAM file (.BAM)
# {output.bam.temp}: temperory file
# {input.sam}: input SAM file (.SAM)
```

- Create index for BAM, __samtools__ is recommended

```bash
samtools index {input.bam} {output.bai}
# {input.bam}: Input sorted BAM file (.BAM)
# {output.bai}: Output index file (.BAM.BAI)
```

- Filter all unmapped reads, __samtools__ is recommended

```bash
samtools view -h -@ {CORES} -bF 4 {input.bam_sorted} -o {output.bam_mapped}
# {CORES}: core number used
# {input.bam_sorted}: Input sorted BAM file (.BAM)
# {output.bam_mapped}: Output BAM file without unmapped reads (.BAM)
```

- Sort by reads name, __samtools__ is recommended
	- Realignment requires input BAM file input with name sorted

```bash
samtools sort -@ {CORES} -m {} -O BAM -n -o {output.bam_name_sorted} -T {output.bam_name_sorted.temp} {input.bam_mapped}
# {CORES}: core number used
# {MEM}: Memory used per core (e.g. 2G)
# {output.bam_name_sorted}: output name sorted BAM file (.BAM)
# {output.bam_name_sorted.temp}: temperory file
# {input.bam_mapped}: input BAM file without unmapped reads (.BAM)
```

- Realign, __scripts here__ is required
	- See scripts details at __Part 3__

```bash
python {Realign_script} --fast -t {CORES} -ms 4.8 \
-x {Reference} -i {input.bam_name_sorted} -o {output.bam_realigned} -f {output.bam_filtered}
# {Realign_script}: use realignment_forward.py for reads aligned to the forward sequences, use realignment_reverse.py for reads aligned to the reverse sequences
# {CORES}: core number used
# {Reference}: Reference used (.fa, same reference as the hisat2 index created from)
# {input.bam_name_sorted}: Output name sorted BAM file (.BAM)
# {output.bam_realigned}: Result BAM file of the realignment (.BAM)
# {output.bam_filtered}: BAM file contains all filtered reads (.BAM)
```


#### 1.2.2 Mapping without realignment

