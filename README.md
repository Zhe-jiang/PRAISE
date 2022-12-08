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
Takara Library is suitable for quantitative purposes of RNA with regular length (e.g. mRNA sequencing). \
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
Takara Library is suitable for semi-quantitative purposes of RNA with regular length (e.g. mRNA sequencing). From our [data](https://www.biorxiv.org/content/10.1101/2022.10.25.513650v1), the KAPA library result is pretty consistant with the Takara library result in whole transcriptome pseudouridine sequencing using PRAISE. However, if you want to get precise quantitative result, we still recomend Takara library since it has UMI for deduplication. \

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
