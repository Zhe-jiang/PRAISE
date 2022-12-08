# PRAISE: A quantitative pseudoridine sequencing method
This is the bioinformatics guide and scripts for PRAISE, a quantitative pseudouridine sequencing method. Our paper is now available at [BioRxiv](https://www.biorxiv.org/content/10.1101/2022.10.25.513650v1). \
For all the codes shown in this document, the parameters are the same as our published data. You can change suitable parameters for yourself.

## 1 Overview of the bioinformatics pipeline of PRAISE
This part is a detailed bioinformatics guide for our quantitative pseudoridine sequencing method: PRAISE.

### 1.1 Prepare your sequences
Here, we provide optimized sequences pre-processing for three different library: eCLIP, KAPA, and Takara (A modified Takara library method, see experimental details in our [Paper](https://www.biorxiv.org/content/10.1101/2022.10.25.513650v1)).

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

#### 1.1.2 KAPA Library


#### 1.1.3 eCLIP Library
