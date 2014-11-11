hm: hyunmin's contribution
==

mpileup.sh:
---
mpileup parser 
 * mpileupToFreq: mpileup parser

```
  usage: 
  ## import functions
  . mpileup.sh
  ## test 
  mpileupToFreq test
  ## run 
  samtools mpileup <bam> | mpileupToFreq -
```
  
* out4Tassa: mpileup output formatter
```
  mpileupToFreq <mpileup_file> | out4Tassa
```

bedfriends.sh : 
---
converters and parsers for bed and its friends (sam,gff,...) 

 
  
  ```
  * bedSeq: extract genomic seqences for each bed interval
  ```
  ## usage:
  bedSeq <bed> <fasta|fasta_dir> <left_flank> <right_flank> <strand_specific>
  
  ## example(flanking 10bp):
  # strand_specific=0
chr22	51062323	51062325	t:1;t:1	51062323,51062324	-	GAGACTCCGT,CT,CAAAAAAAAA
  # strand_specific=1
chr22	51062323	51062325	t:1;t:1	51062323,51062324	-	TTTTTTTTTG,AG,ACGGAGTCTC


  ```
Others:
---
  * nimbleGenGff_to_bam: convert nimblegen GFF to bam file (the log2 ratio is stored in the 4th column)
	```
	# example 
	. hm.sh; nimbleGenGff_to_bam data/nimblegen_chipchip_sample0.gff data/scer3_chrom.size > data/nimblegen_chipchip_sample0.bam
	```

  * count_543bins.sh : make 543 bins, which can be converted to a table format using dcast in R
	```
	# example
	count_543bins.sh -a data/scer3_sgdGene.bed -b data/nimblegen_chipchip_sample0.bam -c 4
	```



