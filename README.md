hm: hyunmin's contribution
==

mpileup parser 
 * mpileup_to_freq: mpileup parser

```
  usage: 
  mpileup_to_freq.sh [options] <bam> [<bam>..]

# OUTPUT example:
chr1    935397  G   A:203,G:229 A:217,T:1,G:360 A:150,G:93  A:74,G:67
chr1    935670  C   C:59,G:125  C:78,G:159  C:43,G:23   C:29,G:18
chr1    940808  C   +:A:5,+:a:1,C:7 +:A:4,+:a:4,C:18    +:A:11,+:a:4,C:57   +:A:11,+:a:8,C:72
..

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



