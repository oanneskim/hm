hm: hyunmin's contribution
==



bedfriehds.sh : 
---
converters and parsers for bed and its friends (sam,gff,...) 

  * mpileupToFreq: mpileup parser
 
  ```
  
  usage: 
  ## import functions
  . bedfriends.sh
  ## test 
  mpileupToFreq test
  ## run 
  samtools mpileup <bam> | mpileupToFreq -
  
  ```
  * bedSeq: extract genomic seqences for each bed interval
  ```
  ## usage:
  bedSeq <bed> <fasta|fasta_dir> <left_flank> <right_flank> <strand_specific>
  
  ## examplar output1 (flanking 10bp):
  chr22   16168803    16168860    gene	0   +   TAAAGCATTT,GATAAAATTCAACATTTTATTGTAAAAACCCTTCAAAAAACAGGGTATAGAAGAAAC,ATACCTCAAC
  
  ## example2 (with strand_specific = 1)
  chr22   16168803    16168860    gene	0   +	TAAAGCATTT,GATAAAATTCAACATTTTATTGTAAAAACCCTTCAAAAAACAGGGTATAGAAGAAAC,ATACCTCAAC
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



