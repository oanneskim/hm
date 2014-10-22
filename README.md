hm: hyunmin's contribution
==

bedfriehds.sh : converters and parsers for bed and its friends (sam,gff,...) 
---

 ```
  
  usage: 
  ## import functions
  . bedfriends.sh
  ## run an example
  mpileupToFreq test
  ## run with a data
  samtools mpileup <bam> | mpileupToFreq -
  
 ```

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



