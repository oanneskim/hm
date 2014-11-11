#!/bin/env bash 

mpileupToFreq(){ usage="
	usage: $FUNCNAME <mpileup>  
"
if [ $# -ne 1 ]; then echo $usage; return; fi

testData='
1	888659	T	14	c$CCCCcccCcCccc	GCECEFGEEBGEE@
1	1268847	T	25	,G...,..GG.....G..GGG,.g^],	GE8BEF==B?EGCB=BA@DFG=EBA
1	1421531	C	27	aAaaAAAAAaaaaAaAAAaAAaaAaa^]a	GEDGEBGEECGCFF-DGEBGGB?FEA:
1	1888193	C	17	.,,.,.,,,,,.,.,..	DGFEGEGBEEEGAGEGG
1	2283117	C	24	.,.....................^].	?GE#GFGEFBEGGGFGGGGGGGDG
1	2938697	T	22	CgGGggCGgGgGGggggggggG	#G=:CG#5E@E#:#A##?@C#E
1	3742257	C	11	,....,.....	GGFGFEGEDDG
1	3743391	C	28	,$t,t,.t,tttt,t,t,t.t,,t,,t,t	G=GEGAGGGFGGEEGGGFGBGGGGCEEE
1	3753813	T	11	.$..........	BGEGEEGFGEF
1	3756074	C	12	............	GEGG==GFGDDE
20	1234	A	20	,,,,,..,.+4CACC.-4CACC....,.,,.^~.	==<<<<<<<<<<<::<;2<<
chr22	16155514	N	1	T+1T	@
chr1	17722	A	52	,,.,,.,,,,,,,,,.,,..........,..,,,.,,,,,...,,..,.,,,	@@DF>DHJJJJGBGJBJE?BBBFH;@IJEJJBDDJ?DCDCJIIDDJJDJDCD
chr1	17723	G	52	,$,$.,,.,,,,,,,,,.,,..........,..,,,.c,,,,...,,..,.,,,	>>BFBDHJJJJIIGJBIIBDDDFFBEJIFJJBDDJ8DCDDIJICDJJDJDCD
chr1	17724	G	50	.,,.,,,,,,,,,.,,..........,..,,,.,,,,,...,,..,.,,,	DFADHJJJJHIIJ@IEDDDDFFDHIIFJFACDJDDCDDJJJ@DJJDJD@D
chr1	17725	A	50	.$,,.,,,,,,,,,.,,..........,..,,,.,,,,,...,,..,.,,,	>F1DHJJJJGHGI?GGD?DDDFDCCGFHD@CDIDDADDJJJCDJJDJD:D
chr1	17726	A	48	,,.,,,,,,,,,.,,..........,..,,,.,,,,,...,,..,.,,	C?DHHJJJFGJJ8GGDBDEEFFAHGFHC;CDFDD>DDJIJDDJJDIDD
chr1	17727	A	47	,,.,,,,,,,,,.,,..........,..,,.,,,,,...,,..,.,,	C8DFHJJJFGIJ>HC>DDCEEEHHHFHADDCDD>DDJJIACJJCJDD
chr1	17728	G	47	,$.,,,,,,,,,.,,..........,..,,,.,,,,,...,,..,.,,	BAFHHJJ?EHJ?GGABDCEECBHGFHH>DDEDD@BDIJG8?JJBIBD
chr1	17729	T	47	.,,,,,,,,,.,,..........,..,,,.,,,,,...,,..,.,,,	AFHHHIFGFI4CG>A:C@A@EEBHEH5DD7DDDDBHHF1@IIBHD4B
chr1	17730	C	46	.,a,aaa,,,.,,..A.AA...A,..a,,.,,,,,...,,..,.,,	CFFHHHAFIJ:<B:ACCCECEH@HFC5DD=DDDDDIIHCDJJDIDD
chr1	17731	C	47	.$,,,,,,,,,.,,..........,..,,,.,,,,,...,,..,.,,,	BFFHHHF?GJ:2G>DCDCECBHCHFB>DDCDCDDDIJICDJJDJD@D	
'

cmd='
	use strict;
	my $offset=33; 
	#my @types= (split / /,"A C G T a c g t \$+ \$- \^+ \^-");
	#print "chrom\tbase0pos\t",join( "\t",@types),"\tinsertion\tdeletion\n";
	while(<STDIN>){
		chomp; next if ($_ eq "");
		my ($chrom, $base1pos,$refseq, $n, $S, $Q) = split /\t/,$_;
		if($S eq "*" || $n < 1){ next;}	

		my %F = (); ## frequencies

		## handle . and , 
		my $tmp = uc $refseq; $S=~s/\./$tmp/g;
		$tmp = lc $refseq; $S=~s/\,/$tmp/g;

		## handle insertion deletion
		my %h=(); while($S=~/([-|+])(\d+)/g){ $h{$1}{$2} = 1; } ## regular expression cannot read the previous pattern
		foreach my $k1 (keys %h){
		my $k1_old = $k1;
		foreach my $k2 (keys %{$h{$k1}}){
			if($k1 eq "+"){ $k1 = "\\+";} ## fixed bug
			while($S =~/$k1$k2([ACGTNacgtn]{$k2})/g){
				$F{"$k1_old:$1"}++;
			}
			$S =~s/$k1$k2[ACGTNacgtn]{$k2}//g;
		}}

		## handle others
		my @A = split //,$S;
		my @B = split //,$Q;
		for(my $i=0; $i <= $#A; $i++){
			my $s= $A[$i];
			my $q= ord($B[$i]) - $offset;
			my $p = 1-exp(-$q/10)/exp(10);
			my $w = int($p+0.5);
			if($i < $#A && $A[$i+1] eq "\$"){
				$s = $s."\$"; $i++; ## additional shifting
			}elsif($s eq "\^"){ $s .= $A[$i+2]; $i +=2; }
			$F{$s} += $w;
		}

		## output 1
		print "$chrom","\t",$base1pos-1,"\t",uc $refseq,"\t",join (",", map{"$_:$F{$_}"} keys %F),"\n";
		## output 2
		
	}
'
if [ $1 == "test" ]; then
	echo "$testData" >&2
	echo "$testData" | perl -e "$cmd";
else
	cat $1 | perl -e "$cmd"
fi
}

out4Tassa(){
	perl -e ' use strict;
	my @w = ("A","C","G","T");
	my @w2 = ();
	foreach my $i (@w){
	foreach my $j (@w){
		push @w2,$i.$j;
	}}
	print "chrom\tpos\t",join("\t",@w2),"\n";
	while(<>){ chomp;
		my ($chrom,$pos,$rs,$tmp) = split /\t/,$_;
		$rs=uc $rs;
		my %res = ();
		foreach my $e (split /,/,$tmp){
			my ($k,$v) = split /:/,$e;
			$k=~s/\$//g; $k=~s/\^.//g; $k= uc $k;
			if($k=~ /[ACGT]/){ $res{$rs.$k} += $v; }
		}
		print $chrom,"\t",$pos;
		foreach my $e (@w2){
			my $v=0;
			$v = $res{$e} if defined $res{$e};
			print "\t",$v
		}
		print "\n";
	}
	'
}
