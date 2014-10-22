mpileupToFreq(){ usage="
    usage: $FUNCNAME <mpileup>  
"
if [ $# -ne 1 ]; then echo $usage; return; fi

testData='
1   888659  T   14  c$CCCCcccCcCccc GCECEFGEEBGEE@
1   1268847 T   25  ,G...,..GG.....G..GGG,.g^], GE8BEF==B?EGCB=BA@DFG=EBA
1   1421531 C   27  aAaaAAAAAaaaaAaAAAaAAaaAaa^]a   GEDGEBGEECGCFF-DGEBGGB?FEA:
1   1888193 C   17  .,,.,.,,,,,.,.,..   DGFEGEGBEEEGAGEGG
1   2283117 C   24  .,.....................^].  ?GE#GFGEFBEGGGFGGGGGGGDG
1   2938697 T   22  CgGGggCGgGgGGggggggggG  #G=:CG#5E@E#:#A##?@C#E
1   3742257 C   11  ,....,..... GGFGFEGEDDG
1   3743391 C   28  ,$t,t,.t,tttt,t,t,t.t,,t,,t,t   G=GEGAGGGFGGEEGGGFGBGGGGCEEE
1   3753813 T   11  .$..........    BGEGEEGFGEF
1   3756074 C   12  ............    GEGG==GFGDDE
20  1234    A   20  ,,,,,..,.-4CACC.-4CACC....,.,,.^~.  ==<<<<<<<<<<<::<;2<<
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
        foreach my $k2 (keys %{$h{$k1}}){
            if($k1 eq "+"){ $k1 = "\\+";} ## fixed bug
            while($S =~/$k1$k2([ACGTNacgtn]{$k2})/g){
                $F{"$k1:$1"}++;
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

        ## output
        print "$chrom","\t",$base1pos-1,"\t",join (",", map{"$_:$F{$_}"} keys %F),"\n";
    }
'

if [ $1 == "test" ]; then
    echo "$testData" >&2
    echo "$testData" | perl -e "$cmd";
else
    cat $1 | perl -e "$cmd"
fi
}


bedSeq(){ usage="
    usage: $FUNCNAME <bed> <fasta|dir> <left_flank> <right_flank> <strand_specific>
     <strand_specific>: [1|0]
"
if [ $# -ne 5 ];then echo "$usage"; return; fi
    BED=$1;FA=$2;LEFT=$3;RIGHT=$4;STRAND=$5;

    cmd='use strict;
    sub revComp{
        my ($seq) = @_;
        $seq =~ tr/ACGTacgt/TGCAtgca/;
        return join("",reverse(split //,$seq));
    }
    my %seq=();
    my $chrom="";
    open (F, "FA") or die "$!";
    while(<F>){ chomp;
        if($_=~/>([\w|\d]+)/){ $chrom=$1; next;}
        $seq{$chrom} .= $_;
    } close(F);
    while(<STDIN>){ chomp; my @a=split /\t/,$_;
        my $sseq = substr($seq{$a[0]},$a[1],$a[2]-$a[1]);   
        my $sseq_left;my $sseq_right;
        if(STRAND==0 || $a[5] eq "+"){
            $sseq_left = substr($seq{$a[0]},$a[1]-LEFT,LEFT) if( LEFT > 0); 
            $sseq_right = substr($seq{$a[0]},$a[2],RIGHT) if(RIGHT > 0);    
        }else{
            $sseq = revComp($sseq);
            $sseq_left = revComp(substr($seq{$a[0]},$a[2],LEFT)) if( LEFT > 0); 
            $sseq_right = revComp(substr($seq{$a[0]},$a[1]-RIGHT,RIGHT)) if(RIGHT > 0); 
        }
        print $_,"\t",uc $sseq_left,",",uc $sseq,",",uc $sseq_right,"\n";
    }
    '
    cmd=${cmd//LEFT/$LEFT};
    cmd=${cmd//RIGHT/$RIGHT};
    cmd=${cmd//STRAND/$STRAND};
    FA=${FA%\/}
    if [ -f $FA ]; then
        cmd=${cmd//BED/$BED};
        cmd=${cmd//FA/$FA};
        cat $BED | perl -e "$cmd"
    elif [ -d $FA ]; then
        tmp_bed=`makeTemp`; cat $BED > $tmp_bed;
        chroms=(`getChroms $tmp_bed`)
        for chrom in ${chroms[@]};do
            fa=$FA/$chrom.fa
            if [ -f $fa ]; then
                cmd1=$cmd;
                cmd1=${cmd1//BED/$tmp_bed};
                cmd1=${cmd1//FA/$fa};
                awk -v CHROM=$chrom '$1==CHROM' $tmp_bed \
                | perl -e "$cmd1"
            fi
        done
    fi
}
