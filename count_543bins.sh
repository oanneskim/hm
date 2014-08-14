#!/bin/bash
FILEA=
FILEB=
FIVE="-500,500,4"
FOUR="500,0,4"
THREE="0,1000,4"
CL=-1;
usage(){
	echo "
	usage : $0 [options]
	 required:
	 -a <bed> : target interval
	 -b <bam> : read, this program takes the center

	 optional:	
	 -5 <num>,<num>,<num> : relative positions of up/downstream, number of bins 
			from the 5' boundary of the target; default $FIVE
	 -3 <num>,<num>,<num> : .. from the 3' boundary of the target; default: $FOUR
	 -4 <num>,<num>,<num> : .. from the body boundary of the target; default: $THREE
	 -c <num> : column to sum in the read file (default not-used)
		Note) use -c 4 for bam converted from chipchip gff
	"
}

while getopts "hc:a:b:5:4:3:" opt; do
	case $opt in
		h) usage; exit 1;;
		a) FILEA=$OPTARG;;
		b) FILEB=$OPTARG;;
		5) FIVE=$OPTARG;;
		3) THREE=$OPTARG;;
		4) FOUR=$OPTARG;;
		c) CL=$OPTARG;;
		f) FLANK=$OPTARG;;
		?) usage; exit ;;
	esac
done
shift $(($OPTIND-1))

## import hm modules
. hm.sh

read_ext(){
	awk -v OFS="\t" -v FIVE=$1 -v THREE=$2  '{ split(FIVE,l,",");split(THREE,r,","); s=$2+l[1]; if(s<0){s=0;} print $1,s,$3+r[2],$0;} '
}

parse_543(){
	cmd='
import numpy, sys;
five=[ FIVE ]; four=[ FOUR ]; three=[ THREE ];
for line in sys.stdin:
	chrom,start,end,name,score,strand,starts,ends,scores = line.rstrip().split("\t");
	start = int(start); end = int(end);
	starts=[ int(x) for x in starts.split(",") ];
	ends=[ int(x) for x in ends.split(",") ];
	centers =[ (x[0]+x[1])/2 for x in zip(starts,ends)];
	length=end-start;

	if five[1] > length or four[0] > length or - four[1] > length or  - three[0] > length:
		continue;

	## calculate relative bin positions
	bins5,step5 = numpy.linspace(five[0],five[1],five[2]+1,endpoint=True,retstep=True);	
	bins4,step4 = numpy.linspace(four[0],four[1]+length,four[2]+1,endpoint=True,retstep=True);	
	bins3,step3 = numpy.linspace(length+three[0],length+three[1],three[2]+1,endpoint=True,retstep=True);	

	## calculate relative positions and scores
	x= map(lambda x: end - x -1, centers);
	y= map(lambda x: float(x),scores.split(","));
	xy=zip(x,y); xy.sort(key=lambda xx: xx[0]);
	x,y = zip(*xy);

	## histogram
	bins = [bins5,bins4,bins3];
	pos= (5,4,3);
	for i in range(0,len(bins)):
		h,b=numpy.histogram(x,bins=bins[i],weights=y);
		for j in range(0,h.size):
			if strand == "+":
				s=int(round(b[j])) + start;
				e=int(round(b[j+1])) + start;
			else:
				s=end - int(round(b[j+1]));
				e=end - int(round(b[j]));
			print "\t".join( map(str, (chrom, s,e,name,".".join(map(str,(pos[i],j))),strand,h[j] )  ));
	';
	cmd=${cmd/FIVE/$1};
	cmd=${cmd/FOUR/$2};
	cmd=${cmd/THREE/$3};
	tmp=`make_temp`
	echo "$cmd" > $tmp;
	python $tmp
}

## number of fields
LA=`head -n 1 $FILEA | awk '{print NF;}'`
CHROM=( `cat $FILEA | read_ext $FIVE $THREE | groupBy -g 1 -c 2,3 -o min,max` )
tmp=`make_temp`
for (( i=0; i<${#CHROM[@]}; i+=3)){
	C=${CHROM[$i]};S=${CHROM[$i+1]};E=${CHROM[$i+2]};
	echo "$C" >&2
	## extract target in a segment 
	echo -e "$C\t$S\t$E" | intersectBed -a $FILEA -b stdin -wa -u | read_ext $FIVE $THREE > $tmp;
	## collect scores
	reg=`echo -e "$C:$S-$E"`;
	samtools view -b $FILEB $reg | bamToBed \
	| intersectBed -a $tmp -b stdin -wa -wb \
	| cut -f 4- | awk -v OFS="\t" -v LA=$LA '{
		print $1,$2,$3,$4,$5,$6,$((LA+2)),$((LA+3)),$((LA+4));  
	}'\
	| groupBy -g 1,2,3,4,5,6 -c 7,8,9 -o collapse,collapse,collapse \
	| parse_543 $FIVE $FOUR $THREE
}
