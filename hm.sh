#!/bin/bash
usage(){
    echo "
	Author: Hyunmin Kim (hyun.kim@ucdenver.edu)
	Descript: this is hyunmin's humble programs

    usage : $0 [options]
     -t : test 
	 -d : debug
     -h : help 
    "
}



DEBUG="no";
while getopts "hd" opt; do
    case $opt in
        h) usage; exit 1;;
        d) DEBUG="yes";;
        ?) usage; exit ;;
    esac
done
shift $(($OPTIND-1))


## hyunmin's most frequent chore tools
make_temp(){
    mktemp 2>/dev/null || mktemp -t ${0##*/}
}

require(){
	prog=$1; 
	if [ !`which $prog` ];then 
		echo "$1 program is ready" >&2; 
	else
		echo "$1 program is not found" >&2; 
		exit -1;
	fi	
}

## some converters
nimbleGenGff_to_bam(){
	require samtools
	require bedToBam
	cat $1 | awk -v OFS="\t" -v BIN=$BIN '{ if(NF > 5){
		print $1,$4,$5,$6,0,$7; ## store score at the name field
	}}' | sort -k 1,1 -k 2,2n | bedToBam -i stdin -g $2 
}

if [ $DEBUG == "yes" ];then
	nimbleGenGff_to_bam data/nimblegen_chipchip_sample0.gff data/scer3_chrom.size > data/nimblegen_chipchip_sample0.bam
fi



