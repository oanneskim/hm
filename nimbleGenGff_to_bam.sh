
nimbleGenGff_to_bed(){
        awk -v OFS="\t" -v BIN=$BIN '{ if(NF > 5){
                print $1,$4,$5,$6,0,$7; ## store score at the name field
        }}'
}

cat $1 | nimbleGenGff_to_bed | bedToBam -i stdin -g $2 | samtools sort -o -

