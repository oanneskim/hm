#!/bin/env bash 

quote(){
    perl -ne 'chomp; my @a = map{ "\"$_\"" } split /,/,$_; print join ",",@a; ';
}

test_ei(){
    ## INPUT: comma separated control and treatment EI file ( bed6 + exclusion + inclusion)
    ## OUTPUT: bed6 + logFC + pvalue 
    a=`echo $1 | quote`
    b=`echo $2 | quote`

    rcmd='
    fa=c(FILEA)
    fb=c(FILEB)
    #fa=c("Events/Wt1/exons_jc.bed","Events/Wt2/exons_jc.bed")  
    #fb=c("Events/C41/exons_jc.bed","Events/C42/exons_jc.bed")  
    out="OUT"

    group=factor(c(rep(1,length(fa)),rep(2,length(fb))));
    D=NULL;
    i=1;
    for( f in c(fa,fb)){
        tt=read.table(f,header=F);
        colnames(tt)=c("chr","start","end","name","score","strand",paste(i,c("exc","inc"),sep="."))
        if(is.null(D)){ D=tt;
        }else{ D=merge(D,tt,by=1:6,all=T); }
        i=i+1;
    }
    D[is.na(D)]=0;
    ix=apply(D[,7:ncol(D)], 1, min) > 0 & apply(D[,7:ncol(D)],1,max) > 10
    d=D[ix,]
    #d[d==0]=0.5
    m=cbind(d[,grep("exc",colnames(d))],d[,grep("inc",colnames(d))])

    library(edgeR)
    #j=ncol(m)/2; #y=DGEList(counts=m[,1:j]+m[,(j+1):ncol(m)],group=group)
    y=DGEList(counts=m,group=rep(group,2))
    y=calcNormFactors(y);

    event.this=factor(rep(1:2,each=length(group)));
    group.this=factor(rep(group,2));
    H1 <- model.matrix(~ event.this + group.this + event.this:group.this, data=y$samples )
    H0 <- model.matrix(~ event.this + group.this )
    coef <- (ncol(H0)+1):ncol(H1)
    #y=estimateCommonDisp(y)
    #y=estimateTagwiseDisp(y, trend="movingave")
    y = estimateGLMCommonDisp(y,H1);
    y = estimateGLMTrendedDisp(y,H1);
    y = estimateGLMTagwiseDisp(y,H1);

    fit=glmFit(m, H1, y$tagwise.dispersion,offset=0,prior.count=0)
    llh=glmLRT(fit,coef=coef)

    ex.h0=apply( m[,group.this == 1 & event.this == 1], 1, sum);
    in.h0=apply( m[,group.this == 1 & event.this == 2], 1, sum);

    res=data.frame(d[,1:6], logIR=log( in.h0/ ex.h0), logFC=llh$table$logFC, pval=llh$table$PValue)
    ## chrom start end logFC pval
    write.table(res, out, col.names=T,row.names=F,sep="\t",quote=F);
    '

    tmp=`mktemp`
    rcmd=${rcmd/FILEA/$a}
    rcmd=${rcmd/FILEB/$b}
    rcmd=${rcmd/OUT/$tmp}
    echo "$rcmd" | R --no-save >&2
    cat $tmp
}
