args = commandArgs(trailingOnly=TRUE)
Transcript_FPKM_file = args[1]
Transcript_file = args[2]

Transcript_FPKM = read.delim(Transcript_FPKM_file,header=T,stringsAsFactors=F,row.names=1)[,c("gene_id","FPKM.RPKM")]
Transcript = read.delim(Transcript_file,stringsAsFactors=F,header=F)[,c(1,4)]
Transcript = Transcript[!grepl("chrKI|^GL", Transcript[,1]), 2]

Gene_FPKM=Transcript_FPKM[Transcript,]
FPKM=Gene_FPKM[!duplicated(Gene_FPKM[,1]),"FPKM.RPKM"]

densMode <- function(x){
        v = density(x)
        peaks_x = v$x[which(diff(sign(diff(v$y)))==-2)]
        peaks_y = v$y[which(diff(sign(diff(v$y)))==-2)]
        plot(v, xlab="FPKM", main="")
        text(x=peaks_x[1],y=peaks_y[1],labels=paste0('x=', round(peaks_x[1],3), "\nquantile=" ,100*round((1-mean(x>peaks_x[1])),3),"%"))
        text(x=peaks_x[2],y=peaks_y[2],labels=paste0('x=', round(peaks_x[2],3), "\nquantile=" ,100*round((1-mean(x>peaks_x[2])),3),"%"))        

        Trough_x = v$x[which(diff(sign(diff(v$y)))==2)]
        Trough_y = v$y[which(diff(sign(diff(v$y)))==2)]
        text(x=Trough_x[1],y=Trough_y[1],labels=paste0('x=', round(Trough_x[1],3), "\nquantile=" ,100*round((1-mean(x>Trough_x[1])),3),"%"))
        
        return(peaks_x[1])
}

pdf_file = paste0(Transcript_FPKM_file,".FPKM_cutoff.pdf")
pdf(pdf_file)
x = FPKM[FPKM>0 & FPKM<0.5]
cutoff = densMode(x)
dev.off()

print(round(cutoff,3))

