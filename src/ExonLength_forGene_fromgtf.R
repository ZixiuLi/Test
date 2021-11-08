library(GenomicFeatures)
args = commandArgs(trailingOnly=TRUE)
gtf_input=args[1]
gtf_output=args[2]

txdb = makeTxDbFromGFF(gtf_input,format="gtf")
exons_gene = exonsBy(txdb,by="gene")
exons_gene_lens=lapply(exons_gene,function(x){sum(width(reduce(x)))})

gene_name=names(exons_gene_lens)
exons_gene_lens_final=unlist(exons_gene_lens)
names(exons_gene_lens_final)=gene_name
write.table(exons_gene_lens_final,gtf_output,sep="\t",quote=F,col.names=F)
