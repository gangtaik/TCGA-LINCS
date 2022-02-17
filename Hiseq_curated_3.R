getwd()
file.loc="/Users/hotcoldbrew/Desktop/OneDrive/TCGA-LINK/"
type.list=dir(path=file.loc)
type.list=type.list[type.list!="Ref"]


for (c.type in type.list){
  #c.type="DLBC"
  library(stringr)
  library(dplyr)
  file2.n=paste0("HiSeqV2_",c.type,"_curated_2.tsv")
  print (c.type)
  setwd(paste0(file.loc,c.type))
  
  cur.df.2=read.table(file2.n,stringsAsFactors = FALSE,header=FALSE,sep="\t",
                      check.names = FALSE, comment.char="")
  colnames(cur.df.2)=paste0("V",seq(1,ncol(cur.df.2)))
  rownames(cur.df.2)=seq(1,nrow(cur.df.2))
  
  gene.id=cur.df.2[,c(1:2)]
  T.df=cbind(gene.id,cur.df.2[,str_detect(cur.df.2[1,],"^T")])
  N.df=cbind(gene.id,cur.df.2[,str_detect(cur.df.2[1,],"^N")])
  R.df=cbind(gene.id,cur.df.2[,str_detect(cur.df.2[1,],"^R")])
  M.df=cbind(gene.id,cur.df.2[,str_detect(cur.df.2[1,],"^M")])
  
  file3.n=paste0("TCGA_RNA-seq_",c.type)
  
  if (ncol(T.df) >=3){
    write.table(T.df,file=paste0(file3.n,"_T.tsv"),sep="\t",quote=FALSE,col.names = FALSE,
                row.names = FALSE)
  }else{
    cat (c.type,"HiSeqV2.gz do not include Tumor sample\n",sep=" ")
  }
  if (ncol(N.df) >=3){
    write.table(N.df,file=paste0(file3.n,"_N.tsv"),sep="\t",quote=FALSE,col.names = FALSE,
                row.names = FALSE)
  }else{
    cat (c.type,"HiSeqV2.gz do not include Normal sample\n",sep=" ")
  }
  if (ncol(R.df) >=3){
    write.table(R.df,file=paste0(file3.n,"_R.tsv"),sep="\t",quote=FALSE,col.names = FALSE,
                row.names = FALSE)
  }else{
    cat (c.type,"HiSeqV2.gz do not include Recurrent sample\n",sep=" ")
  }
  if (ncol(M.df) >=3){
    write.table(M.df,file=paste0(file3.n,"_M.tsv"),sep="\t",quote=FALSE,col.names = FALSE,
                row.names = FALSE)
  }else{
    cat (c.type,"HiSeqV2.gz do not include Metastated sample\n",sep=" ")
  }
  
}
print ("Step4 is Finished")
