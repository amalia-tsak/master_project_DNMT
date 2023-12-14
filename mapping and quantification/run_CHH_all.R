#!/hpc/local/CentOS7/uu_epigenetics/software/miniconda3/envs/quasr/bin/Rscript

print("starting analysis in R")

#
library(BSgenome.Mmusculus.UCSC.mm9)
library(QuasR)
# assign nr of clusters
cluObj=makeCluster(8)


#create a vector holding the names of all the sample names 

sample_names<-c("GSM4594635","GSM1382253","GSM1382254","GSM2533056","GSM4837324","GSM1382255","GSM1382256","GSM1545748","GSM748786","GSM748787")
#create a vector with all the chromosome names
a_vec=c()
for (i in 1:19){
a_vec<-append(a_vec,print(paste("chr",i, sep="")))
}
a_vec<-append(a_vec,"chrX")
a_vec<-append(a_vec,"chrY")
#loop over samples and then over chromosomes to extract chh methylation
for (j in 1:10){
	proj<-readRDS(file=paste("/hpc/shared/uu_epigenetics/Amalia/wgbs/mapped_data/proj_",sample_names[j],".rds",sep = ""))
	for(i in 1:21){
		chromosome<-readRDS(file=paste("/hpc/shared/uu_epigenetics/Amalia/wgbs/mapped_data/",a_vec[i],".rds",sep = ""))
		tr<- qMeth(proj, query= chromosome, mode="allC", collapseBySample=TRUE)
		saveRDS(tr, file=paste("/hpc/shared/uu_epigenetics/Amalia/wgbs/mapped_data/chh_files/",sample_names[j],"_",a_vec[i],".rds",sep = ""))
		}
}
sessionInfo()

