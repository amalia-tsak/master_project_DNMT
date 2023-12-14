#!/hpc/local/CentOS7/uu_epigenetics/software/miniconda3/envs/quasr/bin/Rscript

print("starting analysis in R")
sampleFile_am <- as.data.frame(matrix(NA,2,2))
names(sampleFile_am) <- c("FileName","SampleName")
sampleFile_am$FileName<-c("/hpc/shared/uu_epigenetics/Amalia/wgbs/trimgalore_data/GSM2533056_chopped_1_clumped_trimmed.fq.gz","/hpc/shared/uu_epigenetics/Amalia/wgbs/trimgalore_data/GSM2533056_chopped_2_clumped_trimmed.fq.gz")
sampleFile_am$SampleName<-"GSM2533056"
write.table(sampleFile_am, sep="\t", quote=F, row.names=F, "samples_GSM2533056.tab")

#
library(BSgenome.Mmusculus.UCSC.mm9)
library(QuasR)
# assign nr of clusters
cluObj=makeCluster(16)

#print(getwd())
proj=qAlign(sampleFile="/hpc/shared/uu_epigenetics/Amalia/wgbs/scripts/samples_GSM2533056.tab",
       genome="BSgenome.Mmusculus.UCSC.mm9",
       aligner="Rbowtie",
       bisulfite="undir",
       projectName="TKO_DNMT",
       alignmentsDir="/hpc/shared/uu_epigenetics/Amalia/wgbs/mapped_data",
       clObj=cluObj,
       cacheDir = "/hpc/shared/uu_epigenetics/Amalia/wgbs/temp_data")

qQCReport(proj, "QC_report_GSM2533056.pdf")
saveRDS(proj,"/hpc/shared/uu_epigenetics/Amalia/wgbs/mapped_data/proj_GSM2533056.rds")

## extract stranded CpG data

GSM2533056_add_mCpG <- qMeth(proj, mode="CpG", collapseBySample=TRUE)
saveRDS(GSM2533056_add_mCpG, "/hpc/shared/uu_epigenetics/Amalia/wgbs/mapped_data/GSM2533056_add_mCpG.rds")


sessionInfo()
