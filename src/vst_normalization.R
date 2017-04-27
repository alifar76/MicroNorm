#Rscript vst_normalization.R high_vs_low_otu_table.txt high_low_mapfile.txt Treatment deseq_normalized_otu_table.txt

start.time <- Sys.time()
require(DESeq2)
require(DESeq)

dispersion_normalization <- function(otutable,mapfile,metavariable,outputname){
	MYdata <- read.table(otutable,header = T, sep = "\t", check.names = T, row.names =1, comment.char= "", skip =1)
	# Ignore # at beginning of line. Skip first line which says "Converted from biom"
	# Also for MYdata variable, check.names = F, if the names of samples are numbers (such as 1, 2, 3, etc)
	mapfile <- read.table(mapfile,header = T, sep = "\t", check.names = T, row.names =1, comment.char= "")
	allsamples = length(colnames(MYdata))-1
	MYdata2 <- as.matrix(MYdata[,1:allsamples])
	MYdata2 <- subset(MYdata2,select=rownames(mapfile))
	cds = newCountDataSet(MYdata2, conditions=as.vector(mapfile[,metavariable]))
	# Estimate size factor for each sample (http://seqanswers.com/forums/showpost.php?p=16468&postcount=13)
	cds = estimateSizeFactors(cds)
	# Variable to store size factor
	sizef <- sizeFactors(cds)
	# Normalization done by dividing count of each OTU in a given sample by the sample size factor
	cts = counts(cds, normalized=TRUE)
	cts = cbind(rownames(cts),cts,as.vector(MYdata$taxonomy))
	finalcols <- c("#OTU ID",colnames(cts)[2:(length(colnames(cts))-1)],"taxonomy") 
	colnames(cts) <- finalcols
	suppressWarnings(write.table("# Constructed from biom file",file=outputname,sep="",append=TRUE,row.names=FALSE,col.names=FALSE,quote=FALSE))
	suppressWarnings(write.table(as.matrix(cts), file = outputname,sep = "\t",append=TRUE,row.names=FALSE,quote=FALSE))	
}

# 4 arguments to function: OTU table, Mapping file, Meta-data variable name,Output table name

argv <- commandArgs(TRUE)
dispersion_normalization(argv[1],argv[2],argv[3],argv[4])
print (Sys.time() - start.time)
