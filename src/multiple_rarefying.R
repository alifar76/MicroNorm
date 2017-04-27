##############################################################################################################
######  This is a script to calculate a representative rarefied OTU table from an unrarefied OTU table, ######
######  an alterative to single rarefied tables that stabilizes the random sampling and results         ######
######  in a rarefied table that may be more consistent with the original data than a single            ######  	     
######  rarefied table. Briefly, many single-rarefied OTU tables are calculated, and the distance       ######
######  between the subject-specific rarefied vectors is calculated.  The rarefied vector that is       ######
######  the minimum average (or median) distance from itself to all other rarefied vectors              ######
######  is considered the most representative for that subject and built into the new                   ######
######  rarefied table.  This work is currently in progress by HFHS investigators                       ######
######  (Sitarik A, Levin A, Havstad S) and UCSF investigators (Fujimura K, Lynch S, Faruqi A).         ######
##############################################################################################################


#################################

#Rscript multiple_rarefying.R high_vs_low_otu_table.txt multiple_rarefying_results.txt 3 canberra

start.time <- Sys.time()
library(GUniFrac)

Rarefy <- function (otu.tab, depth = min(rowSums(otu.tab))) 
{
    otu.tab <- as.matrix(otu.tab)
    ind <- (rowSums(otu.tab) < depth)
    sam.discard <- rownames(otu.tab)[ind]
    otu.tab <- otu.tab[!ind, ]
    rarefy <- function(x, depth) {
        y <- sample(rep(1:length(x), x), depth)
        y.tab <- table(y)
        z <- numeric(length(x))
        z[as.numeric(names(y.tab))] <- y.tab
        z
    }
    otu.tab.rff <- t(apply(otu.tab, 1, rarefy, depth))
    rownames(otu.tab.rff) <- rownames(otu.tab)
    colnames(otu.tab.rff) <- colnames(otu.tab)
    return(list(otu.tab.rff = otu.tab.rff, discard = sam.discard))
}


reprare <- function(rawtab=otu_tab_t, ntables=100, distmethod="euclidean", 
summarymeasure=mean, seedstart=500, verbose=TRUE) {

    raretabs = list()
    for (z in 1:ntables) {
        if (verbose==TRUE) {
            print(paste("calculating rarefaction table number", z, sep=" "))
        }
        set.seed(seedstart + z)
        raretabs[[z]] = Rarefy(rawtab, depth = min(rowSums(rawtab)))[[1]]
    }
    raretabsa = array(unlist(raretabs), dim = c(nrow(rawtab), ncol(rawtab), ntables))
    final_tab = c()
    for (y in 1:nrow(rawtab)) {
    if (verbose==TRUE) {
        print(paste("determining rep rarefaction for subject number", y, sep=" "))
    }
    distmat = as.matrix(vegdist(t(raretabsa[y,,]), method=distmethod)) # distance across reps for subject y
    distsummary = apply(distmat, 2, summarymeasure) 
    whichbestrep = which(distsummary == min(distsummary))[1]  # the best rep is the one with the minimum average/median distance to all other reps. (in case of ties, just select the first)
    bestrep = raretabsa[y,,whichbestrep] # select that rep only for subject y
    final_tab = rbind(final_tab, bestrep) # build that rep for subject y into final table
    }
    rownames(final_tab) = rownames(rawtab)
    colnames(final_tab) = colnames(rawtab)
    return(final_tab)
}


inout_file <- function(inputfile,outputfile,ntables,distm){
	MYdata <- read.table(inputfile,header = T, sep = "\t", check.names = T, row.names =1, comment.char= "", skip =1)     # Ignore # at beginning of line. Skip first line which says "Converted from biom"
	allsamples = length(colnames(MYdata))-1    
	MYdata.1<-MYdata[,1:allsamples]
	MYarray<-as.data.frame(t(MYdata.1))
	result <- reprare(rawtab=MYarray,ntables,distmethod=distm)			#ntables=3 just for testing. Default is 100. distmethod=canberra and default is Euclidean.
	taxonomy <- as.vector(MYdata[,"taxonomy"])
	final <- cbind(t(result),taxonomy)
	final <- cbind(rownames(final),final)
	colnames(final) <- c("#OTU ID",colnames(final)[2:ncol(final)])
	write.table("# Constructed from biom file",file = outputfile, append=TRUE,quote = FALSE,sep = "\t",row.names = FALSE,col.names=FALSE)
	suppressWarnings(write.table(final, file = outputfile,row.names = FALSE,append = TRUE ,quote = FALSE, sep = "\t"))
}


# 4 arguments needed for the script: Non-rarefied OTU table, Output file name, # of iterations to perform, Distance metric to use

argv <- commandArgs(TRUE)
inout_file(argv[1],argv[2],as.numeric(argv[3]),argv[4])
print (Sys.time() - start.time)
