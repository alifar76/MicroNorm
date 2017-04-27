# MicroNorm

Background
------

This contains scripts to normalize an OTU table (generated via 16S sequencing) on methods other than the singly rarefied approach.

Required R packages
------

- [DESeq](http://bioconductor.org/packages/release/bioc/html/DESeq.html)
- [DESeq2](http://bioconductor.org/packages/release/bioc/html/DESeq2.html)

How to use
------

There are two scripts in the src folder. They are called:

- ```multiple_rarefying.R```
- ```vst_normalization.R```

The script ```multiple_rarefying.R``` is a part of on ongoing project to find alternatives to the [single_rarefaction.py] (http://qiime.org/scripts/single_rarefaction.html) as implemented in QIIME. This script rarefies the OTU table multiple times and then selects the sample which has the minimum average distance to all of its corresponding multiply rarefied sample.

The ```vst_normalization.R``` is based on the DESeq package based RNA-Seq normalization method, as proposed by [McMurdie et al](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3974642/).

The two scripts are to be run via command line using the Rscript command. 


Running the vst_normalization.R script
------

```vst_normalization.R``` script can be run in the following way:

```Rscript vst_normalization.R high_vs_low_otu_table.txt high_low_mapfile.txt Treatment deseq_normalized_otu_table.txt```

As seen from the command, the script takes in 4 arguments. They are as follows:

1) OTU table generated via QIIME (which is called **high_vs_low_otu_table.txt** in the above example)

2) QIIME compatible mapping file (which is called **high_low_mapfile.txt** in the above example)

3) Column name of the category being compared as labelled in the mapping file (which is called **Treatment** in the above example)

4) Output file contaiing result (which is called **deseq_normalized_otu_table.txt** in the above example)

Please ensure that all the 4 arguments are provided, in the correct order and format. Otherwise, the script will crash and cause problems.


Running the multiple_rarefying.R script
------

This work is currently in progress by HFHS investigators (**Sitarik A, Levin A, Havstad S**) and UCSF investigators (**Fujimura K, Lynch S, Faruqi A**).

```multiple_rarefying.R``` script can be run in the following way:

```Rscript multiple_rarefying.R high_vs_low_otu_table.txt multiple_rarefying_results.txt 3 canberra```

As seen from the command, the script takes in 4 arguments. They are as follows:

1) OTU table generated via QIIME (which is called **high_vs_low_otu_table.txt** in the above example)

2) Output file contaiing result (which is called **multiple_rarefying_results.txt** in the above example)

3) Number of iterations to perform for rarefy the OTU table (which is **3** in the above example)

4) Distance metric to use to calculate distance between various iterations of the sample (which is called **canberra** in the above example)

Please ensure that all the 4 arguments are provided, in the correct order and format. Otherwise, the script will crash and cause problems.
