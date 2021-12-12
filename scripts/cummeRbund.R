
# pre-installs. run only once.
# install.packages("BiocManager”)
# BiocManager::install("cummeRbund")

library(cummeRbund)

cuff_data<-readCufflinks("cuffdiff_output")

# density plot
csDensity(genes(cuff_data))

# volcano plot 1
csVolcano(genes(cuff_data), 'q1', 'q2')

# volcano plot 2
csVolcano(genes(cuff_data), 'q1', 'q2', alpha=0.05, 
          showSignificant=TRUE, features=FALSE, xlimits=c(-8,8))

# scatter plot
csScatter(genes(cuff_data),'q1','q2')

# bar plot
sig_gene<-getGene(cuff_data, "SNZ1")
expressionBarplot(sig_gene)
