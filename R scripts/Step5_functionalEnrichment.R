# Introduction to this script ----
# for the purposes of this script we'll want several data objects generated in previous scripts, including:
# 1) your normalized filtered expression data, in the form of a data matrix with symbols as rownames.
# 2) your study design file
# 3) your contrast matrix that lays out the pairwise comparisons you're interested in testing
# 4) Individual signatures or 'collections' of signatures to test for enrichment in your data.
# These signatures can be downloaded from gene signature databases such as MSigDB
# Signatures can also be custom made based on your interests.
# Signatures can also be pulled from R/Bioconductor as described below

# Sys.unsetenv("R_LIBS_USER")
# .libPaths()
# .libPaths(paste(getwd(), "RLibrary", sep="/"))
# setRepositories()

# install.packages("gplots")
# install.packages("GSEABase")
# install.packages("Biobase")
# install.packages("GSVA")
# install.packages("gprofiler2")
# install.packages("clusterProfiler")
# install.packages("msigdbr")
# install.packages("enrichplot")
# install.packages("qusage")
# install.packages("patchwork")

# Load packages ----
library(tidyverse)
library(limma)
library(gplots) #for heatmaps
library(DT) #interactive and searchable tables of our GSEA results
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
library(enrichplot) # great for making the standard GSEA enrichment plots
library(qusage) # Quantitative Set Analysis for Gene Expression
library(heatmaply)


# Carry out GO enrichment using gProfiler2 ----
# use topTable result to pick the top genes for carrying out a Gene Ontology (GO) enrichment analysis
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=50, sort.by="logFC")
# use the 'gost' function from the gprofiler2 package to run GO enrichment analysis, the list of available species for the "organism" parameters is here: https://biit.cs.ut.ee/gprofiler/page/organism-list 
gost.res <- gost(rownames(myTopHits), organism = "hsapiens", correction_method = "fdr")
# produce an interactive manhattan plot of enriched GO terms
gostplot(gost.res, interactive = T, capped = F) #set interactive=FALSE to get plot for publications
mygostplot <- gostplot(gost.res, interactive = F, capped = F)
# produce a publication quality static manhattan plot with specific GO terms highlighted
# rerun the above gostplot function with 'interactive=F' and save to an object 'mygostplot'
publish_gostplot(
  mygostplot, #your static gostplot from above
  highlight_terms = c("GO:0034987"),
  filename = NULL,
  width = NA,
  height = NA)

#you can also generate a table of your gost results
# publish_gosttable(
#   gost.res,
#   highlight_terms = NULL,
#   use_colors = TRUE,
#   show_columns = c("source", "term_name", "term_size", "intersection_size"),
#   filename = NULL,
#   ggplot=TRUE)
# now repeat the above steps using only genes from a single module from the step 6 script, by using `rownames(myModule)`
# what is value in breaking up DEGs into modules for functional enrichment analysis?


# Competitive GSEA using CAMERA----
# for competitive tests the null hypothesis is that genes in the set are, at most, as often differentially expressed as genes outside the set
# first let's create a few signatures to test in our enrichment analysis
mySig <- rownames(myTopHits) #if your own data is from mice, wrap this in 'toupper()' to make gene symbols all caps
mySig2 <- sample((rownames(v.DEGList.filtered.norm$E)), size = 50, replace = FALSE)
collection <- list(real = mySig, fake = mySig2)
# now test for enrichment using CAMERA
camera.res <- camera(v.DEGList.filtered.norm$E, collection, design, contrast.matrix[,1]) 
camera.df <- as_tibble(camera.res, rownames = "setName")
camera.df

# Self-contained GSEA using ROAST----
# remember that for self-contained the null hypothesis is that no genes in the set are differentially expressed
mroast(v.DEGList.filtered.norm$E, collection, design, contrast=1) #mroast adjusts for multiple testing

# now repeat with an actual gene set collection
# camera requires collections to be presented as a list, rather than a tibble, so we must read in our signatures using the 'getGmt' function
broadSet.C2.ALL <- getGmt("c2.cp.kegg.v2023.1.Hs.symbols.gmt", geneIdType=SymbolIdentifier())
# Option for plants: goto http://structuralbiology.cau.edu.cn/PlantGSEA/download.php
# Download the gene set of interest for the species of interest (e.g. Ara_KEGG.txt)
broadSet.C2.ALL <- getGmt("Ara_KEGG.txt",geneIdType=SymbolIdentifier())
broadSet.C2.ALL <- geneIds(broadSet.C2.ALL)

#extract as a list
broadSet.C2.ALL <- geneIds(broadSet.C2.ALL)
camera.res <- camera(v.DEGList.filtered.norm$E, broadSet.C2.ALL, design, contrast.matrix[,1]) 
camera.df <- as_tibble(camera.res, rownames = "setName")
camera.df

# filter based on FDR and display as interactive table
camera.df <- filter(camera.df, FDR<=0.01)

datatable(camera.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Signatures enriched in leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2,4,5), digits=2)

#as before, add a variable that maps up/down regulated pathways with phenotype
camera.df <- camera.df %>%
  mutate(phenotype = case_when(
    Direction == "Up" ~ "disease",
    Direction == "Down" ~ "healthy"))

#easy to filter this list based on names of signatures using 'str_detect'
#here is an example of filtering to return anything that has 'CD8' or 'CYTOTOX' in the name of the signature
camera.df.sub <- camera.df %>%
  dplyr::filter(str_detect(setName, "CD8|CYTOTOX"))

# graph camera results as bubble chart 
ggplot(camera.df[1:25,], aes(x=phenotype, y=setName)) + 
  geom_point(aes(size=NGenes, color = Direction, alpha=-log10(FDR))) +
  theme_bw()

# Single sample GSEA using the GSVA package----
# the GSVA package offers a different way of approaching functional enrichment analysis.  
# A few comments about the approach:
# In contrast to most GSE methods, GSVA performs a change in coordinate systems,
# transforming the data from a gene by sample matrix to a gene set (signature) by sample matrix. 
# this allows for the evaluation of pathway enrichment for each sample.
# the method is both non-parametric and unsupervised
# bypasses the conventional approach of explicitly modeling phenotypes within enrichment scoring algorithms. 
# focus is therefore placed on the RELATIVE enrichment of pathways across the sample space rather than the absolute enrichment with respect to a phenotype. 
# however, with data with a moderate to small sample size (< 30), other GSE methods that explicitly include the phenotype in their model are more likely to provide greater statistical power to detect functional enrichment.

# be aware that if you choose a large MsigDB file here, this step may take a while
GSVA.res.C2CP <- gsva(v.DEGList.filtered.norm$E, #your data
                      broadSet.C2.ALL, #signatures
                      min.sz=3, max.sz=500, #criteria for filtering gene sets
                      mx.diff=FALSE,
                      method="gsva") #options for method are "gsva", ssgsea', "zscore" or "plage"

# Apply linear model to GSVA result
# now using Limma to find significantly enriched gene sets in the same way you did to find diffGenes
# this means you'll be using topTable, decideTests, etc
# note that you need to reference your design and contrast matrix here
fit.C2CP <- lmFit(GSVA.res.C2CP, design)
ebFit.C2CP <- eBayes(fit.C2CP)

# use topTable and decideTests functions to identify the differentially enriched gene sets
topPaths.C2CP <- topTable(ebFit.C2CP, adjust ="BH", coef=1, number=50, sort.by="logFC")
res.C2CP <- decideTests(ebFit.C2CP, method="global", adjust.method="BH", p.value=0.05, lfc=0)
# the summary of the decideTests result shows how many sets were enriched in induced and repressed genes in all sample types
summary(res.C2CP)

# pull out the gene sets that are differentially enriched between groups
diffSets.C2CP <- GSVA.res.C2CP[res.C2CP[,1] !=0,]
head(diffSets.C2CP)
dim(diffSets.C2CP)



# Create a heatmap of differentially expressed genes ----

heatmaply(diffSets.C2CP, 
          #dendrogram = "row",
          xlab = "Samples", ylab = "KEGG pathways", 
          main = "Responsive KEGG pathways in cutaneous leishmaniasis",
          scale = "column",
          margins = c(60,100,40,20),
          grid_color = "white",
          grid_width = 0.0000001,
          titleX = T,
          titleY = T,
          hide_colorbar = TRUE,
          branches_lwd = 0.1,
          label_names = c("Gene", "Sample:", "Value"),
          fontsize_row = 5, fontsize_col = 5,
          labCol = colnames(diffSets.C2CP),
          labRow = rownames(diffSets.C2CP),
          heatmap_layers = theme(axis.line=element_blank())
)










# the essentials ----
library(tidyverse)
library(limma)
library(gplots) #for heatmaps
library(DT) #interactive and searchable tables of our GSEA results
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
library(enrichplot) # great for making the standard GSEA enrichment plots
gost.res_up <- gost(rownames(myModule_up), organism = "hsapiens", correction_method = "fdr")
gostplot(gost.res_up, interactive = T, capped = T) #set interactive=FALSE to get plot for publications
gost.res_down <- gost(rownames(myModule_down), organism = "hsapiens", correction_method = "fdr")
gostplot(gost.res_down, interactive = T, capped = T) #set interactive=FALSE to get plot for publications

hs_gsea_c2 <- msigdbr(species = "Homo sapiens", # change depending on species your data came from
                      category = "C2") %>% # choose your msigdb collection of interest
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature 

# Now that you have your msigdb collections ready, prepare your data
# grab the dataframe you made in step3 script
# Pull out just the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
mydata.df.sub <- dplyr::select(mydata.df, geneID, LogFC)
mydata.gsea <- mydata.df.sub$LogFC
names(mydata.gsea) <- as.character(mydata.df.sub$geneID)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)

# run GSEA using the 'GSEA' function from clusterProfiler
myGSEA.res <- GSEA(mydata.gsea, TERM2GENE=hs_gsea_c2, verbose=FALSE)
myGSEA.df <- as_tibble(myGSEA.res@result)

# view results as an interactive table
datatable(myGSEA.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Signatures enriched in leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(3:10), digits=2)
# create enrichment plots using the enrichplot package
gseaplot2(myGSEA.res, 
          geneSetID = 47, #can choose multiple signatures to overlay in this plot
          pvalue_table = FALSE, #can set this to FALSE for a cleaner plot
          title = myGSEA.res$Description[47]) #can also turn off this title

# add a variable to this result that matches enrichment direction with phenotype
myGSEA.df <- myGSEA.df %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "disease",
    NES < 0 ~ "healthy"))

# create 'bubble plot' to summarize y signatures across x phenotypes
ggplot(myGSEA.df[1:20,], aes(x=phenotype, y=ID)) + 
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()



------------------------------------------------------------------------
# Step 4 Gene Ontology (GO) enrichment using gProfiler2 and Heatmapping

# installing packages and loading all required libraries 

BiocManager::install("GSEABase")
BiocManager::install("GSVA")
BiocManager::install("clusterProfiler")
BiocManager::install("msigdbr")
BiocManager::install("enrichplot")

library(tidyverse)
library(limma)
library(gplots) #for heatmaps
library(DT) #interactive and searchable tables of our GSEA results
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
library(enrichplot) # great for making the standard GSEA enrichment plots
myTop100Hits <- topTable(ebFit, adjust ="BH", coef=1, number=100, sort.by="logFC") # pick the top 100 genes by logFC value for carrying out GO enrichment analysis
gost.res <- gost(rownames(myTop100Hits), organism = "slycopersicum", correction_method = "fdr") # run GO enrichment analysis
gostplot(gost.res, interactive = T, capped = F) # interactive manhattan plot of enriched GO terms


#save gost plot and table
publish_gostplot(
  gostplot(gost.res, interactive = F, capped = F), # static gostplot
  highlight_terms = c("GO:0005975", "GO:0048046", "GO:0030145", "GO:0003824", "GO:0045735"), # highlight top 7 lowest p-values
  filename = "gostplot100.png",
  width = NA,
  height = NA)

publish_gosttable(
  gost.res$result[order(gost.res$result$p_value, decreasing = F),],
  highlight_terms = NULL,
  use_colors = TRUE,
  show_columns = c("source", "term_name", "term_size", "intersection_size"),
  filename = "gosttable100.pdf",
  ggplot=TRUE)

# create tissue modules based on pearson correlation
clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") # cluster rows by pearson correlation
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete") # cluster columns by pearson correlation
module.assign <- cutree(clustRows, k=2) 

modulePick <- 2 # pick the drought module
myModule_dr <- diffGenes[names(module.assign[module.assign %in% modulePick]),]
length(rownames(myModule_dr))

gost.res.dr <- gost(rownames(myModule_dr), organism = "slycopersicum", correction_method = "fdr") # run GO enrichment analysis
gostplot(gost.res.dr, interactive = T, capped = F) # interactive Manhattan plot of enriched GO terms

# save gost plot and table
publish_gostplot(
  gostplot(gost.res.dr, interactive = F, capped = F), # static gostplot
  highlight_terms = c("GO:0004866", "GO:0061135", "GO:0061134", "GO:0030414", "GO:0004867"), # highlight top 5
  filename = "gostplot_drought.png",
  width = NA,
  height = NA)

publish_gosttable(
  gost.res.dr$result[order(gost.res.dr$result$p_value, decreasing = F),],
  highlight_terms = NULL,
  use_colors = TRUE,
  show_columns = c("source", "term_name", "term_size", "intersection_size"),
  filename = "gosttable_drought.pdf",
  ggplot=TRUE)


modulePick <- 1 # pick the control module
myModule_healthy <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
length(rownames(myModule_healthy))

gost.res.healthy <- gost(rownames(myModule_healthy), organism = "slycopersicum", ordered_query = T, correction_method = "fdr") # run GO enrichment analysis
gostplot(gost.res.healthy, interactive = T, capped = F) # interactive Manhattan plot of enriched GO terms

# save gost plot and table
publish_gostplot(
  gostplot(gost.res.healthy, interactive = F, capped = F), # static gostplot
  highlight_terms = c("GO:0009611", "GO:0034214", "GO:0004177", "GO:0034214", "GO:0070006"), # highlight top 5
  filename = "gostplot_healthy.png",
  width = NA,
  height = NA)

publish_gosttable(
  gost.res.stem$result[order(gost.res.stem$result$p_value, decreasing = F),],
  highlight_terms = NULL,
  use_colors = TRUE,
  show_columns = c("source", "term_name", "term_size", "intersection_size"),
  filename = "gosttable_stem.pdf",
  ggplot= T)

# ----------- Heatmaps -----------------------------------------------------------------------------------------------
myheatcolors <- rev(brewer.pal(name="Spectral", n=11)) # define a color palette
module.color <- rainbow(length(unique(module.assign)), start=0.3, end=0.2) 
module.color <- module.color[as.vector(module.assign)] # assign module colors

# generate a heatmap for all DEGs
heatmap.2(diffGenes, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=myheatcolors, scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20))

# generate a heatmap for the drought species
modulePick <- 2 
hrsub_dr <- hclust(as.dist(1-cor(t(myModule_dr), method="pearson")), method="complete") 
RowSideColors = module.color[module.assign %in% modulePick]

heatmap.2(myModule_dr, 
          Rowv = as.dendrogram(hrsub_dr), 
          Colv = FALSE,  # Set Colv explicitly to FALSE
          labRow = NA,
          col = myheatcolors, scale = "row", 
          density.info = "none", trace = "none", 
          RowSideColors = module.color[module.assign %in% modulePick], 
          margins = c(8, 20))


# generate a heatmap for the healthy species
modulePick <- 1
hrsub_healthy <- hclust(as.dist(1-cor(t(myModule_healthy), method="pearson")), method="complete") 
heatmap.2(myModule_healthy, 
          Rowv=as.dendrogram(hrsub_healthy), 
          Colv=NA, 
          labRow = NA,
          col=myheatcolors, scale="row", 
          density.info="none", trace="none", 
          RowSideColors=module.color[module.assign%in%modulePick], margins=c(8,20))
