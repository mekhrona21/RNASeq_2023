# Step 3 Data Analysis

# installing packages and loading all required libraries 

install.packages("plotly")
install.packages("gt")
install.packages("gplots")
install.packages("gprofiler2")


library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science
library(DT) # interactive tables
library(plotly) # interactive plots
library(gt) # A layered 'grammar of tables'
library(limma) # package for differential gene expression using linear modeling
library(edgeR) # package for differential expression analysis
library(gplots) 
library(gprofiler2) # tools for accessing the GO enrichment results using g:Profiler web resources
library(RColorBrewer)


# Performing Hierarchical clustering 
distance_euc <- dist(t(log2.cpm.filtered.norm), method = "euclidean")  # compute distances based on the normalized and filtered data
clusters_euc <- hclust(distance_euc, method = "complete") # cluster the data based on distance
c_euc <- plot(clusters_euc, labels=sampleLabels) # plot the dendrogram

distance_max <- dist(t(log2.cpm.filtered.norm), method = "maximum") # compute distances based on the normalized and filtered data
clusters_max <- hclust(distance_max, method = "complete") # cluster the data based on distance
c_max <- plot(clusters_max, labels=sampleLabels) # plot the dendrogram

# Performing Principle Component Analysis (PCA)
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
pc.var <- pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1) 
pca.res.df <- as_tibble(pca.res$x)
pca.plot <- ggplot(pca.res.df) +
  aes(x = PC1, y = PC2, label = sampleLabels, color = as.character(group)) +
  geom_point(size = 4) +
  stat_ellipse() +
  xlab(paste0("PC1 (", pc.per[1], "%)")) +
  ylab(paste0("PC2 (", pc.per[2], "%)")) +
  labs(
    title = "PCA plot",
    caption = paste0("produced on ", Sys.time())
  ) +
  coord_fixed() +
  theme_bw()

ggplotly(pca.plot) 

# Gene expression averages

mydata.df <- mutate(log2.cpm.filtered.norm.df,
                    healthy.AVG = (SRR7652567 + SRR7652566 + SRR7652565)/3, 
                    drought.AVG = (SRR7652569 + SRR7652568 + SRR7652563)/3,
                    #making columns comparing each of the averages above that you're interested in
                    LogFC = (drought.AVG - healthy.AVG)) %>% #note that this is the first time you've seen the 'pipe' operator
  mutate_if(is.numeric, round, 2)

write_tsv(mydata.df[,c(1,6:8)], "avg.tsv") #saving table as TSV

datatable(mydata.df[,c(1,6:8)], # view as a table
          extensions = c('KeyTable', "FixedHeader"), 
          filter = 'top',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         dom = "Blfrtip", 
                         buttons = c("copy", "csv", "excel"),
                         lengthMenu = c("10", "25", "50", "100")))

myplot <- ggplot(mydata.df) +
  aes(x=healthy.AVG, y=drought.AVG, 
      text = paste("Symbol:", geneID)) +
  geom_point(shape=16, size=1) +
  ggtitle("healthy vs. drought") +
  theme_bw()

ggplotly(myplot)

# Identification of differentially expressed genes (DEGs) 

group <- factor(targets$group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE)
fit <- lmFit(v.DEGList.filtered.norm, design)
contrast.matrix <- makeContrasts(infection = drought - healthy,
                                 levels=design)

fits <- contrasts.fit(fit, contrast.matrix)
ebFit <- eBayes(fits)
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")

vplot <- ggplot(myTopHits) + 
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", rownames(myTopHits))) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", linewidth=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="green", linewidth=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour= "red", linewidth=1) +
  annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="green") +
  annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="red") +
  labs(
    title = "Volcano plot",
    subtitle = "Solanum lycopersicum",
    caption = paste0("produced on ", Sys.time())
  ) +
  theme_bw()


ggplotly(vplot) #creating volcano plot

results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.01, lfc=7)
# take a look at what the results of decideTests looks like
head(results)
summary(results)
vennDiagram(results, include="both")



results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=1) # output -1 or 1: t-statistic for gene is classified as significant (0 = not significant)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,] # E: numeric matrix of normalized expression values on the log2 scale
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")
write_tsv(diffGenes.df,"DiffGenes.tsv")

datatable(diffGenes.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: DEGs in Solanum lycopersicum',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:6), digits=2)

# Gene Ontology enrichment analysis using gProfiler2

myTop100Hits <- topTable(ebFit, adjust ="BH", coef=1, number=100, sort.by="logFC") # pick the top 100 genes by logFC value for carrying out GO enrichment analysis
gost.res <- gost(rownames(myTop100Hits), organism = "slycopersicum", correction_method = "fdr") # run GO enrichment analysis
gostplot(gost.res, interactive = T, capped = F) # interactive manhattan plot of enriched GO terms

publish_gosttable(
  gost.res$result[order(gost.res$result$p_value, decreasing = F),],
  highlight_terms = NULL,
  use_colors = TRUE,
  show_columns = c("source", "term_name", "term_size", "intersection_size"),
  filename = "gosttable100.pdf",
  ggplot=TRUE)

publish_gostplot(
  gostplot(gost.res, interactive = F, capped = F), # static gostplot
  highlight_terms = c("GO:0019346", "GO:0050667", "GO:0009092", "GO:0008238", "GO:0008237"), # highlight top 5 lowest p-values
  filename = "gostplot100.png",
  width = NA,
  height = NA)

# create tissue modules based on pearson correlation
clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") # cluster rows by pearson correlation
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete") # cluster columns by pearson correlation
module.assign <- cutree(clustRows, k=2) 

modulePick <- 2 # pick the upregulated healthy module
myModule_up_healthy <- diffGenes[names(module.assign[module.assign %in% modulePick]),]
length(rownames(myModule_up_healthy))
gost.res.healthy <- gost(rownames(myModule_up_healthy), organism = "slycopersicum", correction_method = "fdr") # run GO enrichment analysis
gostplot(gost.res.healthy, interactive = T, capped = F)
gostplot(gost.res.healthy, interactive = F, capped = F) # interactive Manhattan plot of enriched GO terms

# save gost plot and table

publish_gosttable(
  gost.res.healthy$result[order(gost.res.healthy$result$p_value, decreasing = F),],
  highlight_terms = NULL,
  use_colors = TRUE,
  show_columns = c("source", "term_name", "term_size", "intersection_size"),
  filename = "gosttable_healthy.pdf",
  ggplot=TRUE)

publish_gostplot(
  gostplot(gost.res.healthy, interactive = F, capped = F), # static gostplot
  highlight_terms = c("GO:0004866", "GO:0061135", "GO:0061134", "GO:0030414", "GO:0009611"), # highlight top 5
  filename = "gostplot_healthy.png",
  width = NA,
  height = NA)

modulePick <- 1 # pick the drought module/downregulated genes
myModule_down_drought <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
length(rownames(myModule_down_drought))

gost.res.drought <- gost(rownames(myModule_down_drought), organism = "slycopersicum", ordered_query = T, correction_method = "fdr") # run GO enrichment analysis
gostplot(gost.res.drought, interactive = T, capped = F) 
gostplot(gost.res.drought, interactive = F, capped = F) # interactive Manhattan plot of enriched GO terms

# save gost plot and table


publish_gosttable(
  gost.res.drought$result[order(gost.res.drought$result$p_value, decreasing = F),],
  highlight_terms = NULL,
  use_colors = TRUE,
  show_columns = c("source", "term_name", "term_size", "intersection_size"),
  filename = "gosttable_drought.pdf",
  ggplot= T)

publish_gostplot(gostplot(gost.res.drought, interactive = F, capped = F), # static gostplot
                 highlight_terms = c("GO:0016762", "GO:0005372", "GO:0015250", "GO:0031298", "GO:0031490"), # lowest 5 p-values
                 filename = "gostplot_drought.png",
                 width = NA,
                 height = NA)

# Heatmapping
rev(brewer.pal(name="RdBu", n=11)) # define a color palette
module.color <- rainbow(length(unique(module.assign)), start=0.3, end=0.2) 
module.color <- module.color[as.vector(module.assign)] # assign module colors

# generate a heatmap for all DEGs
clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete") #cluster rows by pearson correlation
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete")
module.assign <- cutree(clustRows, k=2)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
heatmap.2(diffGenes, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=myheatcolors, scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20))

#Heatmapping for healthy genes 
modulePick <- 2 # module for healthy samples
myModule_up_healthy <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
hrsub_up <- hclust(as.dist(1-cor(t(myModule_up_healthy), method="pearson")), method="complete") 

heatmap.2(myModule_up_healthy, 
          Rowv=as.dendrogram(hrsub_up), 
          Colv=NA, 
          labRow = NA,
          col=myheatcolors, scale="row", 
          density.info="none", trace="none", 
          RowSideColors=module.color[module.assign%in%modulePick], margins=c(8,20))

# Heatmapping downregulated genes (drought samples)

modulePick <- 1 # module for drought samples
myModule_down_drought <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
hrsub_down <- hclust(as.dist(1-cor(t(myModule_down_drought), method="pearson")), method="complete") 

heatmap.2(myModule_down_drought, 
          Rowv=as.dendrogram(hrsub_down), 
          Colv=NA, 
          labRow = NA,
          col=myheatcolors, scale="row", 
          density.info="none", trace="none", 
          RowSideColors=module.color[module.assign%in%modulePick], margins=c(8,20))
