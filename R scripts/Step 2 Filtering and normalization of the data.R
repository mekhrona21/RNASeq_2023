# Step 2 Filtering and normalization of the data

myTPM <- Txi_gene$abundance #abundances in transcripts per million
myTPM.stats <- transform(myTPM, 
                         SD=rowSds(myTPM), 
                         AVG=rowMeans(myTPM),
                         MED=rowMedians(myTPM)) # compute standard statistics for the abundance data
ggplot(myTPM.stats) +
  aes(x = SD, y = MED) +
  geom_smooth(method=lm) +
  geom_hex(show.legend = T, bins=25) +
  labs(y="Median", x = "Standard deviation",
       title="Transcripts per million (TPM)",
       subtitle="unfiltered, non-normalized data",
       caption="DIYtranscriptomics - Spring 2020") +
  theme_bw()

# visualization of counts
sampleLabels <- targets$sample
myDGEList <- DGEList(Txi_gene$counts) # creates a digital gene expression list object
log2.cpm <- cpm(myDGEList, log=TRUE) # use 'cpm' function to get counts per million
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID") # formatting
colnames(log2.cpm.df) <- c("geneID", sampleLabels)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, 
                                  cols = SRR7652567:SRR7652563, 
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # new variable storing the data
p1 <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "samples",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

# filtering
table(rowSums(myDGEList$counts==0)==6) # genes not expressed in all samples
cpm <- cpm(myDGEList) # convert to counts per million
keepers <- rowSums(cpm>1) >= 3 # removes genes with a total count below 1 CPM over any 3 samples
# 3 was chosen as the number of samples in the smallest group of comparison (3) 
myDGEList.filtered <- myDGEList[keepers,] 
log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE) # convert to cpm and log2 scale
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID") # format
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df,
                                           cols = SRR7652567:SRR7652563, 
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)

p2 <- ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "samples",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()


# normalization
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM") # calculate normalization factors
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE) # convert to cpm and log2 scale
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID") # format
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df,
                                                cols = SRR7652567:SRR7652563, 
                                                names_to = "samples", # name of that new variable
                                                values_to = "expression") # name of new variable storing the data

p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "samples",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

plot_grid(p1, p2, p3) # plot all 3 violin plots
