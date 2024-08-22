# 2024 fabio-affaticati
# noct_hypoxia.R


### Research questions:
### 1)  DEGs correlated with T0_OSA (categorical) 
### 2)  DEGs correlated with T0_oAHI (continuous)


library(stringr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(edgeR)
library(dplyr)
library(gridExtra)
library(ggstatsplot)
library(DESeq2)
library(ggpubr)
library(enrichplot)
library(org.Hs.eg.db)
library(clusterProfiler)
library(EnhancedVolcano)


source("src/utils/tools.R")


dataDir = paste0(getwd(), "/data/")
plotsDir = paste0(getwd(), "/plots/")
processedDir = paste0(dataDir, "processed_data/")

create_directories(c(plotsDir, processedDir))



### Load sequencing data (the last column is empty) and the gene reference file
run1 <- read.table(paste0(dataDir,'readcount/readcounts_210602.txt'), sep = '\t', header=T, row.names=1)[,-161]
run2 <- read.table(paste0(dataDir,'readcount/readcounts_210603.txt'), sep = '\t', header=T, row.names=1)[,-161]
gene <- read.delim(paste0(dataDir,'readcount/genes.txt'),header=T,sep="\t",row.names=1)


# Merge read count tables of the two runs
countdata <- merge(x=run1,y=run2,by=0,all=TRUE)
rownames(countdata) <- countdata[,1]
countdata <- countdata[,-1]
countdata[is.na(countdata)] <- 0

# Extract more compact sample ids
name <-  str_extract(colnames(countdata),"\\d+_T\\d")
name_run1 <- str_extract(colnames(run1),"\\d+_T\\d")
name_run2 <- str_extract(colnames(run2),"\\d+_T\\d")

# Sum transcript counts per sample over the lanes
tcountdata <- transpose(countdata)
tcountdata$name <- name
tcountdata_merge <- aggregate(. ~ name, FUN=sum, tcountdata)
samples <- tcountdata_merge[,1]
tcountdata_merge <- tcountdata_merge[,-1]
countdata_merge <- transpose(tcountdata_merge)
colnames(countdata_merge) <- samples
rownames(countdata_merge) <- rownames(countdata)


# Drop first 3 rows of countdata for unwanted non-gene rows
countdata_merge <- countdata_merge[4:nrow(countdata_merge),]

# Find and drop hemoglobin genes
hb_rows <- gene[grep('hemoglobin', gene$description, ignore.case = TRUE), ]
hb_rows <- hb_rows[grepl(paste0("^", 'HB'), hb_rows$name, ignore.case = TRUE), ]
common_rows <- intersect(rownames(countdata_merge), rownames(hb_rows))
countdata_merge <- countdata_merge[!(rownames(countdata_merge) %in% common_rows), ]


plot_library_sizes(countdata_merge, plotsDir)

# Returns filtered countdata 
countdata_merge_filtered <- density_plot_logcpm(countdata_merge, plotsDir)




############################################################ Save filtered readcounts ; can be used as a checkpoint in the analysis
write.csv(countdata_merge_filtered, paste0(processedDir, "merged_reads.csv"), row.names = TRUE)
countdata_merge_filtered <- read.csv(paste0(processedDir, "merged_reads.csv"), row.names=1)




colnames(countdata_merge_filtered) <- sub("^X", "", colnames(countdata_merge_filtered))
gene<-read.delim(paste0(dataDir, 'readcount/genes.txt') ,header=T,sep="\t",row.names=1)


description <- colnames(countdata_merge_filtered)
donor <- str_extract(colnames(countdata_merge_filtered),"^\\d+")
time <- str_extract(description,"T\\d+$")


# Load metadata
sampledata <- read.table(file= paste0(dataDir, 'samples.csv'), header=T,sep="\t",row.names=1)
meta<-read.table(file= paste0(dataDir, 'meta.txt'), header=T,sep="\t")
meta$ID <-as.factor(meta$ID)
meta$Sex <-as.factor(meta$Sex)

# Create df for differential expression analysis
coldata<-merge(x=sampledata,by.x=1,y=meta,by.y="ID",sort=FALSE)
rownames(coldata) <- rownames(sampledata)
coldata$BMI_group <- as.factor(coldata$BMI_group)
coldata$T0_Puberty <- as.factor(coldata$T0_Puberty)
coldata$T0_OSA <- as.factor(coldata$T0_OSA)
coldata$T0_BMI_z <- (coldata$T0_BMI - mean(coldata$T0_BMI)) / sd(coldata$T0_BMI)
coldata$T0_Age_z <- (coldata$T0_Age - mean(coldata$T0_Age)) / sd(coldata$T0_Age)


### keep T0 data
dataset1 <- countdata_merge_filtered[,coldata$time == 'T0']
coldata1 <- coldata[coldata$time == 'T0',]


coldata1$T0_ODI # The ODI is the number of times per hour of sleep that your blood oxygen level drops
                # by a certain degree from baseline.


outlier_boxplot_detect(dataset1, paste0(plotsDir, "gene_exp_beforeoutliers.png"), "LogCPM expression before outliers removal")

############ Drop low counts samples and outliers, 3035, 3041
drops <- !coldata1$donor==3035 & !coldata1$donor==3041
dataset1 <- dataset1[,!coldata1$donor==3035 & !coldata1$donor==3041]
coldata1 <- coldata1[drops,]

outlier_boxplot_detect(dataset1, paste0(plotsDir, "gene_exp_aftereoutliers.png"), "LogCPM expression after outliers removal")


pca_df = plot_pca_batcheffect(dataset1, coldata1, plotsDir)
coldata1$batch <- pca_df$batch


plot_library_sizes_test(dataset1, pca_df, plotsDir)


correlation_tests(coldata1, c('T0_BMI', 'T0_Age', 'T0_Waist', 'T0_Glucose', 'T0_Insulin'), plotsDir)


##################
#     DESeq2     # Genes related to OSA
##################

coldata1$T0_BMI_scaled <- (coldata1$T0_BMI -  mean(coldata1$T0_BMI, na.rm = TRUE)) /  sd(coldata1$T0_BMI, na.rm = TRUE)

dds <- DESeqDataSetFromMatrix(countData = dataset1, 
                              colData = coldata1, design = ~ batch + Sex + T0_Age_z + T0_BMI_scaled + T0_OSA)

dds <- DESeq(dds, minReplicatesForReplace=Inf, parallel = TRUE, fitType = "local")

res <- results(dds, name = "T0_OSA_1_vs_0",  cooksCutoff=FALSE, independentFiltering=FALSE)

res$padj[is.na(res$padj)] <- 1
siggenes <- merge(x = data.frame(res[res$padj < 0.05,]),by.x=0,y=gene,by.y=0,all=F)
rownames(siggenes) <- siggenes$Row.names
write.csv(siggenes[order(siggenes$log2FoldChange),], file= paste0(processedDir, 'OSAsiggenes.csv'),quote=FALSE,row.names = FALSE)

### Define background genes as all the ones in the reference file
background <- merge(x = data.frame(res),by.x=0,y=gene,by.y=0,all=F)

go_enrichment(siggenes, background, "OSA_dotplot", plotsDir)

m <- log2(counts(dds,normalized=TRUE) + 1)
#### Negative correlated genes
neg_gene_correlated <- siggenes[siggenes$log2FoldChange < 0,]
neg_gene_correlated <- neg_gene_correlated[order(neg_gene_correlated$log2FoldChange, decreasing = FALSE), ]
transposed_genes <- m[rownames(neg_gene_correlated),]

##### take top 4 most differentially expressed genes
rownames(transposed_genes) <- neg_gene_correlated$name
transposed_genes <- as.data.frame(t(transposed_genes[1:4,]))

transposed_genes$T0_OSA <- coldata1$T0_OSA
transposed_genes$Sex <- coldata1$Sex
transposed_genes$batch <- coldata1$batch
transposed_genes$T0_oAHI <- coldata1$T0_oAHI
melted_df <- reshape2::melt(transposed_genes,  id.vars = c('T0_OSA', 'T0_oAHI', 'Sex', 'batch'), variable.name = 'series')

melted_df$OSA <- ifelse(as.numeric(as.character(melted_df$T0_OSA)), 'OSA', 'Not OSA')

p <- ggplot(melted_df, aes(x = OSA, y = value)) +
  geom_jitter(width = 0.2) +
  geom_boxplot(outlier.shape = NA, varwidth = TRUE, alpha = 0.2, width = 0.15) +
  facet_wrap(~ series, scales = 'free_y') +
  labs(title = "Negative Log Fold change genes with OSA",
       x = "OSA",
       y = "Log2 Normalised Counts") + 
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), 
                     label.y = max(melted_df$value, na.rm = TRUE)-.2,
                     label.x.npc = 'center',
                     size = 5) +
  stat_compare_means(method = "wilcox.test", label = "p.format", 
                     label.y = max(melted_df$value, na.rm = TRUE) + .5,
                     label.x.npc = 0.42,
                     size = 5)

ggsave(paste0(plotsDir, "neg_OSA_corr_genes.png"), p, scale = 1, width = 8, height = 10)

### Volcano plot
matching_rows <- intersect(rownames(res), rownames(gene))
volcano_data <- res[match(matching_rows, rownames(res)), ]
rownames(volcano_data) <- gene[matching_rows, "name"]

p <- EnhancedVolcano(volcano_data,
                lab = rownames(volcano_data),
                title = "",
                pCutoffCol = "padj",
                pCutoff = 0.05,
                FCcutoff = 0.5,
                legendPosition = "right",
                subtitle = "OSA DEG testing",
                #shape = c(6, 6, 19),
                colAlpha = 0.5,
                x = 'log2FoldChange',
                y = 'pvalue') +
        geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
        annotate("text", x = 2, y = 9, label = "OSA up", angle = 0, vjust = 1.5, color = "black", size = 5, fontface = 2) +
        annotate("segment", x = 1.5, xend = 3, y = 8, yend = 8,
        arrow = arrow(type = "closed", length = unit(0.3, "cm")), 
        color = "black", size = 1) +
        annotate("text", x = -2, y = 9, label = "OSA down", angle = 0, vjust = 1.5, color = "black", size = 5, fontface = 2) +
        annotate("segment", x = -1.5, xend = -3., y = 8, yend = 8,
           arrow = arrow(type = "closed", length = unit(0.3, "cm")), 
           color = "black", size = 1)

ggsave(paste0(plotsDir, "volcanoplot_OSA.png"), p, scale = 1, width = 12, height = 8)
#####################################################################





##################
#     DESeq2     # T0 genes related to oAHI
##################

# Center and scale the variable
coldata1$T0_BMI_scaled <- (coldata1$T0_BMI -  mean(coldata1$T0_BMI, na.rm = TRUE)) /  sd(coldata1$T0_BMI, na.rm = TRUE)
coldata1$T0_oAHI_scaled <- (coldata1$T0_oAHI - mean(coldata1$T0_oAHI, na.rm = TRUE)) /  sd(coldata1$T0_oAHI, na.rm = TRUE)

dds <- DESeqDataSetFromMatrix(countData = dataset1, 
                              colData = coldata1, design = ~ batch + Sex + T0_Age_z + T0_BMI_scaled + T0_oAHI_scaled)


dds <- DESeq(dds, minReplicatesForReplace=Inf, parallel = TRUE, fitType = "local")
res <- results(dds, cooksCutoff=FALSE, independentFiltering=FALSE)

res$padj[is.na(res$padj)] <- 1
siggenes <- merge(x = data.frame(res[res$padj < 0.05,]),by.x=0,y=gene,by.y=0,all=F)
rownames(siggenes) <- siggenes$Row.names

tosave<- siggenes

tosave <- tosave %>%
  filter(abs(log2FoldChange) > 0.5) %>%
  dplyr::select(log2FoldChange, padj, name, description, genebank)

write.csv(tosave[order(tosave$log2FoldChange),], file=paste0(processedDir ,"oAHIsiggenes.csv"),quote=FALSE,row.names = FALSE)


matching_rows <- intersect(rownames(res), rownames(gene))
volcano_data <- res[match(matching_rows, rownames(res)), ]
rownames(volcano_data) <- gene[matching_rows, "name"]


p <- EnhancedVolcano(volcano_data,
                     lab = rownames(volcano_data),
                     title = "",
                     pCutoffCol = "padj",
                     pCutoff = 0.05,
                     FCcutoff = 0.5,
                     legendPosition = "right",
                     subtitle = "oAHI DEG testing",
                     #shape = c(6, 6, 19, 16),
                     colAlpha = 0.5,
                     xlim = c(-3, 1.5),
                     ylim = c(0, 12),
                     x = 'log2FoldChange',
                     y = 'pvalue')+
        geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.5) +
        annotate("text", x = 1, y = 10, label = "AHI up", angle = 0, vjust = 1.5, color = "black", size = 5, fontface = 2) +
        annotate("segment", x = .7, xend = 1.5, y = 9, yend = 9,
                 arrow = arrow(type = "closed", length = unit(0.3, "cm")), 
                 color = "black", size = 1) +
        annotate("text", x = -1, y = 10, label = "AHI down", angle = 0, vjust = 1.5, color = "black", size = 5, fontface = 2) +
        annotate("segment", x = -.7, xend = -1.5, y = 9, yend = 9,
                 arrow = arrow(type = "closed", length = unit(0.3, "cm")), 
                 color = "black", size = 1)

ggsave(paste0(plotsDir, "volcanoplot_continuous.png"),p, scale = 1, width = 12, height = 8)


go_enrichment(siggenes, background, "OSA_cont_dotplot", plotsDir)



m <- log2(counts(dds,normalized=TRUE) + 1)
#### Negative correlated genes
neg_gene_correlated <- siggenes[siggenes$log2FoldChange < 0,]
neg_gene_correlated <- neg_gene_correlated[order(neg_gene_correlated$log2FoldChange, decreasing = FALSE), ]
transposed_genes <- m[rownames(neg_gene_correlated),]

##### Only DUSP21
transposed_genes <- data.frame(DUSP21 = transposed_genes[1,])

transposed_genes$T0_oAHI <- coldata1$T0_oAHI
transposed_genes$T0_oAHI_scaled <- coldata1$T0_oAHI_scaled
transposed_genes$Sex <- coldata1$Sex
transposed_genes$batch <- coldata1$batch
transposed_genes$T0_OSA <- coldata1$T0_OSA
transposed_genes$Date.library.prep <- coldata1$Date.library.prep
melted_df <- reshape2::melt(transposed_genes,  id.vars = c('T0_oAHI', 'T0_OSA', 'T0_oAHI_scaled', 'Sex', 'batch', 'Date.library.prep'), variable.name = 'series')

# Perform Spearman correlation test
cor_test <- cor.test(melted_df$T0_oAHI, melted_df$value, method = "spearman", exact = FALSE)

# Extract correlation coefficient and p-value
cor_coef <- cor_test$estimate
p_value <- cor_test$p.value

# Format the results for annotation
cor_text <- paste0("Spearman rho: ", round(cor_coef, 2))
p_value_text <- paste0("p-value: ", format(p_value, scientific = FALSE, digits = 4))

melted_df$OSA <- ifelse(as.numeric(as.character(melted_df$T0_OSA)), 'OSA', 'Not OSA')
p <- ggplot(melted_df, aes(x = T0_oAHI, y = value, shape = OSA)) +
  geom_point(size = 3) +
  #geom_boxplot() +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ series, scales = 'free_y') +
  #geom_smooth(method = 'glm', aes(color = NULL)) +
  labs(title = "DUSP21 correlation to AHI",
       x = "oAHI",
       y = "Log2 Normalised Counts") + 
  annotate("text", x = 8, y = 9, label = cor_text,size = 4, color = "black") +
  annotate("text", x = 8, y = 8.6, label = p_value_text, size = 4, color = "black")+
  geom_vline(xintercept = 2, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = 2.1, y = 5, label = "Obstructive Sleep Apnea threshold", angle = 90, vjust = 1.5, color = "red", size = 3)

ggsave(paste0(plotsDir, "neg_oAHI_corr_genes.png"), p, scale = 1, width = 6, height = 6)
