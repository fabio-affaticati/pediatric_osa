### Research questions:
### 1)  DEGs correlated with T0_OSA (categorical) 
### 2)  DEGs correlated with T0_oAHI (continuous)


# Install codetools if not already installed
install.packages("codetools")
library(codetools)

# Analyze a function or the entire script
used_globals <- findGlobals(your_function, merge = FALSE)$functions

# Print used globals
print(used_globals)











require(edgeR)
require(ggplot2)
library(reshape2)
library(gridExtra)
require(Biobase)
require(tidyverse)
library(magrittr)
library(tidyr)
library(dplyr)
require(ggrepel)
library(pheatmap)
library(DESeq2)
library(RColorBrewer)
require(data.table)
library(ggpmisc)
#BiocManager::install(c("clusterProfiler", "DOSE", "enrichplot"))
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(ggpubr)
library(ggstatsplot)


workingDir = "/Users/fabioaffaticati/Desktop/Work/pediatric_diet/"
setwd(workingDir) 



outlier_boxplot_detect <- function(data, plot_filename, plot_title) {
  
  y.all <- DGEList(counts=data, genes = rownames(data))
  y.all <- calcNormFactors(y.all)
  
  lcpm <- cpm(y.all, log=TRUE, normalized.lib.sizes = TRUE)
  lcpm_melted <-reshape2::melt(lcpm)
  
  lcpm_melted <- lcpm_melted %>%
    rename(
      Sample_name = Var2,
      LogCPM = value,
    )
  
  lcpm_melted_names <- as.data.frame(do.call(rbind, strsplit(as.character(lcpm_melted$Sample_name), "_")))
  lcpm_melted$Sample <- paste(lcpm_melted_names$V1)
  
  p <- ggplot(lcpm_melted, aes(x=Sample, y=LogCPM, fill=Sample))+
    geom_boxplot()+
    theme_bw(base_size = 15) +
    theme(legend.position = "none", text = element_text(family = "JetMono Brains"),
          plot.title = element_text(family = "JetMono Brains")) +
    labs(title = plot_title) +
    theme(
      axis.text.x = element_text(colour="black",size=12,angle=90,hjust=0.6,vjust=.7,face="plain",family = "JetMono Brains"),
    )
  ggsave(plot_filename, p, scale = 1, width = 8, height = 6)
}






density_plotlogcpm <- function(data, plot_filename, plot_title){
  
  y.all <- DGEList(counts=data, genes = rownames(data))
  y.all <- calcNormFactors(y.all)
  y.all <- cpm(y.all, log=TRUE, normalized.lib.sizes = TRUE)
  
  logCPM_before <-reshape2::melt(y.all)
  
  
  
  ########## Before filtering
  logCPM_before <- logCPM_before %>%
    rename(
      Sample_name = Var2,
      LogCPM = value,
    )
  
  lcpm_melted_names <- as.data.frame(do.call(rbind, strsplit(as.character(logCPM_before$Sample_name), "_")))
  logCPM_before$Sample <- paste(lcpm_melted_names$V1)
  
  
  labels <- logCPM_before %>%
    group_by(Sample) %>%
    summarize(x = density(LogCPM)$x[which.max(density(LogCPM)$y)],
              y = max(density(LogCPM)$y)) %>%
    arrange(desc(y))
  
  
  
  plot_before <- ggplot(logCPM_before, aes(x = LogCPM, color=Sample)) +
    geom_density(alpha = 0.5) +
    geom_text(data = labels[1:3, ], aes(x = x, y = y, label = Sample),
              size = 5, vjust = -0.5) + 
    labs(title = "Density Plot Before Filtering",
         x = "logCPM",
         y = "Density") +
    theme_minimal(base_size = 15) +
    theme(legend.position = "none",
          axis.text.x = element_text(colour="black",size=12,vjust=1)) +
    scale_x_continuous(limits = c(-5, 12))

  # Drop genes with low reads
  y.all <- as.data.frame(y.all)
  y.all <- y.all[apply(y.all,1,sum) >1,] # choose LogCPM threshold
  
  # Find common row names
  common_genes <- intersect(rownames(data), rownames(y.all))
  # Subset df1 to keep only the rows with common names
  data <- data[common_genes, ]
  
  
  
  
  ########## After filtering
  logCPM_after <-reshape2::melt(y.all)
  
  logCPM_after <- logCPM_after %>%
    rename(
      Sample_name = variable,
      LogCPM = value,
    )
  
  lcpm_melted_names <- as.data.frame(do.call(rbind, strsplit(as.character(logCPM_after$Sample_name), "_")))
  logCPM_after$Sample <- paste(lcpm_melted_names$V1)
  
  labels <- logCPM_after %>%
    group_by(Sample) %>%
    summarize(x = density(LogCPM)$x[which.max(density(LogCPM)$y)],
              y = max(density(LogCPM)$y))  %>%
    arrange(desc(y))
  
  # Create density plot for logCPM after filtering
  plot_after <- ggplot(logCPM_after, aes(x = LogCPM, color=Sample)) +
    geom_density(alpha = 0.5) +
    geom_text(data = labels[1:3, ], aes(x = x, y = y, label = Sample),
              size = 5, vjust = -0.5) + 
    labs(title = "Density Plot After Filtering",
         x = "logCPM",
         y = "Density") +
    theme_minimal(base_size = 15) +
    theme(legend.position = "none",
          axis.text.x = element_text(colour="black",size=12,vjust=1)) +
    scale_x_continuous(limits = c(-5, 12))
  
  
  # Arrange the two plots side by side
  p <- grid.arrange(plot_before, plot_after, ncol = 2)
  ggsave(plot_filename, p, scale = 1, width = 16, height = 7)
  
  return(data)
  
}



# Last column is empty
run1<-read.table('pediatric_bmi/readcount/readcounts_210602.txt', sep = '\t', header=T, row.names=1)[,-161]
run2<-read.table('pediatric_bmi/readcount/readcounts_210603.txt', sep = '\t', header=T, row.names=1)[,-161]
gene<-read.delim('pediatric_bmi/readcount/genes.txt',header=T,sep="\t",row.names=1)


# Merge read count tables
countdata <- merge(x=run1,y=run2,by=0,all=TRUE)
rownames(countdata) <- countdata[,1]
countdata <- countdata[,-1]
countdata[is.na(countdata)] <- 0

name <-  str_extract(colnames(countdata),"\\d+_T\\d")
name_run1 <- str_extract(colnames(run1),"\\d+_T\\d")
name_run2 <- str_extract(colnames(run2),"\\d+_T\\d")

tcountdata <- transpose(countdata)
tcountdata$name <- name
tcountdata_merge <- aggregate(. ~ name, FUN=sum,tcountdata)
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


# Plot library sizes
column_sums <- sapply(countdata_merge, sum)
# Create a new dataframe with column names and their sums
sum_df <- data.frame(
  Samples = names(column_sums),
  Library.size = column_sums
)

sum_df$Timepoint <- as.data.frame(do.call(rbind, strsplit(as.character(sum_df$Sample), "_")))$V2
sum_df <- sum_df[sum_df$Timepoint == 'T0',]
sum_df$croppped <- sub("_.*", "", sum_df$Samples)


p <- ggplot(sum_df, aes(x = croppped, y = Library.size)) + #, fill = Timepoint)) +
  geom_bar(stat = "identity") +
  labs(title = "Library sizes of the cohort",
       x = "Sample ID",
       y = "Library size") +
  theme_bw(base_size = 15) +
  scale_y_continuous(labels = scales::comma) +
  theme(plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.text.x = element_text(colour="black",size=12,angle=45,hjust=.5,vjust=.7,face="plain",family = "serif"),
  )

ggsave('results/librarysizes.png', p, scale = 2, width = 14, height = 8)
p




countdata_merge <- density_plotlogcpm(countdata_merge, "results/densityplots.png", 'title')



# Save readcounts
write.csv(countdata_merge, "merged_reads.csv", row.names = TRUE)
countdata_merge <- read.csv("merged_reads.csv", row.names=1)




colnames(countdata_merge) <- sub("^X", "", colnames(countdata_merge))
gene<-read.delim('pediatric_bmi/readcount/genes.txt',header=T,sep="\t",row.names=1)


description <- colnames(countdata_merge)
donor <- str_extract(colnames(countdata_merge),"^\\d+")
time <- str_extract(description,"T\\d+$")


# Load metadata
sampledata <- read.table(file='pediatric_bmi/analysis/samples.csv',header=T,sep="\t",row.names=1)

meta<-read.table(file="pediatric_bmi/meta.txt",header=T,sep="\t")
meta$ID <-as.factor(meta$ID)
meta$Sex <-as.factor(meta$Sex)

coldata<-merge(x=sampledata,by.x=1,y=meta,by.y="ID",sort=FALSE)
rownames(coldata) <- rownames(sampledata)


coldata$BMI_group <- as.factor(coldata$BMI_group)
coldata$T0_Puberty <- as.factor(coldata$T0_Puberty)
coldata$T0_OSA <- as.factor(coldata$T0_OSA)
coldata$T0_BMI_z <- (coldata$T0_BMI - mean(coldata$T0_BMI)) / sd(coldata$T0_BMI)
coldata$T0_Age_z <- (coldata$T0_Age - mean(coldata$T0_Age)) / sd(coldata$T0_Age)



### keep T0 data
dataset1 <- countdata_merge[,coldata$time == 'T0']
coldata1 <- coldata[coldata$time == 'T0',]


coldata1$T0_ODI # The ODI is the number of times per hour of sleep that your blood oxygen level drops
                # by a certain degree from baseline.


  





outlier_boxplot_detect(dataset1, "results/gene_exp_beforeoutliers.png", "LogCPM expression before outliers removal")

############ Drop low counts samples and outliers, 3035, 3029, 3041, 3042
#drops <- !coldata$donor==3035 & !coldata$donor==3029  & !coldata$donor==3041 & !coldata$donor==3042
#countdata_merge <- countdata_merge[,!coldata$donor==3035 & !coldata$donor==3029  & !coldata$donor==3041 & !coldata$donor==3042]
#coldata <- coldata[drops,]

drops <- !coldata1$donor==3035 & !coldata1$donor==3041
dataset1 <- dataset1[,!coldata1$donor==3035 & !coldata1$donor==3041]
coldata1 <- coldata1[drops,]

outlier_boxplot_detect(dataset1, "results/gene_exp_aftereoutliers.png", "LogCPM expression after outliers removal")





y.all <- DGEList(counts=dataset1, genes = rownames(dataset1))
y.all <- calcNormFactors(y.all)
lcpm <- cpm(y.all, log=TRUE, normalized.lib.sizes = TRUE)



# Run PCA
pca <- prcomp(t(lcpm), center = TRUE, scale. = TRUE)
# Combine PCA coordinates with the metadata from the DGEList
to_plot <- data.frame(pca$x)
head(to_plot)
# Calculate how many % of total variance is explained by each principal component
percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
# We focus here on PC1 and PC2
use.pcs <- c(1,2)


to_plot$date <- coldata1$Date.library.prep
to_plot$time <- coldata1$time
to_plot$batch <- as.factor(ifelse(to_plot$PC1 > 20, "batch1", "batch2"))
to_plot$donor <- coldata1$donor
coldata1$batch <- to_plot$batch



labs <- paste0(paste0("PC", use.pcs, " - "), paste0("Var.expl = ", round(percentVar[use.pcs], 2), "%"))

p <- ggplot(to_plot, aes(x = PC1, y = PC2, color = date, label = donor)) +
  geom_point(size = 3, alpha = 0.7)  +
  geom_text_repel(aes(color = date), nudge_x = 0.2, size = 3, segment.color = NA) +
  geom_vline(xintercept = 20, linetype = "dashed", color = "red", linewidth = 1) +
  #scale_color_brewer(palette = "Dark2") +
  theme_bw(base_size = 15) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5)) +
  labs(
    title = "PCA projection of LogCPMs",
    color = "Date",
    x = labs[1],
    y =labs[2]
  ) +  annotate("text", x = 20, y = -40, label = "Batch effect", angle = 90, vjust = 1.5, color = "red", size = 5)
ggsave('results/pca.png', p, scale = 1, width = 10, height = 8)
p



### Plot library sizes test
sum_df <- sum_df[rownames(sum_df) %in% rownames(to_plot), ]
sum_df$batch <- to_plot$batch
p <- ggbetweenstats(data = sum_df,
                    x = batch,
                    y = Library.size,
                    ylab = 'Library Size',
                    xlab = 'Batch',
                    title = 'Statistical comparison of library sizes by Batch',
                    #pairwise.display = 'all',
                    pairwise.display = 's',
                    p.adjust.method = 'BH',
                    digits = 4L,
                    ggtheme = theme_bw(base_size = 12) +
                      theme(
                        text = element_text(size = 12),  # Adjust text size
                        axis.text = element_text(size = 15),  # Adjust axis text size
                        axis.title = element_text(size = 1)  # Adjust axis title size
                      ),
                    boxplot.args = list(width = 0.4, alpha = 0.2, color = "gray"),
                    points.args = list(size = 50, alpha = 1),
                    violin.args = list(width = 0, alpha = 0, color = "lightgray"),
                    plot.type = "boxplot", type = "non-parametric",
                    geom_signif_args = list(textsize = 30),
                    #geom_signif(list(map_signif_level = TRUE)),
                    centrality.point.args = list(size = 3,color = "#8a0f00"),
                    centrality.label.args = list(size = 3)) +
  scale_y_continuous(labels = scales::comma)
ggsave("results/librarysizes_test.png", p, scale = 1, width = 8, height = 8)
p







cor_test <- cor.test(coldata1$T0_oAHI, coldata1$T0_BMI, method = "spearman", exact = FALSE)
cor_coef <- cor_test$estimate
p_value <- cor_test$p.value
cor_text <- paste0("Spearman rho: ", round(cor_coef, 2))
p_value_text <- paste0("p-value: ", format(p_value, scientific = FALSE, digits = 4))


p <- ggplot(coldata1, aes(x=T0_BMI, y=T0_oAHI, shape=T0_OSA, size = 2))+
  geom_point() +
  theme_bw(base_size = 15) +
  labs(x = "BMI",
       y = "AHI") +
  theme(legend.position = "none",
        axis.text.x = element_text(colour="black",size=12,angle=45,hjust=0.6,vjust=.7,face="plain",),
  ) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = 27, y = 12, label = cor_text,size = 6, color = "black") +
  annotate("text", x = 27, y = 11, label = p_value_text, size = 6, color = "black")+
  annotate("text", x = 45, y = 2.5, label = "Obstructive Sleep Apnea threshold", angle = 0, color = "red", size = 4)
p
ggsave('results/BMI_AHI.png', p, scale = 1, width = 8, height = 6)


cor_test <- cor.test(coldata1$T0_oAHI, coldata1$T0_Age, method = "spearman", exact = FALSE)
cor_coef <- cor_test$estimate
p_value <- cor_test$p.value
cor_text <- paste0("Spearman rho: ", round(cor_coef, 2))
p_value_text <- paste0("p-value: ", format(p_value, scientific = FALSE, digits = 4))


p <- ggplot(coldata1, aes(x=T0_Age, y=T0_oAHI, shape=T0_OSA, size = 2))+
  geom_point() +
  theme_bw(base_size = 15) +
  labs(x = "Age",
       y = "AHI") +
  theme(legend.position = "none",
        axis.text.x = element_text(colour="black",size=12,angle=45,hjust=0.6,vjust=.7,face="plain",),
  ) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = 10, y = 10, label = cor_text,size = 6, color = "black") +
  annotate("text", x = 10, y = 9, label = p_value_text, size = 6, color = "black")+
  annotate("text", x = 16.3, y = 2.5, label = "Obstructive Sleep Apnea threshold", angle = 0, color = "red", size = 4)
p
ggsave('results/Age_AHI.png', p, scale = 1, width = 10, height = 6)


gender_data <- coldata1
gender_data$Sex <- ifelse(gender_data$Sex == 1, 'Male', 'Female')
p <- ggplot(gender_data, aes(x = Sex, y = T0_oAHI)) +
  geom_jitter(width = 0.2, aes(shape = T0_OSA)) +
  geom_boxplot(outlier.shape = NA, varwidth = TRUE, alpha = 0.2, width = 0.15) +
  #facet_wrap(~ series, scales = 'free_y') +
  #geom_smooth(method = 'glm', aes(color = NULL)) +
  labs(title = "",
       x = "Sex",
       y = "AHI") + 
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), 
                     label.y = max(coldata1$T0_oAHI, na.rm = TRUE) -.8, # Adjust label position if needed
                     label.x.npc = 'center', # Center the p-value label
                     size = 6) + # Adjust size if needed
  stat_compare_means(method = "wilcox.test", label = "p.format", 
                     label.y = max(coldata1$T0_oAHI, na.rm = TRUE), # Adjust label position if needed
                     label.x.npc = 0.42, # Center the p-value label
                     size = 5) + # Adjust size if needed
  geom_hline(yintercept = 2, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = 1.5, y = 2.5, label = "Obstructive Sleep Apnea threshold", angle = 0, color = "red", size = 4)

p
ggsave(filename = "results/Sex_AHI.png",p, scale = 1, width = 5, height = 5)



cor_test <- cor.test(coldata1$T0_oAHI, coldata1$T0_Waist, method = "spearman", exact = FALSE)
cor_coef <- cor_test$estimate
p_value <- cor_test$p.value
cor_text <- paste0("Spearman rho: ", round(cor_coef, 2))
p_value_text <- paste0("p-value: ", format(p_value, scientific = FALSE, digits = 4))


p <- ggplot(coldata1, aes(x=T0_Waist, y=T0_oAHI, shape=T0_OSA, size = 2))+
  geom_point() +
  theme_bw(base_size = 15) +
  labs(x = "Waist circumference",
       y = "AHI") +
  theme(legend.position = "none",
        axis.text.x = element_text(colour="black",size=12,angle=45,hjust=0.6,vjust=.7,face="plain",),
  ) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = 90, y = 10, label = cor_text,size = 6, color = "black") +
  annotate("text", x = 90, y = 9, label = p_value_text, size = 6, color = "black")+
  annotate("text", x = 105, y = 2.5, label = "Obstructive Sleep Apnea threshold", angle = 0, color = "red", size = 4)
p
ggsave('results/Waist_AHI.png', p, scale = 1, width = 10, height = 6)


cor_test <- cor.test(coldata1$T0_oAHI, coldata1$T0_Glucose, method = "spearman", exact = FALSE)
cor_coef <- cor_test$estimate
p_value <- cor_test$p.value
cor_text <- paste0("Spearman rho: ", round(cor_coef, 2))
p_value_text <- paste0("p-value: ", format(p_value, scientific = FALSE, digits = 4))


p <- ggplot(coldata1, aes(x=T0_Glucose, y=T0_oAHI, shape=T0_OSA, size = 2))+
  geom_point() +
  theme_bw(base_size = 15) +
  labs(x = "Glucose level",
       y = "AHI") +
  theme(legend.position = "none",
        axis.text.x = element_text(colour="black",size=12,angle=45,hjust=0.6,vjust=.7,face="plain",),
  ) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = 90, y = 10, label = cor_text,size = 6, color = "black") +
  annotate("text", x = 90, y = 9, label = p_value_text, size = 6, color = "black")+
  annotate("text", x = 100, y = 2.5, label = "Obstructive Sleep Apnea threshold", angle = 0, color = "red", size = 4)
p
ggsave('results/Glucose_AHI.png', p, scale = 1, width = 10, height = 6)



cor_test <- cor.test(coldata1$T0_oAHI, coldata1$T0_Insulin, method = "spearman", exact = FALSE)
cor_coef <- cor_test$estimate
p_value <- cor_test$p.value
cor_text <- paste0("Spearman rho: ", round(cor_coef, 2))
p_value_text <- paste0("p-value: ", format(p_value, scientific = FALSE, digits = 4))


p <- ggplot(coldata1, aes(x=T0_Insulin, y=T0_oAHI, shape=T0_OSA, size = 2))+
  geom_point() +
  theme_bw(base_size = 15) +
  labs(x = "Insulin level",
       y = "AHI") +
  theme(legend.position = "none",
        axis.text.x = element_text(colour="black",size=12,angle=45,hjust=0.6,vjust=.7,face="plain",),
  ) +
  geom_hline(yintercept = 2, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = 90, y = 10, label = cor_text,size = 6, color = "black") +
  annotate("text", x = 90, y = 9, label = p_value_text, size = 6, color = "black")+
  annotate("text", x = 305, y = 2.5, label = "Obstructive Sleep Apnea threshold", angle = 0, color = "red", size = 4)
p
ggsave('results/Insulin_AHI.png', p, scale = 1, width = 10, height = 6)







##################
#     DESeq2     # T0 genes related to T0_OSA
##################

coldata1$T0_BMI_scaled <- (coldata1$T0_BMI -  mean(coldata1$T0_BMI, na.rm = TRUE)) /  sd(coldata1$T0_BMI, na.rm = TRUE)


dds <- DESeqDataSetFromMatrix(countData = dataset1, 
                              colData = coldata1, design = ~ batch + Sex + T0_Age_z + T0_BMI_scaled + T0_OSA)
#dds <- DESeq(dds, parallel = TRUE, fitType = "local")
#res<-results(dds, name = "T0_OSA_1_vs_0")

dds <- DESeq(dds, minReplicatesForReplace=Inf, parallel = TRUE, fitType = "local")

res <- results(dds, name = "T0_OSA_1_vs_0",  cooksCutoff=FALSE, independentFiltering=FALSE)

res$padj[is.na(res$padj)] <- 1
siggenes <- merge(x = data.frame(res[res$padj < 0.05,]),by.x=0,y=gene,by.y=0,all=F)
rownames(siggenes) <- siggenes$Row.names
siggenes
write.csv(siggenes[order(siggenes$log2FoldChange),],file="results/T0_OSAsiggenes.csv",quote=FALSE,row.names = FALSE)


background <- merge(x = data.frame(res),by.x=0,y=gene,by.y=0,all=F)

#######################
#    GO enrichment    #
#######################

pos_ego <- enrichGO(
  gene = siggenes[siggenes$log2FoldChange > 0,]$name,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  universe = background$name,
  ont = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE,
)


neg_ego <- enrichGO(
  gene = siggenes[siggenes$log2FoldChange < 0,]$name,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  universe = background$name,
  ont = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE,
)

p <- enrichplot::dotplot(pos_ego, showCategory=10, title="GO positive LogFoldChange (increase with T0_OSA == 1)") + 
  facet_grid(ONTOLOGY ~ ., scales="free") + theme_bw(base_size = 12)
ggsave(filename = "results/pos_T0_OSA_1vs0_dotplot.png",p, scale = 1, width = 8, height = 6)
p

p <- enrichplot::dotplot(neg_ego, showCategory=10, title="GO negative LogFoldChange (decrease with T0_OSA == 1)") +
  facet_grid(ONTOLOGY ~ ., scales="free") + theme_bw(base_size = 12)
ggsave(filename = "results/neg_T0_OSA_1vs0_dotplot.png",p, scale = 1, width = 8, height = 8)
p



m <- log2(counts(dds,normalized=TRUE) + 1)
#### Negative correlated genes
neg_gene_correlated <- siggenes[siggenes$log2FoldChange < 0,]
neg_gene_correlated <- neg_gene_correlated[order(neg_gene_correlated$log2FoldChange, decreasing = FALSE), ]
transposed_genes <- m[rownames(neg_gene_correlated),]

##### Normal execution
rownames(transposed_genes) <- neg_gene_correlated$name
transposed_genes <- as.data.frame(t(transposed_genes[1:4,]))

##### If only DUSP21
transposed_genes <- data.frame(DUSP21 = transposed_genes)



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
  #geom_smooth(method = 'glm', aes(color = NULL)) +
  labs(title = "Negative Log Fold change genes with OSA",
       x = "OSA",
       y = "Log2 Normalised Counts") + 
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), 
                     label.y = max(melted_df$value, na.rm = TRUE)-.2, # Adjust label position if needed
                     label.x.npc = 'center', # Center the p-value label
                     size = 5) + # Adjust size if needed
  stat_compare_means(method = "wilcox.test", label = "p.format", 
                     label.y = max(melted_df$value, na.rm = TRUE) + .5, # Adjust label position if needed
                     label.x.npc = 0.42, # Center the p-value label
                     size = 5) # Adjust size if needed
ggsave(filename = "results/neg_T0_OSA_corr_genes.png",p, scale = 1, width = 8, height = 10)
p

p <- ggstatsplot::ggbetweenstats(data = melted_df,
                            x = T0_OSA, y=value,
                            ylab = 'Log2CPMs',
                            title = 'DUSP21 expression',
                            pairwise.display = 's',
                            ggtheme = theme_bw(),
                            boxplot.args = list(width = 0.3, alpha = 0.2, color = "gray"),
                            points.args = list(size = 30, alpha = 1),
                            violin.args = list(width = 0, alpha = 0, color = "lightgray"),
                            plot.type = "boxplot", type = "nonparametric",
                            geom_signif(list(map_signif_level = TRUE)),
                            centrality.point.args = list(size = 1))
ggsave(filename = "results/DUSP21.png",p, scale = 1, width = 6, height = 8)
p


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


ggsave(filename = "results/volcanoplot.png", p, scale = 1, width = 12, height = 8)
p

 #####################################################################





################
#     DESeq2     # T0 genes related to T0_oAHI
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

write.csv(tosave[order(tosave$log2FoldChange),],file="results/T0_oAHIsiggenes.csv",quote=FALSE,row.names = FALSE)


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
ggsave(filename = "results/volcanoplot_continuous.png",p, scale = 1, width = 12, height = 8)
p



#######################
#    GO enrichment    #
#######################

pos_ego <- enrichGO(
  gene = siggenes[siggenes$log2FoldChange >= 0.5,]$name,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  universe = background$name,
  ont = "ALL",  # Biological Process ontology
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE,
)


neg_ego <- enrichGO(
  gene = siggenes[siggenes$log2FoldChange <= -0.5,]$name,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  universe = background$name,
  ont = "ALL",  # Biological Process ontology
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE,
)

p <- enrichplot::dotplot(pos_ego, showCategory=10, title="GO positive LogFoldChange (increase normalised T0_oAHI)") + 
  facet_grid(ONTOLOGY ~ ., scales="free") + theme_bw(base_size = 12)
ggsave(filename = "results/pos_T0_OSA_cont_dotplot.png",p, scale = 1, width = 8, height = 8)
p

p <- enrichplot::dotplot(neg_ego, showCategory=10, title="GO negative LogFoldChange (decrease normalised T0_oAHI)") + 
  facet_grid(ONTOLOGY ~ ., scales="free") + theme_bw(base_size = 12)
ggsave(filename = "results/neg_T0_OSA_cont_dotplot.png",p, scale = 1, width = 8, height = 8)
p






m <- log2(counts(dds,normalized=TRUE) + 1)
#### Negative correlated genes
neg_gene_correlated <- siggenes[siggenes$log2FoldChange < 0,]
neg_gene_correlated <- neg_gene_correlated[order(neg_gene_correlated$log2FoldChange, decreasing = FALSE), ]
transposed_genes <- m[rownames(neg_gene_correlated),]

##### Normal execution
rownames(transposed_genes) <- neg_gene_correlated$name
transposed_genes <- as.data.frame(t(transposed_genes[1:4,]))

##### If only DUSP21
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
  labs(title = "Negatively enrichment for oAHI",
       x = "oAHI",
       y = "Log2 Normalised Counts") + 
  annotate("text", x = 8, y = 9, label = cor_text,size = 4, color = "black") +
  annotate("text", x = 8, y = 8.6, label = p_value_text, size = 4, color = "black")+
  geom_vline(xintercept = 2, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = 2.1, y = 5, label = "Obstructive Sleep Apnea threshold", angle = 90, vjust = 1.5, color = "red", size = 3)
ggsave(filename = "results/neg_T0_oAHI_corr_genes.png", p, scale = 1, width = 6, height = 6)
p
 



p <- ggplot(melted_df, aes(x = T0_oAHI, y = value)) +
  geom_point(size = 5) +
  #geom_boxplot() +
  theme_bw(base_size = 20) +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),   ) +
  facet_grid(~ Date.library.prep, scales="free") +
  #facet_wrap(~ series, scales = 'free_y') +
  #geom_smooth(method = 'glm', aes(color = NULL)) +
  xlim(0, 16) +        
  ylim(0, 10) +
  labs(title = "Negative correlation to AHI",
       x = "T0_oAHI",
       y = "Log2 Normalised Counts") + 
  geom_vline(xintercept = 2, linetype = "dashed", color = "red", linewidth = 2) +
  annotate("text", x = 2.1, y = 5, label = "Obstructive Sleep Apnea threshold", angle = 90, vjust = 1.5, color = "red", size = 6)
ggsave(filename = "results/neg_T0_oAHI_corr_genes_facet.png", p, scale = 1, width = 25, height = 8)
p
