# Function to create directories if they don't exist
create_directories <- function(paths) {
  for (dir_path in paths) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
      cat("Directory created at:", dir_path, "\n")
    } else {
      cat("Directory already exists at:", dir_path, "\n")
    }
  }
}




# Function to plot the histogram plot of the library sizes per sample
plot_library_sizes <- function(countdata, plotsDir) {

  # Sum all the transcript counts together per sample
  column_sums <- sapply(countdata, sum)

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
    scale_y_continuous(labels = scales::comma) +
    theme_bw(base_size = 15) +
    theme(plot.title = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(),  
      axis.text.x = element_text(colour="black",size=12,angle=45,hjust=.5,vjust=.7,face="plain",family = "serif"),
    )

  ggsave(paste0(plotsDir, 'librarysizes.png'), p, scale = 2, width = 14, height = 8)

  return (sum_df)

}





# Function to plot the log counts per million distribution before and after low expression filtering
density_plot_logcpm <- function(data, plotsDir){
  
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
  ggsave(paste0(plotsDir, 'densityplots.png'), p, scale = 1, width = 16, height = 7)
  
  return(data)

}




# Function to plot the transcript distribution on a sample basis for further outlier inspection
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
      axis.text.x = element_text(colour="black",size=12,angle=90,hjust=0.6,vjust=.7,face="plain", family = "serif")
    )

  ggsave(plot_filename, p, scale = 1, width = 8, height = 6)

}




# Function to plot the batch effect separation using Principal Component Analysis (PCA)
plot_pca_batcheffect <- function(countdata, coldata, plotsDir) { 


  # Calculate Log counts per million with edgeR
  y.all <- DGEList(counts=countdata, genes = rownames(countdata))
  y.all <- calcNormFactors(y.all)
  lcpm <- cpm(y.all, log=TRUE, normalized.lib.sizes = TRUE)


  # Run PCA
  pca <- prcomp(t(lcpm), center = TRUE)
  # Combine PCA coordinates with the metadata from the DGEList
  to_plot <- data.frame(pca$x)
  # Calculate how many % of total variance is explained by each principal component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )*100
  # We focus here on PC1 and PC2
  use.pcs <- c(1,2)


  to_plot$date <- coldata$Date.library.prep
  to_plot$time <- coldata$time
  to_plot$batch <- as.factor(ifelse(to_plot$PC1 > 20, "batch1", "batch2"))
  to_plot$donor <- coldata$donor



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
      y =labs[2]) +
    annotate("text", x = 20, y = -40, label = "Batch effect", angle = 90, vjust = 1.5, color = "red", size = 5)


  ggsave(paste0(plotsDir, 'pca_batcheffect.png') , p, scale = 1, width = 10, height = 8)

  return (to_plot)

}




# Function to plot test if there is a difference in sequencing depth between the batch groups
plot_library_sizes_test <- function(countdata, pca_df, plotsDir) {

  # Sum all the transcript counts together per sample
  column_sums <- sapply(countdata, sum)

  # Create a new dataframe with column names and their sums
  sum_df <- data.frame(
    Samples = names(column_sums),
    Library.size = column_sums
  )

  sum_df$Timepoint <- as.data.frame(do.call(rbind, strsplit(as.character(sum_df$Sample), "_")))$V2
  sum_df <- sum_df[sum_df$Timepoint == 'T0',]
  sum_df$croppped <- sub("_.*", "", sum_df$Samples)


  sum_df <- sum_df[rownames(sum_df) %in% rownames(pca_df), ]
  sum_df$batch <- pca_df$batch

  p <- ggbetweenstats(data = sum_df,
                      x = batch,
                      y = Library.size,
                      ylab = 'Library Size',
                      xlab = 'Batch',
                      title = 'Statistical comparison of library sizes by Batch',
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
                      centrality.point.args = list(size = 3,color = "#8a0f00"),
                      centrality.label.args = list(size = 3)) +
                      scale_y_continuous(labels = scales::comma)

  ggsave(paste0(plotsDir, "librarysizes_test.png"), p, scale = 1, width = 8, height = 8)


}




# Function to test correlations between oAHI and metadata variables
correlation_tests <- function(coldata, variables_totest, plotsDir) {


  for (variable in variables_totest) {


    height <- 6
    width <- 10
    text_x <- 90
    text_y <- 10
    threshold_x <- 105
    threshold_y <- 2.5

    if (variable == 'T0_BMI') {
      width <- 8
      text_x <- 27
      text_y <- 12
      threshold_x <- 45
    } else if (variable == 'T0_Age') {
      text_x <- 10
      threshold_x <- 16.3
    } else if (variable == 'T0_Glucose') {
      threshold_x <- 100
    } else if (variable == 'T0_Insulin') {
      threshold_x <- 305
    }


    cor_test <- cor.test(coldata$T0_oAHI, coldata[[variable]], method = "spearman", exact = FALSE)
    cor_coef <- cor_test$estimate
    p_value <- cor_test$p.value
    cor_text <- paste0("Spearman rho: ", round(cor_coef, 2))
    p_value_text <- paste0("p-value: ", format(p_value, scientific = FALSE, digits = 4))


    p <- ggplot(coldata, aes(x=.data[[variable]], y=T0_oAHI, shape=T0_OSA, size = 2))+
      geom_point() +
      theme_bw(base_size = 15) +
      labs(x = sub("^T0_", "", variable),
           y = "AHI") +
      theme(legend.position = "none",
            axis.text.x = element_text(colour="black",size=12,angle=45,hjust=0.6,vjust=.7,face="plain",),
      ) +
      geom_hline(yintercept = 2, linetype = "dashed", color = "red", linewidth = 1) +
      annotate("text", x = text_x, y = text_y, label = cor_text,size = 6, color = "black") +
      annotate("text", x = text_x, y = (text_y - 1), label = p_value_text, size = 6, color = "black")+
      annotate("text", x = threshold_x, y = threshold_y, label = "Obstructive Sleep Apnea threshold", angle = 0, color = "red", size = 4)
    p
    ggsave(paste0(plotsDir, sub("^T0_", "", variable), "_AHI.png"), p, scale = 1, width = width, height = height)

  }


  gender_data <- coldata
  gender_data$Sex <- ifelse(gender_data$Sex == 1, 'Male', 'Female')
  p <- ggplot(gender_data, aes(x = Sex, y = T0_oAHI)) +
    geom_jitter(width = 0.2, aes(shape = T0_OSA)) +
    geom_boxplot(outlier.shape = NA, varwidth = TRUE, alpha = 0.2, width = 0.15) +
    labs(title = "",
         x = "Sex",
         y = "AHI") + 
    theme_bw(base_size = 12) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5)) +
    stat_compare_means(method = "wilcox.test", aes(label = ..p.signif..), 
                       label.y = max(coldata1$T0_oAHI, na.rm = TRUE) -.8, #
                       label.x.npc = 'center', 
                       size = 6) + # Adjust size if needed
    stat_compare_means(method = "wilcox.test", label = "p.format", 
                       label.y = max(coldata1$T0_oAHI, na.rm = TRUE), 
                       label.x.npc = 0.42,
                       size = 5) +
    geom_hline(yintercept = 2, linetype = "dashed", color = "red", linewidth = 1) +
    annotate("text", x = 1.5, y = 2.5, label = "Obstructive Sleep Apnea threshold", angle = 0, color = "red", size = 4)

  
  ggsave(paste0(plotsDir, "Sex_AHI.png"), p, scale = 1, width = 5, height = 5)


}




# Draw GO enrichment plots (if no enriched terms, the error will be caught)
go_enrichment <- function(siggenes, background, plotsname, plotsDir){

  tryCatch(
  {
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
    p <- enrichplot::dotplot(pos_ego, showCategory=10, title="GO positive LogFoldChange (increase with OSA)") + 
                     facet_grid(ONTOLOGY ~ ., scales="free") + theme_bw(base_size = 12)

    ggsave(paste0(plotsDir, "pos_", plotsname, ".png"), p, scale = 1, width = 8, height = 6)
  },
  error = function(e) {
    message("An error occurred: ", e$message)
  }
  )


  tryCatch(
  {
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


  p <- enrichplot::dotplot(neg_ego, showCategory=10, title="GO negative LogFoldChange (decrease with OSA)") +
    facet_grid(ONTOLOGY ~ ., scales="free") + theme_bw(base_size = 12)

  ggsave(paste0(plotsDir, "neg_", plotsname, ".png"), p, scale = 1, width = 8, height = 8)
  },
  error = function(e) {
    message("An error occurred: ", e$message)
  }
  )

}