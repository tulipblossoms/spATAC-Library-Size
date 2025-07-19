#read in packages
library(data.table)
library(dplyr)
library(ggplot2)
library(ggvenn)
library(pheatmap)
library(progress)
library(Seurat)
library(Signac)

#read in and assign clusters for human brain
cluster_assignments<-read.csv("cluster/hippocampus.csv")
cluster_assignments<-cluster_assignments[match(metadata$AACGTGATAACGTGAT,cluster_assignments$X),]
cluster1_idx <- which(cluster_assignments$ATAC_clusters !="C4")
cluster2_idx <- which(cluster_assignments$ATAC_clusters == "C4")

#assign clusters for Slide-tags
cluster1_idx<-which(metadata$cluster_value==1)
cluster2_idx<-which(metadata$cluster_value==2)

#read in human peak matrix & align spots together
peak_matrix1<-readRDS("spatial_cor_data/human_brain_5kb_promoter_matrix.rds")
colnames(peak_matrix1)<-substr(colnames(peak_matrix1),1,nchar(colnames(peak_matrix1))-2)
peak_matrix1<-peak_matrix1[,match(metadata$AACGTGATAACGTGAT,colnames(peak_matrix1))]
spot_sums<-spot_sums[match(metadata$AACGTGATAACGTGAT,spot_sums$cleaned_spots),]

#pre-filter promoters/peaks
#10% non-zero expression#
inds<-which(rowSums(peak_matrix1 != 0) >= 250)

#threshold for library size cutoff
threshold<-quantile(rowSums(peak_matrix1),0.25)
second_inds <- which(rowSums(peak_matrix1) > threshold)
keep<-intersect(inds,second_inds)

#subsection matrix
peak_matrix<-peak_matrix1[keep,]
div_lib_size_matrix1<-sweep(peak_matrix1,2,spot_sums$Freq,"/")
div_lib_size_matrix<-div_lib_size_matrix1[keep,]
seurat_object <- CreateSeuratObject(counts = peak_matrix1, assay = "peaks")
seurat_object <- RunTFIDF(seurat_object, assay = "peaks",method=1) #default method
tfidf_matrix1 <- GetAssayData(seurat_object, layer = "data", assay = "peaks")
tfidf_matrix1<-as.matrix(tfidf_matrix1)
tfidf_matrix<-tfidf_matrix1[keep,]

##determine differential peaks##
#Wilcoxon function
wilcox_func <- function(group1, group2) {
  test <- wilcox.test(as.numeric(group1), as.numeric(group2))
  fdr_val <- p.adjust(test$p.value, method = "fdr")
  alpha <- 0.05
  return(fdr_val<alpha)
}

##determine whether or not peaks are differential between the matrices##
#normal matrix
normal_counter <- sapply(seq_len(nrow(peak_matrix)), function(i) {
  spot_strand<-peak_matrix[i,]
  g1 <- as.numeric(spot_strand[cluster1_idx])
  g2 <- as.numeric(spot_strand[cluster2_idx])
  
  result <- tryCatch({
    wilcox_func(as.numeric(g1), as.numeric(g2))
  }, error = function(e) {
    NA
  })
  return(result)
})

#divide by library size matrix
residual_counter <- sapply(seq_len(nrow(div_lib_size_matrix)), function(i) {
  spot_strand <- div_lib_size_matrix[i, ]
  g1 <- as.numeric(spot_strand[cluster1_idx])
  g2 <- as.numeric(spot_strand[cluster2_idx])
  
  result <- tryCatch({
    wilcox_func(g1, g2)
  }, error = function(e) {
    NA 
  })
  return(result)
})

#TF-IDF matrix
tfidf_counter <- sapply(seq_len(nrow(tfidf_matrix)), function(i) {
  spot_strand <- tfidf_matrix[i, ]
  g1 <- as.numeric(spot_strand[cluster1_idx])
  g2 <- as.numeric(spot_strand[cluster2_idx])
  
  result <- tryCatch({
    wilcox_func(g1, g2)
  }, error = function(e) {
    NA
  })
  return(result)
})


#calculate log2 fold changes between clusters
log2_fc_normal <- apply(peak_matrix, 1, function(row) {
  mean1 <- mean(row[cluster1_idx], na.rm = TRUE)
  mean2 <- mean(row[cluster2_idx], na.rm = TRUE)
  return(log2((mean2 + 1e-6) / (mean1 + 1e-6)))
})


log2_fc_residual <- apply(div_lib_size_matrix, 1, function(row) {
  mean1 <- mean(row[cluster1_idx], na.rm = TRUE)
  mean2 <- mean(row[cluster2_idx], na.rm = TRUE)
  return(log2((mean2 + 1e-6) / (mean1 + 1e-6)))
})


log2_fc_tfidf <- apply(tfidf_matrix, 1, function(row) {
  mean2 <- mean(row[cluster2_idx], na.rm = TRUE)
  mean1 <- mean(row[cluster1_idx], na.rm = TRUE)
  return(log2((mean2 + 1e-6) / (mean1 + 1e-6)))
  })

#assign significantly up-regulated, significantly down-regulated, or no significant difference code#
normal_result <- rep(NA_integer_, length(normal_counter))
normal_result[!normal_counter] <- 0
remaining <- is.na(normal_result)
normal_result[remaining & log2_fc_normal > 0] <- 2
normal_result[remaining & log2_fc_normal <= 0] <- 1


residual_result <- rep(NA_integer_, length(residual_counter))
residual_result[!residual_counter] <- 0
remaining <- is.na(residual_result)
residual_result[remaining & log2_fc_residual > 0] <- 2
residual_result[remaining & log2_fc_residual <= 0] <- 1


tfidf_result <- rep(NA_integer_, length(tfidf_counter))
tfidf_result[!tfidf_counter] <- 0
remaining <- is.na(tfidf_result)
tfidf_result[remaining & log2_fc_tfidf > 0] <- 2
tfidf_result[remaining & log2_fc_tfidf <= 0] <- 1


###code to generate figure 4###

#venn diagrams#
ggvenn(
  list(
    "Raw" = which(normal_result!=0),
    "Norm" = which(residual_result!=0)
  ),
  fill_color = c("skyblue", "lightpink"),
  stroke_size = 0.5,
  text_size = 16,
  set_name_size = 16
)

ggvenn(
  list(
    "Raw" = which(normal_result!=0),
    "TFIDF" = which(tfidf_result!=0)
  ),
  fill_color = c("skyblue", "lightpink"),
  stroke_size = 0.5,
  text_size = 16,
  set_name_size = 16
)

## % difference stacked barplots ##
comparison <- c("Norm", "TFIDF")
changed <- c(length(which(normal_result!=0&residual_result!=0&normal_result!=residual_result)), length(which(normal_result!=0&tfidf_result!=0&normal_result!=tfidf_result)))       # Number of spots where direction changed
not_changed <- c(length(which(normal_result!=0&residual_result!=0&normal_result==residual_result)), length(which(normal_result!=0&tfidf_result!=0&normal_result==tfidf_result)))   # Number where direction stayed the same

# Convert to percentages
total <- changed + not_changed
pct_changed <- changed / total * 100
pct_not_changed <- not_changed / total * 100

# Long format for ggplot2
df <- data.frame(
  Comparison = rep(comparison, each = 2),
  Direction  = rep(c("Changed", "Not Changed"), 2),
  Percent    = c(pct_changed[1], pct_not_changed[1], pct_changed[2], pct_not_changed[2]),
  Count      = c(changed[1], not_changed[1], changed[2], not_changed[2])
)
df$Direction <- factor(df$Direction, levels = c("Not Changed", "Changed"))
df$adjust_vjust <- ifelse(df$Direction == "Not Changed",
                          0.1,   # Move this label upward
                          1.05)   # Default for others
ggplot(df, aes(x = Comparison, y = Percent, fill = Direction)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = Count,vjust = adjust_vjust), position = position_stack(vjust = 0.5), size = 16) +
  scale_fill_manual(values = c("Changed" = "tomato", "Not Changed" = "steelblue")) +
  labs(y = "Percent", x = "", fill = "Direction", title = "") +
  theme_minimal() +
  theme(
    text = element_text(size = 60),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 56),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none"
  )


##Before vs after normalization specific example comparisons##

#human hippocampus example in figure 4#
lib_matrix <- matrix(0, nrow = 50, ncol = 50)

#match spot sums with metadata
spot_sums<-spot_sums[match(metadata$ACAAGCTAAACGTGAT,spot_sums$cleaned_spots),]
# Assign values based on row and column indices
for (i in seq_len(nrow(metadata))) {
  r <- as.numeric(metadata$row[i])+1
  c <- as.numeric(metadata$col[i])+1
  lib_matrix[r,c] <- peak_matrix[1042,i] #may also change to div_lib_size_matrix or tfidf_matrix to view
}
lib_matrix<-lib_matrix[,ncol(lib_matrix):1]
pheatmap(lib_matrix,cluster_cols=FALSE,cluster_rows=FALSE,fontsize=32)


#Slide-tags example in figure 4#
metadata$orig<-as.numeric(peak_matrix[783,])
metadata$norm<-as.numeric(div_lib_size_matrix[783,])
metadata$tfidf<-as.numeric(tfidf_matrix[783,])
ggplot(metadata, aes(x = X, y = Y, color = orig, size = orig)) +
  geom_point(shape = 16) +
  scale_color_gradientn(
    colors = c("#E8E8E8", "#FEEEE1", "#FFD799", "#7B68EE", "#5E3C99"),
    values = scales::rescale(quantile(metadata$orig, probs = c(0, 0.885, 0.95, 0.99, 1))),
    limits = range(quantile(metadata$orig, probs = c(0, 0.885, 0.95, 0.99, 1))),
    oob = scales::squish
  ) +
  scale_size_continuous(
    range = c(2, 5),
    trans = "sqrt"
  ) +
  labs(color = "") +
  guides(size = "none") + 
  theme_void() +
  theme(
    legend.key.height = unit(2, "cm") 
  )+
  theme(
    legend.title = element_text(size = 36),
    legend.text = element_text(size = 36)
  )

