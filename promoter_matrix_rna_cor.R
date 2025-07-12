#Read in packages
library(data.table)
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(progress)
library(rtracklayer)
library(Seurat)
library(Signac)

###make promoter matrices for either human or mouse
raw_mouse_df <- import("C:/Users/kelly/R/spatial_atac_software/data/mm10.refGene.gtf.gz") ##downloaded from UCSC Genome Browser Database
raw_mouse_df<-as.data.frame(raw_mouse_df)

upstream<-5000 ##may also change to 2000 or adjust as needed
downstream<-5000
promoters <- raw_mouse_df %>%
  mutate(
    tss = ifelse(strand == "+", start, end),
    promoter_start = ifelse(strand == "+", tss - upstream, tss - downstream),
    promoter_end = ifelse(strand == "+", tss + downstream, tss + upstream),
    promoter_start = pmax(promoter_start, 0)  # ensure non-negative
  ) %>%
  dplyr::select(chrom = seqnames, start = promoter_start, end = promoter_end, gene = gene_name, strand)

promoters$chrom<-as.character(promoters$chrom)
promoters<-promoters[which(nchar(promoters$chrom)<=5),]
colnames(promoters)[1:3]<-c("seqnames","start","end")

#Orders bed data peak list regions
##input: data frame of the bed data
order_regions<-function(data_df){
  data_df<-as.data.frame(data_df)
  data_df$seqnames<-as.character(data_df$seqnames)
  data_df$seqnames <- ifelse(grepl("^chr", data_df$seqnames), 
                             data_df$seqnames, 
                             paste0("chr", data_df$seqnames))
  data_df$seqnames<-as.character(data_df$seqnames)
  #convert x and y chromosomes to 23 and 24
  x_chrs<-which(grepl("X",data_df$seqnames))
  data_df$seqnames[x_chrs]<-"chr23"
  y_chrs<-which(grepl("Y",data_df$seqnames))
  data_df$seqnames[y_chrs]<-"chr24"
  remove_nums<-suppressWarnings(which(is.na(as.numeric(substr(data_df$seqnames,4,nchar(data_df$seqnames))))))
  keep_nums<-setdiff(1:nrow(data_df),remove_nums)
  data_df<-data_df[keep_nums,]
  chr_nums<-as.numeric(substr(data_df$seqnames,4,nchar(data_df$seqnames)))
  ordered_data_df <- data_df[order(chr_nums,as.numeric(data_df$start),as.numeric(data_df$end)), ]
  #return ordered data frame
  ordered_data_df<-as.data.frame(ordered_data_df)
  ordered_data_df[,2]<-as.numeric(ordered_data_df[,2])
  ordered_data_df[,3]<-as.numeric(ordered_data_df[,3])
  return(ordered_data_df)
}

ordered_promoters<-order_regions(promoters)
cleaned_promoters<-ordered_promoters[!duplicated(ordered_promoters[,c("seqnames","start","end")]),]
cleaned_promoters<-cleaned_promoters[!duplicated(cleaned_promoters[,c("gene")]),]

##common peaks/genes
common<-intersect(cleaned_promoters$gene,rownames(rna_spatial_exp))
rna_spatial_exp_cleaned<-rna_spatial_exp[which(rownames(rna_spatial_exp)%in%common),]
cleaned_promoters<-cleaned_promoters[which(cleaned_promoters$gene%in%common),]
rna_spatial_exp_cleaned<-rna_spatial_exp_cleaned[match(cleaned_promoters$gene,rownames(rna_spatial_exp_cleaned)),]

#MAKE PROMOTER BY SPOT MATRIX (PEAK MATRIX)
cells<-paste0(as.character(spot_sums$cleaned_spots),"-1")
fragments<-CreateFragmentObject(path="GSM6801813_ME13_50um_fragments.tsv.gz",cells=cells)
gr.peaks <- makeGRangesFromDataFrame(cleaned_promoters)
peak_matrix <- suppressWarnings(FeatureMatrix(
  fragments = fragments,
  features = gr.peaks,
  cells = cells 
))

#filtering + matching rows up
colnames(peak_matrix)<-substr(colnames(peak_matrix),1,nchar(colnames(peak_matrix))-2)
keep<-which(paste0(cleaned_promoters$seqnames,"-",cleaned_promoters$start,"-",cleaned_promoters$end)%in%rownames(peak_matrix))
cleaned_promoters<-cleaned_promoters[keep,]
rna_spatial_exp_cleaned<-rna_spatial_exp_cleaned[keep,]
peak_matrix<-peak_matrix[,match(colnames(rna_spatial_exp_cleaned),colnames(peak_matrix))]
keep<-which(rowSums(peak_matrix)>0&rowSums(rna_spatial_exp_cleaned)>0&!duplicated(rownames(peak_matrix)))
peak_matrix<-peak_matrix[keep,]
rna_spatial_exp_cleaned<-rna_spatial_exp_cleaned[keep,]
cleaned_promoters<-cleaned_promoters[keep,]

##Prepare Raw vs. Normalized + TF-IDF peak matrices & Normalized RNA matrix
spot_sums<-spot_sums[match(colnames(peak_matrix),spot_sums$cleaned_spots),]
rna_libsizes<-rna_libsizes[match(spot_sums$cleaned_spots,names(rna_libsizes))]
div_lib_size_matrix<-sweep(peak_matrix,2,spot_sums$Freq,"/")
div_lib_size_matrix<-as.matrix(div_lib_size_matrix)
rna_spatial_exp_cleaned<-rna_spatial_exp_cleaned[,match(names(rna_libsizes),colnames(rna_spatial_exp_cleaned))]
rna_div_lib_size_matrix<-sweep(rna_spatial_exp_cleaned,2,rna_libsizes,"/")
rna_div_lib_size_matrix<-as.matrix(rna_div_lib_size_matrix)
seurat_object <- CreateSeuratObject(counts = peak_matrix, assay = "peaks")
seurat_object <- RunTFIDF(seurat_object, assay = "peaks",method=1) #default method
tfidf_matrix <- GetAssayData(seurat_object, layer = "data", assay = "peaks")
tfidf_matrix<-as.matrix(tfidf_matrix)

peak_matrix<-as.matrix(peak_matrix)
rna_spatial_exp_cleaned<-as.matrix(rna_spatial_exp_cleaned)

##Calculate row-wise correlations across datasets (gene vs. peak) for spatial-ATAC-RNA-seq and spatial-Mux-seq
rowwise_cor <- function(mat1, mat2) {
  stopifnot(all(dim(mat1) == dim(mat2)))
  
  # Center each row (subtract row mean)
  mat1_centered <- mat1 - rowMeans(mat1, na.rm = TRUE)
  mat2_centered <- mat2 - rowMeans(mat2, na.rm = TRUE)
  
  # Compute numerator: sum of products
  numerator <- rowSums(mat1_centered * mat2_centered, na.rm = TRUE)
  
  # Compute denominator: product of sqrt sum of squares
  denom <- sqrt(rowSums(mat1_centered^2, na.rm = TRUE) * rowSums(mat2_centered^2, na.rm = TRUE))
  
  # Final correlation
  cor_vals <- numerator / denom
  return(cor_vals)
}

#Correlation Calculations
raw_atac_norm_rna<-rowwise_cor(rna_div_lib_size_matrix,peak_matrix)
raw_both<- rowwise_cor(rna_spatial_exp_cleaned,peak_matrix)
norm_both<- rowwise_cor(rna_div_lib_size_matrix,div_lib_size_matrix)
norm_atac_raw_rna<- rowwise_cor(rna_spatial_exp_cleaned,div_lib_size_matrix)
tfidf_atac_raw_rna<-rowwise_cor(tfidf_matrix,rna_spatial_exp_cleaned)
tfidf_atac_norm_rna<-rowwise_cor(tfidf_matrix,rna_div_lib_size_matrix)


#Plot violins (figure 2) for correlation calculations#
correlation_df<-cbind(raw_atac_norm_rna,raw_both,norm_both,norm_atac_raw_rna,tfidf_atac_raw_rna,tfidf_atac_norm_rna)
means<-colMeans(correlation_df,na.rm=TRUE)
df <- data.frame(
  value = c(raw_atac_norm_rna, raw_both, norm_both, norm_atac_raw_rna, tfidf_atac_raw_rna, tfidf_atac_norm_rna),
  group = factor(rep(
    c("raw_atac_norm_rna", "raw_atac_raw_rna", "norm_atac_norm_rna",
      "norm_atac_raw_rna", "tfidf_atac_raw_rna", "tfidf_atac_norm_rna"),
    each = nrow(peak_matrix)
  ))
)

# Set factor levels for consistent ordering
df$group <- factor(df$group, levels = c(
  "raw_atac_raw_rna", "norm_atac_raw_rna", "tfidf_atac_raw_rna",
  "raw_atac_norm_rna", "norm_atac_norm_rna", "tfidf_atac_norm_rna"
))

# ✅ Split into raw and normalized groups
df_raw <- df %>% filter(grepl("raw_rna", group))
df_norm <- df %>% filter(grepl("norm_rna", group))

df_raw$group<-factor(rep(c("raw_atac","norm_atac","tfidf_atac"),each=nrow(peak_matrix)))
df_norm$group<-factor(rep(c("raw_atac","norm_atac","tfidf_atac"),each=nrow(peak_matrix)))
df_raw$group <- factor(df_raw$group, levels = c(
  "raw_atac", "norm_atac", "tfidf_atac"
))

df_norm$group <- factor(df_norm$group, levels = c(
  "raw_atac", "norm_atac", "tfidf_atac"
))

# Set colors for each subset
raw_colors <- c("raw_atac" = "pink","norm_atac" = "pink","tfidf_atac" = "pink")
norm_colors <- c("raw_atac" = "skyblue","norm_atac" = "skyblue","tfidf_atac" = "skyblue")

# Calculate mean for each group
label_df_raw <- df_raw %>%
  group_by(group) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            ymax = max(value, na.rm = TRUE)) %>%
  mutate(label_y = ymax -0.02)

label_df_norm <- df_norm %>%
  group_by(group) %>%
  summarise(mean_value = mean(value, na.rm = TRUE),
            ymax = max(value, na.rm = TRUE)) %>%
  mutate(label_y = ymax - 0.03)

# Raw RNA Violin Plots
ggplot(df_raw, aes(x = group, y = value, fill = group)) +
  geom_violin(trim = FALSE) +
  geom_point(data = label_df_raw, aes(x = group, y = mean_value), colour = "black", size = 2, inherit.aes = FALSE) +
  geom_label(
    data = label_df_raw,
    aes(x = group, y = label_y, label = sprintf("%.3f", mean_value)),
    size = 15,
    fill = "white",
    label.size = 0.2,
    alpha = 0.8,
    position = position_nudge(x = 0.05),
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = raw_colors) +
  theme_minimal() +
  labs(title = "", x = "", y = "Correlation") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 45),
    axis.title = element_text(size = 50),
    plot.title = element_text(size = 30),
    axis.text = element_text(size = 45),
    panel.grid = element_blank(),
    legend.position = "none"
  )

#Normalized RNA Violin Plots
ggplot(df_norm, aes(x = group, y = value, fill = group)) +
  geom_violin(trim = FALSE) +
  geom_point(data = label_df_norm, aes(x = group, y = mean_value), colour = "black", size = 2, inherit.aes = FALSE) +
  geom_label(
    data = label_df_norm,
    aes(x = group, y = label_y, label = sprintf("%.3f", mean_value)),
    size = 15,
    fill = "white",
    label.size = 0.2,
    alpha = 0.8,
    position = position_nudge(x = 0.05),
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = norm_colors) +
  theme_minimal() +
  labs(title = "", x = "", y = "Correlation") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text = element_text(size = 45),
    axis.title = element_text(size = 50),
    plot.title = element_text(size = 30),
    axis.text = element_text(size = 45),
    panel.grid = element_blank(),
    legend.position = "none")

#save matrices for future analyses
saveRDS(peak_matrix,"peak_matrix.rds")
saveRDS(rna_spatial_exp_cleaned,"rna_cleaned_matrix.rds")
