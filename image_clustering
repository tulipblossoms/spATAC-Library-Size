##read in packages
library(data.table)
library(dplyr)
library(EBImage)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(readr)

##########spatial-ATAC-RNA-sequencing & spatial-MUX-sequencing datasets analysis##########
#read in image (see example on mouse embryo below)
img<-readImage("ME13_50um_50bc/ME13_50um_spatial/tissue_lowres_image.png")
pixel_values <- channel(img,"gray")

#smooth and cluster high resolution image
sliding_window_function <- function(slide_k, cluster_num,valid_idx=NA) {
  # Smoothing with mean filter
  kernel <- matrix(1, nrow = slide_k, ncol = slide_k)
  kernel <- kernel / sum(kernel)
  smoothed_img <- filter2(img, kernel)
  smoothed_matrix <- as.array(smoothed_img[,,1])
  pixel_vector <- as.vector(smoothed_matrix)
  set.seed(1)
  # Flatten and apply k-means clustering
  if(!is.na(valid_idx)){
    km <- kmeans(pixel_vector[valid_idx], centers = cluster_num)
    cluster_vector <- rep(0, length(pixel_vector))
    cluster_vector[valid_idx] <- km$cluster
    cluster_matrix <- matrix(cluster_vector, nrow = nrow(smoothed_matrix), ncol = ncol(smoothed_matrix))
    rgb_array <- array(0, dim = c(nrow(cluster_matrix), ncol(cluster_matrix), 3))
  }else{
    km <- kmeans(pixel_vector, centers = cluster_num)
    cluster_matrix <- matrix(km$cluster, nrow = nrow(smoothed_matrix), ncol = ncol(smoothed_matrix))
  return(list(cluster_matrix=cluster_matrix))
  }
}

img_cluster<-sliding_window_function(51,5)
cluster_matrix <-img_cluster$cluster_matrix

# Get original dimensions
nrow_img <- nrow(cluster_matrix)
ncol_img <- ncol(cluster_matrix)

# Initialize RGB image
rgb_array <- array(0, dim = c(nrow_img, ncol_img, 3))

# Fill in RGB channels
for (i in 1:5) {
  mask <- cluster_matrix == i
  rgb_array[,,1][mask] <- rgb_colors[i, 1]  # Red
  rgb_array[,,2][mask] <- rgb_colors[i, 2]  # Green
  rgb_array[,,3][mask] <- rgb_colors[i, 3]  # Blue
}

# Convert to EBImage RGB object
color_img <- Image(rgb_array, colormode = "Color")

# Display the image
display(color_img, method = "raster")


#use for 5 cluster analyses
cluster_colors <- c("#000000","red","blue","#006400","purple","orange")

#use for 10 cluster analyses
cluster_colors<-c("black","red","blue","darkgreen","#00CED1","orange","#A65628","#F781BF", "#999999","#66C2A5","#FFD700","purple")
custom_colors <- t(col2rgb(cluster_colors)) / 255

background_remove_cluster<-sliding_window_function(51,5,which(cluster_matrix!=4))

# Fill image
for (i in 0:5) {
  mask <- cluster_matrix == i
  rgb_array[,,1][mask] <- custom_colors[i + 1, 1]
  rgb_array[,,2][mask] <- custom_colors[i + 1, 2]
  rgb_array[,,3][mask] <- custom_colors[i + 1, 3]
}

# Display
display(Image(rgb_array, colormode = "Color"))

#downsample image function
downsample_img<-function(large_dim, k,cluster_matrix){
  block_size <- large_dim / k 
  downsampled_matrix <- matrix(0, nrow = k, ncol = k)
  # Loop over each block
  for (i in 1:k) {
    for (j in 1:k) {
      row_start <- floor((i - 1) * block_size) + 1
      row_end <- floor(i * block_size)
      col_start <- floor((j - 1) * block_size) + 1
      col_end <- floor(j * block_size)
      
      block <- cluster_matrix[row_start:row_end, col_start:col_end]
      downsampled_matrix[i, j] <- get_mode(as.vector(block))
    }
  }
  return(downsampled_matrix)
}

#get mode function
get_mode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#example on mouse embryo data
downsampled_matrix<-downsample_img(1715,50,cluster_matrix)
transposed_img_matrix<-t(downsampled_matrix)

#mouse brain only
downsampled_matrix<-downsampled_matrix[,ncol(downsampled_matrix):1]

breaks <- 0:6
pheatmap(
  transposed_img_matrix,
  color =cluster_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  legend = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  main = "",
  fontsize=32
)

clusters<-as.vector((downsampled_matrix))

###READ IN ATAC + RNA DATA + CALCULATE LIBRARY SIZE###

metadata<-read.csv("spatial/tissue_positions_list.csv")
metadata<-rbind(colnames(metadata),metadata)
metadata[1,2:6]<-c(1,0,9,11,1768)
colnames(metadata)[3:4]<-c("row","col")

#RNA
rna_spatial_exp <- fread("ME13_50um_50bc/GSM6799937_ME13_50um_matrix_merge.tsv/ME13_50um_matrix_merge.tsv")
rna_spatial_exp<-as.data.frame(rna_spatial_exp)
rownames(rna_spatial_exp)<-rna_spatial_exp$V1
rna_spatial_exp<-rna_spatial_exp[,2:ncol(rna_spatial_exp)]

#ATAC
atac_spatial_exp<-fread("GSM8189706_ME13_50um_3_ATAC.fragments.tsv/GSM8189706_ME13_50um_3_ATAC.fragments.tsv")
colnames(atac_spatial_exp)[1:4]<-c("seqnames","start","end","spot")
atac_spatial_exp$spot<-sub("*-.","",atac_spatial_exp$spot)
cleaned_spots<-sub("*-.","",atac_spatial_exp$spot)

#ATAC Library Size
spots<-unique(cleaned_spots)
spot_sums<-table(cleaned_spots)
spot_sums<-as.data.frame(spot_sums)

##lib size pheatmap
lib_matrix <- matrix(NA, nrow = 50, ncol = 50)
spot_sums<-spot_sums[match(metadata$ACAAGCTAAACGTGAT,spot_sums$cleaned_spots),]

# Assign values based on row and column indices
for (i in seq_len(nrow(metadata))) {
  r <- as.numeric(metadata$row[i])+1
  c <- as.numeric(metadata$col[i])+1
  lib_matrix[r,c] <- spot_sums$Freq[i]
}

lib_matrix<-lib_matrix[,ncol(lib_matrix):1]
pheatmap(lib_matrix,cluster_cols=FALSE,cluster_rows=FALSE,fontsize=32,main=" ")

##SPECIFICALLY FOR THE MOUSE BRAIN 20um
pheatmap(t(log2(lib_matrix))[,ncol(lib_matrix):1],cluster_cols=FALSE,cluster_rows=FALSE)

##violin + t-test
libsizes<-as.vector(t(lib_matrix))
clusters<-as.vector(downsampled_matrix)
df<-data.frame(clusters=clusters,libsizes=libsizes)
df$clusters<-as.factor(df$clusters)
df<-df[which(df$clusters!=0),]


my_comparisons <- list(
  c("1", "2"), c("1", "3"), c("1","4"), c("1","5"), c("2", "3"), c("2","4"), c("2","5"), c("3","4"), c("3","5"), c("4","5"))
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns")).

ggplot(df, aes(x = clusters, y = log2(libsizes+1), fill = clusters)) +
  geom_violin(trim = FALSE, color = "black") +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_fill_brewer(palette = "Set1") +
  stat_compare_means(comparisons = my_comparisons, method = "t.test",size=8,symnum.args=symnum.args) +  # pairwise comparisons
  theme_minimal() +
  labs(
    title="",
    x = "Cluster",
    y = "log2 Library Size"
  )+theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+theme(
    text = element_text(size = 16),          
    axis.title = element_text(size = 30),    
    axis.text = element_text(size = 26),     
    plot.title = element_text(size = 18, hjust = 0.5),
    strip.text = element_text(size = 14)     
  )+theme(legend.position = "none")

##GENERATE STACKED BARPLOT
df$quantile <- cut(df$libsizes,breaks = quantile(df$libsizes, probs = seq(0, 1, 0.1), na.rm = TRUE),include.lowest = TRUE,labels = paste0("Q", 1:10))

plot_df <- df %>%
  group_by(quantile, clusters) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(quantile) %>%
  mutate(percent = count / sum(count) * 100)

ggplot(plot_df, aes(x = quantile, y = percent, fill = clusters)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = c(
    "1" = "red",
    "2" = "blue",
    "3" = "#4CAF50",
    "4" = "purple",
    "5" = "orange"
  )) +
  theme_minimal() +
  labs(
    title = " ",
    x = "Library Size Quantile",
    y = "Percent",
    fill = "Cluster"
  )+ theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+theme(
    text = element_text(size = 26), 
    axis.title = element_text(size = 30), 
    axis.text = element_text(size = 26),  
    plot.title = element_text(size = 18, hjust = 0.5), 
    strip.text = element_text(size = 14)
  )+ theme(legend.position = "none")


#Calculate Correlations with RNA Library Sizes

#Read in RNA dataset#
#Calculate library sizes + match with ATAC datasets
rna_libsizes<-colSums(rna_spatial_exp)
names(rna_libsizes)<-colnames(rna_spatial_exp)
spot_sums1<-spot_sums[which(spot_sums$cleaned_spots%in%names(rna_libsizes)),]
rna_libsizes1<-rna_libsizes[match(spot_sums1$cleaned_spots,names(rna_libsizes))]

#Calculate + Plot Correlation
plot_df <- data.frame(
  atac = log2(as.numeric(spot_sums1$Freq) + 1),
  rna = log2(as.numeric(rna_libsizes1$freq) + 1)
)

# Fit model between ATAC & RNA library size
fit <- lm(atac ~ rna, data = plot_df)

# Compute correlation
correlation <- cor(plot_df$atac, plot_df$rna)

# Plot
ggplot(plot_df, aes(x = rna, y = atac)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1.5) +
  labs(
    title = paste0("Cor = ", round(correlation, 2)),
    x = "Log2 RNA Library Size",
    y = "Log2 ATAC Library Size"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 42, hjust = 0.4),
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 34),panel.grid = element_blank() ) # 👈 removes both major and minor grid lines
  
##########Pipeline for Slide-tags Dataset##########
#read in image
img<-readImage("C:/Users/kelly/R/spatial_atac_software/data/SPATE/Data/slide_tags/raw/tissue_lowres_img.png")
pixel_values <- channel(img,"gray")
cluster_matrix<-matrix(NA,443,421)
img_height <- nrow(cluster_matrix)  
img_width  <- ncol(cluster_matrix)

#read in metadata and assign clusters
#read in metadata
metadata<-read.csv("C:/Users/kelly/R/spatial_atac_software/data/SPATE/Data/slide_tags/raw/HumanMelanomaMultiome_spatial.csv")
metadata<-metadata[2:nrow(metadata),]
metadata$NAME<-substr(metadata$NAME,1,nchar(metadata$NAME)-2)
colnames(metadata)[1]<-"name"

# Step 2: Rescale spot coordinates to image indices
# Normalize x and y to [1, width] and [1, height]
metadata$X<-as.numeric(metadata$X)
metadata$Y<-as.numeric(metadata$Y)
x_scaled <- (metadata$X - min(metadata$X)) / (max(metadata$X) - min(metadata$X)) * (img_width - 1) + 1
y_scaled <- (metadata$Y - min(metadata$Y)) / (max(metadata$Y) - min(metadata$Y)) * (img_height - 1) + 1

# Step 3: Round to nearest pixel index
metadata$col_img <- round(x_scaled)
metadata$row_img <- round(y_scaled)

center_row <- nrow(cluster_matrix) / 2
center_col <- ncol(cluster_matrix) / 2

angle_deg <- 38
slope <- tan(angle_deg * pi / 180)
intercept <- center_row - slope * center_col

# Initialize cluster matrix
cluster_matrix <- matrix(0, nrow = 443, ncol = 421)

for (i in 1:nrow(cluster_matrix)) {
  for (j in 1:ncol(cluster_matrix)) {
    if (i < slope * j + intercept) {
      cluster_matrix[i, j] <- 1  # upper side of 38° diagonal from BL to TR
    } else {
      cluster_matrix[i, j] <- 2
    }
  }
}

metadata$cluster_value <- mapply(
  function(r, c) cluster_matrix[r, c],
  metadata$row_img,
  metadata$col_img
)


#read in ATAC-seq dataset
atac_df<-read.csv(gzfile("HumanMelanomaMultiome_atac.csv.gz"), row.names = 1)
spot_sums<-data.frame(name=colnames(atac_df),Freq=colSums(atac_df))
spot_sums$name<-substr(spot_sums$name,1,nchar(spot_sums$name)-2)
colnames(atac_df)<-substr(colnames(atac_df),1,nchar(colnames(atac_df))-2)
char_before_dash <- regexpr("-", rownames(atac_df)) - 1
keep<-which(char_before_dash<=5)
atac_df<-atac_df[keep,]
keep<-which(rowSums(atac_df)>0)
atac_df<-atac_df[keep,]

#Visualize library size
ggplot(metadata, aes(x = as.numeric(X), y = as.numeric(Y), fill = as.numeric(libsizes))) +
  geom_tile(width = 100, height = 100) +
  coord_fixed() +
  scale_fill_manual(color = colorRampPalette(c("navy", "white", "firebrick3"))(5)) +  
  xlim(0, 6000) +
  ylim(0, 6000) +
  theme_void() +
  guides(fill = guide_legend(title = NULL))+theme(
    legend.text = element_text(size = 45)
  )

#Correlation with RNA
#read in RNA
library(Matrix)
rna_matrix_raw<-readMM("slide_tags/raw/HumanMelanomaMultiome_rna/matrix.mtx/matrix.mtx")
files <- list.files(path = "slide_tags/raw/HumanMelanomaMultiome_rna/barcodes.tsv/", full.names = TRUE)
files<-substr(files,nchar(files)-17,nchar(files))
files<-substr(files,3,nchar(files))
colnames(rna_matrix_raw)<-files
rna_libsizes_raw<-data.frame(name=colnames(rna_matrix_raw),freq=colSums(rna_matrix_raw))
rna_libsizes1<-rna_libsizes_raw[which(rownames(rna_libsizes_raw)%in%spot_sums$name),]
rna_libsizes1<-rna_libsizes1[match(spot_sums$name,rownames(rna_libsizes1)),]
#####PROCEED WITH spatial-ATAC-RNA-seq & spatial-MUX-seq code to produce correlation between ATAC and RNA libsizes#####
