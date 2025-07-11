#Read in packages
library(aricode)
library(BASS)
library(data.table)
library(pheatmap)
library(ggplot2)

#Clean metadata & assign ground truth clusters
spot_sums<-spot_sums[match(metadata$ACAAGCTAAACGTGAT,spot_sums$cleaned_spots),]

residual_matrix<-matrix(NA,50,50)
for (i in seq_len(nrow(metadata))) {
  r <- as.numeric(metadata$row[i])+1
  c <- as.numeric(metadata$col[i])+1
  residual_matrix[r,c] <- i
}

residual_matrix<-residual_matrix[,ncol(residual_matrix):1]
inds_order<-as.vector(t(residual_matrix))
df_mapping<-data.table(spot_names=spot_sums$cleaned_spots[inds_order],cluster=clusters)
df_mapping<-df_mapping[match(metadata$ACAAGCTAAACGTGAT,df_mapping$spot_names),]

#filtering out background
metadata$cluster_value<-df_mapping$cluster
metadata<-metadata[which(metadata$cluster_value!=0),]
spot_sums<-spot_sums[which(spot_sums$cleaned_spots%in%metadata$spot_id),]

#making spot metadata
metadata <- metadata[colnames(peak_matrix), ]
metadata$pixelx<-as.numeric(metadata$pixelx)
metadata$pixely<-as.numeric(metadata$pixely)
metadata_xy<-metadata[,cbind("pixelx","pixely")]
metadata_xy$pixelx<-as.numeric(metadata_xy$pixelx)
metadata_xy$pixely<-as.numeric(metadata_xy$pixely)
xy <- list(metadata_xy)


#preparing data (raw, norm, tfidf)
promoter_matrix_log<-log2(peak_matrix+1) #log transform
promoter_matrix_normalized<-sweep(peak_matrix,2,spot_sums$Freq,"/")
promoter_matrix_normalized_log<-promoter_matrix_normalized*10000
promoter_matrix_normalized_log<-log2(promoter_matrix_normalized_log+1) #log transform normalized peaks

seurat_object <- CreateSeuratObject(counts = promoter_matrix, assay = "peaks")
seurat_object <- RunTFIDF(seurat_object, assay = "peaks",method=1) #default method
tfidf_matrix <- GetAssayData(seurat_object, layer = "data", assay = "peaks")
tfidf_matrix<-as.matrix(tfidf_matrix)

#preparing + running bass on raw counts
cnts <- list(promoter_matrix_log)
BASS <- createBASSObject(cnts, xy, C = 15, R = 5, beta_method = "SW",burnin=2000,k=6,tol=0.0005,step_size=0.15) #may also change R=2 for slidetags and R=10 for k=10 clusters
BASS <- BASS.preprocess(BASS, doLogNormalize = FALSE, geneSelect = "hvgs",doPCA=TRUE,nPC=30)
BASS <- BASS.run(BASS)
BASS <- BASS.postprocess(BASS,adjustLS=FALSE)
metadata$orig<-BASS@results$z[[1]] #save results

#preparing + running bass on norm counts
cnts<-list(promoter_matrix_normalized_log)
BASS <- createBASSObject(cnts, xy, C = 15, R = 5, beta_method = "SW",burnin=2000,k=6,tol=0.0005,step_size=0.15) #may also change R=2 for slidetags and R=10 for k=10 clusters
BASS <- BASS.preprocess(BASS, doLogNormalize = FALSE, geneSelect = "hvgs",doPCA=TRUE,nPC=30)
BASS <- BASS.run(BASS)
BASS <- BASS.postprocess(BASS,adjustLS=FALSE)
metadata$norm<-BASS@results$z[[1]] #save results

#preparing + running bass on TF-IDF counts
cnts<-list(tfidf_matrix)
BASS <- createBASSObject(cnts, xy, C = 15, R = 5, beta_method = "SW",burnin=2000,k=6,tol=0.0005,step_size=0.15) #may also change R=2 for slidetags and R=10 for k=10 clusters
BASS <- BASS.preprocess(BASS, doLogNormalize = FALSE, geneSelect = "hvgs",doPCA=TRUE,nPC=30)
BASS <- BASS.run(BASS)
BASS <- BASS.postprocess(BASS,adjustLS=FALSE)
metadata$tfidf<-BASS@results$z[[1]] #save results

#make spatial domain heatmap
spatial_domains<-matrix(NA,nrow=50,ncol=50)
for (i in seq_len(nrow(metadata))) {
  r <- as.numeric(metadata$x[i])+1
  c <- as.numeric(metadata$y[i])+1
  spatial_domains[r,c] <- metadata$orig[i] #may also replace with norm or tfidf
}
spatial_domains<-(spatial_domains)[,ncol(spatial_domains):1]
pheatmap(spatial_domains,cluster_rows=FALSE,cluster_cols=FALSE,color=c("red","blue","#32CD32","purple","orange"),main=paste0("ARI: ",round(ARI(metadata$orig, metadata$cluster_value),4)),fontsize=32)


#specificially for slidetags
ggplot(metadata, aes(x = as.numeric(X), y = as.numeric(Y), color = as.factor(orig))) + #may also replace orig with norm and tfidf
  geom_point(size = 1.5) +
  coord_fixed() +
  scale_color_manual(values = c("red", "blue")) +
  theme_minimal() +
  labs(
    title = paste0("ARI: ", round(ARI(metadata$cluster_value, metadata$orig), 4)), 
    x = "X",
    y = "Y") +
  xlim(0, 6000) +
  ylim(0, 6000) +
  theme(
    plot.title = element_text(size = 40, face = "bold",hjust=0.5),
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 36),
    axis.text.x = element_text(angle = 0),
    legend.position = "none",
    panel.grid = element_blank())
