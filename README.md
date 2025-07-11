# spATAC-Library-Size
Library size has been typically regarded as a technical confounder in spATAC-seq. However, analyses indicated that library size is correlated with biology and that traditional normalization methods do not improve/may worsen the performance of downstream analyses tasks like spatial domain detection, and change biological trends in differential promoter activity.


# Description of each File 
## image_clustering.R 
1. Reads in the image and processes it via a smoothing function before clustering based on grayscale intensity.
2. Downsample image size to correspond to spatial library size.
3. Compares library sizes between clusters and also calculates the library size correlation between ATAC and RNA library size.
4. Includes code to produce Figure 1.

## promoter_matrix_rna_cor.R
1. Creates a list of promoter regions and generates a peak by spot matrix.
2. Aligns promoter regions with RNA gene expression matrix. and calculates row-wise correlation.
3. Creates a divide by library size normalized and TF-IDF corrected matrix for the peak by spot table and a normalized matrix for the RNA matrix.
4. Calculate row-wise correlation between promoters and RNA gene expression.
5. Includes code to produce Figure 2.

## spatial_domain_detection.R
1. Read in the metadata and assign clusters produced in image_clustering.R file.
2. Prepare original, normalized, and TF-IDF matrix.
3. Run BASS on the datasets and store the results.
4. Calculate Adjusted Rand Index (ARI) and visualize identified spatial domains.
5. Includes code to produce Figure 3.

## differential_promoters.R
1. Assign clusters to compare differential promoter activity (i.e. GCL vs non-GCL in Human hippocampus and the 2 tumor clusters in the Slide-tags dataset).
2. Filters out low-read promoters.
3. Determine which promoters are differential promoters using a Wilcoxon test and false discovery rate. Determine which differential promoters are shared between raw vs. normalized conditions.
4. Calculate log2 fold change for each peak between the different clusters. Determine % of differential promoters raw vs. normalized conditions whose trends shift (i.e. significantly upregulated to significantly downregulated). 
5. Includes code to produce Figure 4.
