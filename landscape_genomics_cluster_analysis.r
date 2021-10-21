###########################################
### Landscape genomics cluster analysis ###
###########################################

# ### pre-processing 2018 collection genotype data mapped on the Lope Byrne V1 reference genome
# dat = read.csv("CLETH_MERGED_ALL_MAF0.001_DEPTH50_ALLELEFREQ.csv", F)
# dat_names = read.csv("CLETH_MERGED_ALL_bam.list", F)
# X = dat[,4:ncol(dat)]
# vec_names = gsub(".bam", "", matrix(unlist(strsplit(dat_names$V1, "/")), ncol=6, byrow=TRUE)[,6])
# vec_idx_bool = !grepl("Pool", vec_names)
# X = X[, vec_idx_bool]
# out = cbind(dat[,1:3], X)
# colnames(out) = c("SCAFFOLD", "POSITION", "REFERENCE_ALLELE", vec_names[vec_idx_bool])
# write.csv(out, file="Lolium_2018_collection-ByrneV1_reference-ALLELEFREQ.csv", quote=FALSE, row.names=FALSE)

### Load allele frequencies
df_freq = read.csv("Lolium_2018_collection-ByrneV1_reference-ALLELEFREQ.csv")

### Load population ID
df_ID = read.csv("Lolium_collection.csv")

### Merge the alele frequency dataframe with the population ID dataframe
df_loc = data.frame(POPULATION=gsub("_", "", df_ID$Accession_Name), 
                    LOCATION=df_ID$Location,
                    NORTH=df_ID$Coordinates_N,
                    EAST=df_ID$Coordinates_E,
                    CROP=df_ID$Crop)
df_geno = data.frame(colnames(df_freq)[4:ncol(df_freq)], t(df_freq[, 4:ncol(df_freq)]))
colnames(df_geno) = c("POPULATION", paste(df_freq[, 1], df_freq[, 2], df_freq[, 3], sep="_"))
df_main = merge(df_loc, df_geno, by="POPULATION")

### Remove outlier population which is most likely *Lolium loliaceum* based on plant morphology
df_main = df_main[!grepl("ACC21", df_main$POPULATION), ]

### Extract the allele frequency matrix
mat_freq = df_main[, 6:ncol(df_main)]

### Compute the pairwise Euclidean distance between population based on allele frequencies
mat_distance = dist(mat_freq, method="euclidean")

### Cluster using average distances and plot the dendrogram
obj_cluster = hclust(mat_distance, method="average")
plot(obj_cluster, hang=-1)

### Setup the matrix of allele frequency-based pairwise distances
mat_dG = as.matrix(mat_distance)

### Setup the matrix of spatial pairwise distances
fun_haversine_distance = function(vec_degree_coordinates_1, vec_degree_coordinates_2, radius_km=6371){
    # ############################
    # ### TEST
    # vec_degree_coordinates_1 = c(-35.004448, 147.464267)
    # vec_degree_coordinates_2 = c(-38.584017, 146.013588)
    # radius_km = 6371
    # geosphere::distHaversine(rev(vec_degree_coordinates_1), rev(vec_degree_coordinates_2))
    # ############################
    lat_1 = vec_degree_coordinates_1[1] *pi/180
    lon_1 = vec_degree_coordinates_1[2] *pi/180
    lat_2 = vec_degree_coordinates_2[1] *pi/180
    lon_2 = vec_degree_coordinates_2[2] *pi/180
    a = sin( (lat_2 - lat_1)/2 )^2
    b = cos(lat_1) * cos(lat_2) * sin( (lon_2 - lon_1)/2 )^2
    d = 2 * radius_km * asin( sqrt(a + b) )
    return(d)
}

mat_dL = matrix(NA, nrow=nrow(df_main), ncol=nrow(df_main))
for (i in 1:nrow(df_main)){
    vec_degree_coordinates_1 = c(df_main$NORTH[i], df_main$EAST[i])
    for (j in 1:nrow(df_main)){
        vec_degree_coordinates_2 = c(df_main$NORTH[j], df_main$EAST[j])
        mat_dL[i, j] = fun_haversine_distance(vec_degree_coordinates_1, vec_degree_coordinates_2)
    }
}


### regress genetic distance against geographic distance
vec_genetic_distances = mat_dG[lower.tri(mat_dG)]
vec_geographic_distances = mat_dL[lower.tri(mat_dL)]

svg("Prelim_distance_regression.svg", width=8, height=10)
layout(matrix(c(1,2,3,3,3,3), ncol=2, byrow=TRUE))
image(mat_dG, main="Genetic distance\n(allele frequency-based Euclidean distances)")
image(scale(mat_dL), main="Geographical distance\n(Haversine distances)")
plot(vec_genetic_distances ~ vec_geographic_distances,
    pch=19, col=rgb(0.2, 0.2, 0.2, alpha=0.5),
    xlab="Geographical distance (km)", ylab="Genetic distance (|allele frequency|)")
grid()
mod = lm(vec_genetic_distances ~ vec_geographic_distances)
R2adj = round(summary(mod)$adj.r.squared * 100, 2)
slope = round(mod$coef[2], 4)
abline(mod, col="red")
legend("topleft", legend=c(paste0("R2 adjusted = ", R2adj, "%"), paste0("Slope = ", slope)))
dev.off()
