library(readxl)
library(dplyr)
library(tibble)
library(vegan)
library(ggplot2)
# Read file
path <-"C:/Users/Dell/OneDrive/Desktop/Pea Research Pruess Lab/Drought microbiome/Bidya 11-4-25 Core files/raw seq data"
setwd(path)
otu_raw <- read_excel("v3-v499.xlsx")
# Read metadata
meta <- read_excel("Bidya metadata 10-21-25.xlsx")

# abundance table is already numeric-only
otu <- as.data.frame(otu_raw)
otu[] <- lapply(otu, as.numeric)
otu[is.na(otu)] <- 0

# transpose so rows = samples, columns = taxa
otu_t <- t(otu)

# alpha diversity
alpha_div <- data.frame(
  Sample = rownames(otu_t),
  Shannon = diversity(otu_t, index = "shannon"),
  Simpson = diversity(otu_t, index = "simpson"),
  Richness = specnumber(otu_t)
)

# clean metadata
meta <- as.data.frame(meta)
colnames(meta) <- trimws(colnames(meta))
meta$OTUID <- trimws(meta$OTUID)
alpha_div$Sample <- trimws(alpha_div$Sample)

# merge using Sample from alpha_div and OTUID from metadata
alpha_meta <- merge(alpha_div, meta, by.x = "Sample", by.y = "OTUID")

# check
head(alpha_meta)

#####write treatments in order
alpha_meta$Treatment <- factor(alpha_meta$Treatment,
                               levels = c(
                                 "ACyt",
                                 "ACyt Drought",
                                 "AEx",
                                 "AEx Drought",
                                 "A",
                                 "A Drought",
                                 "Cyt",
                                 "Cyt Drought",
                                 "Ex",
                                 "Ex Drought",
                                 "NB",
                                 "NB Drought"
                               )
)
###create shannon boxplot
ggplot(alpha_meta, aes(x = Treatment, y = Shannon)) +
  geom_boxplot(
    fill = "white",        # remove color inside box
    color = "black",
    linewidth= 0.6, # keep box outline
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 2,
    color = "black"
  ) +
  scale_x_discrete(drop = FALSE) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold",size=10, color="black"),
    axis.text.y = element_text(face = "bold",size=10,color="black"),
    axis.title.y = element_text(face = "bold"),   # bold y label
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    legend.position = "none"
  ) +
  labs(x= NULL,
    title = "Shannon Diversity Across Treatments",
    y = "Shannon Index"
  )
#####Simpson and Richness plots:
  
ggplot(alpha_meta, aes(x = Treatment, y = Simpson)) +
  geom_boxplot(
    fill = "white",      
    color = "black",
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 2,
    color = "black"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "Simpson Diversity Across Treatments",
    x = "Treatment",
    y = "Simpson Index"
  )
######For richness
ggplot(alpha_meta, aes(x = Treatment, y = Richness)) +
  geom_boxplot(
    fill = "white",      
    color = "black",
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.12,
    size = 2,
    color = "black"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "Richness Across Treatments",
    x = "Treatment",
    y = "Observed Richness"
  )

###### For statistics across all treatments:
  
kruskal.test(Shannon ~ Treatment, data = alpha_meta)
pairwise.wilcox.test(alpha_meta$Richness, alpha_meta$Treatment, p.adjust.method = "BH")

#############################
###PCA ONLY
library(factoextra)
library(ggplot2)
library(ggrepel)

fviz_pca_ind(
  pca_res,
  repel = TRUE,          
  col.ind = "blue",      # treatment points
  pointshape = 16,
  pointsize = 3,
  labelsize = 4
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(face = "bold", size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "none"
  )

###biplot
install.packages("readxl")
install.packages("factoextra")
install.packages("ggplot2")

library(readxl)
library(factoextra)
library(ggplot2)
RA <- read_excel("Relativeabundancecosine.xlsx")
mat <- as.data.frame(RA[ , -1])
rownames(mat) <- RA[[1]]
mat_tr <- t(mat)
treatment_order <- c(
  "ACyt", "ACyt Drought",
  "AEx", "AEx Drought",
  "A", "A Drought",
  "Cytidine", "Cytidine Drought",
  "Ex", "Ex Drought",
  "NB", "NB Drought"
)
pca_res <- prcomp(mat_tr, scale. = TRUE)
fviz_pca_biplot(
  pca_res,
  repel = TRUE,
  col.var = "red",      # phylum arrows
  col.ind = "blue",     # treatment points
  arrowsize = 0.8,
  label = "all"
)

###Show treatment names clearly
fviz_pca_biplot(
  pca_res,
  repel = TRUE,
  geom.ind = c("point", "text"),
  col.ind = "blue",
  col.var = "red",
  pointsize = 4,
  arrowsize = 0.8
) +
  theme_bw()


######################
library(factoextra)
library(ggplot2)
library(ggrepel)

fviz_pca_biplot(
  pca_res,
  repel = TRUE,                 # avoids overlapping labels
  col.ind = "blue",             # sample/treatment points
  col.var = "red",              # taxa/phylum arrows
  pointshape = 16,
  pointsize = 3,
  labelsize = 4,
  arrowsize = 0.8
) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(face = "bold", size = 10),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    legend.position = "none"
  )

##################################Cosinus Similarity
# Install if needed:
# install.packages(c("readxl","vegan","pheatmap"))

library(readxl)

path <-"C:/Users/Dell/OneDrive/Desktop/Pea Research Pruess Lab/Drought microbiome/Bidya 11-4-25 Core files/raw seq data"
setwd(path)
df<- read_excel("Relativeabundancecosine.xlsx")
install.packages("lsa")
install.packages("pheatmap")
library(lsa)
library(pheatmap)
library(readxl)
library(pheatmap)


####################
###RA cosine
library(readxl)
library(FactoMineR)
library(pheatmap)
cosRA <- read_excel("Relativeabundancecosine.xlsx")
mat <- as.data.frame(cosRA[,-1])
rownames(mat) <- cosRA[[1]]
ca_res <- CA(mat, graph = FALSE)
phylum_co <- ca_res$row$coord[, 1:2]
treatment_co <- ca_res$col$coord[, 1:2]
cos_sim <- function(a, b) {
  sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))
}
cosine_matrix <- matrix(NA,
                        nrow = nrow(phylum_co),
                        ncol = nrow(treatment_co))

rownames(cosine_matrix) <- rownames(phylum_co)
colnames(cosine_matrix) <- rownames(treatment_co)

for (i in 1:nrow(phylum_co)) {
  for (j in 1:nrow(treatment_co)) {
    cosine_matrix[i, j] <- cos_sim(phylum_co[i, ], treatment_co[j, ])
  }
}
treatment_order <- c(
  "ACyt","ACyt Drought",
  "AEx","AEx Drought",
  "A","A Drought",
  "Cyt","Cyt Drought",
  "Ex","Ex Drought",
  "NB","NB Drought"
)
cosine_matrix <- cosine_matrix[, treatment_order]
cosine_matrix
pheatmap(cosine_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         axis.title.x = element_text(face = "bold"),
         main = "Phylum vs Treatment Cosine Similarity")
############


