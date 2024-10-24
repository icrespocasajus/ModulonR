
# ModulonR

ModulonR identifies “modulons”—clusters of TFs exhibiting coordinated
activity across different cell states, which are further characterized
based on their ability to discriminate a particular cell state from the
rest. Next, by leveraging exploration of the inferred gene regulatory
network (GRN), ModulonR predicts the set of TFs that might be essential
for the commitment/stability of a particular cell state. ModulonR takes,
as input, the transcription factor (TF) activity matrix with single-cell
resolution and a GRN. The Analysis consists of the following three
steps: 1) Modulon identification through hierarchical clustering of TFs
by their activity across cell states. The resulting clustering tree is
iteratively pruned at different heights (k clustering parameter), and
the resulting cluster-derived signatures are calculated for each cell.
Subsequently, a multivariate discriminant analysis (OPLSDA) for each
pre-known T cell state is performed, considering all the cluster
signatures together. A global discriminant score (GDS) is then
calculated for each k to summarize the discriminant power for all cells
into a single numerical value. The k with the best GDS determines the
modulons; 2) modulon selection. The modulon signatures are calculated
for each cell and used as input for a cell state discriminant analysis.
The top discriminant modulon of a given state is selected; 3) modulon
perturbation analysis. Finding combinations of an arbitrary number of KO
target TFs for optimal disruption of the modulon activity and the T cell
states that depend on it. To this end, both the strength of the TF-gene
regulations (GRN weights) and the discriminancy (assumed as relative
importance) of the regulated gene for specific T cell states are
considered. Combinations of modulon elements and are ranked using the
weighted coverage score (WCS), where high values represent combinations
of TFs that strongly regulate many genes with high discriminancy within
the cell state to be disrupted.

## Installation

You can install the development version of ModulonR from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("icrespocasajus/ModulonR")
```

## Example

This is a basic example which shows you how to use ModulonR on a mouse
TILs reference dataset by Andreatta, M. et al. \[1\], including 5 CD8+ T
cell functional states.

### Load the package:

``` r
library(ModulonR)
```

### Load other required packages:

``` r
library(AUCell)
library(ggplot2)
library(pheatmap)
library(dplyr)
```

The UMAP projection of the input mouse TILs reference dataset showing
the 5 functional states in different colors:
<center>
<img src="./man/figures/README-UMAP_TILs_Dataset.png" width="75%" />
</center>
The violin plots show the expression of important marker genes across
functional states:
<center>
<img src="./man/figures/README-Markers_TILs_Dataset.png" width="25%" />
</center>

A regulon activity matrix and the gene regulatory network (GRN) were
derived from the scRNASeq data TILs reference dataset using SCENIC
\[2\].

### Load input data:

``` r
# Regulon activity matrix
regulon.activity.matrix = readRDS(file="./data-raw/Regulon_Activity_Matrix_TILs.Rds")

# Cell annotation
annotation = readRDS(file="./data-raw/Annotation_TILs.Rds")

# Regulon signature
regulons = readRDS(file="./data-raw/Regulons.TILS.w.Tox.Rds")

# Gene expression matrix
gene.expression.matrix = readRDS(file="./data-raw/Gene_Expression_Matrix_TILs.Rds")

# TF-target weights in Tex
Tex.GENIE3.links = readRDS(file="./data-raw/1.4_GENIE3_linkList_CD8_Tex.w.Tox.Rds")
```

### STEP 1: Modulon Identification

### Identify modulons:

``` r
ModulonIdent.results = ModulonIdent(
  data=regulon.activity.matrix,
  annotation=annotation,
  BackgroundClasses=NULL,
  QueryClasses=NULL,
  k.range = c(2:10))
```

### Explore results: plot the Global Discriminant Score (GDS)

``` r
GDS = ModulonIdent.results[['GDS']]

# Identify the best clustering resolution
Best.k = GDS[which.max(GDS$avg_max_weight),'NumCluster']

# Plot the GDS for all clustering resolutions (k)
p <- ggplot2::ggplot(
  GDS,aes(x = NumCluster, y = avg_max_weight, fill = avg_max_weight)) +
  coord_cartesian(ylim=c(0.70,1)) + 
  geom_col(color='blue') +
  scale_x_continuous(breaks = seq(1, 10, 5)) +
  theme_light() +
  theme(legend.position = "right") +
  ggtitle("Cell State Discriminancy") +
  xlab("Number of clusters (k)") +
  ylab("Global Discriminant Score (GDS)") +
  labs(fill = "GDS") +

# Add vertical line at x = 5
  geom_vline(xintercept = Best.k, linetype = "dashed", color = "red", size = 1) +
  
  # Add label for the vertical line
  annotate("text", x = Best.k, y = 0.95, label = paste0("Optimal k = ",Best.k), color = "red", vjust = -0.5)

plot(p+theme(legend.position = "none"))
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

### Explore results: visualize modulons

``` r
# Plot a heatmap
expMat.AUC.TILs.aggregated = aggregate(t(regulon.activity.matrix),by=list(annotation),FUN= mean)

rownames(expMat.AUC.TILs.aggregated)=expMat.AUC.TILs.aggregated[,1]
expMat.AUC.TILs.aggregated=expMat.AUC.TILs.aggregated[,-1]
expMat.AUC.TILs.aggregated=t(expMat.AUC.TILs.aggregated)

samples=c('CD8_NaiveLike','CD8_EarlyActiv','CD8_EffectorMemory','CD8_Tpex','CD8_Tex')
phm.input = expMat.AUC.TILs.aggregated[,samples]

ann_colors = list(
  State = c(
    
    `CD8_Tex`=alpha("#A1D9C7",alpha = 1),
    `CD8_Tpex`=alpha("#B3DD6F",alpha = 1),
    `CD8_EffectorMemory`=alpha("#E19175",alpha = 1),
    `CD8_EarlyActiv`=alpha("#B6A8D7",alpha = 1),
    `CD8_NaiveLike`=alpha("#C060D4",alpha = 1)
    
  ),
  #Condition = c('Tumor infiltrating' = 'black'),
  Modulon = c(
    "1"="#fe9288",
    "2"="#05c7ff",
    "3"="#dab702",
    "4"="#fe82f8",
    "5"="#00da8d",
    "6"="#a7aeff",
    "7"="#03ddc6")
  )

annotation.c = data.frame(State = colnames(phm.input))
rownames(annotation.c)=colnames(phm.input)

# Create a data frame with the modulon annotation
modulons = ModulonIdent.results[['Modulons']]
names(modulons) = as.character(strsplit2(names(modulons),'[.]')[,2])


modulon_annotation <- do.call(rbind, lapply(names(modulons), function(name) {
  data.frame(Modulon = name, TF = modulons[[name]])
}))
rownames(modulon_annotation) <- modulon_annotation$TF
modulon_annotation$TF <- NULL
modulon_annotation = modulon_annotation[rownames(phm.input),,drop=F]
annotation.r = modulon_annotation


phm.output.AUC.TILs=pheatmap::pheatmap(phm.input,
                             main = "Regulon activity",
                             cluster_rows = T,
                             cluster_cols = F,
                             fontsize_row = 0.5,
                             #gaps_col = 4,
                             annotation_col = annotation.c,
                             annotation_colors = ann_colors ,
                             annotation_row  = annotation.r,
                             annotation_names_row=T,
                             annotation_names_col=T,
                             cutree_rows = Best.k,
                             show_rownames = F,
                             cellwidth = 12,
                             cellheight = 0.5,
                             fontsize_col = 12,
                             scale='row'
)
```

<img src="man/figures/README-unnamed-chunk-11-1.png" width="100%" />

## STEP 2: Modulon Selection

### Modulon discriminant analysis of the CD8_Tex target state:

``` r
TargetState = c('CD8_Tex')

ModulonSelect.results = ModulonSelect(
  data=regulon.activity.matrix,
  modulons=modulons,
  annotation=annotation,
  BackgroundClasses=NULL,
  TargetState = TargetState)
#> [1] "Running discriminant analysis for CD8_Tex"
ModulonSelect.results[["Selected_Modulon"]][[TargetState]]
#> [1] "3"
modulons[[ModulonSelect.results[["Selected_Modulon"]][[TargetState]]]]
#>  [1] "Ybx1"    "Crem"    "Nr3c1"   "Rora"    "Fli1"    "Gata3"   "Cebpd"  
#>  [8] "Foxn2"   "Tbx21"   "Hmgb2"   "Zfp433"  "Sap30"   "E4f1"    "Eomes"  
#> [15] "Runx3"   "Runx2"   "Gtf2f1"  "Snai3"   "Srebf2"  "Maf"     "Rarb"   
#> [22] "Rara"    "Nono"    "Zfp821"  "Srf"     "Zfp959"  "Zfp961"  "Stat5a" 
#> [29] "Tbp"     "Uqcrb"   "Trp73"   "Tead1"   "Gfi1"    "E2f6"    "Mxd3"   
#> [36] "Nfe2"    "Vezf1"   "Zmiz1"   "Batf"    "Gabpb1"  "Smarca4" "Borcs8" 
#> [43] "Psmd12"  "Nup133"  "Foxc1"   "Gm14295"
```

### Explore results: plot modulon discriminancy for the CD8_Tex target state

``` r
# Plot the discriminant score of all modulons for a given T cell state

input.tmp = ModulonSelect.results[["Modulon_DA"]][[TargetState]]
input.tmp$label = input.tmp$feature
input.tmp = input.tmp[order(input.tmp$weightStarMN,decreasing = T),]

input.tmp$order = factor(c(1:nrow(input.tmp)))

custom_colors  = c(
  "1" ="#fe9288",
  "2" ="#05c7ff",
  "3" ="#dab702",
  "4" ="#fe82f8",
  "5" ="#00da8d",
  "6" ="#a7aeff",
  "7" = "#03ddc6")

p <- ggplot2::ggplot(
  input.tmp,aes(x = order, y = weightStarMN, fill = label)) +
  coord_cartesian(ylim=c(-1,1)) + 
  geom_col(color='black') +
  theme_light() +
 theme(legend.position = "none",
        axis.text.x = element_blank(),  # Remove x-axis tick labels
        axis.ticks.x = element_blank()) +  # Optionally remove x-axis ticks
  ggtitle(paste0(TargetState," State Modulon Discriminancy")) +
  xlab("Modulons") +
  ylab("Discriminant Score") +
  labs(fill = "Modulon") +
  # Add specific colors
  scale_fill_manual(values = custom_colors) +

  # Add labels on top or bottom of bars based on their value
  geom_text(aes(label = label, 
                vjust = ifelse(weightStarMN >= 0, -0.5, 1.5)),
            color = "black")

plot(p)
```

<img src="man/figures/README-unnamed-chunk-14-1.png" width="100%" />

### STEP 3: Modulon Perturbation

### Rank combinations by expected impact on the target state:

``` r
ModulonPert.results = ModulonPert(
  Regulons = regulons,
  Modulons = modulons,
  ExpMat = gene.expression.matrix,
  annotation= annotation,
  BackgroundClasses = NULL,
  TargetState='CD8_Tex',
  TargetModulon=ModulonSelect.results[["Selected_Modulon"]][["CD8_Tex"]],
  CombSize=3,
  Weights=Tex.GENIE3.links
  )
```

### Explore results: see top 5 combinations of 3 KOs

``` r
head(ModulonPert.results[["Combinations"]])[c(1:5),c(2:5)]
#>                  Element_1 Element_2 Element_3      WCS
#> Combination_1262      Crem     Tbx21     Zmiz1 1.554944
#> Combination_1263      Crem     Tbx21      Batf 1.501018
#> Combination_1079      Crem      Fli1     Tbx21 1.472272
#> Combination_1240      Crem     Tbx21     Runx2 1.440801
#> Combination_1243      Crem     Tbx21    Srebf2 1.437606
```

## References

1.  Andreatta, M. et al. “Interpretation of T cell states from
    single-cell transcriptomics data using reference atlases.” Nature
    communications vol. 12,1 2965. 20 May. 2021,
    <doi:10.1038/s41467-021-23324-4>.

2.  Aibar, S. et al. “SCENIC: single-cell regulatory network inference
    and clustering.” Nature methods 14.11 (2017): 1083-1086.

# Authors

<img src="./man/figures/README-isaaccrespo.jpg" width="10%" />

Isaac Crespo, phD  
Senior Computational Scientist  
CHUV/UNIL \| George Coukos group  
Ludwig Institute for Cancer Research \| Lausanne Branch  
AGORA, Bugnon 25A, 1005 Lausanne (Switzerland), 4th floor, Room 026  
<isaaccrespo@hotmail.com>

<img src="./man/figures/README-anarodriguez.png" width="10%" />

Ana Rodriguez Sanchez-Archidona, phD  
Bioinformatician  
CHUV/UNIL \| George Coukos group  
Ludwig Institute for Cancer Research \| Lausanne Branch  
AGORA, Bugnon 25A, 1005 Lausanne (Switzerland), 4th floor, Room 190  
<arsanchezarchidona@gmail.com>

<img src="./man/figures/README-remypetremand.jpg" width="10%" />

Remy Petremand  
PhD student  
CHUV/UNIL \| Alexandre Harari group  
Ludwig Institute for Cancer Research \| Lausanne Branch  
AGORA, Bugnon 25A, 1005 Lausanne (Switzerland), 4th floor, Room 026  
<remy.petremand@chuv.ch>
