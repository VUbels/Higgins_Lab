################
# Parameter File
################

#    This file contains all the parameters necessary to run the full semi-automated annotation
#    script provided in the Higgins Lab github repository in addition to being an instruction
#    manual on how to properly run the entire pipeline. The pipeline contains X scripts that the
#    user should run back to back. Each script provides the user with an interim object that can
#    be loaded into the next script. Between each script the output (usually in the form of multiple
#    expression dimplots to determine the labels for the prepared scRNA datset) of that script
#    should be used to inform the user on how to set up the parameters for the subsequent script

###############################################################
# 1: Initial filtering, removing doublets, removing ambient RNA
###############################################################


libraries <- c(
    "Seurat", "SeuratDisk", "knitr", "devtools", "cowplot", "data.table", "dplyr",
    "future", "ggplot2", "ggrepel", "gridExtra", "Matrix", "patchwork",
    "readxl", "scales", "stringr", "tidyverse", "writexl", "DoubletFinder",
    "celda", "harmony", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db",
    "Rmagic", "viridis", "phateR", "ggthemes", "magick",
    "tidyr", "ggrastr", "renv", "extrafont", "ape", "openxlsx"
  )


#########################
# Global Parameters 
#########################

#Project and dataset name
project_name = "Greenleaf23_aLSI"
dataset = "Greenleaf23"

#The directory containing all raw files, currently only works for datasets comprising barcodes/features/matrix files

directory_mtx = "D:/scRNA_datasets/Greenleaf_scRNA_dataset/"

#The directory in which all output will be generated

output_directory_base = "D:/scRNA_output/"
output_directory = paste0("D:/scRNA_output/", project_name, "/")

#The directory containing the R scripts and sample_cmap.rds file

script_directory = "C:/Users/UVict/Desktop/Higgins_Lab/Automated_Annotation/"

#Only use if your know how many cores your CPU has, otherwise comment out.
#Depending on OS you might want to just set nThreads = 1 to avoid potential errors from parallelization

nThreads = 16
plan("multicore", workers = nThreads)

#Harmony parameters, these do not require changing and should not be touched

harmonize <- c() 
covariates <- c() 

##########################
# Parameters for Script_01 
##########################

#Sets the cutoff point for mitochondrial RNA

Mito_Cutoff = 30

#Sets initial minimal RNA features required per cell whilst loading dataset

minFeatures = 200

#Sets the minimal mRNA strand count required from each cell

minCounts = 500

#LSI parameters for initial clustering purposes, these do not require changing

nVarGenes <- 4000
nPCs <- 1:25
resolution <- c(0.2, 0.4, 0.8)
recluster_resolution <- c(0.1, 0.3, 0.6)

###############################################################
# 2: scRNA clustering and annotation
###############################################################

###########################
# Parameters for Scripts_02
###########################

#Set clusters to drop based on discrepancy in disease status dimplot results from Script_01. This value is set to 17 only as an
#example as the dataset used in the example file requires it. If no cluster drop is required set idents_drop to blank 

#idents_drop= c("17")

###########################
# Parameters for Script_03 
###########################

#    In script_03 the data is ready to be divided into various subgroups based on the generated UMAP plots, see Split_By_Disease.pfd
#    DotPlot_Markers.pdf, and Broad_UMAPs and Specific_CellType_UMAP folders generated in script_2. This library can be tailored 
#    to user desires as long as each element is connected in both the subClusterGroups list and the clustNames list.
#    This example file divides all calculated clusters into 5 broad clusters; Lymphoid, Keratinocytes, Fibroblasts, Endothelial,
#    and Other. The dataset used in this example parameter file resulted in 23 individual cell clusters.
#    As such, each of the 23 cell clusters were divided into the 5 broad cluster groups. Ensure that each individual cell cluster
#    retains the same structure, namely r (for RNA, ATAQ sequence data uses the 'a' identification) + celltype + number of
#    specific cluster to be added to the broad cluster. Eaxmple: cluster 8 comprises fibroblasts and is numerically the second
#    cluster to be added to the fibroblast broad cluster group. Therefore "8" = "rFb2". This structure is non-negotiable and must
#    be used to summarize subgroups into broad cluster groups as the script relies on the second and third character in rFb2 for
#    identification.

#Determine which broad cluster groups exist in the dataset and set rules for subclusters. All cluster names in clustNames must be
#attributed to a larger family in subClusterGroups list.

subClusterGroups <- list(
  "Lymphoid" = c("Tc"), 
  "Keratinocytes" = c("Kc"),
  "Fibroblasts" = c("Fb"),
  "Endothelial" = c("Ed", "Le"),
  "Other" = c("Ma", "Me", "Mu", "Im", "Dc", "Oh")
  ) %>% invertList()

clustNames <- list(
    "0" = "rFb1", 
    "1" = "rKc1",
    "2" = "rKc2",
    "3" = "rMu1",
    "4" = "rFb2",
    "5" = "rEd1",
    "6" = "rKc3",
    "7" = "rKc4",
    "8" = "rKc5",
    "9" = "rKc6", 
    "10" = "rEd2", 
    "11" = "rEd3", 
    "12" = "rKc7", 
    "13" = "rTc1", 
    "14" = "rMe1", 
    "15" = "rKc7", 
    "16" = "rKc8", 
    "17" = "rKc9", 
    "18" = "rFb3",
    "19" = "rOh1",
    "20" = "rKc10",
    "21" = "rTc2",
    "22" = "rOh2",
    "23" = "rOh3",
    "24" = "rKc11"
)

###########################
# Parameters for Script_04 
###########################

# In script_3 all broad cluster are assigned and divided into individual cell specific sub-objects that are going to be subjected
# to another round of LSI and clustering through MAGIC to ensure optimal cell clusterisation. 

#Subcluster formation, subgroups contains identical names to the subClusterGroups list

subgroups <- c("Lymphoid", "Keratinocytes",  "Fibroblasts", "Endothelial", "Other")

#subClusterTag is used to assing the x-axis labels in the generated DotPlots for further cluster identification

subClusterGroups_Cluster <- list(
  "Lymphoid" = c("rTc"), 
  "Keratinocytes" = c("rKc"),
  "Fibroblasts" = c("rFb"),
  "Endothelial" = c("rEd"),
  "Other" = c("rOh")
  )

#Give a tag to each subgroup that corresponds to the subgroups vector.

subgroup_tag <- c("rTc", "rKc", "rFb", "rEd", "rOh")

#Subclustering parameters are identical across all subgroups, if more groups are added in subClusterGroups simply 
#copy another group into the paramDict with identical structure as each contained subClusterGroup in the paramDict list.
#If a new group is added set nNeighbors to 40 as this is the default strategy.

#Subclustering parameters:

paramDict <- list(
  "Lymphoid" = list(
    "nVarGenes" = 2000,
    "nPCs" = 1:15,
    "lsiRes" = c(0.1, 0.3),
    "nNeighbors" = 35,
    "minDist" = 0.4,
    "pointSize" = 0.75,
    "harmonize" = c(1,2), # Which iterations should be 'harmonized'
    "covariates" = c('sample')
    ),
  "Keratinocytes" = list(
    "nVarGenes" = 2000,
    "nPCs" = 1:15,
    "lsiRes" = c(0.2, 0.4),
    "nNeighbors" = 35,
    "minDist" = 0.4,
    "pointSize" = 0.75,
    "harmonize" = c(1,2),
    "covariates" = c('sample')
    ),
  "Fibroblasts" = list(
    "nVarGenes" = 2000,
    "nPCs" = 1:15,
    "lsiRes" = c(0.1, 0.3),
    "nNeighbors" = 35,
    "minDist" = 0.4,
    "pointSize" = 1.00,
    "harmonize" = c(1,2),
    "covariates" = c('sample')
    ),
  "Endothelial" = list(
    "nVarGenes" = 2000,
    "nPCs" = 1:15,
    "lsiRes" = c(0.1, 0.3),
    "nNeighbors" = 40,
    "minDist" = 0.35,
    "pointSize" = 1.00,
    "harmonize" = c(1,2),
    "covariates" = c('sample')
    ),
  "Other" = list(
    "nVarGenes" = 2000,
    "nPCs" = 1:15,
    "lsiRes" = c(0.1, 0.3),
    "nNeighbors" = 40,
    "minDist" = 0.35,
    "pointSize" = 1.00,
    "harmonize" = c(1,2),
    "covariates" = c('sample')
    )
  )

###########################
# Parameters for Script_05 
###########################

# Script_04 has generated many UMAP expression plots as well as various DotPlots, using these the specific celltype can be infered.
# In the following Script_05 we are going to override the subtags with the actual cell tags and rejoin all the individual cluster
# seurat objects back into a single useable file for further downstream processing such as velocity, intercellular communication,
# intracellular communication analysis etc.

subCluster_Annotation <- list(
    
    #Endothelial
    "Endothelial" = c("rEd1"),
    "Vascular.Ed" = c("rEd3"), 
    "Stromal" = c("rEd2"), 
   
    #Fibroblasts
    "D.Sheath" = c("rFb1"),
    "D.Papilla" = c("rFb2"),
    "D.Sheath" = c("rFb3"),
    "D.Sheath" = c("rFb4"),

    #Keratinocytes
    "HM.Matrix" = c("rKc1_1"),
    "HM.Proliferating" = c("rKc1_2"),
    "HM.Cortex" = c("rKc1_3"), 
    "HM.Cuticle" =c("rKc1_4"),
    "HM.Matrix" = c("rKc1_5"),
    "ORS.Suprabasal" = c("rKc2_1"),
    "ORS.Companion" = c("rKc2_2"),
    "ORS.Basal" = c("rKc3"),
    "Bulge" = c("rKc5"),
    "IRS.Huxley" = c("rKc6"),
    "IRS.Henle" = c("rKc7"),
    "Differentiating.KCs" = c("rKc9"),

    #Lymphoids
    "Th17.Tc" = c("rTc1_1"),
    "Monocytes.Tc" = c("rTc1_2"),
    "Macrophage" = c("rTc1_3"),
    "Myeloid.Tc" = c("rTc1_4"),
    "CD8NK.Tc" = c("rTc1_5"),
    "Dendritic" = c("rTc3"), 

    #Other
    "Myofibroblast-like" = c("rOh1"),
    "Myofibroblasts" = c("rOh2"),
    "Melanocytes" = c("rOh3"),
    "Mast-cells" = c("rOh4"),
    "Glandular" = c("rOh5"),
    "S.Gland" = c("rOh6")

)  %>% invertList()

