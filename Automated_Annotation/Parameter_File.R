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

#########################
# Global Parameters 
#########################

#Project and dataset name
project_name = "Shim"
dataset = "Shim"

#The directory containing all raw files, currently only works for datasets comprising barcodes/features/matrix files

directory_mtx = "C:/Users/VictorUPC/Documents/CellChat/Shim_sRNA_dataset/"

#The directory in which all output will be generated

output_directory_base = "C:/Users/VictorUPC/Documents/scRNA_output/"
output_directory = paste0("C:/Users/VictorUPC/Documents/scRNA_output/", project_name, "/")

#The directory containing the R scripts and sample_cmap.rds file

script_directory = "C:/Users/VictorUPC/Desktop/Higgins_Lab-main/Automated_Annotation/"

#Sets tags for various conditions that might be available within each dataset. For example, the Green example dataset contains data
#from alopecia areata (AA), punch biopsies (C_PB), and peripheral surgical tissue (C_SD) patient samples.
 
tags = c('Hair_Follicle')


#Only use if your know how many cores your CPU has, otherwise comment out

nThreads = 8
plan("multicore", workers = nThreads)

#Harmony parameters, these do not require changing and should not be touched

harmonize <- c() 
covariates <- c() 

##########################
# Parameters for Script_01 
##########################

#Sets the cutoff point for mitochondrial RNA

Mito_Cutoff = 20

#Sets initial minimal RNA features required per cell whilst loading dataset

minFeatures = 200

#Sets the minimal RNA features required per cell

NFeaturesRNA = 1000

#Sets the minimal mRNA strand count required from each cell

minCounts = 1000

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
    "2" = "rTc1",
    "3" = "rKc2",
    "4" = "rTc2",
    "5" = "rEd1",
    "6" = "rMu1",
    "7" = "rKc3",
    "8" = "rFb2",
    "9" = "rKc4", 
    "10" = "rKc5", 
    "11" = "rIm1", 
    "12" = "rTc3", 
    "13" = "rIm2", 
    "14" = "rMu2", 
    "15" = "rEd2", 
    "16" = "rMe1", 
    "17" = "rLe1", 
    "18" = "rDc1",
    "19" = "rMa1",
    "20" = "rIm3",
    "21" = "rFb3",
    "22" = "rFb4"
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
    "lsiRes" = c(0.1, 0.3),
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
    "Vas.Endo" = c("rEd1"),
    "Vas.Endo" = c("rEd2"), 
    "Vas.Endo" = c("rEd3"), 
    "Lymph.Endo" = c("rEd4"), 
    "Angiogenic" = c("rEd5"),

    #Fibroblasts
    "D.Sheath" = c("rFb2"),
    "D.Papilla" = c("rFb6"),
    "D.Fibroblasts" = c("rFb1"),
    "D.Fibroblasts" = c("rFb3"),
    "D.Fibroblasts" = c("rFb4"),
    "D.Fibroblasts" = c("rFb5"),

    #Keratinocytes
    "Basal.Kc" = c("rKc1"),
    "Spinous.Kc" =c("rKc2"),
    "Spinous.Kc" =c("rKc3"), 
    "Spinous.Kc" =c("rKc4"),
    "Proliferating.Kc" = c("rKc5"),
    "Inf.Segment" = c("rKc6"),
    "Outer Root Sheath" = c("rKc7"),
    "Isthmus" = c("rKc8"),
    "Sebaceous" = c("rKc9"),
    "Eccrine" = c("rKc10"),

    #Lymphoids
    "CD4.Tc" = c("rTc1"), 
    "CD4.Tc" = c("rTc2"), 
    "CD8.Tc" = c("rTc3"), 
    "Treg" = c("rTc4"), 
    "NK" = c("rTc5"),
    "Cyc.Tc" = c("rTc6"),

    #Other
    "D.Sheath" = c("rOh11"),
    "M2.Mac" = c("rOh8"),
    "TREM2.Mac" = c("rOh5"),
    "CLEC9.DC" = c("rOh7"),
    "cDC2" = c("rOh1"),
    "Bulge" = c("rOh2"),
    "Dermal.Mac" = c("rOh3"),
    "Hair Matrix" = c("rOh6"),
    "Bulge_2" = c("rOh4"),
    "M1.Mac" = c("rOh9"),
    "IRF4.cells" = c("rOh10")


)  %>% invertList()

