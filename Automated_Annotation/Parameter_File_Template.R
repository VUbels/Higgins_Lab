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

#Minimal libraries required to run the entire annotation pipeline for human data
libraries <- c(
    "Seurat", "knitr", "devtools", "cowplot", "data.table", "dplyr",
    "future", "ggplot2", "ggrepel", "gridExtra", "Matrix", "patchwork",
    "readxl", "scales", "stringr", "tidyverse", "writexl", "DoubletFinder",
    "celda", "harmony", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db",
    "Rmagic", "viridis", "phateR", "ggthemes", "magick", "reticulate",
    "tidyr", "ggrastr", "renv", "extrafont", "ape", "openxlsx"
  )

#########################
# Global Parameters 
#########################

#Project and dataset name
project_name = "ACG_2024"
dataset = "ACG_2024"
directory_mtx = "D:/scRNA_ACG_2024_Input"

#The directory containing all raw files, currently only works for datasets comprising barcodes/features/matrix files

#directory_mtx = paste0(directory, "/", project_name)

#The directory in which all output will be generated

output_directory_base = "D:/scRNA_ACG_2024_Output"
output_directory = paste0(output_directory_base, "/", project_name,"/")

#The directory containing the R scripts and sample_cmap.rds file

script_directory = "C:/Users/VictorUPC/Desktop/Higgins_Lab-main/Automated_Annotation/"

# Sets tags for various conditions that might be available within each dataset. The conditions must be specified according to:
# pattern <- paste0("_(", paste(tags, collapse="|"), ")(?=\\d)")
# As an example, GSM6532925_C_PB3, GSM6532920_AA4 all fill this criteria such that the tag for GSM6532920_AA4 is 'AA' etc.
# If no pattern is found or dictated disease state will be reverted to orig ident.

tags = c('A', 'C', 'T')


#Set CPU availability, on windows simply open task manager > performance > check Logical processors

nThreads = 8
plan("multicore", workers = nThreads)

#Harmony parameters, set harmonize to 1,2 if batch correction is required. IE: harmonize <- c(1,2) 

harmonize <- c() 
covariates <- c() 

##########################
# Parameters for Script_01 
##########################
#Indicate whether or not to apply batch correction
Batch_Correction = FALSE

#Sets the cutoff point for mitochondrial RNA
Mito_Cutoff = 50

#Sets initial minimal RNA features required per cell whilst loading dataset
minFeatures = 150

#Sets the minimal RNA features required per cell
NFeaturesRNA = 150

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
#example as the dataset used in the example file requires it. If no cluster drop is comment out idents_drop with #

#idents_drop <- c("17")

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
  "Keratinocytes" = c("Kc", "Irs", "Ors"),
  "Fibroblasts" = c("Fb"),
  "Endothelial" = c("Ed", "Le"),
  "Other" = c("Ma", "Me", "Mu", "Im", "Dc", "Oh")
  ) %>% invertList()

clustNames <- list(
    "0" = "rFb1", 
    "1" = "rMa1",
    "2" = "rKc1",
    "3" = "rMu1",
    "4" = "rFb2",
    "5" = "rEd1",
    "6" = "rKc2",
    "7" = "rIrs1",
    "8" = "rMa2",
    "9" = "rKc3", 
    "10" = "rEd2", 
    "11" = "rEd3", 
    "12" = "rKc4", 
    "13" = "rTc1", 
    "14" = "rMe1", 
    "15" = "rIrs2", 
    "16" = "rKc5", 
    "17" = "rMa3", 
    "18" = "rFb3",
    "19" = "rIm1",
    "20" = "rKc6",
    "21" = "rTc2",
    "22" = "rOh1",
    "23" = "rOh2",
    "24" = "rTc3"
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
    "BTNL9_Vas.End" = c("rEd1"),
    "VWF_Vas.End" = c("rEd2"), 
    "Angionegic.End" = c("rEd3"), 
    "Lymph.End" = c("rEd4"),

    #Fibroblasts
    "D.Sheath" = c("rFb2"),
    "D.Papilla_2" = c("rFb6"),
    "D.Sheath" = c("rFb1"),
    "D.Papilla_1" = c("rFb3"),
    "D.Papilla_2" = c("rFb4"),
    "D.Sheath" = c("rFb5"),

    #Keratinocytes
    "ORS_Basal" = c("rKc1"),
    "ORS_Suprabasal" =c("rKc2"),
    "Bulge" =c("rKc3"), 
    "IRS_Cuticle_Huxley" =c("rKc4"),
    "IRS_Henle" = c("rKc5"),
    "Companion_Layer" = c("rKc6"),

    #Lymphoids
    "Macrophage.M1" = c("rTc1"), 
    "Macrophage.M2" = c("rTc2"), 
    "Lymphocyte" = c("rTc3"), 
    "CD8.Tc" = c("rTc4"), 
    "TFH.Tc" = c("rTc5"),
    "Immune.Keratinocytes" = c("rTc6"),
    "Langerhans.cells" = c("rTc7"),

    #Other
    "Melanocytes.2" = c("rOh11"),
    "Mast.cells" = c("rOh8"),
    "Proliferating.HM" = c("rOh5"),
    "Melanocytes.1" = c("rOh7"),
    "Cuticle" = c("rOh1"),
    "MFB-like cell" = c("rOh2"),
    "HM" = c("rOh3"),
    "Proliferating.HM" = c("rOh6"),
    "APM" = c("rOh4"), #Arrector Pili Muscle
    "Secretory.Duct" = c("rOh9"),
    "Secretory.Gl" = c("rOh10")


)  %>% invertList()


# Visualize similar clusters together by setting their main catagory. Requires identical labels between
# subCluster_annotation and BigCluster_Annotation. 

BigCluster_Annotation <- list(
    
    #Endothelial
    "ICAM1_Endothelial" = c("Endothelial"), #VFW/ICAM1/SELE
    "Vascular_Endothelial" = c("Endothelial"), #VWF/CCDC3

    #Fibroblasts
    "Dermal_Sheath" = c("Dermal_Sheath"), #ACTA2/SOX5/CD36
    "Dermal_Endothelial" = c("Endothelial"), #CD74/PECAM/PLVAP Initiatlly missclassfied to fibroblasts due to VIM
    "Dermal_Papilla" = c("Dermal_Papilla"), #CORIN/VCAN/CD44
    "Myofibroblast_like_Cells" = c("MFBs"), #RGS5/DSP
    "Myofibroblasts" = c("MFBs"), #IL6/SORBS2/TAGLN

    #Epithelial
    "Bulge_SC_Permanent" = c("Bulge_SCs"), #VAV3/KTR15/CD34 No CD200 expression
    "Isthmus_Sebocytes" = c("Isthmus"), #KRT14/SOX9/NFATC1/PRBM1(BLIMP1)
    "Bulge_SC_Maintenance" = c("Bulge_SCs"), #KRT15/CD200/TGFB2/IL31A Contains immune specific genes
    "Suprabasal_ORS" = c("Suprabasal_ORS"), #KRT6B/KRT16/KRT17 Additionally contains Companion Layer KRT75
    "Basal_ORS" = c("Basal_ORS"),  #KRT5/KRT14/S100A2
    "Isthmus" = c("Isthmus"), #KRT15/TRPC6/SEMA5A/LGR6 Mechanotransductive properties
    "Granular_Layer" = c("Granular_Layer"), #KRT1/KRTDAP/S100A8
    "Commmited_Bulge_SC" = c("Bulge_SCs"), #DIO2/KRT15/FGF14/CXCL14
    "Growth_Stimulating_KCs" = c("Proliferating_KCs"), #MDM2/MYOB5/TMCC3
    "Insulin-Responsive_Cells" = c("Infindibulum"), #IRS2/-SON Also has high mitochondrial activity
    "Lower_Spinous_Layer" = c("Spinous_Layer"), #IVL/SLC10A6/
    "TACs" = c("Bulge_SCs"), #KRT8/CDKN1C/FGL2 Perhaps lineage commited stem cells a better name?
    "Upper_Spinous_Layer" = c("Spinous_Layer"), #KRT1/KRT10/LYD6
    "Basal_Layer" = c("Basal_Layer"), #KRT14/SERPINB2/COL17A1/SOX6/PCDH7
    "Terminal_Control_Basal_KCs" = c("Interfollecular_Epidermis"), #MCU/NABP1/KRT15
    "Cuticle" = c("IRS"), #KRT32
    "Proliferating_KCs" = c("Matrix"), #KRT35/KRT85/LEF1
    "Matrix" = c("Matrix"), #MKI67
    "IRS" = c("Matrix"), #KRT25/KRT28
    "Germ_Cells" = c("Matrix"), #ENSG00000286533/CCNG2/INHBA

    #Immune
    "Dentritic_Cells" = c("Dendritic_Cells"), #PTPRC/CD74/HLA-DRA
    "Innate_Lymphoid_Cells" = c("Lymphoid_Cells"), #CDC14A/CD69/PTPRC/CD3D-
    "Gamma_Delta_T_Cells" = c("Lymphoid_Cells"), #PTPRC/CD3D/CD48/IL32

    #Other
    "Melanocytes_1" = c("Melanocytes"), #MLANA/CD200
    "Sebaceous_Gland" = c("Sebaceous_Gland"), #SCGB2A2/KRT15/WWC1
    "Melanocyte_2" = c("Melanocytes") #SOX5/BCL2/LRMDA/RAB38

) 