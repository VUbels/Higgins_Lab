### Markers for broad cell clusters ###

featureSets_broad <- list(
    "Keratinocytes" = c("KRT5", "KRT10", "KRT14", "KRT15", "KRT6A", "KRT6B"),
    "Fibroblasts" = c("THY1", "COL1A1", "COL1A2"), 
    "T_cells" = c("CD3D", "CD8A", "CD4", "FOXP3", "IKZF2", "CCL5"), # PTPRC = CD45 
    "B_cells" = c("CD19", "MS4A1", "MZB1"), # MS41A = CD20, MZB1 = marginal zone B cell specific protein
    "APCs" = c("CD86", "CD74", "CCR7", "CD163"), # Monocyte lineage (FCGR3A = CD16, FCGR1A = CD64, CD74 = HLA-DR antigens-associated invariant chain)
    "Melanocytes" = c("MITF", "TYR", "SOX10", "MLANA"), # Melanocyte markers
    "Endothlial" = c("VWF", "PECAM1", "SELE"), # Endothlial cells (PECAM1 = CD31), SELE = selectin E (found in cytokine stimulated endothelial cells)
    "Lymphatic" = c("ACKR2", "FLT4", "LYVE1"),  # Lymphatic endothelial (FLT4 = VEGFR-3)
    "Muscle" = c("TPM1", "TAGLN", "MYL9"), # TPM = tropomyosin, TAGLN = transgelin (involved crosslinking actin in smooth muscle), MYL = Myosin light chain
    "Mast_cells" = c("KIT", "FCER1A", "IL1RL1", "TPSB2"), # KIT = CD117, ENPP3 = CD203c, FCER1A = IgE receptor alpha, TPSB2 is a beta tryptase, which are supposed to be common in mast cells
    "Langerhans_cells" = c("ITGAX", "CD1A", "CLEC1A", "CD207"), # CD11c = ITGAX, Langerin = CD207, CLEC4K
    "HF_surface_markers" = c("GJB6", "ITGB8", "CD200", "FZD7")
)

### Markers for specific cell types ###

featureSets_specific <- list(
    "Basal_epithelial" = c("KRT15", "KRT5", "COL17A1", "KRT6", "KRT14"),
    "Spinous" = c("KRT1"),
    "HF_keratinocytes" = c("KRT75", "SOX9", "LHX2","ITGB8", "KRT16", "KRT17", "RUNX3", "KRT35", "LEF1", "GATA3", "TCHH", "KRT85"),
    "Glandular" = c("KRT7"),
    "T_cells" = c("CD3D", "CD8A", "CD4", "FOXP3", "IKZF2", "IFNG"),
    "B_cells" = c("MS4A1"), # MS41A = CD20, MZB1 = marginal zone B cell specific protein
    "M1_macs" = c("CCL20", "CD80", "CD86"),
    "M2a_macs" = c("CD163", "TGFB2"),
    "TREM2_macs" = c("TREM2", "OSM"),
    "FOLR2_macs" = c("FOLR2"),
    "CD1a1c_DCs" = c("CD1A", "CD1C", "ITGAX", "ITGAM"), # SIRPA = CD172a
    "CD1a141_DCs" = c("CLEC9A", "XCR1"), # THBD = CD141 = BDCA-3 (thrombomodulin)
    "Mast_cells" = c("KIT", "TPSB2"), # KIT = CD117, ENPP3 = CD203c, FCER1A = IgE receptor alpha, TPSB2 is a beta tryptase, which are supposed to be common in mast cells
    "Melanocytes" = c("MITF", "SOX10", "MLANA"), # Melanocyte markers
    "Endothlial" = c("VWF", "PECAM1", "SELE"), # Endothlial cells (PECAM1 = CD31), SELE = selectin E (found in cytokine stimulated endothelial cells)
    "Lymphatic" = c("FLT4", "LYVE1", "CCL21"),  # Lymphatic endothelial (FLT4 = VEGFR-3)
    "Angiogenic" = c("SEMA3G"),
    "Muscle" = c("TPM1", "TAGLN"), # TPM = tropomyosin, TAGLN = transgelin (involved crosslinking actin in smooth muscle), MYL = Myosin light chain
    "Fibroblasts" = c("THY1", "COL1A1"),
    "Dermal_sheath" = c("SOX2", "COL11A1", "DCN"), # Dermal Sheath? Hard to find clearly defined markers...
    "Papillary_dermis" = c("COL6A5", "APCDD1"), # PMID: 29391249
    "Reticular_dermis" = c("CD36"), # CD36 seems more specific for the 'muscle 2' cluster... Myofibroblast?
    "Dermal_Papilla" = c("BMP7", "HHIP", "PTCH1", "SOX18", 'THBS2', 'CCN2', 'RSPO2', 'RSPO3', "TWIST1", "TWIST2", "VIM", "RNX3"),
    "cycling" = c("MKI67", "CDK1", "TOP2A")
)

featureSets_Overview <- list(
    "Keratin_overview" = c("KRT1", "KRT2", "KRT3", "KRT4", "KRT5", "KRT6A", "KRT6B", "KRT6C", "KRT7", "KRT8", "KRT9", "KRT10", "KRT12", "KRT13", "KRT14", "KRT15", "KRT16", "KRT17", "KRT18", "KRT19", "KRT20", "KRT21", "KRT22", "KRT23", "KRT24", "KRT25", "KRT26", "KRT27", "KRT28", "KRT31", "KRT32",
    "KRT33A", "KRT33B", "KRT34", "KRT35", "KRT36", "KRT37", "KRT38", "KRT39", "KRT40", "KRT71", "KRT72",
    "KRT73", "KRT74", "KRT75", "KRT76", "KRT77", "KRT78", "KRT79", "KRT80", "KRT81", "KRT82", "KRT83", "KRT84", "KRT85", "KRT86", "SOX9", "TCHH"),
    "Fibroblast_overview" = c("TWIST1", "TWIST2", "PDGFRA", "ACTA2", "LEF1", "CORIN", "VIM", "DCN"),
    "Endothelial_cells" = c("PECAM1", "VWF"),
    "Melanocytes_cells" = c("MLANA", "MITF")


)


### Genes to be plotted for each subcluster after subclustering ###

featureSets_Im <- list(
    # T-cell subtypes:
    "Tregs" = c("FOXP3", "CD4", "IL2RA", "IKZF2"), # IL2RA = CD25
    "Th1_cells" = c("CCR1", "CCR5", "CXCR3", "TNF", "LTA", "TBX21"), # TBX21 = T-bet, TNF = TNF alpha, LTA = TNF-beta
    "Th17_cells" = c("RORA", "CCL20", "BATF", "IL1RL1", "IL6R", "IL17A", "IL17F", "IL21", "IL22"),
    "Th22_cells" = c("AHR", "CCR4", "CCR6", "CCR10", "IL6R", "IL10", "IL13", "IL22"),
    "TFH_cells" = c("IL21", "CXCL13", "IL17A", "IL17F", "BCL6"),
    "Memory_Tcells" = c("CD44", "IL7R", "CCR7", "BCL6"), 
    "CD8NK_Tcells" = c("CD8A", "KLRK1", "KLRD1", "IFNG", "CCL5", "GZMA", "GZMB"), # GZMK = Granzyme K, SELL = CD62L
    "Proliferating" = c("MKI67", "CDK1", "TOP2A"),
    "Tissue_res" = c("CCR7", "ITGAE", "SELL", "KLRG1",  # ITGAE = CD103; SELL = CD62L = L-selectin
        "CCR4", "CCR8", "CCR10",  "SELPLG") # SELPLG = CLA
)

featureSets_Ep <- list(
    "Matrix_Broad" = c("MSX2", "KRT31", "GATA3", "TCHH", "KRT35", "KRT81", "KRT83", "KRT85", "KRT86", "MKI67", "HOXC13", "LEF1"),
    "hHFSC" = c("KRT15", "COL17A1", "CD200", "LGR5", "CD34", "SOX9", "LHX2", "NFATC1", "FST"), # CD200 only expressed at top of bulge, not lower part of bulge
    "mORS" = c("KRT14", "KRT15", "KRT10"),
    "IBL" = c("KRT14", "KRT6C", "KRT75", "LHX2"),
    "SG_Precursor" = c("KRT15", "BLIMP1"),
    "IRS" = c("KRT71", "KRT73", "KRT74"),
    "CL" = c("KRT75"),
    "ICU" = c("KRT71", "KRT72"),
    "CU" = c("KRT82"),
    "Cortex" = c("KRT37", "KRT38"),
    "Spinous" = c("KRT1", "KRT10"),
    "Inner_Bulge_Telogen" = c("TIMP3", "PAI2"),
    "ORS_Differentiating" = c("PIA2"),
    "Infundibulum" = c("S100A8", "S100A9", "IVL"),
    "Isthmus" = c("PLET1", "LRIG1"),

    # "HFSCs" = c("SOX9", "LHX2", "NFATC1", "TCF3", # Key HFSC TFs
    # "ITGA6", "CD200", "FRZB", "IL31RA", "IL13RA1", "OSMR", # Other bulge markers (IL31RA pairs w/ OSMR)
    # "CD34", "CDH3", "LGR5", "LGR6", "RUNX1" # Hair germ markers
    # ), # CDKN2A = P16
    "Basal_epithelial" = c("KRT15", "KRT14", "KRT5", "COL17A1", "TNFRSF12A", "FN14"),
    "Granular" = c("DSC1", "KRT2", "IVL", "TGM3"),
    "RUNX3_high" = c("RUNX1", "RUNX2", "RUNX3", "KRT23", "KRT18"),
    "HairGerm" = c("CD34", "CDH3", "LGR5", "CDKN2A", "RUNX1", "PCDH7"), # Hair germ marker
    "TFs" = c("SOX9", "LHX2", "NFATC1", "TCF3"), # Key HFSC TFs # Sebaceous
    "Glandular" = c("SCGB2A2", "SCGB1D2", "KRT7", "KRT8", "KRT19", "AQP5"),
    "Proliferating" = c("MKI67", "CDK1", "TOP2A")
)

featureSets_Fb <- list(
    "Fibroblasts" = c("THY1", "COL1A1", "COL1A2", "COL3A1", "DCN", "MGP", "COL6A2", "CEBPB", "APOD", "CFD"),
    "HF_associated" = c("APCDD1", "VCAN", "CORIN", "PTGDS", "SOX2", "COL11A1", "DCN", "ALPL"), # Dermal Sheath
    "ImmuneRecruiting" = c("CXCL1", "CXCL2", "CXCL14", "CD44"),
    "Papillary_dermis" = c("COL6A5", "APCDD1", "HSPB3", "WIF1", "ENTPD1"),
    "Reticular_dermis" = c("CD36"),
    "Dermal_Papilla" = c("WNT5A", "BMP4", "BMP7", "HHIP", "PTCH1", "ITGA9", "SOX18", "RUNX1", "RUNX3", "ALX4", "TNFSF12", "VIM", "RSPO2", "RSPO3"),
    "Dermal_Cup" = c("CD133", "SOX2", "ACTA2", "LEF1", "RSPO4", "PAX1", "RSPO2", "ITGA5"), #LEF1 Should be not expressed
    "Dermal_Sheath_Inner" = c("COL4A1", "ITGA8", "ITGA5"),
    "Dermal_Sheath_Outer" = c("LAM", "ITGA8", "ITGA5")
)

featureSets_Ed <- list(
    "Endothlial" = c("VWF", "PECAM1", "SELE", "FLT1"), # Endothlial cells (PECAM1 = CD31), SELE = selectin E (found in cytokine stimulated endothelial cells), FLT1 = VEGFR1
    "Lymphatic" = c("ACKR2", "FLT4", "LYVE1", "CCL21", "TFF3", "APOD"),  # Lymphatic endothelial (FLT4 = VEGFR-3)
    "Angiogenic" = c("SEMA3G", "FBLN5", "NEBL", "CXCL12", "PCSK5", "SPECC1"),
    "cluster_markers" = c("SEMA3G", "FBLN5", "CXCL12", "TXNIP", "ZFP36", "FOS", "SOCS3", 
        "CSF3", "IL6", "MPZL2", "FKBP11", "RGCC", "RBP7", "BTNL9")
)

featureSets_Oh <- list(
    "Mast_cells" = c("KIT", "ENPP3", "FCER1A", "IL1RL1", "TPSB2"), # KIT = CD117, ENPP3 = CD203c, FCER1A = IgE receptor alpha, TPSB2 is a beta tryptase
    "Macrophages" = c("CD163", "LGMN", "FCGR2A", "C1QB", "C5AR1", "MAFB", "FOLR2"),
    "M1_macs" = c("CCL20", "CXCL3", "IL1B", "IL6", "IL12A", "IFNG", "TNF", "CD163"), # CD163 should be NEGATIVE
    "M2a_macs" = c("CD163", "CD200", "IRF4", "TGFB1", "TGFB2", "CCL2", "STAT6"),
    "TREM2_macs" = c("TREM2", "C3", "FCGBP", "FCGR3A", "OSM", "APOE"),
    "Langerhans_cells" = c("ITGAX", "CD1A", "CLEC1A", "CD207", "EPCAM"), # CD11c = ITGAX, Langerin = CD207, CLEC4K
    "pDC" = c("CCR7", "PTPRC", "CD209", "CLEC4C"), # PTPRC = CD45RA; IFNA1 = interferon alpha
    "moDC" = c("CD14", "CD1A", "CD1C", "ITGAX", "ITGAM", "SIRPA"), # SIRPA = CD172a
    "cDC1" = c("BTLA", "ITGAE", "CD1A", "ITGAM", "CLEC9A", "XCR1", "THBD"), # THBD = CD141 = BDCA-3 (thrombomodulin)
    "cDC2" = c("CD14", "CD163", "CLEC10A", "NOTCH2", "ITGAM", "SIRPA", "CX3CR1", "CD1C", "CD2"), # THBD = CD141 = BDCA-3 (thrombomodulin)
    "TCR_macs" = c("CD3D", "TRAC", "TRBC1", "SPOCK2", "CD14", "CD2"),
    "Basal" = c("KRT14", "KRT5", "KRT15", "COL17A1"), # Basal epithelia
    "HFSCs" = c("ITGA6", "ITGB1", "CD200", "LGR5","LHX2", "FRZB", "FZD1", "FZD5", "FZD10",  "IL31RA", "OSMR"), # HFSCs
    "HairGerm" = c("CD34", "CDH3", "LGR5", "CDKN2A", "RUNX1"), # Hair germ markers
    "Matrix" = c("KRT81", "KRT83", "HOXC13", "LEF1"), # Matrix hair keratins/genes
    "Sheath" = c("KRT71", "KRT75"), # IRS/ORS keratins / genes
    "TFs" = c("SOX9", "LHX2", "NFATC1", "TCF3") # Key HFSC TFs
    
)

featureSetsList <- list(featureSets_Im, featureSets_Ep, featureSets_Fb, featureSets_Ed, featureSets_Oh)

plot_genes_Im <- c(
    "CCL5", "CD8A", "GZMK", "IFNG", "TNF",
    "XCL1", "GNLY", "NKG7", "KLRK1", "PTPRC",
    "NR3C1", "IL7R", "CD3D", "CD69", "CDC14A",
    "RUNX3", "JUN", "CD4", "CD69", "CD48", "IL32",
    "IKZF2", "IL2RA", "FOXP3", "HLA-DRA",
    "MKI67", "TOP2A"
)

plot_genes_Ep <- c(
  "VAV3", "KTR15", "CD34", "KRT14", "SOX9", "NFATC1", "PRBM1", "KRT23", "MKI67", 
  "KRT15", "CD200", "TGFB2", "IL31A", "KRT6B", "KRT16", "KRT17", "KRT75", "TOP2A",
  "KRT5", "S100A2", "TRPC6", "SEMA5A", "LGR6", "KRTDAP", "S100A8", "DIO2", "ITGB8",
  "FGF14", "CXCL14", "MDM2", "MYOB5", "TMCC3", "IRS2", "SON", "IVL", "RUNX3",
  "SLC10A6", "KRT8", "CDKN1C", "FGL2", "KRT10", "LYD6", "SERPINB2", "KRT7", "AQP5",
  "COL17A1", "SOX6", "PCDH7", "MCU", "NABP1", "KRT25", "KRT28", "LHX2", "KRT19",
  "KRT35", "KRT85", "LEF1", "ANXA2", "ENSG00000286533", "CCNG2", "INHBA", "ITGA6"
)

plot_genes_Fb <- c(
    "THY1", "COL1A1", "COL1A2", 
    "CXCL1", "CXCL2", # Fb_1
    "CCL19", "CXCL12",
    "APCDD1", "COL18A1", # rVe3
    "WISP2", "SCN7A", "CCL21", "PDGFD",
    "COL11A1", "EDNRA", "SOX2",
    "HHIP", "PTCH1", "WNT5A",
    "TNFSF12", "TNFRSF12A",
    "ACTA2", "SOX5", "CD36", "CD74", "PECAM", "PLVAP",
    "CORIN", "VCAN", "CD44", "RGS5", "DSP", "IL6", "SORBS2", "TAGLN"
)

plot_genes_Ed <- c(
    "VWF", "PECAM1", "SELE", "ACKR1", "IL6", "SOD2",
    "BTNL9", "RGCC", 
    "PRCP", "TXNIP", 
    "FLT4", "LYVE1", "CCL21",
    "SEMA3G", "HEY1", "NEBL", "CXCL12"
)

plot_genes_Oh <- c(
    "IL15", "CCR7", "CCL19", "CCL17",
    "CD3D", "TRAC", "MLANA", "WWC1", "SOX5",
    "FOLR2", "C1QA", "CD163", "BCL2", "LRMDA", "RAB38",
    "CXCL2", "CCL20",
    "TREM2", "OSM", "CD14",
    "FCER1A", "CLEC10A", "CD1C",
    "CLEC9A", "XCR1", 
    "IFNG", "CD207",
    "JCHAIN", "KRT15", "KRT14", "COL17A1", "ITGA6",
    "KRT1", "KRT10",
    "SOX9", "LHX2",
    "RUNX3", "KRT23",
    "KRT75", 
    "KRT7", "KRT19", "AQP5",
    "ITGB8", "CD200", "HLA-A",
    "MKI67", "TOP2A"
    
)

plot_genesList <- list(plot_genes_Im, plot_genes_Ep, plot_genes_Fb, plot_genes_Ed, plot_genes_Oh)

clustOrder_Im <- c(
    "CD4.Tc",
    "CD8.Tc",
    "Treg",
    "NK",
    "Cyc.Tc"
)

clustOrder_Ep <- c(
    "Basal.Kc",
    "Spinous.Kc",
    "Proliferating.Kc",
    "Inf.Segment",
    "Outer Root Sheath",
    "Isthmus",
    "Sebaceous",
    "Eccrine" 
)

clustOrder_Fb <- c(
    "D.Fibroblasts", # "Fb_1", # CXCL1,2,3
    "D.Sheath", # "Fb_2", # CCL19, CXCL12
    "D.Papilla" # "Fb_3", # APCDD1, COL18A1, F13A1
)

clustOrder_Ed <- c(
    "Vas.Endo", 
    "Lymph.Endo",
    "Angiogenic"
)

clustOrder_Oh <- c(
    "M2.Mac",
    "TREM2.Mac",
    "CLEC9.DC",
    "cDC2",
    "Bulge",
    "Dermal.Mac",
    "Hair Matrix",
    "Bulge_2",
    "M1.Mac", 
    "IRF4.cells"
)

clustOrderlist <- list(clustOrder_Im, clustOrder_Ep, clustOrder_Fb, clustOrder_Ed, clustOrder_Oh)

### Marker sets for 

# Labels for clusters

BroadClust <- list(
    "Tc" = "Lymphoid",
    "My" = "Myeloid",
    "Ma" = "Mast",
    "Kc" = "Keratinocytes",
    "Fb" = "Fibroblasts",
    "Ve" = "Vascular",
    "Le" = "Lymphatic",
    "Mu" = "Muscle", 
    "Me" = "Melanocytes",
    "Bc" = "Plasma",
    "Other" = "Other"
)


s.genes <- c("MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", 
"GINS2", "MCM6", "CDCA7", "DTL", "PRIM1", "UHRF1", "MLF1IP", 
"HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", 
"SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2", "RAD51", "RRM2", 
"CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", 
"USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8")

g2m.genes <- c("HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", 
"NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", 
"FAM64A", "SMC4", "CCNB2", "CKAP2L", "CKAP2", "AURKB", "BUB1", 
"KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", 
"HN1", "CDC20", "TTK", "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", 
"DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA", 
"PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", 
"GAS2L3", "CBX5", "CENPA")


