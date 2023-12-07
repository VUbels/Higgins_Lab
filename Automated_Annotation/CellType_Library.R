### Markers for broad cell clusters ###

featureSets_broad <- list(
    "Keratinocytes" = c("KRT5", "KRT10", "KRT14", "KRT15", "KRT6A", "KRT6B"), 
    "Fibroblasts" = c("THY1", "COL1A1", "COL1A2"), 
    "T_cells" = c("CD3D", "CD8A", "CD4", "FOXP3", "IKZF2", "CCL5"), # PTPRC = CD45 
    "B_cells" = c("CD19", "MS4A1", "MZB1"), # MS41A = CD20, MZB1 = marginal zone B cell specific protein
    "APCs" = c("CD86", "CD74", "CCR7", "CD163"), # Monocyte lineage (FCGR3A = CD16, FCGR1A = CD64, CD74 = HLA-DR antigens-associated invariant chain)
    "Melanocytes" = c("MITF", "TYR", "SOX10", "MLANA"), # Melanocyte markers
    "Endothlial" = c("VWF", "PECAM1", "SELE"), # Endothlial cells (PECAM1 = CD31), SELE = selectin E (found in cytokine stimulated endothelial cells)
    "Lymphatic" = c("ACKR2", "FLT4", "LYVE1"),  # Lymphatic endothelial (FLT4 = VEGFR-3))
    "Muscle" = c("TPM1", "TAGLN", "MYL9"), # TPM = tropomyosin, TAGLN = transgelin (involved crosslinking actin in smooth muscle), MYL = Myosin light chain
    "Mast_cells" = c("KIT", "FCER1A", "IL1RL1", "TPSB2"), # KIT = CD117, ENPP3 = CD203c, FCER1A = IgE receptor alpha, TPSB2 is a beta tryptase, which are supposed to be common in mast cells
    "Langerhans_cells" = c("ITGAX", "CD1A", "CLEC1A", "CD207"), # CD11c = ITGAX, Langerin = CD207, CLEC4K
    "HF_surface_markers" = c("GJB6", "ITGB8", "CD200", "FZD7")
)

### Markers for specific cell types ###

featureSets_specific <- list(
    "Basal_epithelia" = c("KRT15", "KRT5", "COL17A1"),
    "Spinous" = c("KRT1"),
    "HF_keratinocytes" = c("KRT75", "SOX9", "LHX2","ITGB8", "KRT16", "KRT17", "RUNX3"),
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
    "Lymphatic" = c("FLT4", "LYVE1", "CCL21"),  # Lymphatic endothelial (FLT4 = VEGFR-3))
    "Angiogenic" = c("SEMA3G"),
    "Muscle" = c("TPM1", "TAGLN"), # TPM = tropomyosin, TAGLN = transgelin (involved crosslinking actin in smooth muscle), MYL = Myosin light chain
    "Fibroblasts" = c("THY1", "COL1A1"),
    "Dermal_sheath" = c("SOX2", "COL11A1"), # Dermal Sheath? Hard to find clearly defined markers...
    "Papillary_dermis" = c("COL6A5", "APCDD1"), # PMID: 29391249
    "Reticular_dermis" = c("CD36"), # CD36 seems more specific for the 'muscle 2' cluster... Myofibroblast?
    "Dermal_Papilla" = c("BMP7", "HHIP", "PTCH1", "SOX18", 'THBS2'),
    "cycling" = c("MKI67", "CDK1", "TOP2A")
)


### Genes to be plotted for each subcluster after subclustering ###

featureSets_Tc <- list(
    # T-cell subtypes:
    "Tregs" = c("FOXP3", "CD4", "IL2RA", "IKZF2"), # IL2RA = CD25
    "Th1_cells" = c("CCR1", "CCR5", "CXCR3", "TNF", "LTA", "TBX21"), # TBX21 = T-bet, TNF = TNF alpha, LTA = TNF-beta
    "Th17_cells" = c("RORA", "CCL20", "BATF", "IL1RL1", "IL6R", "IL17A", "IL17F", "IL21", "IL22"),
    "Th22_cells" = c("AHR", "CCR4", "CCR6", "CCR10", "IL6R", "IL10", "IL13", "IL22"),
    "TFH_cells" = c("IL21", "CXCL13", "IL17A", "IL17F"),
    "Memory_Tcells" = c("CD44", "IL7R", "CCR7", "BCL6"), 
    "CD8NK_Tcells" = c("CD8A", "KLRK1", "KLRD1", "IFNG", "CCL5", "GZMA", "GZMB"), # GZMK = Granzyme K, SELL = CD62L
    "Proliferating" = c("MKI67", "CDK1", "TOP2A"),
    "Tissue_res" = c("CCR7", "ITGAE", "SELL", "KLRG1",  # ITGAE = CD103; SELL = CD62L = L-selectin
        "CCR4", "CCR8", "CCR10",  "SELPLG") # SELPLG = CLA
)

featureSets_Kc <- list(
    "HFSCs" = c(
    "SOX9", "LHX2", "NFATC1", "TCF3", # Key HFSC TFs
    "ITGA6", "CD200", "FRZB", "IL31RA", "IL13RA1", "OSMR", # Other bulge markers (IL31RA pairs w/ OSMR)
    "CD34", "CDH3", "LGR5", "LGR6", "RUNX1" # Hair germ markers
    ), # CDKN2A = P16
    "Basal_epithelia" = c("KRT15", "KRT14", "KRT5", "COL17A1"),
    "Spinous" = c("KRT1", "KRT10"),
    "Granular" = c("DSC1", "KRT2", "IVL", "TGM3"),
    "RUNX3_high" = c("RUNX1", "RUNX2", "RUNX3", "KRT23", "KRT18"),
    "HF_keratinocytes" = c(
      "KRT81", "KRT82", "LEF1", # Matrix hair keratins/genes
      "KRT75", # ORS keratin
      "SOX9", "LEF1", "HOXC13", # HF TFs
      "MGST1", "COL14A1", "CD200",
      "ELOVL3" # Sebaceous 
      ),
    "Glandular" = c("SCGB2A2", "SCGB1D2", "KRT7", "KRT8", "KRT19", "AQP5"),
    "Proliferating" = c("MKI67", "CDK1", "TOP2A")
)

featureSets_Fb <- list(
    "Fibroblasts" = c("THY1", "COL1A1", "COL1A2", "COL3A1", "DCN", "MGP", "COL6A2", "CEBPB", "APOD", "CFD"),
    "HF_associated" = c("APCDD1", "VCAN", "CORIN", "PTGDS", "SOX2", "COL11A1"), # Dermal Sheath
    "ImmuneRecruiting" = c("CXCL1", "CXCL2", "CXCL14", "CD44"),
    "Papillary_dermis" = c("COL6A5", "APCDD1", "HSPB3", "WIF1", "ENTPD1"), # PMID: 29391249
    "Reticular_dermis" = c("CD36"),
    "Dermal_Papilla" = c("WNT5A", "BMP4", "BMP7", "HHIP", "PTCH1", "SOX18", "RUNX1", "RUNX3", "ALX4")
)

featureSets_Ed <- list(
    "Endothlial" = c("VWF", "PECAM1", "SELE", "FLT1"), # Endothlial cells (PECAM1 = CD31), SELE = selectin E (found in cytokine stimulated endothelial cells), FLT1 = VEGFR1
    "Lymphatic" = c("ACKR2", "FLT4", "LYVE1", "CCL21", "TFF3", "APOD"),  # Lymphatic endothelial (FLT4 = VEGFR-3))
    "Angiogenic" = c("SEMA3G", "FBLN5", "NEBL", "CXCL12", "PCSK5", "SPECC1"),
    "cluster_markers" = c("SEMA3G", "FBLN5", "CXCL12", "TXNIP", "ZFP36", "FOS", "SOCS3", 
        "CSF3", "IL6", "MPZL2", "FKBP11", "RGCC", "RBP7", "BTNL9")
)

featureSets_Oh <- list(
    # APC subtypes:
    #"Mast_cells" = c("KIT", "ENPP3", "FCER1A", "IL1RL1", "TPSB2"), # KIT = CD117, ENPP3 = CD203c, FCER1A = IgE receptor alpha, TPSB2 is a beta tryptase
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

featureSetsList <- list(featureSets_Tc, featureSets_Kc, featureSets_Fb, featureSets_Ed, featureSets_Oh)

plot_genes_Tc <- c(
    "CCL5", "CD8A", "GZMK", "IFNG", "TNF",
    "XCL1", "GNLY", "NKG7", "KLRK1",
    "NR3C1", "IL7R",
    "RUNX3", "JUN", "CD4", "CD69",
    "IKZF2", "IL2RA", "FOXP3",
    "MKI67", "TOP2A"
)

plot_genes_Kc <- c(
    "KRT15", "KRT14", "COL17A1", "ITGA6",
    "KRT1", "KRT10",
    "SOX9", "LHX2",
    "RUNX3", "KRT23",
    "KRT75",
    "KRT7", "KRT19", "AQP5",
    "ITGB8", "CD200", "HLA-A",
    "MKI67", "TOP2A"
)

plot_genes_Fb <- c(
    "THY1", "COL1A1", "COL1A2", 
    "CXCL1", "CXCL2", # Fb_1
    "CCL19", "CXCL12",
    "APCDD1", "COL18A1", # rVe3
    "WISP2", "SCN7A", "CCL21", "PDGFD",
    "COL11A1", "EDNRA", "SOX2",
    "HHIP", "PTCH1", "WNT5A"
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
    "CD3D", "TRAC",
    "FOLR2", "C1QA", "CD163",
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

plot_genesList <- list(plot_genes_Tc, plot_genes_Kc, plot_genes_Fb, plot_genes_Ed, plot_genes_Oh)

clustOrder_Tc <- c(
    "CD4.Tc",
    "CD8.Tc",
    "Treg",
    "NK",
    "Cyc.Tc"
)

clustOrder_Kc <- c(
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

clustOrderlist <- list(clustOrder_Tc, clustOrder_Kc, clustOrder_Fb, clustOrder_Ed, clustOrder_Oh)

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

rna.NamedClust <- list(
    "rFb1" = "D.Fib", # Papillary/Reticular dermal fibroblasts
    "rFb2" = "D.Sheath",
    "rTc1" = "Tregs",
    "rTc2" = "CD4.Tc",
    "rTc3" = "CD8.Tc",
    "rMy1" = "DCs_1",
    "rMy2" = "Macs_1",
    "rMy3" = "CLEC9a.DC",
    "rMy4" = "M1.macs",
    "rMa1" = "Mast", # Mast cells
    "rKc1" = "Spinous.Kc_1",
    "rKc2" = "Spinous.Kc_2",
    "rKc3" = "HF.Kc_1",
    "rKc4" = "Basal.Kc_1",
    "rKc5" = "Basal.Kc_2",
    "rMu1" = "Muscle",
    "rMu2" = "Pericytes",
    "rVe1" = "Vas.Endo",
    "rLe1" = "Lymph.Endo",
    "rMe1" = "Melanocytes_1",
    "rMe2" = "Melanocytes_2",
    "rBc1" = "Plasma"
)

rna.FineClust <- list(
    # Lymphoid / T-cells
    "rIm1" = "CD4.Tc", # T-helper: NR3C1, RORA, IL7R, CREM, etc. 
    "rIm2" = "CD4.Tc", # T-helper: JUN, FOS, HSP, CD69 etc.
    "rIm3" = "CD8.Tc",  # Cytotoxic T-cells: CCL4, CCL5, CD8A, GZMK, IFNG  
    "rIm4" = "Treg", # Regulatory T cells: IKZF2, IL2RA, CTLA4, FOXP3 
    "rIm5" = "NK", # Natural Killer cells: XCL1, XCL2, GNLY, NKG7, KLRD1
    "rIm6" = "Cyc.Tc", # Cycling T cells: MKI67, TOP2A,
    
    "CD4.Tc" = "CD4.Tc",
    "CD8.Tc" = "CD8.Tc",
    "Treg" = "Treg",
    "NK" = "NK",
    "Cyc.Tc" = "Cyc.Tc",
    
    # Myeloid
    "rMy1" = "cDC2", # CD1c, CLEC10a (conventional DCs - type 2)
    "rMy2" = "M2.macs_1", # C1Qa/b/c, FOLR2, CD14, CD163, (CCL13)
    "rMy3" = "M2.macs_2", # CXCL2, CXCL3, (CCL20, S100A8/9) 
    "rMy4" = "CLEC9a.DC", # CLEC9a, CLEC4C, XCR1 https://www.frontiersin.org/articles/10.3389/fimmu.2014.00239/full
    "rMy5" = "M1.macs", # IL15, IL32, CCR7 (CCL19, CCL17)
    "rMy6" = "TCR.macs", # CD3, TCR gene positive macrophages
    "rMy7" = "TREM2.macs", # TREM2
    "rMy8" = "Plasma.contam", # (MS4A1, IGHM, CD79A; likely contaminating plasma cells / doublets)
    
    # Keratinocytes
    "rKc1" = "Basal.Kc",
    "rKc2" = "Spinous.Kc",
    "rKc3" = "Spinous.Kc", 
    "rKc4" = "Spinous.Kc", # Borderline HF (SOX9 / KRT75 weak)
    "rKc5" = "Proliferating.Kc", # Cdk1/Top2a/Mki67
    "rKc6" = "Inf.Segment", # Lhx2, LGR5 high (Inferior segment)
    "rKc7" = "Outer Root Sheath", #Krt23/Sox9/ELOVL/RUNX3
    "rKc8" = "Isthmus", # CD200, AR, PPARG, KRT7/CXCL14 (Isthmus/Sebaceous)
    "rKc9" = "Sebaceous", # KRT17, KRT18, KRT23, CDH11, S100P, RUNX3
    "rKc10" = "Eccrine", # AQP5, KRT7
    
    "Basal.Kc" = "Basal.Kc",
    "Spinous.Kc" = "Spinous.Kc",
    "Proliferating.Kc" = "Proliferating.Kc",
    "Inf.Segment" = "Inf.Segment",
    "Outer Root Sheath" = "Outer Root Sheath",
    "Isthmus" = "Isthmus",
    "Sebaceous" = "Sebaceous",
    "Eccrine" = "Eccrine", 
    
    # Fibroblasts
    "rFb1" = "D.Fibroblasts", # CXCL1,2,3
    "rFb2" = "D.Sheath", # COL11A1, EDNRA
    "rFb3" = "D.Fibroblasts", # CCL19, CXCL12
    "rFb4" = "D.Fibroblasts", # APCDD1, COL18A1, F13A1
    "rFb5" = "D.Fibroblasts", # WISP2, AOX1, ARFGEF3
    "rFb6" = "D.Papilla", # HHIP, PTCH1, etc.
    
    "D.Fibroblasts" = "D.Fibroblasts",
    "D.Sheath" = "D.Sheath",
    "D.Papilla" = "D.Papilla",
    
    # Endothelial
    "rEd1" = "Vas.Endo", #
    "rEd2" = "Vas.Endo", #
    "rEd3" = "Vas.Endo", #
    "rEd4" = "Lymph.Endo", #
    "rEd5" = "Angiogenic", #
    
     "Vas.Endo" = "Vas.Endo",
    "Lymph.Endo" = "Lymph.Endo",
    "Angiogenic" = "Angiogenic",

     # Other
    "rOh11" = "D.Sheath",
    "rOh8" = "M2.Mac",
    "rOh5" = "TREM2.Mac",
    "rOh7" = "CLEC9.DC",
    "rOh1" = "cDC2",
    "rOh2" = "Bulge",
    "rOh3" = "Dermal.Mac",
    "rOh6" = "Hair.Matrix",
    "rOh4" = "Bulge.2",
    "rOh9" = "M1.Mac",
    "rOh10" = "IRF4.cells",
 
    "M2.Mac" = "M2.Mac",
    "TREM2.Mac" = "TREM2.Mac",
    "CLEC9.DC" = "CLEC9.DC",
    "cDC2" = "cDC2",
    "Bulge" = "Bulge",
    "Dermal.Mac" = "Dermal.Mac",
    "Hair Matrix" = "Hair.Matrix",
    "Bulge_2" = "Bulge.2",
    "M1.Mac" = "M1.Mac", 
    "IRF4.cells" = "IRF4.cells"
)

atac.NamedClust <- list(
    "aFb1" = "D.Fib", # Papillary/Reticular dermal fibroblasts
    "aFb2" = "D.Sheath",
    "aTc1" = "CD4.Tc",
    "aTc2" = "CD8.Tc",
    "aTc3" = "Tregs",
    "aMy1" = "DCs_1",
    "aMy2" = "Macs_1",
    "aMy3" = "CLEC9a.DC",
    "aKc1" = "Basal.Kc_1",
    "aKc2" = "Spinous.Kc_1",
    "aKc3" = "Spinous.Kc_2",
    "aKc4" = "HF.Kc_1",
    "aKc5" = "HF.Kc_2",
    "aKc6" = "HF.Kc_3",
    "aKc7" = "HF.Kc_4",
    "aMu1" = "Muscle",
    "aMu2" = "Pericytes",
    "aVe1" = "Vas.Endo_1",
    "aVe2" = "Vas.Endo_2",
    "aLe1" = "Lymph.Endo",
    "aMe1" = "Melanocytes",
    "aBc1" = "Plasma"
)

atac.FineClust <- list(
    # Lymphoid / T-cells
    "aTc1" = "CD4.Tc_1", # T-helper: NR3C1, RORA, IL7R, CREM, etc. 
    "aTc2" = "CD4.Tc_2", # T-helper: JUN, FOS, HSP, CD69 etc.
    "aTc3" = "CD8.Tc",  # Cytotoxic T-cells: CCL4, CCL5, CD8A, GZMK, IFNG  (Also NK cells absorbed here)
    "aTc4" = "Treg", # Regulatory T cells: IKZF2, IL2RA, CTLA4, FOXP3 
    "aTc5" = "CD4.Tc_3",  
    # Myeloid
    "aMy1" = "M2.macs_1",
    "aMy2" = "cDC2_1", # CD1c, CLEC10a (conventional DCs - type 2)
    "aMy3" = "M2.macs_2", # (between M2 macs and DCs)
    "aMy4" = "M2.macs_3", # CXCL8
    "aMy5" = "CLEC9a.DC", # CLEC9a, CLEC4C, XCR1 https://www.frontiersin.org/articles/10.3389/fimmu.2014.00239/full
    "aMy6" = "M1.macs", # IL15, IL32, CCR7 (CCL19, CCL17)
    "aMy7" = "cDC2_2",
    # Keratinocytes
    "aKc1" = "Basal.Kc_1",
    "aKc2" = "Spinous.Kc_2",
    "aKc3" = "Spinous.Kc_1",
    "aKc4" = "Infundibulum", # SOX9, DKK3
    "aKc5" = "Inf.Segment_1", # Lhx2, LGR5 high
    "aKc6" = "Sebaceous", # RUNX3, KRT23, KRT18
    "aKc7" = "Inf.Segment_2", # Lhx2, LGR5 high
    "aKc8" = "Isthmus", # CD200, AR high
    "aKc9" = "Matrix", 
    "aKc10" = "Eccrine", # AQP5
    "aKc11" = "Unknown_1", # (Suspected doublet)
    # Fibroblasts
    "aFb1" = "D.Fib_1", 
    "aFb2" = "D.Fib_2", 
    "aFb3" = "D.Sheath", # COL11A1
    "aFb4" = "D.Fib_3", # NECAB1, SCN7A (rFb5)
    "aFb5" = "D.Fib_4",
    "aFb6" = "D.Papilla", # HHIP, PTCH1, etc.
    # Endothelial
    "aVe1" = "Vas.Endo_1",
    "aVe2" = "Vas.Endo_2",
    "aVe3" = "Vas.Endo_3",
    "aVe4" = "Unknown_2",
    "aLe1" = "Lymph.Endo",
    # Non-subclustered
    "aMu1" = "Muscle",
    "aMu2" = "Pericytes",
    "aMe1" = "Melanocytes",
    "aBc1" = "Plasma",
    "Other" = "Other",
    # Hair Folicle Subclustered
    "aHF1" = "Sheath_1",
    "aHF2" = "Sheath_2",
    "aHF3" = "Migrating",
    "aHF4" = "HG",
    "aHF5" = "Matrix",
    "aHF6" = "HFSCs"
)

