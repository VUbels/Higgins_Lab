# Higgins_Lab

Automated scRNA analysis and annotation pipelines. Each folder contains an individual pipeline with dedicated Renvironments to ensure compatability independent of the user.

-Automated_Annotation: Pseudo-automated annotation of Hair/Skin scRNA datasets using Adjusted Iterative Latent Semantic Idexing (LSI). Requires MTX, H5, 10x as input.

-Microarray_Analysis: Simple automatated microarray preprocessing and analysis using Oligo and GAGE. Requires Microarray (CEL) or batch RNA seq as input.

-CellChat_Single: infers intercellular communication. Runs a standard CellChat pipeline on a single annotated scRNA dataset. Requires Seurat object as input.

-CellChat_Comparison: infers intercellular communication. Runs a comparison CellChat pipeline between multiple scRNA datasets. Requires Seurat object(s) as input.

-SIGNET pipeline: Infers intracellular communication. Runs a SIGNET pipeline to establish transcriptional network using MLP training model. Requires raw expression matrix as input.

-scATAC_Pipeline: Uses ArchR and Cicero to generate GRNs for the purpose of CellOracle pertubartion analysis. Can modify ArchR to annotate scATAC data easily. 

**Important prerequisites:**

-Download and install R (4.3.2+): https://www.r-project.org/

-Download and install RTools (4.3+): https://cran.r-project.org/bin/windows/Rtools/

-Download and install RStudio (3.3+): https://posit.co/download/rstudio-desktop/

-Download and install Python (3.11+): https://www.python.org/downloads/

-Download and install Microsoft C++ Build Tools: https://visualstudio.microsoft.com/visual-cpp-build-tools/

(In the Visual Studio Installer make sure to highlight the Desktop development with C++ and from Individual Components the Windows 11 or Windows 10 SDK)
