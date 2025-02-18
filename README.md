# Higgins-Lab - Automated scRNA/scATAC/RNA/Microarray analysis pipelines.

### Each folder contains an individual pipeline with dedicated **REnv/Conda** (renv or yaml) environments to ensure compatability.
|Pipeline|Language|Data Type|OS|
|---|---|---|---|
|Automated_Annotation|R/Py(Reticulate)|scRNA|![Windows](./Assets/icon_windows.png)![Linux](./Assets/icon_linux.png)![iOS](./Assets/icon_apple.png)
|Microarray_Analysis|R|Microarray/Bulk RNA|![Windows](./Assets/icon_windows.png)![Linux](./Assets/icon_linux.png)![iOS](./Assets/icon_apple.png)
|CellChat_Single|R/Py(Reticulate)|scRNA|![Windows](./Assets/icon_windows.png)![Linux](./Assets/icon_linux.png)![iOS](./Assets/icon_apple.png)
|CellChat_Comparison|R/Py(Reticulate)|scRNA|![Windows](./Assets/icon_windows.png)![Linux](./Assets/icon_linux.png)![iOS](./Assets/icon_apple.png)
|SIGNET_Pipeline|Py|scRNA|![Windows](./Assets/icon_windows.png)![Linux](./Assets/icon_linux.png)![iOS](./Assets/icon_apple.png)
|Alignment_Pipeline|Py/R|RNA/scRNA|![Linux](./Assets/icon_linux.png)
|scATAC_Pipeline|R|scATAC|![Linux](./Assets/icon_linux.png)
|CellOracle_Pipeline|Py|scRNA/scATAC|![Linux](./Assets/icon_linux.png)
|CellRank_Pipeline|Py|scRNA|![Linux](./Assets/icon_linux.png)

> [!TIP]
**Pipeline quick summarization:**

* Alignment_Pipeline: Automated Kallisto|Bustools sequence alignment. Contains R based Droplet removal for droplet-based sequence methods. Requires FASTQ files as input.

* Automated_Annotation: Pseudo-automated annotation of Hair/Skin scRNA datasets using Adjusted Iterative Latent Semantic Idexing (LSI). Requires MTX, H5, 10x as input.

* Microarray_Analysis: Simple automated Microarray/RNA preprocessing and analysis using Oligo and GAGE. Requires Microarray (CEL) or batch RNA seq as input.

* CellChat_Single: infers intercellular communication. Runs a standard CellChat pipeline on a single annotated scRNA dataset. Requires Seurat object as input.

* CellChat_Comparison: infers intercellular communication. Runs a comparison CellChat pipeline between multiple scRNA datasets. Requires Seurat object(s) as input.

* SIGNET_Pipeline: Infers intracellular communication. Runs a SIGNET pipeline to establish transcriptional network using MLP training model. Requires raw expression matrix as input.

* scATAC_Pipeline: Uses ArchR and Cicero to generate GRNs for the purpose of CellOracle pertubartion analysis. Can modify ArchR to annotate scATAC data easily.

* CellOracle_Pipeline: Uses CellOracle and customizable GRNs from scATAC_Pipeline to perform gene perturbation analysis on scRNA dataset to track lineage changes.

* CellRank_Pipeline: Uses Palentir pseudotime and CellRank for Unified fate mapping in multiview single-cell data


> [!IMPORTANT]
>**Important minimum prerequisites:**

Download and install R (4.3.2+): https://www.r-project.org/

Download and install RTools (4.3+): https://cran.r-project.org/bin/windows/Rtools/

Download and install RStudio (3.3+): https://posit.co/download/rstudio-desktop/

Download and install Python (3.11+): https://www.python.org/downloads/

Download and install Microsoft C++ Build Tools: https://visualstudio.microsoft.com/visual-cpp-build-tools/

(In the Visual Studio Installer make sure to highlight the Desktop development with C++ and from Individual Components the Windows 11 or Windows 10 SDK)

For iOS users:

Download and install XQuartz https://www.xquartz.org/

Download and install XCode https://developer.apple.com/xcode/
