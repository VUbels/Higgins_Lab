# Higgins_Lab

Automated scRNA analysis and annotation pipelines. Each folder contains an individual pipeline with dedicated Renvironments to ensure compatability independent of the user.

-Automated_Annotation: Pseudo-automated annotation of Hair/Skin scRNA datasets using Adjusted Iterative Latent Semantic Idexing (LSI). Requires MTX, H5, 10x as input.

-CellChat_Single: infers intercellular communication. Runs a standard CellChat pipeline on a single annotated scRNA dataset. Requires Seurat object as input.

-CellChat_Comparison: infers intercellular comumunication. Runs a comparison CellChat pipeline between multiple scRNA datasets. Requires Seurat object(s) as input.

-SIGNET pipeline: Infers intracellular communication. Runs a SIGNET pipeline to establish transcriptional network using MLP training model. Requires raw expression matrix as input.
