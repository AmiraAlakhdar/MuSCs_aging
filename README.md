This code was used to analyse the scRNA data for the paper titled:

# Apoptotic processes direct the matrix-induced fibrogenic conversion of aged muscle stem cells (in preperation)

The stiffening of the extracellular matrix (ECM) with age hinders muscle regeneration by causing intrinsic muscle stem cell (MuSC) dysfunction through a poorly understood mechanism. Here, young and aged MuSCs were seeded onto substrates engineered to mimic a soft and stiff ECM microenvironment to study those molecular changes using single-cell RNA-sequencing (scRNA). Data revealed the presence of a branching point within the trajectory leading to the emergence of an age-related fibroblastic population characterized by activation of the TNF-related apoptosis-inducing ligand (TRAIL) pathway, which was significantly activated in aged cells cultured on stiff substrates. Next, using the collagen cross-linking inhibitor β-aminopropionitrile (BAPN) in vivo, we elucidated stiffness changes on TRAIL downstream apoptotic targets (caspase 8 and caspase 3) using immunostaining. TRAIL activity was significantly inhibited by BAPN in aged animals, indicating a complex mechanism of age-related declines in muscle function through inflammatory and apoptotic mediators.
***
## What to expect
Here, I included the following R scripts: <br />
1- Loading data and quality control (loading_count_data_and_quality_control.R) <br />
2- Data normalization and visualization (visualizing_cell_cycle_effect_and_normalization.R) <br />
3- Data integration (data_integration.R) <br />
4- Performing diffderential expression between the four samples (differential_expression_between_samples.R) <br />
5- Clustering and cell type annotation (clustering_and_cell_type_annotation.R) <br />
6- Investigating transcription factors activities differences between cell types (TFs_activities_between_cell_types_dorothea.R) <br />
7- Investigating transcription factors activities differences between the four samples (TFs_activities_between_samples_dorothea.R) <br />
8- Trajectory analysis (TA_with_monocle.R) <br />
9- Investigating pathways activities differences between cell types (investigating_pathways_in_annotated_cell_types_progney.R) <br />
10- Investigating pathways activities differences between the four samples (investigating_pathways_in_samples_progney.R) <br />
11- scWGCNA analysis between MuSCs in the late differentiation stage and_ fibroblasts (scwgcna_analysis_between_late_differentiation_and_fibro.R) <br />
12- Investigating transcription factors and opathways activities in fibroblasts emerging from different samples (fibroblasts_analysis.R) <br />

-----
## Contributors:
This work was done by Amira A. Alakhdar, Sruthi Sivakumar, Rylee M. Kopchak, Allison N. Hunter, Fabrisia Ambrosio, Newell R. Washburn



