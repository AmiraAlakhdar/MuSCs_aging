This code was used to analyse the scRNA data for the paper titled:

# Age‐Related ECM Stiffness Mediates TRAIL Activation in Muscle Stem Cell Differentiation

## Summary
The stiffening of the extracellular matrix (ECM) with age hinders muscle regeneration by causing intrinsic muscle stem cell (MuSC) dysfunction through a poorly understood mechanism. Here, young and aged MuSCs were seeded onto substrates engineered to mimic a soft and stiff ECM microenvironment to study those molecular changes using single-cell RNA-sequencing (scRNA). Data revealed the presence of a branching point within the trajectory leading to the emergence of an age-related fibroblastic population characterized by activation of the TNF-related apoptosis-inducing ligand (TRAIL) pathway, which was significantly activated in aged cells cultured on stiff substrates. Next, using the collagen cross-linking inhibitor β-aminopropionitrile (BAPN) in vivo, we elucidated stiffness changes on TRAIL downstream apoptotic targets (caspase 8 and caspase 3) using immunostaining. TRAIL activity was significantly inhibited by BAPN in aged animals, indicating a complex mechanism of age-related declines in muscle function through inflammatory and apoptotic mediators.
***
## What to expect

<mark>Marked text</mark>
Here, I included the following R scripts: <br />
1- Loading data and quality control (<mark>`loading_count_data_and_quality_control.R`</mark>) <br />
2- Data normalization and visualization (<mark>`visualizing_cell_cycle_effect_and_normalization.R`</mark>) <br />
3- Data integration (<mark>`data_integration.R`</mark>) <br />
4- Performing diffderential expression between the four samples (<mark>`differential_expression_between_samples.R`</mark>) <br />
5- Clustering and cell type annotation (<mark>`clustering_and_cell_type_annotation.R`</mark>) <br />
6- Investigating transcription factors activities differences between cell types (<mark>`TFs_activities_between_cell_types_dorothea.R`</mark>) <br />
7- Investigating transcription factors activities differences between the four samples (<mark>`TFs_activities_between_samples_dorothea.R`</mark>) <br />
8- Trajectory analysis (<mark>`TA_with_monocle.R`</mark>) <br />
9- Investigating pathways activities differences between cell types (<mark>`investigating_pathways_in_annotated_cell_types_progney.R`</mark>) <br />
10- Investigating pathways activities differences between the four samples (<mark>`investigating_pathways_in_samples_progney.R`</mark>) <br />
11- scWGCNA analysis between MuSCs in the late differentiation stage and_ fibroblasts (<mark>`scwgcna_analysis_between_late_differentiation_and_fibro.R`</mark>) <br />
12- Investigating transcription factors and opathways activities in fibroblasts emerging from different samples (<mark>`fibroblasts_analysis.R`</mark>) <br />

-----
## Contributors:
This work was done by: <br />
Amira A. Alakhdar, Sruthi Sivakumar, Rylee M. Kopchak, Allison N. Hunter, Fabrisia Ambrosio, Newell R. Washburn



