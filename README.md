# dynDysRegNet
This repository belongs to the Guided Research Project report 'Studying dynamic gene regulatory networks in single cells with DysRegNet' by Anna Ketteler (anna.ketteler@tum.de). For details on the method, please consult the report. Do not hesitate to contact me for questions about the code if answers cannot be taken from the report.

The project uses the tool DysRegNet:
DysRegNet: Patient-specific and confounder-aware dysregulated network inference
Olga Lazareva, Zakaria Louadi, Johannes Kersting, Jan Baumbach, David B. Blumenthal, Markus List
bioRxiv 2022.04.29.490015; doi: https://doi.org/10.1101/2022.04.29.490015

...and diffusion pseudo-time:
Haghverdi, L., Büttner, M., Wolf, F. et al. Diffusion pseudotime robustly reconstructs lineage branching. Nat Methods 13, 845–848 (2016). https://doi.org/10.1038/nmeth.3971

...and scANVI:
A Python library for probabilistic analysis of single-cell omics data
Adam Gayoso, Romain Lopez, Galen Xing, Pierre Boyeau, Valeh Valiollah Pour Amiri, Justin Hong, Katherine Wu, Michael Jayasuriya, Edouard Mehlman, Maxime Langevin, Yining Liu, Jules Samaran, Gabriel Misrachi, Achille Nazaret, Oscar Clivio, Chenling Xu, Tal Ashuach, Mariano Gabitto, Mohammad Lotfollahi, Valentine Svensson, Eduardo da Veiga Beltrame, Vitalii Kleshchevnikov, Carlos Talavera-López, Lior Pachter, Fabian J. Theis, Aaron Streets, Michael I. Jordan, Jeffrey Regier & Nir Yosef
Nature Biotechnology 2022 Feb 07. doi: 10.1038/s41587-021-01206-w

...and also includes harmony in the code, which is not used by default:
Korsunsky, I., Millard, N., Fan, J. et al. Fast, sensitive and accurate integration of single-cell data with Harmony. Nat Methods 16, 1289–1296 (2019). https://doi.org/10.1038/s41592-019-0619-0

...and GRNBoost2:
Thomas Moerman, Sara Aibar Santos, Carmen Bravo González-Blas, Jaak Simm, Yves Moreau, Jan Aerts, Stein Aerts, GRNBoost2 and Arboreto: efficient and scalable inference of gene regulatory networks, Bioinformatics, Volume 35, Issue 12, June 2019, Pages 2159–2161, https://doi.org/10.1093/bioinformatics/bty916

...mostly in the scanpy framework:
Wolf, F., Angerer, P. & Theis, F. SCANPY: large-scale single-cell gene expression data analysis. Genome Biol 19, 15 (2018). https://doi.org/10.1186/s13059-017-1382-0

### How to run the code:
- for the full pipeline, you need the scripts 'runner_dataset_sub.py', 'infer_template_network.py', and 'dysregnet_script.py'. To identify which metacell type (ann_level_2) are needed for template network inference, you further need to analyze which metacell types have the lowest mean DPT values. For this purpose (and visualizations), setting up a notebook is great.
- first, adjust the data paths in the three scripts. For now, we have one data_path and one out_dir which are two different directories. The paths are set according to the directory structure of the workstation of Markus' chair.
- after putting the data (e.g. HLCA core in h5ad format, and the list of human known transcription factors) into the data_dir, run the runner_dataset_sub.py to obtain the metacells data in different versions:
        hlca_core_metacell_<dataset_name>.h5ad: contains the metacells from the single dataset with dataset_name. The 'unclustered metacell' is already removed.
        metacell_all.h5ad: contains the merged metacells from all datasets. Here, also the uninformative metacells according to Shannon entropy are removed (by default, the ann_level_2 annotation is used for this filtering, but you can change this. Just be careful that on finer annotation levels, we inherit annotations from coarser annotation levels, as in ann_level_3-5, there are 'None' annotations. The entropy filtering will then be performed according to the inherited annotation levels.)
        metacells_all_batch_corrected.h5ad: contains the batch corrected embedding, probably in the adata.obsm section of your anndata object, depending on the correction methods you use. Watch out: scANVI, for example, stores only the genes it used for correction in the metacells_all_batch_corrected.h5ad. Therefore, we do not use it for template network inference and dysregnet.
        metacell_all_dpt.h5ad: contains the DPT value annotation in the adata.obsannotations.
- so now, you have the plain metacells data (metacells_all.h5ad), and you have metacells_all_dpt.h5ad, which contains the DPT values inferred from a batch corrected embedding (of course, you technically don't have to use a batch corrected embedding to infer pseudo-time). Next, you have to identify which metacell types concerning the ann_level_2 annotation you want to use as controls, i.e., as 'immature' metacell cohort, because these you use for template network inference. Analyze this in a separate notebook. Then, put the according ann_level_2 categories in the infer_template_network.py in the network_metacell_types list variable and start the script.
- Lastly, specify in dysregnet_script.py which covaraites you want to include to the analysis, and run the script.

### Ideas for future studies:
- regarding the tests where we include DPT pseudo time as a covariate, and we obtain much more dysregulated edges than without the covariate: either, the models become better in predicting the controls (overfitting), or, they indeed become worse in predicting the cases. To distinguish, check the RMSE of the models; from this, you can see, whether the predictions on the controls become better (--> overfitting) or not (--> no overfitting).
- do not only check the number of dysregulated edges per metacell, but also computethe overlap of dysregulated edges between the metacells --> are always the same edges dysregulated, or lots of differences? On edges that are frequently dysregulated, we could perform gene set enrichment analysis.
- Use multiple control/case splitsand 'rotate' through the dataset to see how dysregulations change through comparisons of different metacell cohorts. I.e., partition the metacells into controls and cases at different DPT thresholds.
- Test batch correction methods that apply batch correction directly to the count data, and apply correction as early as possible (on the raw cells of the HLCA core).
- incorporate pseudo-time branchings into the analysis (differentiate not only by absolute pseudo-time value in the final analysis, but also differentiate between branches that a metacell is located on).
- run the code on different datasets, e.g. parts of the HLCA full. Adjustments will be necessary, because the annotations have slightly different names.
