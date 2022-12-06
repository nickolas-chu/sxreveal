# sxreveal
A toolkit for demultiplexing single nuclei RNA sequencing samples by sex 
Create a Seurat object from your data set, following tutorial steps until you have clusters (https://satijalab.org/seurat/articles/pbmc3k_tutorial.html#run-non-linear-dimensional-reduction-umap-tsne-)
Make sure the default assay is set to RNA and the active identities are set to the cluster numbers.

Use:

probabilities <- femaleProb(Seuratobject). 

probabilities_with_sex_assignment <- sexAssign(probabilities).

Seuratobject_with _sex_assignemnts <- quickAdd(probabilities_with_sex_assignment, seuratobject).

Probability cutoffs for male and female can be specified in sexAssign. Default is >= 0.8 female, <= 0.2 male. Remaining cells are labled 'soup'
