# Coalescence project 
Script used in the coalescence project:
> When the river meets the sea: bacterial community dynamics after coalescence.
> Xiu Jia, Torsten Schubert, Rick Beeloo, Aristeidis Litos, Swapnil Doijad, Pim van Helvoort, Theodor Sperlea, Matthias Labrenz, and Bas E. Dutilh

## Pipeline to analyze raw reads of full length 16S RNA gene
We used the Nanopore FLO-MIN114 flowcell (SQK-LSK114 kit) sequenced the 16S library. 
Please click here for [**raw reads basecalling, demultiplexing and taxonomic assignment scripts**](https://github.com/Jia-Xiu/coalescence_project/tree/main/16S_analysis).


## Downstream community analysis
### Table cleaning
For table cleanning, check this script [feature_table_clean_MPA_family_genus.Rmd](https://github.com/Jia-Xiu/coalescence_project/blob/main/com_analysis_scripts/feature_table_clean_MPA_family_genus.Rmd), which includes
- Filter-out non bacteria reads.
- Removing taxa with a total reads less than 3 and that occur in fewer than 3 samples.
- Save table at family and genus level.
- CLR transformation.
- Combine "feature/OTU" table, taxonomy table, and metadata as phyloseq objects.


### Rarefaction curves
To check the sequencing depth, I generated rarefaction curves by `vegan` pacakge in R, see script [rarefaction_family_genus.Rmd](https://github.com/Jia-Xiu/coalescence_project/blob/main/com_analysis_scripts/rarefaction_family_genus.Rmd)
