# Coalescence project 
Script used in the coalescence project:
> When the river meets the sea: bacterial community dynamics after coalescence.
> Xiu Jia, Torsten Schubert, Rick Beeloo, Aristeidis Litos, Swapnil Doijad, Pim van Helvoort, Theodor Sperlea, Matthias Labrenz, and Bas E. Dutilh


## Pipeline to analyze raw reads of full length 16S RNA gene
The 16S amplicon library was sequenced by Nanopore FLO-MIN114 flowcell with SQK-LSK114 kit. 
Please click here for scripts regarding [**raw reads processing**](https://github.com/Jia-Xiu/coalescence_project/tree/main/16S_analysis).
- Basecalling
- Demultiplexing
- Taxonomic assignment 


## Downstream community analysis
### Combine species table, taxonomic information and metadata as pyloseq object
Check this script [feature_table_phyloseq_preparation.Rmd](https://github.com/Jia-Xiu/coalescence_project/blob/main/com_analysis_scripts/feature_table_phyloseq_preparation.Rmd), which includes combine "species/OTU" table, taxonomy table, and metadata as phyloseq objects including the OUT_table in:
- raw reads
- CLR transformation
- relative abundance
> [!TIP]
> To view intermediate data and outputs generated at each step of the analysis, please refer to the accompanying HTML files produced with R Markdown (download and open by your browser, e.g. Chrome).


### Rarefaction curves
To check the sequencing depth, I generated rarefaction curves by `vegan` pacakge in R, see script [rarefaction_curves_species.Rmd](https://github.com/Jia-Xiu/coalescence_project/blob/main/com_analysis_scripts/rarefaction_curves_species.Rmd)

### Differential abundance analysis between freshwater and seawater source communities

I used [ALDEx2](https://www.bioconductor.org/packages/release/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.html) performed the differential abundance analysis between freshwater and seawater source communities, see script:
- [ALDEx2_differential_abundance_analysis_source_coms_species.Rmd](https://github.com/Jia-Xiu/coalescence_project/blob/main/com_analysis_scripts/ALDEx2_differential_abundance_analysis_source_coms_species.Rmd)
- [ALDEx2_differential_abundance_analysis_source_coms_genus.Rmd](https://github.com/Jia-Xiu/coalescence_project/blob/main/com_analysis_scripts/ALDEx2_differential_abundance_analysis_source_coms_genus.Rmd)


### Diversity analysis
To analyse the community diveristy changes after coalescence, see script [div_analysis_species_level_rclr.Rmd](https://github.com/Jia-Xiu/coalescence_project/blob/main/com_analysis_scripts/div_analysis_species_level_rclr.Rmd)

### Co-occurence network analysis
Co-occurence network analysis was performed using [SpeSpeNet](https://utrecht-university.shinyapps.io/SpeSpeNet_v1/) with the following settings: Normalization: CLR, Correlation method: Spearman.\n I did further plotting in R with a Fruchterman–Reingold layout, you can find the script [here](https://github.com/Jia-Xiu/coalescence_project/blob/main/com_analysis_scripts/Network_SpeSpeNet_rawreads.Rmd). 
