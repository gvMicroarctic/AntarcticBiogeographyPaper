# AntarcticBiogeographyPaper

The scripts used to perform all the analyses in our paper, 'Biogeographic survey of soil bacterial communities across Antarctica,' can be found on this GitHub page. The 17 datasets were first analyzed separately. The R code for these separate analyses is available in the folder '/single_datasets_analyses.' The R objects resulting from the analysis of these 17 datasets were then imported into R combined into one phyloseq object (combined_dataset.R). The bioclimatic variables, elevation, and distance from coastline values were extracted using 'extract_bioclim_bioregion_information.R.' The following scripts were used to create the paper's figures and tables:

- 'combined_dataset.R' for tables S5 and S6

- 'alpha_diversity.R' for figures 2, S3, S4, and S5

- 'pcoa.R' for figures 3A-B and S10

- 'comparison_by_bioregion.R' for figures 3C and S8

- 'dbRDA_variation_partitioning.R' for figures 4 and S8, as well as table S8

- 'tangleglam.R' for figures 5 and S9

- 'island_analysis.R' for figures 6A-C

- 'mainland_analysis.R' for figures 6D-E

- 'random_forest_indicator.R' for figures 7 and S13, and table S9

- 'bioclimatic_data_representation.R' for figure S6 and table S7

- 'genera_bioclimatic_correlation.R' for figure S7

RData objects and input files necessary to run the R code can be found in the folder '/data'.

For any questions, feel free to email me at gilda.varliero@wsl.ch.
