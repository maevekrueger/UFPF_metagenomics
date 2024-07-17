# UFPF_metagenomics

## INCLUDES ANALYSIS OF THE FOLLOWING DATASETS: 

- **UFPF dataset** (PD, IBD, Control patient cohorts collected from clinics at the University of Florida)  (https://www.ncbi.nlm.nih.gov/bioproject/1096686)

- **HMP2 IBD dataset** (IBD and non-IBD patients that were part of the Human Microbiome Project 2- Metadata can be downloaded at https://www.ibdmdb.org/results (Date Updated: 2018-08-20) and merged tables used can be obtained at https://www.ibdmdb.org/downloads/html/products_MGX_2017-08-12.html)

- **Wallen PD dataset** (also referred to as Payami, contains PD and control patients that were collected as part of the publication by Wallen et al. 2022 (https://doi.org/10.1038/s41467-022-34667-x) and can be accessed at https://zenodo.org/doi/10.5281/zenodo.7246184 as specified in their publication) 



### Source Data used in our analysis/figure generation is available on Zenodo along with a readme file explaining each table: DOI: https://zenodo.org/doi/10.5281/zenodo.10912505


## ANALYSIS: Detailed methods for how analyses were performed for these datasets is outlined in our manuscript. All scripts utilized are publicly available on this project's Github page. 

### Navigating the UFPF_metagenomics Github repository:

### Bioinformatic Processing of Sequences
- Contains the code used to perform QC on our metagenomic sequences. Raw sequences (with human sequences removed) can be found on NCBI Sequence Read Archive (SRA) under BioProject PRJNA1096686. (https://www.ncbi.nlm.nih.gov/bioproject/1096686) 

### UFPF
- This folder contains files and data associated with the UFPF dataset of PD, IBD, and control samples collected at the University of Florida.
  - **Metadata.rds:** Contains the demographic and metadata of the UFPF samples.
  - **Metadata All Samples.csv:** Excel-friendly version of demographic and metadata of the UFPF samples.
  - **Metaphlan output**
    - **Metaphlan_all_combined.tsv:** Output from running Metaphlan4.0 the first time.
      - **Rel_ab_cleaned.rds:** Cleaned up relative abundance output of Metaphlan_all_combined.tsv used in downstream analysis.
    - **Metaphlan_unknown_all_combined.tsv:** Output from running Metaphlan4.0 a second time, now including the -unclassified-estimation flag.
      - **Rel_ab_unkn_cleaned.rds:** Cleaned up relative abundance output with unknown/unclassified microbes from the output metaphlan_unknown_all_combined.tsv used in downstream analysis.
    - **Counts:** Contains the raw counts calculated using the relative abundance output from Metaphlan and the total read count per sample.
      - **Counts WITHOUT unclassified.rds:** Raw counts calculated from the relative abundance values that did not include the unclassified-estimation when running Metaphlan.
      - **Counts w Unclassified.rds:** Raw counts calculated from the relative abundance values when Metaphlan was run with the unclassified-estimation added.
      - **Genus raw counts.rds:** Filtered Counts w Unclassified.rds to only the genus level.
      - **Phylum raw counts.rds:** Filtered Counts w Unclassified.rds to only the phyla level.
      - **Species raw counts.rds:** Filtered Counts w Unclassified.rds to only the species level.
  - **Humann output**
    - **Pathway_abundance_all.tsv:** Contains the MetaCyc pathway output from running HUMAnN3.5.
    - **Consolidated pathways.RData:** Due to the size of the original output Pathway_abundance_all.tsv, pathways were consolidated and examined at the community level.
  - **ANCOMBC2**
    - **UFPF ancombc2 genus.rds:** Genus-level output from running ANCOMBC2 on taxonomic data.
    - **UFPF ancombc2 pathways.rds:** Pathway output from running ANCOMBC2 on functional (pathway) data.
    - **UFPF ancombc2 species.rds:** Species-level output from running ANCOMBC2 on taxonomic data.
    - **Bias corrected abund genus.rds:** Bias-corrected genus-level abundances extracted from ANCOMBC2 output run at the genus level. Used downstream in the module analysis.
    - **Bias corrected abund pathways.rds:** Bias-corrected pathway abundances extracted from ANCOMBC2 output of pathway data. Used downstream in the module analysis.
    - **Bias corrected abund species.rds:** Bias-corrected species-level abundances extracted from ANCOMBC2 output run at the species level. Used downstream in the module analysis.

### UFPF_scripts
- This folder contains the scripts used to perform our analysis. These scripts rely on files and data present in the UFPF folder.
  - **Demographics metadata table.R:** Generates a metadata table including demographic information and questionnaire results from our UFPF PD, IBD, and control subjects.
  - **Consolidating pathway data.R:** Consolidates pathways from the original HUMAnN3.5 output Pathway_abundance_all.tsv for examination at the community level.
  - **Getting individual counts.R:** Computes raw counts from the relative abundance output from Metaphlan4.0.
  - **Phylum genus species raw counts.R:** Filters raw counts from the relative abundance output from Metaphlan4.0 (done in getting individual counts.R) into separate files for species, genus, and phylum levels.
  - **PCoA Aitchison distances.R:** Computes and visualizes Aitchison distances using species-level count data in a principal coordinate analysis.
  - **Creating taxa phyloseq object + ancombc2.R:** Creates phyloseq objects at the species and genus levels used as input for running ANCOMBC2 differential abundance analysis.
  - **Creating pathway phyloseq object + ancombc2.R:** Creates a pathways phyloseq object used as input for running ANCOMBC2 differential abundance analysis.
  - **Getting bias corrected abundances ancombc2.R:** Extracts bias-corrected abundances from the various ANCOMBC2 runs to utilize the species-level and pathway data in downstream module analysis.
  - **Creating tables of ANCOMBC2 output.R:** Displays ANCOMBC2 output in table format.
  - **Pathway venn diagrams.R:** Creates Venn diagrams to display the PD and IBD-associated enriched and depleted pathways after ANCOMBC2 analysis from the script creating pathway phyloseq object + ancombc2.R.

### HMP2_Payami
- This folder contains files and data associated with the HMP2 IBD dataset and the Wallen PD dataset. Notable files include:

  #### HMP2 IBD-related
  - **Taxonomic_profiles.tsv.gz:** Contains the original HMP2 IBD taxonomic data produced in Metaphlan2.0, downloaded from https://www.ibdmdb.org/results.
  - **Pathabundance.tsv.gz:** Contains the original HMP2 IBD functional (pathway) data produced in HUMAnN, downloaded from https://www.ibdmdb.org/results.
  - **Hmp2_metadata_2018-08-20.csv:** Contains the original HMP2 IBD metadata, downloaded from https://www.ibdmdb.org/results.
  - **IBD Metadata Age Filtered.rds:** Contains the HM2 IBD metadata file utilized in downstream analysis. This was age-filtered to only include subjects that were 40 years and older.
  - **HMP2 Pathways Age Filtered.rds:** Contains pathway output age-filtered to match the samples included in the HMP2 IBD metadata file.
  - **HMP2 Consolidated and Age Filtered Pathways.rds:** Contains pathway output age-filtered and consolidated to perform pathway analysis at the community level.
  - **HMP2 IBD Age Filtered All Levels relab.rds:** Contains taxonomic relative abundance values filtered to only include subjects that match those in the age-filtered metadata file.
  - **HMP2 IBD Age Filtered All Levels Counts.rds:** Contains taxonomic raw counts calculated using the relative abundance values from Metaphlan and total sample read counts, filtered to only include subjects that match those in the age-filtered metadata file.

  #### Wallen PD-related
  - **PAYAMI DATA - Source_Data_24Oct2022.xlsx:** Original source data from Wallen et al. 2022 accessed from their public Zenodo repository https://zenodo.org/doi/10.5281/zenodo.7246184.
  - **Wallen Counts.rds:** Count data from the Wallen PD dataset.
  - **Wallen PD Metadata.rds:** Demographic and metadata from the Wallen PD dataset.

  #### Phyloseq Objects
  - **Wallen PD species phyloseq object.rds:** Phyloseq object created from species-level Metaphlan output to later be used in ANCOMBC2 analysis.
  - **Wallen PD genus phyloseq object.rds:** Phyloseq object created from genus-level Metaphlan output to later be used in ANCOMBC2 analysis.
  - **Wallen PD pathways phyloseq object.rds:** Phyloseq object created from pathway abundance HUMAnN output to later be used in ANCOMBC2 analysis.
  - **HMP2 IBD species phyloseq object.rds:** Phyloseq object created from species-level Metaphlan output to later be used in ANCOMBC2 analysis.
  - **HMP2 IBD genus phyloseq object.rds:** Phyloseq object created from genus-level Metaphlan output to later be used in ANCOMBC2 analysis.
  - **HMP2 IBD pathways phyloseq object.rds:** Phyloseq object created from pathway abundance HUMAnN output to later be used in ANCOMBC2 analysis.

  #### ANCOMBC2
  - Contains output files from running ANCOMBC2 differential abundance analysis in Wallen PD and HMP2 IBD datasets respectively, at the species and genus levels and with pathway abundance data. Bias-corrected abundances were also extracted from these ANCOMBC2 runs to be used in downstream module analysis.

### HMP2_Payami_scripts
- This folder contains the scripts used to perform our analysis. These scripts rely on files and data present in the HMP2_Payami folder.
  - **Demographics tables.R:** Used to generate demographic tables for the HMP2 IBD and Wallen PD datasets.
  - **Cleaning metaphlan relab output Wallen HMP2.R:** Cleaning up the relative abundance output from Metaphlan for the Wallen PD and HMP2 IBD datasets prior to further analyses.
  - **Calculating counts + age filtering HMP2 dataset.R:** Computes counts from the relative abundance Metaphlan2 output and total read counts found in the metadata. The HMP2 IBD dataset was also filtered to only include samples of subjects that were 40 years or older.
  - **CLR transformation on species counts.R:** Performs a CLR transformation on species-level count data for PCoA creation.
  - **PCoA Species CLR Counts.R:** Generates PCoAs for the Wallen PD and HMP2 IBD datasets.
  - **Consolidating HMP2 pathways for phyloseq.R:** Consolidates pathway data for community-level analysis to simplify further analyses.
  - **Creating HMP2 pathway phyloseq object.R:** Creates pathway phyloseq object for later ANCOMBC2 differential abundance analysis of the HMP2 IBD dataset.
  - **Creating Wallen pathway phyloseq object.R:** Creates pathway phyloseq object for later ANCOMBC2 differential abundance analysis of the Wallen PD dataset.
  - **Creating species level phyloseq objects.R:** Creates phyloseq objects of taxa at the species level for downstream ANCOMBC2 differential abundance analysis for HMP2 IBD and Wallen PD datasets.
  - **Creating genus level phyloseq objects.R:** Creates phyloseq objects of taxa at the genus level for downstream ANCOMBC2 differential abundance analysis for HMP2 IBD and Wallen PD datasets.
  - **HMP2 Wallen ancombc2.R:** Performs ANCOMBC2 on the various phyloseq objects created for the Wallen PD and HMP2 IBD datasets.
  - **Taxa and pathway venn diagrams.R:** Visualizes the HMP2 IBD-associated features and Wallen PD-associated features as determined from ANCOMBC2 analysis of these respective datasets using Venn diagrams.
  - **Getting bias corrected log abundances from ancombc2.R:** Extracts the bias-corrected abundances from the ANCOMBC2 analyses of the Wallen PD and HMP2 IBD datasets to be utilized in the module analysis.

### Module Analysis
- The folder contains the scripts used to perform the module analysis utilized to compare all three datasets: UFPF, Wallen, and HMP2.
  - **Taxonomic module analysis.R:** Performs a module analysis of the species-level taxonomic data and relies on the metadata, bias-corrected species abundances from ANCOMBC2, and the ANCOMBC2 output file for each dataset.
