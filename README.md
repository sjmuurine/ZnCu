Contents

This repository contains supplementary code, data-files and scripts from the pig growth promoter study by S. Johanna Muurinen, Jacob Richert, Carmen Wickware, Brian Richert and Timothy Johnson

Included files
1. Sequence clustering shell scripts
    * Scripts that were run on the Purdue University ITaP Research Computing clusters to produce the primary data from Illumina read files
    * A batch file shell script (ZnCu.batch.txt)
        * Sequences are available at xxxx.xxx.xxx
        * SILVA-based bacterial reference alignment and mothur-formatted version of the RDP training set (v.9) can be downloaded at: https://www.mothur.org/wiki/MiSeq_SOP#OTUs
2. Data files: 
    * Rarefied and subsampled OTU table (Phylo.tx.1.pick.1.subsample.shared_correctedID.txt)
    * TSS normalized OTU table (Phylo.tx.1.pick.shared_correctedID.txt)
    * Metadata for OTU tables (ZnCu_metadata.txt)
    * Taxonomy data (Phylo.cons.taxonomy)
    * Wafergen ARG qPCR array results (6 files):
        * IIJohnson_Chip1_CorrectedNames_sorted.txt
        * IIJohnsonChip2_CorrectedNames_sorted.txt
        * IIJohnsonChip3_CorrectedNames_sorted.txt
        * JohnsonChip4_CorrectedNames_sorted.txt
        * IIJohnsonChip5rerun_sorted.txt
        * IIJohnsonChip6_CorrectedNames_sorted.txt
    * Metadata for ARGs (Wafergen_metadata.txt)
    * ARG annotation (Primerset2) (Primerset2_0_tnpAs_fixed.txt)
3. Data analysis scripts and Rdata that use the mentioned data files
    * Script for data analysis in R (Data_analysisARGsAND16SOTUs_script-f.R)
    * R data file that can be loaded into R (Rdata_ARGsAND16SOTUs.RData)
