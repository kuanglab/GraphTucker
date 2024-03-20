## Step 0: Downloading and preprocessing data

In this part of the tutorial, we will walkthrough how to download and preprocess a Visium dataset from the 10x Genomics website. Click [https://www.10xgenomics.com/datasets?query=&page=2&configure%5BhitsPerPage%5D=50&configure%5BmaxValuesPerFacet%5D=1000&refinementList%5Bproduct.name%5D%5B0%5D=Spatial%20Gene%20Expression]{here} to go to the 10x Genomics website to see a list of all "Spatial Gene Expression" datasets. If it's not already selected, on the left panel under "Filter datasets", select the box to the left of "Spatial Gene Expression" under the "Products" section. For this tutorial, we will be using the [https://www.10xgenomics.com/datasets/human-breast-cancer-block-a-section-1-1-standard-1-0-0]{Human Breast Cancer (Block A Section 1)} dataset from Space Ranger v.1.1.0. 

Clicking on the dataset will take you to an overview of the data. Scroll down to the bottom and download two files: "Feature / bracode matrix (filtered)" and "Spatial Imaging data". Move these downloaded files to a dedicated folder for the dataset somewhere in the GraphTucker folder e.g. GraphTucker/data/BRCA1/. Extract both files twice (.gz -> .tar -> folder), and you should have two folders "filtered_feature_bc_matrix" and "spatial". You can delete/remove the .tar and .gz files. Your file/folder structure should look like this and include the following:

	
        . <data-folder>
        ├── ...
        ├── <tissue-folder>
        │   ├── filtered_feature_bc_matrix
        │   │   ├── barcodes.tsv.gz
        │   │   ├── features.tsv.gz
        │   │   └── matrix.mtx.gz
        │   ├── spatial
        │   │   └── tissue_positions_list.csv
        └── ...


Once you've made sure your file structure is correct, we are going to run the R script Visium2Tensor.R in the 'processing' folder to convert the raw data into MATLAB tensor format. You can then run the script by opening the terminal/command prompt inside the GraphTucker folder and running the following command:

Rscript Visium2Tensor.R --input data --output processed_data
            
This will output two processed data files into the processed_data folder, BRCA1.mat and BRCA1_gene.csv. If you have multiple datasets in the data folder, this script will process them all and output each file into the processed_data folder.

Next, we will run the Convert2Sptensor.m script in the 'processing' folder to convert our new .mat data file into a sparse tensor and save a significant amount of space. Open the Convert2Sptensor script  in MATLAB, and change the "path_data_path" and "data_name" variables accordingly so the filepath correctly points to the .mat data file. After running this, you can delete the original file if you want to save space. The data is now ready to be used for running GraphTucker.

We also provide an example MOSTA Stereo-seq mouse embryo (day 9.5) dataset that has been preprocessed into tensor format. We also provide the corresponding gene name and mouse PPI for [BioGRID](https://thebiogrid.org/) that can be used for imputation/graph regularization. This data is provided in the repository in `/processed_data/MOSTA_9.5.zip`. Simply unzip it into its own folder to use in this code.

We will also be providing detailed preprocessing steps on how to convert raw Stereo-seq data into a spatial gene expression tensor (coming soon).

Here is a description of each data file and their general formatting: 

1. MOSTA_9.5_tensor.mat
   - spatial gene expression tensor, stored as a `int64`. Loaded in and used as a `tensor` object.
   - size (n_x, n_y, n_g) before pre-processing: (80, 107, 11153) 
   - size (n_x, n_y, n_g) after pre-processing: (80, 107, 9687) 
3. MOSTA_9.5_gene_names.csv
   - list of gene names corresponding to 11153 genes in the spatial gene expression tensor. This file contains the Ensembl ID and gene symbol for each gene.
3. MUS_PPI.mat
   - Version 4.4.226 Mus musculus PPI network from BioGRID that was used for testing.
   - Contains variables `MUS_GENE` which gives gene symbols and `MUS_BIOGRID` which is the adjacency matrix for the network (size: 11153 x 11153).
   - Can be replaced with any version of the network, provided the .mat file formatting is the same.
4. MOSTA_9.5_val_subs.mat
   - Contains valid spot indices to indcate non-background spots
   - This file is strictly used for visualization purposes, and will be saved when running GraphTucker_tutorial1.m
