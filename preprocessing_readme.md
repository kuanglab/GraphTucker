## Step 0: Downloading and preprocessing data

# Downloading and processing a 10x Genomics Visium dataset

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

Next, we will run the Convert2Sptensor.m script in the 'processing' folder to convert our new .mat data file into a sparse tensor and save a significant amount of space. Open the Convert2Sptensor script  in MATLAB, and change the "path_data_path" and "data_name" variables accordingly so the filepath correctly points to the .mat data file. After running this, you can delete the original file if you want to save space. The data is now ready to be used.

# Downloading and processing a PPI network from BioGRID

Here, we will show how to download a PPI network from BioGRID and process it into a .mat file so that it can be used for GraphTucker. First go to the [BioGRID](https://downloads.thebiogrid.org/BioGRID/), and go to the 'Current-Release' to see the most up-to-date networks. At the time of writing this the most up to date version is 4.4.231. Click on the 'BIOGRID-ORGANISM-4.4.231.tab3.zip' to download it. Note that this zip file will contain networks for many different organisms, but here we will only focus on the Human PPI, which will be named BIOGRID-ORGANISM-Homo_sapiens-4.4.231.tab3.txt (or differently depending on what version you want). We recommend you move this file to wherever you are storing the raw Visium data, for example in 'GraphTucker/data/BRCA1'.

Next we will be running the 'process_ppi.py' script in the 'processing' folder to build the adjacency matrix for this PPI network and get a list of its gene names. Open the 'process_ppi.py' script and ensure that 'PPI_data_path' variable points to the this newly download PPI .txt file, the 'ppi_name' variable is named as you want (the output file will be saved using this name), the 'save_path' variable points to your desired output folder, and 'features_path' points to the features.tsv.gz file in the raw Visium data you will be using. By default, this script will build the PPI matrix based on the intersection of genes in the Visium data and PPI network.

After running this script, it will save the adjacency matrix and PPI gene names together (variables 'HSA_BIOGRID' and 'HSA_UNIQ_BIOGRID_GENE') in a single HSA_ppi.mat file (to wherever you directed the output). The next step is optional but highly recommended to save space, since you will notice this new .mat file is very large.

Open matlab and load the new HSA_PPI.mat. You want to convert the adjacency matrix to a sparse matrix which will significantly reduce the file size, which you can do by calling:

HSA_BIOGRID = sparse(HSA_BIOGRID);

You can then resave these variables using:

save('HSA_PPI.mat', 'HSA_BIOGRID', 'HSA_UNIQ_BIOGRID_GENE')

The PPI will then be overwritten with a much smaller size - note that this will save the new file to whatever directory MATLAB is in when you open it, so use 'pwd' in the MATLAB command line to check where this file was saved, and move it to the corresponding 'processed_data/BRCA1' folder.

Another important note is that if you are using a PPI processed in this way, please use the 'preprocessing/data_prep_alt.m' function for preprocessing the data as used in the 'GraphTucker_tutorial2.m' script, instead of the default 'data_prep_visium()'  (at line 28). Also verify that in 'data_prep_alt.m' the path at line 10 points correctly to your processed PPI.mat file, and the variables at lines 62-64 correspond to the variable names you save in the 'processing_ppi.py' script (e.g. 'HSA_BIOGRID' and 'HSA_UNIQ_BIOGRID_GENE').

# Example preprocessed data

We also provide an example MOSTA Stereo-seq mouse embryo (day 9.5) dataset and Visium Human Breast Cancer dataset that have been preprocessed into tensor format. We also provide the corresponding gene name and mouse/human PPI that can be used for imputation/graph regularization. This data is provided in the repository in `/processed_data/MOSTA_9.5.zip` and `/processed_data/BRCA1.zip`. Simply unzip them into their own folder to use in this code.

We will also be providing detailed preprocessing steps on how to convert raw Stereo-seq data into a spatial gene expression tensor (coming soon).

Here is a description of each data file and their general formatting (same structure for BRCA1): 

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
   - Contains valid spot indices to indicate non-background spots
   - This file is strictly used for visualization purposes, and will be saved when running GraphTucker_tutorial1.m
  
The BRCA1 contains an additional file 'BRCA1_val_indices.mat' that saves the a set of tissue indices that are used for the spatial component analysis, similar to the '_val_subs.mat' file. This file will be automatically generated and saved when the 'data_prep_visium()' function is called.
