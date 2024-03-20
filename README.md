# GraphTucker

# An example of running GraphTucker on a mouse embryo Stereo-seq dataset and visualizing the resulting components

This repository provides code we used to run GraphTucker on a MOSTA Stereo-seq mouse embryo (day 9.5) dataset<sup>1</sup> to obtain an imputed spatial gene expression tensor in Tucker form, and visualize the resulting spatial components. Please review the system and add-on requirements below. 

![Figure 1](./figure1.png)

Figure A4 shows the Tucker decomposition of a spatial gene expression tensor into its components, which consists of a core tensor and three factor matrices (one for each mode). This decomposition is graph-regularized by a Cartesian product graph consisting of two spatial chain graphs and a PPI network as shown in Figure B. Figure C  demonstrates how the spatial component tensor can be constructed from the Tucker components. 

The code provided here demonstrates how to run GraphTucker on a spatial gene expression tensor to obtain graph-regularized Tucker components, and then visualizing the spatial components found.

## System and add-on requirements

Code for this mouse embryo dataset was tested on a machine with Linux (Ubuntu 20.04.6 LTS) with the following specifications:

CPU: Intel® Xeon® E52687W v33.10GHz, 25M Cache
Cores: 20
Memory: 256GB

Note that the necessary memory/RAM for running GraphTucker on the MOSTA_9.5 dataset is significantly lower - we observed ~13GB of memory being used during the algorithm itself, so we recommend at least 16GB of memory.

We use the most up to date MATLAB version R2023b for all implementations. 

MATLAB add-on requirements:
The only add-on needed is the Tensor Toolbox for MATLAB, Version 3.6<sup>2</sup> (https://www.tensortoolbox.org/), which is included in the repository in the `GT_utils` folder.

## Step 0: Downloading and preprocessing data

In this part of the tutorial, we will walkthrough how to download and preprocess a Visium dataset from the 10x Genomics website. Click [https://www.10xgenomics.com/datasets?query=&page=2&configure%5BhitsPerPage%5D=50&configure%5BmaxValuesPerFacet%5D=1000&refinementList%5Bproduct.name%5D%5B0%5D=Spatial%20Gene%20Expression]{here} to go to the 10x Genomics website to see a list of all "Spatial Gene Expression" datasets. If it's not already selected, on the left panel under "Filter datasets", select the box to the left of "Spatial Gene Expression" under the "Products" section. For this tutorial, we will be using the [https://www.10xgenomics.com/datasets/human-breast-cancer-block-a-section-1-1-standard-1-0-0]{Human Breast Cancer (Block A Section 1)} dataset from Space Ranger v.1.1.0. 

Clicking on the dataset will take you to an overview of the data. Scroll down to the bottom and download two files: "Feature / bracode matrix (filtered)" and "Spatial Imaging data". Move these downloaded files to a dedicated folder for the dataset somewhere in the GraphTucker folder e.g. GraphTucker/data/BRCA1/. Extract both files twice (.gz -> .tar -> folder), and you should have two folders "filtered_feature_bc_matrix" and "spatial". You can delete/remove the .tar and .gz files. Your file/folder structure should look like this and include the following:

-- GraphTucker
   -- data
      -- BRCA1 
         -- filtered_feature_bc_matrix
            - barcodes.tsv.gz
            - features.tsv.gz
            - matrix.mtx.gz
         -- spatial
            - tissue_positions_list.csv

Once you've made sure your file structure is correct, we are going to run the R script Visium2Tensor.R to convert our raw data into MATLAB tensor format. First, create a new folder where you want the processed data to be output e.g. GraphTucker/processed_data. You can then run the script by opening the terminal/command prompt inside the GraphTucker folder and running the following command:

Rscript Visium2Tensor.R --input data --output processed_data
            
This will output two processed data files into the processed_data folder, BRCA1.mat and BRCA1_gene.csv. If you have multiple datasets in the data folder, this script will process them all and output each file into the processed_data folder.

Next, we will run the Convert2Sptensor.m script to convert our new .mat data file into a sparse tensor and save a significant amount of space. Open the Convert2Sptensor script in MATLAB, and change the "path_data_path" and "data_name" variables accordingly so the filepath correctly points to the .mat data file. After running this, you can delete the original file if you want to save space. The data is now ready to be used for running GraphTucker.

We also provide an example MOSTA Stereo-seq mouse embryo (day 9.5) dataset that has been preprocessed into tensor format. We also provide the corresponding gene name and mouse PPI for [BioGRID](https://thebiogrid.org/) that can be used for imputation/graph regularization. This data is provided in the repository in `data.zip`. Simply unzip it into its own folder to use in this code.

We will also be providing detailed preprocessing steps on how to convert raw Stereo-seq data into a spatial gene expression tensor (coming soon).

Here is a description of each data file and their general formatting: 

1. MOSTA_9.5_tensor.mat
   - spatial gene expression tensor, stored as a `int64`. Loaded in and used as a `tensor` object.
   - size (n_x, n_y, n_g) before pre-processing: (80, 107, 11153) 
   - size (n_x, n_y, n_g) after pre-processing: (80, 107, 9687) 
3. MOSTA_9.5_gene_names.csv
   - list of gene names corresponding to 11153 genes in the spatial gene expression tensor. This file contains the Ensembl ID and gene symbol for each gene.
5. MUS_PPI.mat
   - Version 4.4.226 Mus musculus PPI network from BioGRID that was used for testing.
   - Contains variables `MUS_GENE` which gives gene symbols and `MUS_BIOGRID` which is the adjacency matrix for the network (size: 11153 x 11153).
   - Can be replaced with any version of the network, provided the .mat file formatting is the same.
  
# Tutorial 1: Spatial component visualization for a Stereo-seq mouse embryo dataset

## Step 1: Running GraphTucker to obtain imputed spatial gene expression tensor and Tucker components

GraphTucker can be run on the mouse embryo dataset by running the `GraphTucker_tutorial1.m` script. The only necessary changes needed to run the script are to change the `work_path` and `data_path` variables to the paths the GraphTucker and MOSTA_9.5 folders were downloaded to, respectfully. A mouse PPI network is included with the data as well.

GraphTucker parameters are defaulted to rank=(64,64,64), lambda=0.1, with 1000 max iterations. Runtime for the MOSTA_9.5 dataset on our machine takes ~6-7 hours. The max iterations can be lowered to significantly reduce runtime, though the  approximation error may suffer as a result.

## Step 2: Visualizing spatial components output by GraphTucker

To visualize the spatial components output by GraphTucker, run `visualize_GraphTucker_spatial_components.m` script. Again, ensure the path variables are correct. This script will automatically save visualizations (as .pngs) of each spatial component into a distinct folder depending on the rank and lambda used. We include the spatial component images for MOSTA_9.5 obtained from running GraphTucker with the default parameters. Here are three example spatial components that show how they are visualized:

**Spatial component 1**
![Spatial component 1](./vis/MOSTA_9.5/rank=64-64-64/lambda=0.1/g_comp=1.png)


**Spatial component 15**
![Spatial component 15](./vis/MOSTA_9.5/rank=64-64-64/lambda=0.1/g_comp=15.png)


**Spatial component 60**
![Spatial component 1](./vis/MOSTA_9.5/rank=64-64-64/lambda=0.1/g_comp=60.png)

# Tutorial 2: Spatial component analysis and clustering on a Human breast cancer dataset

## Step 1: Running GraphTucker to obtain imputed spatial gene expression tensor and Tucker components

Similar to Tutorial 1, GraphTucker can be run on the mouse embryo dataset by running the `GraphTucker_tutorial2.m` script. The only necessary changes needed to run the script are to change the `work_path` and `data_path` variables to the paths the GraphTucker and BRCA1 processed data folders. You will also need a PPI network file - we provide a human PPI network 'HSA_ppi.mat' in the data/ppi_networks folder. Copy or move this file to where the processed BRCA1 data is located before running.

For analysis on this dataset, we recommend running GraphTucker with rank=(30, 30, 20), lambda=1, with 1000 max iterations (these are the default parameters). For other datasets, you may have to tune the rank. Tuning the x/y rank will affect the resolution (higher rank = higher resolution spatial components). Since we are trying to analyze how GraphTucker can capture 20 known regions, we set the the gene rank to 20. This should be adjusted based on the dataset being used, especially if prior knowledge is known about gene modules for the particular tissue. 

After verifying the PPI network file has been placed in the BRCA1 data folder and the proper parameters have been selected, you can run this script. This should take 1-2 hours depending on your computer/server.

## Step 2: Analyzing results: finding spatial component matches and perform clustering

After GraphTucker finishes running, the output decomposition will be saved to 'res/BRCA1'. We also provide an example output file that can be directly used for analysis, which we will use here.

The script we'll be runnning first is 'spatial_component_analysis.m'. If you are using you're own dataset, make sure the path variables are set correctly so the results from Step 1 can be found and loaded in. Otherwise for this tutorial, no changes are needed and it can be run as is. This script requires a regional annotation specific to a dataset - we provide one that can be used for the BRCA1 dataset that's being used in this tutorial.

This script will load in the GraphTucker decomposition computed in Step 1, and find the best spatial component matches to the given region annotations. It also will check a variety of top n% core tensor entries to see whether only taking a select number of top interactions can help find spatial component matches. This is determined by selecting the n% with the lowest avg. ED of its matched spatial components. After the best core tensor sparsity has been computed, we output each matched spatial component with its corresponding AUC and ED to 'vis/DATASET_NAME/spatial_componentmatches' folder.

The second script we run for this analysis is 'spatial_component_clustering.m' for clustering analysis. This script is similar to the previous in terms of preprocessing the output data, but we do not check different n% of core tensor entries and instead default to using all entries since we observed best clustering ARI with this choice. After the decomposition is loaded and the spatial component tensor is constructed, we run k-means clustering on the spatial component tensor with k=num_regions (for this example this is k=20). We then calculate the ARI between this resulting clustering and the tissue annotations as ground truth labels. The original annotated regions are saved as an image with a specific color mapping, and some additional processing is done to find a best color matching between the original region colors and the k-means clustering output. This is to provide a better visualization for identifying which clusterings best correspond to original regions. The original regions and the color mapped clusterings are both saved to the 'vis/BRCA1/clustering/' folder.

## Final notes

All scripts have been commented throughout to provide some guidance at main steps. Regardless, if any questions or issues come up please feel free to leave comments or ask for help.

## References

1. Chen, A., Liao, S., Cheng, M., Ma, K., Wu, L., Lai, Y., ... & Wang, J. (2022). Spatiotemporal transcriptomic atlas of mouse organogenesis using DNA nanoball-patterned arrays. Cell, 185(10), 1777-1792. DOI: 10.1016/j.cell.2022.04.003
2. Brett W. Bader, Tamara G. Kolda and others, Tensor Toolbox for MATLAB, Version 3.6, www.tensortoolbox.org, September 28, 2023. 
3. Stark, C., Breitkreutz, B. J., Reguly, T., Boucher, L., Breitkreutz, A., & Tyers, M. (2006). BioGRID: a general repository for interaction datasets. Nucleic acids research, 34(Database issue), D535–D539. https://doi.org/10.1093/nar/gkj109
