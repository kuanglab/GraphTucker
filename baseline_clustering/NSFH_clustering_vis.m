% addpath('FIST_utils/')
clear all
rng('default')

dataset = 'BRCA1'; % dataset to visualize, either 'hb' or 'brca1' (bcba1)

if strcmp(dataset, 'HB')
    M = readmatrix('res/NSFH_HB_spatial_factors_clustered.csv');
    k = 7;
    % coloring for 7 regions
    cmap = parula(k);
    n = [59 76];  % tissue size for Human Brain
    savepath = ['vis/NSFH_HB_clustering.png'];
elseif strcmp(dataset, 'BRCA1')
    M = readmatrix('res/NSFH_BRCA1_spatial_factors_clustered.csv');
    k = 20; % num regions
    % coloring for 20 regions
    cmap = ['1f77b4', 'aec7e8', 'ff7f0e', 'ffbb78', '2ca02c', '98df8a', 'd62728', 'ff9896', '9467bd', 'c5b0d5', '8c564b', 'c49c94', 'e377c2', 'f7b6d2', '7f7f7f', 'c7c7c7', 'bcbd22', 'dbdb8d', '17becf', '9edae5'];
    cmap = hex2rgb(cmap);
    cmap = map_colormap(cmap);
    n = [60 77]; % tissue size for Breast Cancer (BCBA1/BRCA1)
    savepath = ['vis/NSFH_BRCA1_clustering.png'];
end

anno_vec = M(:,6);
anno_mat = reshape(anno_vec, n);
anno_mat = anno_mat';

clustering_vec = M(:,7);
clustering = reshape(clustering_vec, n);
clustering = clustering';

clustering(anno_mat <= -1) = -1;

num_regions = k;
cost_matrix = zeros(num_regions);
for i = 0:num_regions-1
    for j = 0:num_regions-1
        anno_reg = (anno_mat == j);
        reg = (clustering == i);
        overlap = anno_reg + reg;
        count = length(overlap(overlap == 2));
        if count == 0
            count = 1e-6;
        end
        cost_matrix(i+1,j+1) = count;
    end
end

cost_matrix = 1 ./ cost_matrix;

[assignment, cost_matrix_] = munkres(cost_matrix);

clustering_mapped = clustering;
for i = 0:num_regions-1
    clustering_mapped(clustering == i) = assignment(i+1);
end

img_data = clustering_mapped;
img_data(anno_mat <= -1) = nan;
img_data(clustering_mapped <= -1) = nan;


figure
h = imagesc(img_data);
colormap(cmap)
set(h, 'AlphaData', ~isnan(img_data));
set(gcf, 'Position', [100, 100, 300, 300]);
axis image off
s.FontName = 'TimesNewRoman';
saveas(gcf, savepath);   

% custom mapping function to map correct colors for BRCA1/BCBA1 dataset
function cmap = map_colormap(cmap)
    cmap_cpy = cmap;
    cmap_cpy(12,:) = cmap(1,:);
    cmap_cpy(5,:) = cmap(2,:);
    cmap_cpy(4,:) = cmap(3,:);
    cmap_cpy(18,:) = cmap(4,:);
    cmap_cpy(9,:) = cmap(5,:);
    cmap_cpy(1,:) = cmap(6,:);
    cmap_cpy(8,:) = cmap(7,:);
    cmap_cpy(11,:) = cmap(8,:);
    cmap_cpy(10,:) = cmap(9,:);
    cmap_cpy(6,:) = cmap(10,:);
    cmap_cpy(2,:) = cmap(11,:);
    cmap_cpy(19,:) = cmap(12,:);
    cmap_cpy(20,:) = cmap(13,:);
    cmap_cpy(3,:) = cmap(14,:);
    cmap_cpy(17,:) = cmap(15,:);
    cmap_cpy(16,:) = cmap(16,:);
    cmap_cpy(15,:) = cmap(17,:);
    cmap_cpy(13,:) = cmap(18,:);
    cmap_cpy(14,:) = cmap(19,:);
    cmap_cpy(7,:) = cmap(20,:);
    cmap = cmap_cpy;
end