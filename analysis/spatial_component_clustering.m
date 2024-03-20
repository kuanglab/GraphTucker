clear all
addpath('../GT_utils/')
rng('default')

vis_path = '../vis/BRCA1'; % path to where images will be saved to

% check these variables to see if results will be loaded correctly
dataset = 'BRCA1';
num_regions = 20;
rank_list = [30 30 num_regions];
rank_str = [num2str(rank_list(1)), '-', num2str(rank_list(2)), '-', num2str(rank_list(3))];
lambda = 1;

% load data
load(['../res/', dataset, '/GT_', dataset, '_rank=', rank_str ,'_lambda=', num2str(lambda) ,'_EXAMPLE.mat']);
G = permute(G, [3 2 1]); % permute core tensor for visualization
A_x = A_set{1}{2};
A_y = A_set{1}{1};
n = [size(A_x, 1), size(A_y, 1)];

% coloring for 20 regions
cmap = ['1f77b4', 'aec7e8', 'ff7f0e', 'ffbb78', '2ca02c', '98df8a', 'd62728', ...
    'ff9896', '9467bd', 'c5b0d5', '8c564b', 'c49c94', 'e377c2', 'f7b6d2', ...
    '7f7f7f', 'c7c7c7', 'bcbd22', 'dbdb8d', '17becf', '9edae5'];
cmap = hex2rgb(cmap);

% construct spatial component tensor
X = ttm(G, {A_set{1}{2} A_set{1}{1}}, [1 2]); 

G_xy_mat = tenmat(X, 3, [1 2]); % compute spatial component tensor
G_xy = G_xy_mat.data;
G_xy = G_xy';

% load tissue annotation file for dataset
load(['../aux_data/tissue_annotations/', dataset, '_annotations.mat'])
annotations(:,1:2) = annotations(:,1:2) + 1;

% load subs file which contains indices for valid spots
load(['../processed_data/', dataset, '/', dataset, '_val_indices.mat'])

% arrange annotation onto spots correctly
annomat = zeros(size(A_x, 1), size(A_y, 1)) - 1;
for i = 1:size(annotations, 1)
    annomat(annotations(i,1), annotations(i,2)) = annotations(i,3);
end

% remove boundary spots so annotation matches size of (A_x, A_y)
annomat = annomat(ind1, ind2);

% adjusting annotation for better image (remove background)
anno_img = annomat;
anno_img(annomat == -1) = nan;
anno_img = anno_img * 100;

figure
colormap(cmap)
h = imagesc(anno_img);
set(h, 'AlphaData', ~isnan(anno_img));
set(gcf, 'Position', [100, 100, 300, 300]);
axis image off

img_path = ['../vis/', dataset, '/clustering/ground_truth_clustering.png'];
saveas(gcf, img_path);

% vectorize labels
cluster_label = annomat(:);
label_list = cluster_label(:);

% find non-background indices
val_idx = (label_list ~= -1);

% select non-background spots from spatial component tensor
G_xy = G_xy(val_idx,:);
G_xy = normr(G_xy); % normalize

% select non-background stots from cluster labelings
label_list = label_list(val_idx); %

% run k-means with k=num_regions, 5 replicates
idx = kmeans(G_xy, num_regions, 'Replicates', 5);

vis = zeros([size(A_x, 1), size(A_y, 1)]);
vis(val_idx) = idx;

% calculate ARI
ari = rand_index(idx, label_list, 'adjusted');

% some to do processing to figure out best matched region so cluster colors
% match up to original colorings
clustering = annomat;
clustering(val_idx) = idx - 1;

cost_matrix = zeros(num_regions);
for i = 0:num_regions-1
    for j = 0:num_regions-1
        anno_reg = (annomat == j);
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

% now change clustering labels based on their best matched regions
clustering_mapped = clustering;
for i = 0:num_regions-1
    clustering_mapped(clustering == i) = assignment(i+1);
end

% adjust background of clusterings for better visualization
img_data = clustering_mapped;
img_data(annomat == -1) = nan;

figure
h = imagesc(img_data);
colormap(cmap) 
set(h, 'AlphaData', ~isnan(img_data));
set(gcf, 'Position', [100, 100, 300, 300]);
axis image off

% save clustering results
savepath = ['../vis/', dataset, '/clustering/clustered_spatial_components.png'];
s = subtitle(['ARI=', num2str(round(ari,3))]);
s.FontName = 'TimesNewRoman';
saveas(gcf, savepath);       

close all