% perform k-means on spatial component tensor
rng('default')
dataset = 'HB';

if strcmp(dataset, 'HB')
    rank = 7;
    num_regions = 7;
    tissue = '151673';
    %coloring for 7 regions
    cmap = parula(num_regions);
elseif strcmp(dataset, 'BRCA1')
    rank = 20;
    num_regions = 20;
    tissue = 'BRCA1';
    % coloring for 20 regions
    cmap = ['1f77b4', 'aec7e8', 'ff7f0e', 'ffbb78', '2ca02c', '98df8a', 'd62728', 'ff9896', '9467bd', 'c5b0d5', '8c564b', 'c49c94', 'e377c2', 'f7b6d2', '7f7f7f', 'c7c7c7', 'bcbd22', 'dbdb8d', '17becf', '9edae5'];
    cmap = hex2rgb(cmap);        
end

% 
load(['res/', dataset, '/NMF_', dataset, '_rank=', num2str(rank) ,'.mat']);

load(['data/', dataset, '/', dataset, '_annotations.mat'])
annotations(:,1:2) = annotations(:,1:2) + 1;

load(['data/', dataset, '/', dataset, '_val_indices.mat'])

load(['res/', dataset, '/GT_', tissue,'_rank=10-10-', num2str(rank), '_lambda=0.mat']);

A_x = A_set{1}{2};
A_y = A_set{1}{1};
annomat = zeros(size(A_x, 1), size(A_y, 1)) - 1;
for i = 1:size(annotations, 1)
    annomat(annotations(i,1), annotations(i,2)) = annotations(i,3);
end

annomat = annomat(ind2, ind1);
% imagesc(annomat)
cluster_label = annomat(:);
label_list = cluster_label(:);

val_idx = (label_list ~= -1);

H_mat = H';

nx = size(A_x, 1);
ny = size(A_y, 1);

H_ = zeros(size(H_mat));
counter = 1;
for i = 1:size(A_x, 1)
    for j = 1:size(A_y, 1)
        row = i;
        col = (j - 1)* nx;
        ind = row + col;
        H_(ind,:) = H_mat(counter, :);
        counter = counter + 1;
    end
end

H_mat = H_;
H_mat = H_mat(val_idx,:);
H_mat = normr(H_mat);
label_list = label_list(val_idx);
idx = kmeans(H_mat, rank, 'Replicates', 5);

vis = zeros([size(A_x, 1), size(A_y, 1)]);
vis(val_idx) = idx;
img = imagesc(vis);

ari = rand_index(idx, label_list, 'adjusted');

ari = round(ari , 3);
ari_list(counter) = ari;
        
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

clustering_mapped = clustering;
for i = 0:num_regions-1
    clustering_mapped(clustering == i) = assignment(i+1);
end

img_data = clustering_mapped;
img_data(annomat == -1) = nan;

figure
h = imagesc(img_data);
colormap(cmap)
set(h, 'AlphaData', ~isnan(img_data));
set(gcf, 'Position', [100, 100, 300, 300]);
axis image off

savepath = ['vis/NMF_' dataset, '_clustering.png'];
disp(['ARI=', num2str(round(ari,3))]);
saveas(gcf, savepath);       
% close all