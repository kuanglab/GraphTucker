% addpath('FIST_utils/')
clear all
rng('default')

vis_path = 'vis/BRCA1'; % path to where images will be saved to

dataset = 'BRCA1';
num_regions = 20;
rank_list = [30 30 num_regions];
rank_str = [num2str(rank_list(1)), '-', num2str(rank_list(2)), '-', num2str(rank_list(3))];
lambda = 1;

load(['res/', dataset, '/GT_', dataset, '_rank=', rank_str ,'_lambda=', num2str(lambda) ,'_EXAMPLE.mat']);
G = permute(G, [3 2 1]);
A_x_ = A_set{1}{2};
A_y_ = A_set{1}{1};
n = [size(A_x_, 1), size(A_y_, 1)];

% sort core tensor entries for selecting top K% of entries
[vals, idx] = sort(G.data(:), 'descend');
num_el = prod(rank_list);

pct_list = [1 2 5 10 25 100];
avg_ed_list = zeros(1, length(pct_list));

for pcti = 1:length(pct_list)

pct = pct_list(pcti) / 100;
pct_str = num2str(pct);

G_tmp = G;

if pct ~= 100
    stop_idx = ceil(num_el * pct);
    exclude_idx = idx(stop_idx+1:end);
    G_tmp(exclude_idx) = 0; % remove top pct% of core tensor entries
end

% construct spatial component tensor
X = ttm(G_tmp, {A_set{1}{2} A_set{1}{1}}, [1 2]); 

A_x = A_set{1}{2};
A_y = A_set{1}{1};

G_xy_mat = tenmat(X, 3, [1 2]); % matricize spat. comp tensor
G_xy = G_xy_mat.data;
G_xy = G_xy';

% load tissue annotation file for dataset
load(['data/tissue_annotations/', dataset, '_annotations.mat'])
annotations(:,1:2) = annotations(:,1:2) + 1;

% load subs file which contains indices for valid spots
load(['processed_data/', dataset, '/', dataset, '_val_indices.mat'])

% arrange annotation onto spots correctly
annomat = zeros(size(A_x, 1), size(A_y, 1)) - 1;
for i = 1:size(annotations, 1)
    annomat(annotations(i,1), annotations(i,2)) = annotations(i,3);
end

% remove boundary spots so annotation matches size of (A_x, A_y)
annomat = annomat(ind1, ind2);

aucmat = zeros(rank_list(3), num_regions); % pairwise AUCs
edmat = zeros(rank_list(3), num_regions); % pairwise euclidean distances


for i=1:rank_list(3)
    comp = G_xy(:,i); % get spatial component
    dist_vec = zeros(num_regions,1);
    auc_vec = zeros(num_regions,1);
    for label = 0:num_regions-1

        %remove background spots from component
        comp(annomat(:) == -1) = 0;

        % binarize + vectorize annotation for current region
        anno_vec = annomat;
        anno_vec(anno_vec ~= label) = -1;
        anno_vec(anno_vec == label) = 1;
        anno_vec(anno_vec == -1) = 0;
        anno_vec = anno_vec(:);
        
        % normalize component values
        comp = normr(comp')';

        % get spots not in current region
        inval_idx = find(anno_vec(:) == 0);
        
        % normalize annotation vector
        anno_vec_nn = anno_vec;
        anno_vec = normr(anno_vec')';     
        
        % calc AUC between current spatial component and the region
        [~, ~, ~, AUC] = perfcurve(anno_vec_nn, comp, 1); 
        auc_vec(label+1) = AUC;
        auc_mat(i, label+1) = AUC;

        %calculate ED between current spatial component and the region
        ed_mat(i, label+1) = pdist2(anno_vec', comp');
    end

end

% do a greedy selection for choosing region-spatial-component matches
match_T = greedy_comp_select(round(ed_mat, 3));
best_auc = zeros(size(match_T, 1),1);

% get corresponding AUCs for select region-component matches
for i = 1:size(match_T, 1)
    row = match_T(i,2) + 1;
    col = match_T(i,3);
    best_auc(i) = auc_mat(row, col);
end

match_T = [match_T round(best_auc,3)];

match_T(:,2) = match_T(:,2) + 1;
match_T = sortrows(match_T, 3); 
match_T = [match_T(:,2:4) match_T(:,1)];

% calculate Avg. ED across matches
avg_ed_list(pcti) = mean(match_T(:,4));
end

% choose core tensor percentage with lowest Avg. ED
[~, min_idx] = min(avg_ed_list);
min_pct = pct_list(min_idx);
disp(['Lowest avg. ED using top ', num2str(min_pct), '% of core tensor entries'])

% repeat above with best selected core tensor percentage
% 
pct = min_pct / 100;
pct_str = num2str(pct);

G_tmp = G;

if pct ~= 100
    stop_idx = ceil(num_el * pct);
    exclude_idx = idx(stop_idx+1:end);
    G_tmp(exclude_idx) = 0; 
end

X = ttm(G_tmp, {A_set{1}{2} A_set{1}{1}}, [1 2]);

A_x = A_set{1}{2};
A_y = A_set{1}{1};

G_xy_mat = tenmat(X, 3, [1 2]); % matricize spat. comp tensor
G_xy = G_xy_mat.data;
G_xy = G_xy';

aucmat = zeros(rank_list(3), num_regions); % AUCs
edmat = zeros(rank_list(3), num_regions); % euclidean distances


for i=1:rank_list(3)
    comp = G_xy(:,i); % get spatial component
    for label = 0:num_regions-1

        %remove background spots from component
        comp(annomat(:) == -1) = 0;

        % binarize + vectorize annotation for current region
        anno_vec = annomat;
        anno_vec(anno_vec ~= label) = -1;
        anno_vec(anno_vec == label) = 1;
        anno_vec(anno_vec == -1) = 0;
        anno_vec = anno_vec(:);
        
        % normalize component values
        comp = normr(comp')';

        % get spots not in current region
        inval_idx = find(anno_vec(:) == 0);
        
        % normalize annotation  vector
        anno_vec_nn = anno_vec;
        anno_vec = normr(anno_vec')';     
        
        [~, ~, ~, AUC] = perfcurve(anno_vec_nn, comp, 1);
        auc_mat(i, label+1) = AUC;
        ed_mat(i, label+1) = pdist2(anno_vec', comp');
    end
end

match_T = greedy_comp_select(round(ed_mat, 3));
best_auc = zeros(size(match_T, 1),1);

for i = 1:size(match_T, 1)
    row = match_T(i,2) + 1;
    col = match_T(i,3);
    best_auc(i) = auc_mat(row, col);
end

match_T = [match_T round(best_auc,3)];

match_T(:,2) = match_T(:,2) + 1;
match_T = sortrows(match_T, 3); 
match_T = [match_T(:,2:4) match_T(:,1)];

% save visualizations  of spatial components 
% that match with each of the regions
for i = 1:size(match_T, 1)
    figure
    comp_num = match_T(i,1);

    comp_vis = reshape(G_xy(:,comp_num), n);
    img_data = comp_vis;

    img_data(annomat == -1) = nan;

    h = imagesc(img_data);
    set(h, 'AlphaData', ~isnan(img_data))
    set(gcf, 'Position', [100, 100, 300, 300])
    axis image off
    title_text = ['(', num2str(match_T(i,3), 3), '), (', num2str(match_T(i,4), 3), ')'];
    title(title_text)
    img_path = [vis_path, '/spatial_component_matches/region=', num2str(i), '_comp=', num2str(comp_num), '.png'];
    saveas(gcf, img_path);
end
disp(['Avg AUC of matched spatial components: ', num2str(mean(match_T(:,3)))]);
disp(['Avg ED of matched spatial components: ', num2str(mean(match_T(:,4)))]);

close all

% function for doing greedy selection of region-spatial-component matches
function dist_mat = greedy_comp_select(match_Table)
    dist_mat = zeros(size(match_Table, 2), 3);
    match_T = match_Table; % copy of table we can delete columns from
    comp_counts = zeros(size(match_Table, 1), 1);
    for j = 1:size(match_Table, 1)
        [M, I] = min(match_T, [], "all", 'linear');
        dist_mat(j,1) = M;
        [min_row, min_col] = ind2sub(size(match_T), I);
        dist_mat(j,3) = min_col;
        dist_mat(j,2) = min_row - 1;

        match_T(:,min_col) = inf;
    end
end
