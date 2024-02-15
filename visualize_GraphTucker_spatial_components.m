clear all;clc

tissue_name = 'MOSTA_9.5'; % name of data being used
work_path = 'GraphTucker/'; % path to where GraphTucker folder is located
data_path = ['GraphTucker/data/', tissue_name, '/']; % path to where spatial gene expression data is located
utils_path =[work_path, 'GT_utils/'];
res_path = [work_path, 'res/', tissue_name, '/'];
dst_path = ['GraphTucker/vis/', tissue_name, '/']; % path to where spatial component visualizations should be saved to
addpath(work_path) % add path to ensure scripts can be found correctly
addpath(utils_path) % add path so util files can be found

% set to tested rank and lambda
rank_list = [64 64 64];
lambda = 0.1;

% load GraphTucker components.
% this script loads an example set of GraphTucker components run 
% with rank=(64,64,64), lambda=0.1 with 5000 iterations, denoted by EXAMPLE.
% Change the path string to load a different set of results.
rank_str = [num2str(rank_list(1)), '-', num2str(rank_list(2)), '-', num2str(rank_list(3))];
load([res_path, 'GT_', tissue_name, '_rank=', rank_str ,'_lambda=', num2str(lambda) ,'_EXAMPLE.mat']);

% permute core tensor for n-mode multiplication
G = permute(G, [3 2 1]);

% construct spatial component tensor X
% X is of size (n_x, n_y, r_g)
X = ttm(G, {A_set{1}{2} A_set{1}{1}}, [1 2]);

A_x = A_set{1}{2};
A_y = A_set{1}{1};

% need to load list of valid spots to distinguish from background spots
% for better visualization
load([data_path, tissue_name, '_val_subs.mat']);
val_subs = [T_subs(:,1), T_subs(:,2)];
val_subs = unique(val_subs, 'rows');
bg_mat = zeros(size(A_x,1), size(A_y,1));
val_subs = sub2ind(size(bg_mat), val_subs(:,2), val_subs(:,1));
bg_mat(val_subs) = 1;

% loop through all spatial components
for i=1:rank_list(3)
    i

    % create distinct folders depending on parameters to avoid clutter
    foldername = [dst_path, 'rank=', rank_str , '/'];
    if ~exist(foldername, 'dir')
       mkdir(foldername)
    end
    lam_foldername = [foldername, 'lambda=', num2str(lambda) , '/'];
    if ~exist(lam_foldername, 'dir')
       mkdir(lam_foldername)
    end
    savefolder = lam_foldername;
    
    % grab i-th spatial component
    X_slice = X.data(:,:,i);
    img_data = X_slice;
    
    figure
    % set background spots to nan for visualization purposes
    img_data(bg_mat == 0) = nan;
    h = imagesc(img_data);
    set(h, 'AlphaData', ~isnan(img_data));

    savepath = [savefolder,'/g_comp=', num2str(i), '.png'];
    axis image off
    saveas(gcf, savepath);
    clf
end
close all