clear all;clc

opts.stopcrit = 10^-4; % stopping criteria for GraphTucker convergence
opts.maxiters = 5; % maxiters=5000 used for all experiments
opts.mode = 1; % set to 0 to run GraphTucker with sparse data tensors, or 1 for non-sparse
% we observe faster runtimes when using non-sparse tensors, so we set this
% to 1

tissue_name = 'MOSTA_9.5';
work_path = 'GraphTucker/'; % path to where GraphTucker folder is located
data_path = ['GraphTucker/data/', tissue_name, '/']; % path to where spatial gene expression data is located
utils_path =[work_path, 'GT_utils/'];
res_path = [work_path, 'res/', tissue_name, '/'];
addpath(work_path) % add path to ensure scripts can be found correctly
addpath(utils_path) % add path so util files can be found

% initialize GraphTucker parameters
lambda = 0.1; % graph regularization parameter
tucker_rank = [64, 64, 64]; %core tensor/factor matrix ranks: r_x, r_y, r_g
opts.rank_g = tucker_rank(3);
opts.rank_set = tucker_rank(1:2);
rank_str = [num2str(opts.rank_set(1)), '-', num2str(opts.rank_set(2)), '-', num2str(opts.rank_g)];

data_name = tissue_name;

% prepare data tensor and graphs
[T, W] = data_prep_mosta(data_name, data_path, utils_path);

disp(['Spatial gene expression tensor size:', num2str(size(T))])
disp(['Density: ', num2str(length(T.vals) / prod(size(T)))]);

% list of non-background spots. Needed for spatial component visualization
T_subs = T.subs;
save([data_path, data_name, '_val_subs.mat'], 'T_subs'); 

rng('default');
% create mask tensor. No cross validation so mask is all ones
M = tenones(size(T));

%permute tensor for GraphTucker input
T = permute(T, [3 1 2]);
T_set = {T};
W_set = {{W{1}, W{2}}};
W_g = W{3};

% permute mask tensor to match Y0
M = permute(M, [3 1 2]);
M_set = {M};

opts.lambda = lambda;
opts.tissue_name = data_name;
opts.loss_iters = 50; % calculate training loss every 50 iterations. Can be changed to any number of iterations

T_perm = permute(T, [3 1 2]); % used for calculating training loss

tic
[A_set, A_g, G, train_loss] = GraphTucker(T_set, M_set, W_set, W_g, T_perm, opts);
toc

disp(['Time taken: ', num2str(toc - tic)])

name=[res_path,'GT_', tissue_name,'_rank=', rank_str, '_lambda=', num2str(lambda), '.mat'];
save(name,'A_set', 'A_g', 'G', 'lambda', 'train_loss', 'opts.maxiters', '-v7.3');

