function [A_set, A_g, G, train_loss] = FIST_tucker_loss(T_set, M_set, W_set, W_g, V, opts)
%% Input
% T_set: the set of input tensors for replicates (2D: n_g by n_x by n_y; 3D: n_g by n_x by n_y by n_z)
% M_set: the set of mask tensors for replicates (2D: n_g by n_x by n_y; 3D: n_g by n_x by n_y by n_z)
% W_set: the set of chain graph along spatial coordinates in cell arrays for replicates
% W_g: protein-protein interaction graph shared by replicates
% opts: parameter object
%    - rank_set: the set of tensor Tucker rank along different spatail coordinates
%    - rank_g: tensor Tucker rank along gene coordinate
%    - lambda: graph hyperparameter
%    - mode: 0 denotes sparse FIST for Visium data, 1 denotes dense FIST for ST data
%    - stopcrit: stop criteria
%    - maxiters: maximum iteration
%% Output
% A_set: the set of Tucker spatial component matrices in cell arrays for replicates
% A_g: Tucker gene component matrix shared by replicates
% G: Tucker core tensor shared by replicates


train_loss = [];

eps = 1e-10; % numerical stability

n_reps = length(T_set); % the number of replicates
n_coords = length(W_set{1}); % the number of spatial coordinates
n_ways = ndims(T_set{1}); % the number of ways
n_g = size(W_g, 1); % the number of genes

% Obtain tensor dimensions for replicates
% 2D: n_y by n_x by n_g
% 3D: n_z by n_y by n_x by n_g
n_set = cell(1, n_reps);
for i = 1:n_reps
    n_set{i} = flip(size(T_set{i})); % tensor dimensions
end

% Set tensor ranks
% 2D: n_y by n_x by n_g
% 3D: n_z by n_y by n_x by n_g
ranks = [flip(opts.rank_set), opts.rank_g];

% Convert both input tensors and masks into dense tensors
if opts.mode == 1 
    for i = 1:n_reps
        M_set{i} = tensor(M_set{i});
        T_set{i} = tensor(T_set{i});
    end
end

% Normalize spatial chain graphs for replicates
% 2D: n_y by n_x by n_g
% 3D: n_z by n_y by n_x by n_g
D_set = cell(1, n_reps);
for i = 1:n_reps
    W = flip(W_set{i});
    D = cell(1, n_coords);
    for j = 1:n_coords
        W{j} = W{j} - diag(diag(W{j}));
        d = sum(W{j}, 2);
        d(d ~= 0) = (d(d ~= 0)) .^ -(0.5);
        W{j} = W{j} .* d;
        W{j} = d' .* W{j};
        D{j} = diag(sum(W{j}, 2));
    end
    W_set{i} = W;
    D_set{i} = D;
end

% Normalize protein-protein interaction graph
W_g = W_g - diag(diag(W_g));
d = sum(W_g, 2);
d(d ~= 0) = (d(d ~= 0)) .^ -(0.5);
W_g = W_g .* d;
W_g = d' .* W_g;
D_g = diag(sum(W_g, 2));

% Initialize spatial components for replicates
% and auxiliary variables
% 2D: n_y by n_x by n_g
% 3D: n_z by n_y by n_x by n_g
A_set = cell(1, n_reps);
ATA_set = cell(1, n_reps);
WA_set = cell(1, n_reps);
DA_set = cell(1, n_reps);
ATWA_set = cell(1, n_reps);
ATDA_set = cell(1, n_reps);
for i = 1:n_reps
    A_set{i} = cell(1, n_coords);
    ATA_set{i} = cell(1, n_coords);
    WA_set{i} = cell(1, n_coords);
    DA_set{i} = cell(1, n_coords);
    ATWA_set{i} = cell(1, n_coords);
    ATDA_set{i} = cell(1, n_coords);
end

for i = 1:n_reps
    for j = 1:n_coords
        rng(i);
        A_set{i}{j} = rand(n_set{i}(j), ranks(j));
        ATA_set{i}{j} = A_set{i}{j}' * A_set{i}{j};
        WA_set{i}{j} = W_set{i}{j} * A_set{i}{j};
        DA_set{i}{j} = diag(D_set{i}{j}) .* A_set{i}{j};
        ATWA_set{i}{j} = A_set{i}{j}' * WA_set{i}{j};
        ATDA_set{i}{j} = A_set{i}{j}' * DA_set{i}{j};
    end
end

% Initialize gene component shared by replicates
% and auxiliary variables
rng(0);
A_g = rand(n_g, opts.rank_g);
ATA_g = A_g' * A_g;
WA_g = W_g * A_g;
DA_g = diag(D_g) .* A_g;
ATWA_g = A_g' * WA_g;
ATDA_g = A_g' * DA_g;

% Initialize core tensor shared by replicates
% 2D: n_y by n_x by n_g
% 3D: n_z by n_y by n_x by n_g
rng(0);
G = tensor(rand(ranks), ranks);

% Initialize tensor related auxiliary variables
MT = cell(1, n_reps);
MT_ = cell(1, n_reps);
for i = 1:n_reps
    % Reorder the modes of the tensor
    % 2D: n_y by n_x by n_g
    % 3D: n_z by n_y by n_x by n_g
    M_set{i} = permute(M_set{i}, flip(1:n_ways));
    T_set{i} = permute(T_set{i}, flip(1:n_ways));

    MT{i} = cell(1, n_coords+1);
    MT_{i} = cell(1, n_coords+1);

    % 2D: {A_y, A_x, A_g}
    % 3D: {A_z, A_y, A_x, A_g}
    A = [A_set{i}, A_g];
    for j = 1:(n_coords+1)
        MT{i}{j} = tennmat_fc(M_set{i} .* T_set{i}, j);
        MT_{i}{j} = tennmat_fc(M_set{i} .* full(ttensor(G, A)), j);
    end
end

% training
for iter = 1:opts.maxiters

    % Record previous results
    A_set_old = A_set;
    A_g_old = A_g;
    G_old = G;
    disp('Updating Ax/y')

    % Update spatial component matrix along j-axis for i-th replicate
    for i = 1:n_reps
        for j = 1:n_coords
            A = [A_set{i}, A_g];
            A{j} = eye(ranks(j));
            T_ = full(ttensor(G, A));
    
            % J1 negative
            J1n = MT{i}{j} * transpose(tennmat_fc(T_, j));
    
            % J1 positive
            J1p = MT_{i}{j} * transpose(tennmat_fc(T_, j));
    
    
            G_ = tenvec_fc(G, j+1); % vectorize the core tensor along j+1-axis
            G_ = reshape(G_, prod(ranks)/ranks(j), ranks(j)); % vec-transpose
    
            % J2 negative
            J2n = 0;
            J2p = 0;
    
            if opts.lambda ~= 0
                for k = 1:(n_coords+1)
                    ATA = [ATA_set{i}, ATA_g];
                    ATA{j} = eye(ranks(j));
                    ATWA = [ATWA_set{i}, ATWA_g];
                    ATWA{j} = eye(ranks(j));
                    if k ~= j
                        ATA{k} = ATWA{k};
                        G__ = tenvec_fc(ttm(G, ATA, 1:n_ways), j+1); % vectorize the tensor along j+1-axis (see eq (15) and (16) in slides)
                        G__ = reshape(G__, prod(ranks)/ranks(j), ranks(j)); % vec-transpose
                        J2n = J2n + A_set{i}{j} * G_' * G__;
                    else
                        G__ = tenvec_fc(ttm(G, ATA, 1:n_ways), j+1); % vectorize the tensor along j+1-axis (see eq (15) and (16) in slides)
                        G__ = reshape(G__, prod(ranks)/ranks(j), ranks(j)); % vec-transpose
                        J2n = J2n + WA_set{i}{j} * G_' * G__;
                    end
                end
        
                % J2 positive
                J2p = 0;
                for k = 1:(n_coords+1)
                    ATA = [ATA_set{i}, ATA_g];
                    ATA{j} = eye(ranks(j));
                    ATDA = [ATDA_set{i}, ATDA_g];
                    ATDA{j} = eye(ranks(j));
                    if k ~= j
                        ATA{k} = ATDA{k};
                        G__ = tenvec_fc(ttm(G, ATA, 1:n_ways), j+1); % vectorize the tensor along j+1-axis (see eq (15) and (16) in slides)
                        G__ = reshape(G__, prod(ranks)/ranks(j), ranks(j)); % vec-transpose
                        J2p = J2p + A_set{i}{j} * G_' * G__;
                    else
                        G__ = tenvec_fc(ttm(G, ATA, 1:n_ways), j+1); % vectorize the tensor along j+1-axis (see eq (15) and (16) in slides)
                        G__ = reshape(G__, prod(ranks)/ranks(j), ranks(j)); % vec-transpose
                        J2p = J2p + DA_set{i}{j} * G_' * G__;
                    end
                 end
    	    end

            num = J1n + opts.lambda * J2n + eps;
            denom = J1p + opts.lambda * J2p + eps;
            A_set{i}{j} = A_set{i}{j} .* (num ./ denom);
    
            % Update spatial components related auxiliary variables
            ATA_set{i}{j} = A_set{i}{j}' * A_set{i}{j};
            WA_set{i}{j} = W_set{i}{j} * A_set{i}{j};
            ATWA_set{i}{j} = A_set{i}{j}' * WA_set{i}{j};
            DA_set{i}{j} = diag(D_set{i}{j}) .* A_set{i}{j};
            ATDA_set{i}{j} = A_set{i}{j}' * DA_set{i}{j};
            A = [A_set{i}, A_g];
            T_ = full(ttensor(G, A));
            for k = 1:(n_coords+1)
                MT_{i}{k} = tennmat_fc(M_set{i} .* T_, k);
            end
        end
    end

    % Update gene component matrix
    num = 0;
    denom = 0;
    disp('Updating Ag')
    for i = 1:n_reps

        A = [A_set{i}, A_g];
        A{n_ways} = eye(ranks(n_ways));
        T_ = full(ttensor(G, A));

        % J1 negative
        J1n = MT{i}{n_ways} * transpose(tennmat_fc(T_, n_ways));

        % J1 positive
        J1p = MT_{i}{n_ways} * transpose(tennmat_fc(T_, n_ways));


        G_ = tenvec_fc(G, 1); % vectorize the core tensor along j+1-axis
        G_ = reshape(G_, prod(ranks)/ranks(n_ways), ranks(n_ways)); % vec-transpose

        % J2 negative
        J2n = 0;
	    J2p = 0;
	
	if opts.lambda ~= 0 % we don't need to compute all this if lambda=0
        for k = 1:(n_coords+1)
            ATA = [ATA_set{i}, ATA_g];
            ATA{n_ways} = eye(ranks(n_ways));
            ATWA = [ATWA_set{i}, ATWA_g];
            ATWA{n_ways} = eye(ranks(n_ways));
            if k ~= n_ways
                ATA{k} = ATWA{k};
                G__ = tenvec_fc(ttm(G, ATA, 1:n_ways), 1); % vectorize the tensor along first axis (see eq (11) in slides)
                G__ = reshape(G__, prod(ranks)/ranks(n_ways), ranks(n_ways)); % vec-transpose
                J2n = J2n + A_g * G_' * G__;
            else
                G__ = tenvec_fc(ttm(G, ATA, 1:n_ways), 1); % vectorize the tensor along first axis (see eq (11) in slides)
                G__ = reshape(G__, prod(ranks)/ranks(n_ways), ranks(n_ways)); % vec-transpose
                J2n = J2n + WA_g * G_' * G__;
            end
        end

        % J2 positive
        for k = 1:(n_coords+1)
            ATA = [ATA_set{i}, ATA_g];
            ATA{n_ways} = eye(ranks(n_ways));
            ATDA = [ATDA_set{i}, ATDA_g];
            ATDA{n_ways} = eye(ranks(n_ways));
            if k ~= n_ways
                ATA{k} = ATDA{k};
                G__ = tenvec_fc(ttm(G, ATA, 1:n_ways), 1); % vectorize the tensor along first axis (see eq (11) in slides)
                G__ = reshape(G__, prod(ranks)/ranks(n_ways), ranks(n_ways)); % vec-transpose
                J2p = J2p + A_g * G_' * G__;
            else
                G__ = tenvec_fc(ttm(G, ATA, 1:n_ways), 1); % vectorize the tensor along first axis (see eq (11) in slides)
                G__ = reshape(G__, prod(ranks)/ranks(n_ways), ranks(n_ways)); % vec-transpose
                J2p = J2p + DA_g * G_' * G__;
            end
        end
	end

        num = num + J1n + opts.lambda * J2n;
        denom = denom + J1p + opts.lambda * J2p;

    end

    num = num + eps;
    denom = denom + eps;
    A_g = A_g .* (num ./ denom);

    % Update gene component related auxiliary variables
    ATA_g = A_g' * A_g;
    WA_g = W_g * A_g;
    DA_g = diag(D_g) .* A_g;
    ATWA_g = A_g' * WA_g;
    ATDA_g = A_g' * DA_g;

    for i = 1:n_reps
        A = [A_set{i}, A_g];
        T_ = full(ttensor(G, A));
        for k = 1:(n_coords+1)
            MT_{i}{k} = tennmat_fc(M_set{i} .* T_, k);
        end
    end

    % Update core tensor
    disp('Updating core tensor')
    num = 0;
    denom = 0;

    for i = 1:n_reps
        A = [A_set{i}, A_g];
        T_ = full(ttensor(G, A));
        A = cellfun(@transpose, A, 'UniformOutput', false);
        J1n = ttm(M_set{i} .* T_set{i}, A, 1:n_ways);
        J1p = ttm(M_set{i} .* T_, A, 1:n_ways);

        % J2 negative
        J2n = 0;
	J2p = 0;

	if opts.lambda ~= 0
        for k = 1:(n_coords+1)
            ATA = [ATA_set{i}, ATA_g];
            ATWA = [ATWA_set{i}, ATWA_g];
            ATA{k} = ATWA{k};
            J2n = J2n + ttm(G, ATA, 1:n_ways);
        end

        % J2 positive
        J2p = 0;
        for k = 1:(n_coords+1)
            ATA = [ATA_set{i}, ATA_g];
            ATDA = [ATDA_set{i}, ATDA_g];
            ATA{k} = ATDA{k};
            J2p = J2p + ttm(G, ATA, 1:n_ways);
        end
	end

        num = num + J1n + opts.lambda * J2n;
        denom = denom + J1p + opts.lambda * J2p;

    end

    num = num + eps;
    denom = denom + eps;
    G = G .* (num ./ denom);

    % Update core tensor related auxiliary variables
    for i = 1:n_reps
        A = [A_set{i}, A_g];
        T_ = full(ttensor(G, A));
        for k = 1:(n_coords+1)
            MT_{i}{k} = tennmat_fc(M_set{i} .* T_, k);
        end
    end

    if mod(iter, opts.loss_iters) == 0 % calculate loss every loss_iters iterations
       A = [A_set{i}, A_g];
       Y_curr = tensor(ttensor(G, A)); % y x g
       Y_curr = permute(Y_curr, [3 2 1]); % y x g -> g x y
       
       M_train = permute(M_set{1}, [3 2 1]); % y x g -> g x y

       train_len = numel(M_set{1}.data);

       Y_train = tensor(Y_curr .* M_train);
       V_train = tensor(V .* M_train);

       train_rmse = sqrt(sum((Y_train.data - V_train.data).^2, 'all') / train_len);

       train_loss = [train_loss train_rmse];
       disp('Train loss:')
       disp(train_rmse(end))
    end

    norm_matrices = cell(n_ways, 1);
    % normalize factor matrices with diagonal matrices U,V, and W.
    % First, create U with dimensions: n_y x n_y
    for k = 1:(n_coords)
        % First, create U with dimensions: n_y x n_y
        U = diag(sum(A_set{1}{k}, 1));
        % Multiple A_y by U^-1
        A_set{i}{k} = A_set{1}{k} * pinv(U);
        norm_matrices{k} = U;

    end

    % Lastly, create W with dim: n_g by n_g
    W = diag(sum(A_g, 1));
    % Normalize A_g with W^-1
    A_g = A_g * pinv(W);
    norm_matrices{n_ways} = W;
    % Take n-mode product of core tensor G with U,V, and W
    G = ttm(G, norm_matrices, [1:n_ways]);

    res = compute_res(A_set, A_set_old, A_g, A_g_old, G, G_old, n_reps);
    disp(['GraphTucker training...residual: ', num2str(res)]);
    disp(['GraphTucker training...iteration: ', num2str(iter)]);
    if res < opts.stopcrit
        break;
    end
end

% Change the order of spatial components
% 2D: {A_g, A_x, A_y}
% 3D: {A_g, A_x, A_y, A_z}
for i = 1:n_reps
    A_set{i} = flip(A_set{i});
end
% Reorder the mode of core tensor
% 2D: n_g by n_x by n_y
% 3D: n_g by n_x by n_y by n_z
G = permute(G, flip(1:n_ways));

end

function T_n = tennmat_fc(T, n)
% convert tensor to matrix (mode-n matricization in forward cyclic manner)
n_ways = ndims(T);
T_n = tenmat(T, n, [(n+1):n_ways, (1:n-1)]);
T_n = T_n.data;
end

function t = tenvec_fc(T, n)
% convert tensor to vector (vectorization after mode-n matricization in forward cyclic manner)
n_ways = ndims(T);
t = tenmat(T, [n:n_ways, (1:n-1)]);
t = t.data;
end

function res = compute_res(A_set, A_set_old, A_g, A_g_old, G, G_old, n_reps)
% compute residual over all components and core tensor
res_num = 0;
res_denom = 0;
for i = 1:n_reps
    for j = 1:length(A_set{i})
        res_num = res_num + sum(sum((A_set{i}{j} - A_set_old{i}{j}).^2));
        res_denom = res_denom + sum(sum(A_set_old{i}{j}.^2));
    end
end
res_num = res_num + sum(sum((A_g - A_g_old).^2));
res_denom = res_denom + sum(sum(A_g_old.^2));
res_num = res_num + collapse((G - G_old).^2);
res_denom = res_denom + collapse(G_old.^2);
res = sqrt(res_num / res_denom);
end 
