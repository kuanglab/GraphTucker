
% V: data tensor 
% W: spatial chain graphs and PPI network
function [V,W] = data_prep_mosta(data_name, data_path, utils_path)

% load data tensor
load([data_path, data_name,'_tensor.mat']); 

% load PPI network. This script loads the Mouse PPI (MUS) by default
load([data_path,'MUS_PPI.mat']);

% load gene ids
genes = importdata([data_path,data_name,'_gene.csv']);
genes = genes(2:end);
for i = 1:size(genes)
    tmp = split(genes{i},',');
    genes{i} = tmp{2};
end

cd(utils_path);
V = sptensor(tensor(double(expr_mat)));
V = permute(V, [2 3 1]);
vals = V.vals; 
vals(vals==2) = 1; % set to zeros if UMI counts < 3. 
V = sptensor(V.subs,log(vals),V.size); % log normalization



% eliminate boundaries
Z = double(tenmat(V,1));
ind1 = find(sum(Z,2) ~= 0);
Z = double(tenmat(V,2));
ind2 = find(sum(Z,2) ~= 0);
V  = V(ind1,ind2,:);

% Filter out the low-density genes
Z = double(tenmat(V,3));
Z(Z>0)=1;
total_counts = sum(Z,2);
ind3 = find(total_counts>3); % remove genes that are expressed in less than 4 spots
V  = V(:,:,ind3);
n = [size(V,1),size(V,2),size(V,3)];
genes = genes(ind3);

W = cell(3,1);
% build spatial graphs
for netid = 1:2 
    W{netid} = diag(ones(n(netid)-1,1),-1) + diag(ones(n(netid)-1,1),1);
end

% build PPI network
% This script loads mouse PPI by default
W{3} = zeros(n(3),n(3));

[ia,ib] = ismember(genes,MUS_GENE);
W{3}(ia,ia) = MUS_BIOGRID(ib(ib>0),ib(ib>0));

end
