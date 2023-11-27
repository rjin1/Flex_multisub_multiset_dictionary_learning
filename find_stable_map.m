% Top level code for choosing best DL run from all the given runs
addpath_recurse("GroupICATv4.0c");
path = "experiments\real\K38\results\1";
% path = "experiments\real\K30-138(110-28)-109(87-22)-correct\runs\0.008 0.363 0.16 0.055 0.99";
% path = "experiments\real\K30-138(110-28)-109(87-22)-correct\runs\0.009 0.36 0.15 0.054 0.99";
costs = [];
Zs = {}; 
recons = [];
for i = 1:50
    data_fn = fullfile(path, string(i), 'Z_hat_ks.mat');
    load(data_fn); 
    data_fn = fullfile(path, string(i), 'Z_check_ks.mat');
    load(data_fn);
    data_fn = fullfile(path, string(i), 'costs.mat');
    load(data_fn); 
    %{
    costs = [costs; [i, cost]];
    recons = [recons; [i, cost_recon]];
    %}
    Z_ks = Z_hat_ks + Z_check_ks;
    [S, num_voxels, num_components] = size(Z_ks);
    Z = zeros(S * num_components, num_voxels);
    for k = 1: num_components
        Z(k, :) = Z_ks(1, :, k);
        Z(num_components + k, :) = Z_ks(2, :, k);
        Z(2 * num_components + k, :) = Z_ks(3, :, k);
    end
    Zs{i} = Z;
end 
%{
sorted_cost = sortrows(costs, 2);
sorted_recons = sortrows(recons, 2);
fprintf("Run with min cost is %d\n", sorted_cost(1, 1));
fprintf("Run with min reconstruction error is %d\n", sorted_recons(1, 1));
%}
[RunCorrT, BestZ, BestZIndex, IndexAssign, Tmaps] = DL_bestRunSelection(Zs);
 
fprintf("BestZIndex = %d\n", BestZIndex);
%% 
clear; close all;