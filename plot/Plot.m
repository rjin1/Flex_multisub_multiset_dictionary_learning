filePath = ".\most stable run - 34\34";
plot_col = [10, 10, 10];
map_plot_ind = [1,2,3,4,5,6,7,9,10,11]; % 1,2,3,4,5,6,7,9,10,11,13 %15,18,19, 20,23,24,25,27,28,29,30
ind_start = 0;
disp(length(map_plot_ind));
% filePath = "..\data\1 0.008 0.35 0.3 0.05 0.99 30 50 0";
files = ["D", "D_hat", "Z_check_ks", "Z_hat_ks", "costs"];
for i = 1: size(files, 2)
    fn = fullfile(filePath, files(i) + ".mat");
    load(fn); 
end 

D = D(:, [map_plot_ind, map_plot_ind+30, map_plot_ind+60]);
D_hat = D_hat(:, [map_plot_ind, map_plot_ind+30, map_plot_ind+60]);
Z_check_ks = Z_check_ks(:,:,map_plot_ind);
Z_hat_ks = Z_hat_ks(:,:,map_plot_ind);
K = length(map_plot_ind);

S = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Construct Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z1 = zeros(K, num_voxels);
Z2 = zeros(K, num_voxels);
Z3 = zeros(K, num_voxels);
Z_hat1 = zeros(K, num_voxels);
Z_hat2 = zeros(K, num_voxels);
Z_hat3 = zeros(K, num_voxels);
Z_check1 = zeros(K, num_voxels);
Z_check2 = zeros(K, num_voxels);
Z_check3 = zeros(K, num_voxels);
Z_ks = Z_hat_ks + Z_check_ks;
for k = 1: K
    Z1(k, :) = Z_ks(1, :, k); 
    Z2(k, :) = Z_ks(2, :, k); 
    Z3(k, :) = Z_ks(3, :, k); 
    Z_hat1(k, :) = Z_hat_ks(1, :, k); 
    Z_hat2(k, :) = Z_hat_ks(2, :, k); 
    Z_hat3(k, :) = Z_hat_ks(3, :, k); 
    Z_check1(k, :) = Z_check_ks(1, :, k); 
    Z_check2(k, :) = Z_check_ks(2, :, k); 
    Z_check3(k, :) = Z_check_ks(3, :, k); 
end
D_tilde = D - D_hat;
D1 = D(:, 1: K);
D2 = D(:, K + 1: 2 * K);
D3 = D(:, 2 * K + 1: 3 * K);

D1_tilde = D_tilde(:, 1: K);
D2_tilde = D_tilde(:, K + 1: 2 * K);
D3_tilde = D_tilde(:, 2 * K + 1: 3 * K);

D1_hat = D_hat(:, 1: K);
D2_hat = D_hat(:, K + 1: 2 * K);
D3_hat = D_hat(:, 2 * K + 1: 3 * K);

%%%%%%%%%%%%%%%%% Compute p Values Of D/D_tilde/D_hat %%%%%%%%%%%%%%%%%%%%%%%
p1 = zeros(1, K);
p2 = zeros(1, K);     
p3 = zeros(1, K);
sele = ones(S * K);
idx = 1;
for k = 1: K     
    [~, p1(k), ~, tstats] = ttest2(D1(1:subjects_range(2), k), D1(subjects_range(2) + 1:subjects_range(3),k), .05);
    sele(idx) = sign(tstats.tstat) * idx;

    [~, p2(k), ~, tstats] = ttest2(D2(1:subjects_range(2), k), D2(subjects_range(2) + 1:subjects_range(3),k), .05);
    sele(idx + K) = sign(tstats.tstat) * (idx + K);
    
    [~, p3(k), ~, tstats] = ttest2(D3(1:subjects_range(2), k), D3(subjects_range(2) + 1:subjects_range(3),k), .05);
    sele(idx + 2 * K) = sign(tstats.tstat) * (idx + 2 * K);
    idx = idx + 1;
end 
pp_D = [p1;p2;p3]; 
for k = 1: K     
    [~, p1(k), ~, ~] = ttest2(D1_tilde(1:subjects_range(2), k), D1_tilde(subjects_range(2) + 1:subjects_range(3),k), .05);
    [~, p2(k), ~, ~] = ttest2(D2_tilde(1:subjects_range(2), k), D2_tilde(subjects_range(2) + 1:subjects_range(3),k), .05);
    [~, p3(k), ~, ~] = ttest2(D3_tilde(1:subjects_range(2), k), D3_tilde(subjects_range(2) + 1:subjects_range(3),k), .05);
end 
pp_D_tilde = [p1; p2; p3];

pp_D_sign = pp_D >= 0.05;
pp_D_tilde_sign = pp_D_tilde >= 0.05;
pp_diff = (pp_D_sign - pp_D_tilde_sign) ~= 0;

%%%%%%%%%%%%%%%%% Compute Correlation of Components %%%%%%%%%%%%%%%%%%%%%%%
coef12 = zeros(K, K);
coef13 = zeros(K, K);
coef23 = zeros(K, K);
for i = 1:K
    for j = 1:K
        tmp = corrcoef(Z1(i, :)', Z2(j, :)');
        coef12(i, j) = tmp(1,2);
        tmp = corrcoef(Z1(i, :)', Z3(j, :)');
        coef13(i, j) = tmp(1,2);
        tmp = corrcoef(Z2(i, :)', Z3(j, :)');
        coef23(i, j) = tmp(1,2);
    end
end 

coef_hat12 = zeros(K, K);
coef_hat13 = zeros(K, K);
coef_hat23 = zeros(K, K);
for i = 1:K
    for j = 1:K
        tmp = corrcoef(Z_hat1(i, :)', Z_hat2(j, :)');
        coef_hat12(i, j) = tmp(1,2);
        tmp = corrcoef(Z_hat1(i, :)', Z_hat3(j, :)');
        coef_hat13(i, j) = tmp(1,2);
        tmp = corrcoef(Z_hat2(i, :)', Z_hat3(j, :)');
        coef_hat23(i, j) = tmp(1,2);
    end
end 

Z_hat_groupNotZeroRate = zeros(1, K);
Z_check_rowNotZero = zeros(S, K);
Z_sparsity = zeros(S, K);
D_hat_columnNonZero = zeros(S, K);
for k = 1: K
    for j = 1: num_voxels
        for s = 1: S
             if 0 ~= Z_ks(s, j, k)
                 Z_sparsity(s, k) = Z_sparsity(s, k) + 1;
             end
        end 
        if norm(Z_hat_ks(:, j, k)) ~= 0
            Z_hat_groupNotZeroRate(k) = Z_hat_groupNotZeroRate(k) + 1;
        end
    end
    for s = 1: S
        Z_check_rowNotZero(s, k) = norm(Z_check_ks(s, :, k)) ~= 0;
        D_hat_columnNonZero(s, k) = norm(D_hat(:, (s - 1) * K + k)) ~= 0;
    end
end
Z_hat_groupNotZeroRate = string(Z_hat_groupNotZeroRate./num_voxels * 100) + "%";
Z_sparsity = string(Z_sparsity./num_voxels* 100) + "%";

%%%%%%%%%%%%%%%%% Plot Costs and Change of blocks %%%%%%%%%%%%%%%%%%%%%%%
names = ["costs", "Delta_D"];
for i = 1: size(names, 2)
     fn = fullfile(filePath, names(i) + ".mat");
     load(fn);
end

figure;
plot(Delta_D(2: end));  
xlabel("Iteration Number t"); 
ylabel("Change of D");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = num_subjects;
load('indices.mat');
addpath('spm12');%spm:statistical parametric mapping
fileName=['test']; %filename for saving the figures
fprintf('- Reading in Anatomical Data\n');
anat = spm_read_vols(spm_vol('nsingle_subj_T1_2_2_3.hdr'));
z_choice = 4:5:46;
monthresh = 2; %set z-score threshold (higher threshold means fewer plotted voxels

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Plot Maps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sf=[Z1; Z2; Z3];

%% Create 4-D fMRI data and read in anatomical map
dat3f = zeros(K * 3, 53*63*46); 
dat3f(:,index_) = Sf(:, (1:length(index_)));
dat3f = reshape(dat3f, 3 * K, 53, 63, 46);%why use reshape here?

for k = 1: length(plot_col)
    figure;
    for s = 1: S
        for j = 1: plot_col(k)
            if k == 1
                i = (s - 1) * K + j; % Plot the i-th component
            else
                i = (s - 1) * K + (k - 1) * plot_col(k-1) + j; % Plot the i-th component
            end
            tmp_map(:,:,:) = dat3f(abs(sele(i)), :, :, z_choice);%x,y,time,slice?
            % m=(sign(sele(i))*Sf/stdN(Sf)); % use if you want to normalize based on all voxels (as opposed to only the plotted ones (see line 43))

            subplot_tight(length(plot_col), plot_col(k), (s - 1) * plot_col(k) + j, [0.1,0.02]);
            make_composite(sign(sele(i))*tmp_map/stdN(tmp_map),anat(:,:,z_choice),monthresh); %plot component
            xticks([]);
            yticks([]);
            pp_D_tmp = pp_D';
            xlabel(['{\it p} = ', num2str(pp_D_tmp(i), 2)], 'FontSize', 20);
            if s == 1
                title(['{\it k} = ', num2str(i + ind_start)], 'FontSize', 30, 'FontWeight','normal');
            end
            % make_composite(tmp,anat(:,:,z_choice),monthresh); %plot t-statistics with probability
            clear tmp_map
        end
    end
end 

nonZeroColumn = [];
for k = 1: size(D_hat, 2)
     if norm(D_hat(:, k)) ~= 0
          nonZeroColumn = [nonZeroColumn, k];
     end
end
figure;
num_plot = 3;
for i = 1: min(num_plot, size(nonZeroColumn, 2))
    subplot(1, num_plot, i);
    cdfplot(abs(D_tilde(:, i)));
    hold on;
    cdfplot(abs(D_hat(:, i)));
end
pp_correct = pp_D_sign .* D_hat_columnNonZero;
fprintf("D_hat correct rate %f\n", sum(pp_correct(:))/numel(pp_correct));
% clear; close all;
%% save files as PDF
% 
% print('-depsc', fileName);
% hgsave(fileName);
% save(fileName);
% eps2pdf([fileName, '.eps'], 'C:\Program Files\gs\gs9.14\bin\gswin64.exe');   %for windows, the path of gs should be different
