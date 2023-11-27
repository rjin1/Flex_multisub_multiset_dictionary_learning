filePath = ".\results_toplot\visual_result.mat";
plot_col = [3,3,3];
load(filePath)
S = 3;
K = 3;
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

            subplot_tight(length(plot_col), plot_col(k), (s - 1) * plot_col(k) + j, [0.075, 0.0001]);
            make_composite(sign(sele(i))*tmp_map/stdN(tmp_map),anat(:,:,z_choice),monthresh); %plot component
            xticks([]);
            yticks([]);
            pp_D_tmp = pp_D';
            xlabel(['{\it p} = ', num2str(pp_D_tmp(i), 2)], 'FontSize', 20);
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
