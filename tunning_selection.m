Files=dir("D:\GitHub\multitask-palm\simulation\3_1_1_1\ranges1\tunning_results");
idx = 3;    
results = []; 

while idx <= length(Files)
    sp = strsplit(Files(idx).name, '_');
    [lambda1, lambda2, lambda3, lambda4, alpha, acc_com_NC, acc_com_NN, recon, cost, seed] = deal(str2double(sp{1}), str2double(sp{2}), str2double(sp{3}), str2double(sp{4}), str2double(sp{5}), str2double(sp{6}), str2double(sp{7}), str2double(sp{8}), str2double(sp{9}), str2double(sp{10}));
    results = [results; [lambda1, lambda2, lambda3,  lambda4, alpha, acc_com_NC, acc_com_NN, recon, cost, seed]];
    idx = idx + 1;
end     

sorted_lambda1 = sortrows(results, 1, 'ascend');
sorted_lambda2 = sortrows(results, 2, 'ascend');
sorted_lambda3 = sortrows(results, 3, 'ascend');
sorted_lambda4 = sortrows(results, 4, 'ascend');
sorted_alpha = sortrows(results, 5, 'ascend');

sorted_com_NC = sortrows(results, 6, 'descend');
sorted_com_NN = sortrows(results, 7, 'descend');
sorted_recon = sortrows(results, 8, 'ascend');
sorted_cost = sortrows(results, 9, 'ascend');

clear; close all;

