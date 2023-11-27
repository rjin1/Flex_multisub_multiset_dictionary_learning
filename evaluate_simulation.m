function evaluate_simulation()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load Results And Synthetic %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    filePath = "D:\GitHub\multitask-palm\data\results_saved\subjNorm\0.08 2 0.17 5.64 0.99 6 200 1\2";
    %filePath = "D:\GitHub\multitask-palm\simulation\3_1_1_1\runs\0.06_1.1_0.1_1.02_0.99_6_200_1\1";
    %filePath = "data\E1\2";

    files = ["D", "D_hat", "Z_check_ks", "Z_hat_ks", "costs"];
    for i = 1: size(files, 2)
        fn = fullfile(filePath, files(i) + ".mat");
        load(fn); 
    end 
    files = ["Synthesized_D1", "Synthesized_D2", "Synthesized_D3", "Synthesized_Z1", "Synthesized_Z2", "Synthesized_Z3"];
    for i = 1: size(files, 2)
        fn = fullfile("data", files(i) + ".mat");
        load(fn); 
    end 
    
    S = 3;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Construct Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Z1 = zeros(K, num_voxels);
    Z2 = zeros(K, num_voxels);
    Z3 = zeros(K, num_voxels);
    Z_hat1 = zeros(K, num_voxels);
    Z_hat2 = zeros(K, num_voxels);
    Z_hat3 = zeros(K, num_voxels);
    D_tilde = D - D_hat;
    Z_ks = Z_hat_ks + Z_check_ks;
    for k = 1: K
        Z1(k, :) = Z_ks(1, :, k); 
        Z2(k, :) = Z_ks(2, :, k); 
        Z3(k, :) = Z_ks(3, :, k); 
        Z_hat1(k, :) = Z_hat_ks(1, :, k); 
        Z_hat2(k, :) = Z_hat_ks(2, :, k); 
        Z_hat3(k, :) = Z_hat_ks(3, :, k); 
    end
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
    for k = 1: K     
        [~, p1(k), ~, ~] = ttest2(D1(1:subjects_range(2), k), D1(subjects_range(2) + 1:subjects_range(3),k), .05);
        [~, p2(k), ~, ~] = ttest2(D2(1:subjects_range(2), k), D2(subjects_range(2) + 1:subjects_range(3),k), .05);
        [~, p3(k), ~, ~] = ttest2(D3(1:subjects_range(2), k), D3(subjects_range(2) + 1:subjects_range(3),k), .05);
    end 
    pp_D = [p1;p2;p3]; 
    
    for k = 1: K     
        [~, p1(k), ~, ~] = ttest2(D1_tilde(1:subjects_range(2), k), D1_tilde(subjects_range(2) + 1:subjects_range(3),k), .05);
        [~, p2(k), ~, ~] = ttest2(D2_tilde(1:subjects_range(2), k), D2_tilde(subjects_range(2) + 1:subjects_range(3),k), .05);
        [~, p3(k), ~, ~] = ttest2(D3_tilde(1:subjects_range(2), k), D3_tilde(subjects_range(2) + 1:subjects_range(3),k), .05);
    end 
    pp_D_tilde = [p1; p2; p3];

    for k = 1: K     
        [~, p1(k), ~, ~] = ttest2(D1_hat(1:subjects_range(2), k), D1_hat(subjects_range(2) + 1:subjects_range(3),k), .05);
        [~, p2(k), ~, ~] = ttest2(D2_hat(1:subjects_range(2), k), D2_hat(subjects_range(2) + 1:subjects_range(3),k), .05);
        [~, p3(k), ~, ~] = ttest2(D3_hat(1:subjects_range(2), k), D3_hat(subjects_range(2) + 1:subjects_range(3),k), .05);
    end 
    pp_D_hat = [p1; p2; p3];
    
    for k = 1: K     
        [~, p1(k), ~, ~] = ttest2(Synthesized_D1(1:subjects_range(2), k), Synthesized_D1(subjects_range(2) + 1:subjects_range(3),k), .05);
        [~, p2(k), ~, ~] = ttest2(Synthesized_D2(1:subjects_range(2), k), Synthesized_D2(subjects_range(2) + 1:subjects_range(3),k), .05);
        [~, p3(k), ~, ~] = ttest2(Synthesized_D3(1:subjects_range(2), k), Synthesized_D3(subjects_range(2) + 1:subjects_range(3),k), .05);
    end 
    pp_D_Syn = [p1; p2; p3];
    
    
    %%%%%%%%%%%%%%%%% Compute Correlation of Components %%%%%%%%%%%%%%%%%%%%%%%
    coef12 = zeros(K, K);
    coef13 = zeros(K, K);
    coef23 = zeros(K, K);
    coef1 = zeros(K, K);
    coef2 = zeros(K, K);
    coef3 = zeros(K, K);
    for i = 1:K
        for j = 1:K
            tmp = corrcoef(Z1(i, :)', Z2(j, :)');
            coef12(i, j) = tmp(1,2);
            tmp = corrcoef(Z1(i, :)', Z3(j, :)');
            coef13(i, j) = tmp(1,2);
            tmp = corrcoef(Z2(i, :)', Z3(j, :)');
            coef23(i, j) = tmp(1,2);
            tmp = corrcoef(Z1(i, :)', Synthesized_Z1(j, :)');
            coef1(i, j) = tmp(1,2);
            tmp = corrcoef(Z2(i, :)', Synthesized_Z2(j, :)');
            coef2(i, j) = tmp(1,2);
            tmp = corrcoef(Z3(i, :)', Synthesized_Z3(j, :)');
            coef3(i, j) = tmp(1,2);
        end
    end 
    
    coef_hat12 = zeros(K, K);
    coef_hat13 = zeros(K, K);
    coef_hat23 = zeros(K, K);
    coef_hat1 = zeros(K, K);
    coef_hat2 = zeros(K, K);
    coef_hat3 = zeros(K, K);
    for i = 1:K
        for j = 1:K
            tmp = corrcoef(Z_hat1(i, :)', Z_hat2(j, :)');
            coef_hat12(i, j) = tmp(1,2);
            tmp = corrcoef(Z_hat1(i, :)', Z_hat3(j, :)');
            coef_hat13(i, j) = tmp(1,2);
            tmp = corrcoef(Z_hat2(i, :)', Z_hat3(j, :)');
            coef_hat23(i, j) = tmp(1,2);
            tmp = corrcoef(Z_hat1(i, :)', Synthesized_Z1(j, :)');
            coef_hat1(i, j) = tmp(1,2);
            tmp = corrcoef(Z_hat2(i, :)', Synthesized_Z2(j, :)');
            coef_hat2(i, j) = tmp(1,2);
            tmp = corrcoef(Z_hat3(i, :)', Synthesized_Z3(j, :)');
            coef_hat3(i, j) = tmp(1,2);
        end
    end 
    
    %%%%%%%%%%%%%%%%% Check if the permutations are right %%%%%%%%%%%%%%%%%%%%%%%
    [maxVal1, maxIdx1] = max(abs(coef1)');
    [maxVal2, maxIdx2] = max(abs(coef2)');
    [maxVal3, maxIdx3] = max(abs(coef3)');
    allRight = true;
    for k = 1: K
        if pp_D(1, k) < 0.05 && pp_D_Syn(1, maxIdx1(k)) >= 0.05 || pp_D(1, k) >= 0.05 && pp_D_Syn(1, maxIdx1(k)) < 0.05
            allRight = false;
            fprintf("Task #1, Component %d not right w.r.t p values !\n", k);
        end
        if pp_D(2, k) < 0.05 && pp_D_Syn(2, maxIdx2(k)) >= 0.05 || pp_D(2, k) >= 0.05 && pp_D_Syn(2, maxIdx2(k)) < 0.05
            allRight = false;
            fprintf("Task #2, Component %d not right w.r.t p values !\n", k);
        end
        if pp_D(3, k) < 0.05 && pp_D_Syn(3, maxIdx3(k)) >= 0.05 || pp_D(3, k) >= 0.05 && pp_D_Syn(3, maxIdx3(k)) < 0.05
            allRight = false;
            fprintf("Task #3, Component %d not right w.r.t p values !\n", k);
        end
    end
    
    K_bM = 0; % number of common components
    [maxVal12, maxIdx12] = max(abs(coef12)');
    [maxVal13, maxIdx13] = max(abs(coef13)');  
    common_com = [];
    correlation_threshold = 0.3;
    for k = 1: K
        [a, b, c] = deal(k, maxIdx12(k), maxIdx13(k));
        if abs(coef12(a, b)) >= correlation_threshold && abs(coef13(a, c)) >= correlation_threshold && abs(coef23(b, c)) < correlation_threshold
            allRight = false;
            fprintf("Component %d not right w.r.t mutual correlation !\n", k);
        end
        if abs(coef12(a, b)) >= correlation_threshold && abs(coef13(a, c)) >= correlation_threshold && abs(coef23(b, c)) > correlation_threshold
            common_com = [common_com,[a; b; c]];
            K_bM = K_bM + 1;
        end
    end
    if allRight
        disp("Permutations are right for all components! ");
    end

    
    %%%%%%%%%%%%%%%%% Check Sparsity %%%%%%%%%%%%%%%%%%%%%%%
    Z_group_sparsity_com = zeros(1, K_bM);
    i = 1;
    for k = 1: K
         [a, b, c] = deal(k, maxIdx12(k), maxIdx13(k));
         if abs(coef12(a, b)) >= correlation_threshold
             z = [Z1(a, :); Z2(b, :); Z3(c, :)]; 
             for j = 1: num_voxels
                 if norm(z(:, j)) == 0
                      Z_group_sparsity_com(i) = Z_group_sparsity_com(i) + 1;
                 end
             end
             Z_group_sparsity_com(i) = Z_group_sparsity_com(i)/num_voxels;
             i = i + 1;
         end
         
    end
    
    Z_check_rowSparsity = zeros(S, K);
    Z_hat_groupSparsity = zeros(1, K);
    D_hat_colNorm = zeros(S, K);
    Z_groupSparsity = zeros(1, K);
    Z_hat_norm = zeros(S, K);
    Z_sparsity = zeros(S, K);
    for k = 1: K
        for j = 1: num_voxels
            for s = 1: S
                 if 0 == Z_check_ks(s, j, k)
                     Z_check_rowSparsity(s, k) = Z_check_rowSparsity(s, k) + 1;
                 end
                 if 0 == Z_ks(s, j, k)
                     Z_sparsity(s, k) = Z_sparsity(s, k) + 1;
                 end
            end 
            if norm(Z_hat_ks(:, j, k)) == 0
                Z_hat_groupSparsity(k) = Z_hat_groupSparsity(k) + 1;
            end
            if norm(Z_ks(:, j, k)) == 0
                Z_groupSparsity(k) = Z_groupSparsity(k) + 1;
            end
        end
        for s = 1: S
            D_hat_colNorm(s, k) = norm(D_hat(:, (s - 1) * K + k)); 
            Z_hat_norm(s, k) = norm(Z_hat_ks(s, :, k));
        end
        Z_hat_k = Z_hat_ks(:, :, k);
        Z_check_k = Z_check_ks(:, :, k);
        Z_k = Z_hat_ks(:, :, k) + Z_check_ks(:, :, k);
    end
    Z_hat_groupSparsity = Z_hat_groupSparsity./num_voxels;
    Z_check_rowSparsity = Z_check_rowSparsity./num_voxels;
    Z_groupSparsity = Z_groupSparsity./num_voxels;
    Z_sparsity = Z_sparsity./num_voxels;
    
    flag_component = true;
    flag_D = true;
    for k = 1: K
        for s = 1: S
            isCommon = any(common_com(s, :) == k);
            isRowEmpty = Z_check_rowSparsity(s, k) == 1;
            if isCommon
                if ~isRowEmpty
                    fprintf("Task #%d, %d-th row of Z_check is not empty, but it's a common component !\n", s, k);
                    flag_component = false;
                end
            else
                if isRowEmpty
                    fprintf("Task #%d, %d-th row of Z_check is empty, but it's a distinct component !\n", s, k);
                    flag_component = false;
                end
            end
            
        end
    end
    for k = 1: K
        for s = 1: S
            isEmptyColmn =  D_hat_colNorm(s, k) == 0;
            isDiscriminative = pp_D(s, k) < 0.05;
            if isDiscriminative
                if ~isEmptyColmn
                    fprintf("Task #%d, %d-th column of D_hat is not empty, but it's a discriminative component !\n", s, k);
                    flag_D = false;
                end
            else
                if isEmptyColmn
                    fprintf("Task #%d, %d-th column of D_hat is empty, but it's not a discriminative component !\n", s, k);
                    flag_D = false;
                end
            end
        end
    end
    if flag_component
        disp("Row sparsity of components all match"); 
    end
    if flag_D
        disp("Column sparsity of dictionaries all match"); 
    end
    disp("Row sparsity of Z_check");
    disp(Z_check_rowSparsity);
    disp("Column norm of D_hat");
    disp(D_hat_colNorm);
    disp("P values of D");
    disp(pp_D);
    
    
    
    %%%%%%%%%%%%%%%%% Plot Costs and Change of blocks %%%%%%%%%%%%%%%%%%%%%%%
    names = ["costs", "Delta_D", "Delta_D_hat", "Delta_Z_hat_ks", "Delta_Z_check_ks"];
    for i = 1: size(names, 2)
         fn = fullfile(filePath, names(i) + ".mat");
         load(fn);
    end
    figure;
    plot(costs(2: end));  
    xlabel("Iteration Number t"); 
    ylabel("cost\_t");
    
    figure;
    plot(Delta_D(2: end));  
    xlabel("Iteration Number t"); 
    ylabel("Change of D");
    
    figure;
    plot(cost_D_tildes(2: end));  
    xlabel("Iteration Number t"); 
    ylabel("Cost of cost\_D\_tildes");
    
    figure;
    plot(cost_D_hats(2: end));  
    xlabel("Iteration Number t"); 
    ylabel("Cost of cost\_D\_hats");
    
    figure;
    plot(cost_recons(2: end));  
    xlabel("Iteration Number t"); 
    ylabel("Cost of cost\_recons");
    
    figure;
    plot(cost_Zs(2: end));  
    xlabel("Iteration Number t"); 
    ylabel("Cost of cost\_Zs");
    
    figure;
    plot(Delta_D_hat(2: end));  
    xlabel("Iteration Number t"); 
    ylabel("Change of D\_hat");
    
    for k = 3:4
        figure;
        plot(Delta_Z_hat_ks(2: end , k));  
        xlabel("Iteration Number t"); 
        ylabel("Change of Z\_hat_" + string(k));
        
        figure;
        plot(Delta_Z_check_ks(2: end, k));  
        xlabel("Iteration Number t"); 
        ylabel("Change of Z\_check_" + string(k));
    end
    
    clear; close all;
    %{
    changeOfCost = [];
    slope = 1e-2;
    interval = 50;
    theshold_break = slope * interval;
    for i = interval + 1: size(costs, 1)
        changeOfCost = [changeOfCost,  abs(costs(i) - costs(i - interval))/interval];
    end
    figure;
    plot(changeOfCost(interval: end));  
    xlabel("Iteration Number t"); 
    ylabel("Change of cost");
    %}
end