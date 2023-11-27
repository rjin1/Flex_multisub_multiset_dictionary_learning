function [] = simulation_generator()
    rng(10);
    K = 6;
    K_bM = 4;
    K_tM = 2;
    M = 247; % total number of subjects
    Y_range = [0, 138, 247];
    N = 10000; % number of column of Z
        
    step = 1;
    Synthesized_D1 = normrnd(0, 1, [M, K]);
    Synthesized_D2 = normrnd(0, 1, [M, K]);
    Synthesized_D3 = normrnd(0, 1, [M, K]);
    Synthesized_D3(:, 4) = normrnd(0, 1, [M, 1]); 
    for i = Y_range(2) + 1: Y_range(3)
        Synthesized_D1(i, 1) = Synthesized_D1(i, 1) + step;
        Synthesized_D2(i, 1) = Synthesized_D2(i, 5) + step;
        Synthesized_D3(i, 1) = Synthesized_D3(i, 1) + step;

        Synthesized_D1(i, 2) = Synthesized_D1(i, 2) + step;
        Synthesized_D2(i, 2) = Synthesized_D2(i, 2) + step;

        Synthesized_D3(i, 3) = Synthesized_D3(i, 3) + step;

        Synthesized_D1(i, 5) = Synthesized_D1(i, 5) + step;
        Synthesized_D2(i, 5) = Synthesized_D2(i, 5) + step;
        Synthesized_D3(i, 5) = Synthesized_D3(i, 5) + step;
    end
    % dkdk<=1
    for k = 1:K
        Synthesized_D1(:, k) = Synthesized_D1(:, k)/max(norm(Synthesized_D1(:, k)), 1);
        Synthesized_D2(:, k) = Synthesized_D2(:, k)/max(norm(Synthesized_D2(:, k)), 1);
        Synthesized_D3(:, k) = Synthesized_D3(:, k)/max(norm(Synthesized_D3(:, k)), 1);
    end 
    
    %%%%%%%%%%%%%% check the discrimination ability of dictionary %%%%%%%%%%%%%%%%%%%%
    p1 = zeros(1, K); 
    p2 = zeros(1, K);
    p3 = zeros(1, K);  
    for k = 1: K
        [~, p1(k), ~, ~] = ttest2(Synthesized_D1(1:Y_range(2),k),Synthesized_D1(Y_range(2) + 1:Y_range(3),k),.05);
        [~, p2(k), ~, ~] = ttest2(Synthesized_D2(1:Y_range(2),k),Synthesized_D2(Y_range(2) + 1:Y_range(3),k),.05);
        [~, p3(k), ~, ~] = ttest2(Synthesized_D3(1:Y_range(2),k),Synthesized_D3(Y_range(2) + 1:Y_range(3),k),.05);
    end
    pp = [p1;p2;p3]; 
    figure; imagesc(Synthesized_D1); xlabel("Profile Index"); ylabel("Image Scalar Of Syhthesized D1");
    figure; imagesc(Synthesized_D2); xlabel("Profile Index"); ylabel("Image Scalar Of Syhthesized D2");
    figure; imagesc(Synthesized_D3); xlabel("Profile Index"); ylabel("Image Scalar Of Syhthesized D3");
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulate the Z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Synthesized_Z1 = normrnd(0, 1, [K, N]);
    Synthesized_Z2 = normrnd(0, 1, [K, N]); 
    Synthesized_Z3 = normrnd(0, 1, [K, N]);
 
    rho = 0.97;
    for i = 1:1:K_bM
        z = normrnd(0, 1, [1, N]);
        w1 = normrnd(0, 1, [1, N]);
        w2 = normrnd(0, 1, [1, N]);
        w3 = normrnd(0, 1, [1, N]);
        Synthesized_Z1(i, :) = sqrt(rho)*z + sqrt(1- rho)*w1;
        Synthesized_Z2(i, :) = sqrt(rho)*z + sqrt(1- rho)*w2;
        Synthesized_Z3(i, :) = sqrt(rho)*z + sqrt(1- rho)*w3;
    end
      

    %%%%%%%% make the components sparse by setting the absolute smallest 50% of
    %%%%%%%% elements to be 0
    sparsity = 0.5;
    for i = 1:K
        threshold = prctile(abs(Synthesized_Z1(i, :)), 100*sparsity);
        for j = 1:N
            Synthesized_Z1(i, j) = max(abs(Synthesized_Z1(i, j)) - threshold, 0)*sign(Synthesized_Z1(i, j));
        end
        
        threshold = prctile(abs(Synthesized_Z2(i, :)), 100*sparsity);
        for j = 1:N
            Synthesized_Z2(i, j) = max(abs(Synthesized_Z2(i, j)) - threshold, 0)*sign(Synthesized_Z2(i, j));
        end
         
        threshold = prctile(abs(Synthesized_Z3(i, :)), 100*sparsity);
        for j = 1:N
            Synthesized_Z3(i, j) = max(abs(Synthesized_Z3(i, j)) - threshold, 0)*sign(Synthesized_Z3(i, j));
        end
    end 
   
    figure; histogram(Synthesized_Z1(1, :));
    figure; histogram(Synthesized_Z2(1, :));
    figure; histogram(Synthesized_Z3(1, :)); 
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Check the correlation of
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% components
    coefSynZ12 = zeros(K, K);
    coefSynZ13 = zeros(K, K);
    coefSynZ23 = zeros(K, K);
    for i = 1:1:K
        for j = 1:1:K
            tmp = corrcoef(Synthesized_Z1(i, :)', Synthesized_Z2(j, :)');
            coefSynZ12(i, j) = tmp(1, 2);
            
            tmp = corrcoef(Synthesized_Z1(i, :)', Synthesized_Z3(j, :)');
            coefSynZ13(i, j) = tmp(1, 2);
            
            tmp = corrcoef(Synthesized_Z2(i, :)', Synthesized_Z3(j, :)');
            coefSynZ23(i, j) = tmp(1, 2);
        end
    end

    sparse_threshold = 0;
    sparsities = zeros(K, 3);
    for i = 1:K
        for j = 1:N
            if abs(Synthesized_Z1(i, j)) == sparse_threshold
                 sparsities(i, 1) = sparsities(i, 1) + 1;
            end
            if abs(Synthesized_Z2(i, j)) == sparse_threshold
                 sparsities(i, 2) = sparsities(i, 2) + 1;
            end
            if abs(Synthesized_Z3(i, j)) == sparse_threshold
                 sparsities(i, 3) = sparsities(i, 3) + 1;
            end
        end
    end
    sparsities = sparsities./N;
    
    group_sparsities = zeros(K, 1);
    for i = 1:K
        for j = 1:N
            if abs(Synthesized_Z1(i, j)) == sparse_threshold && abs(Synthesized_Z2(i, j)) == sparse_threshold && abs(Synthesized_Z3(i, j)) == sparse_threshold
                 group_sparsities(i) = group_sparsities(i) + 1;
            end
        end 
    end
    group_sparsities = group_sparsities./N;
        
    Y = Synthesized_D1*Synthesized_Z1;
    Y_std = std(Y(:));
    gaussian_noise = normrnd(0, sqrt(5)*Y_std/5, [M, N]);
    Y = Y + gaussian_noise;
    fn = fullfile("data", "Synthesized_dataset1.mat");
    save(fn, "Y", "Y_range");
    
    Y = Synthesized_D2*Synthesized_Z2; 
    Y_std = std(Y(:));
    gaussian_noise = normrnd(0, sqrt(5)*Y_std/5, [M, N]);
    Y = Y + gaussian_noise;
    fn = fullfile("data", "Synthesized_dataset2.mat");
    save(fn, "Y", "Y_range");
    
    Y = Synthesized_D3*Synthesized_Z3;
    Y_std = std(Y(:));
    gaussian_noise = normrnd(0, sqrt(5)*Y_std/5, [M, N]);
    Y = Y + gaussian_noise;
    fn = fullfile("data", "Synthesized_dataset3.mat");
    save(fn, "Y", "Y_range");
    
    fn = fullfile("data", "Synthesized_D1.mat");
    save(fn, "Synthesized_D1");
    fn = fullfile("data", "Synthesized_D2.mat");
    save(fn, "Synthesized_D2");
    fn = fullfile("data", "Synthesized_D3.mat");
    save(fn, "Synthesized_D3");
    
    fn = fullfile("data", "Synthesized_Z1.mat");
    save(fn, "Synthesized_Z1");
    fn = fullfile("data", "Synthesized_Z2.mat");
    save(fn, "Synthesized_Z2");
    fn = fullfile("data", "Synthesized_Z3.mat");
    save(fn, "Synthesized_Z3");
    clear; close all;
end