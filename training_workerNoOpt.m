function training_workerNoOpt(E_id)
    addpath('utils');  
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% read parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fid = fopen("data/training_parameters.txt", "r"); 
    pars = fgetl(fid);
    fclose(fid);
    pars_cells = strsplit(pars, " ");
    lambda1 = str2num(pars_cells{1});
    lambda2 = str2num(pars_cells{2});
    lambda3 = str2num(pars_cells{3}); 
    lambda4_coef = str2num(pars_cells{4});
    alpha = str2num(pars_cells{5});
    K = str2num(pars_cells{6});
    numRuns = str2num(pars_cells{7});
    isSimulation = str2num(pars_cells{8}); % 1 for simulation, 0 for real datasets
    [gamma_k, gamma_1, gamma_2] = deal(1, 1, 1); % Fix gamma_hat_k = gamma_check_k = gamma_k = gamma_1 = gamma_2 = 1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load Datasets %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    datasets = [];
    if 0 == isSimulation
        datasets = ["New_AOD3", "New_SIRP", "New_SM"];
    else
        datasets = ["Synthesized_dataset_test_1", "Synthesized_dataset_test_2", "Synthesized_dataset_test_3"];
    end
   %{
    data_fn = fullfile('data', strcat(datasets(1), '.mat'));
    load(data_fn);
    X1 = normc(Y')';
    subjects_range = Y_range;  % subjects_range is same for all three datasets becase they have the same range.
    data_fn = fullfile('data', strcat(datasets(2), '.mat'));
    load(data_fn);
    X2 = normc(Y')';
    data_fn = fullfile('data', strcat(datasets(3), '.mat'));
    load(data_fn);
    X3 = normc(Y')';
     %}
    N_train = [97, 76]; % number of subjects to train per class
    N_vals = [41, 33]; % number of subjects to validate of [HC, SZ].
    
    %%%%%%%%%%%%%%%%%%%%%%% Fetch Datasets and Their Sizes. The Three Datasets Have The Same Sizes %%%%%%%%%%%%%%%%%%%%%%%
    seed = 10000;
    [X1, ~, label_train, ~, idx_train] = pickTrainTest(datasets(1), N_train, N_vals, seed);
    [X2, ~, ~, ~, ~] = pickTrainTest(datasets(2), N_train, N_vals, seed);
    [X3, ~, ~, ~, ~] = pickTrainTest(datasets(3), N_train, N_vals, seed);
    subjects_range = label_to_range(label_train);
    [num_subjects, num_voxels] = size(X1);
    X = [X1', X2', X3']';
    X_prime = [X1, X2, X3];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Define H of Fisher's Discri. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    C = 2; % Number of classes.
    H1 = [];
    for c = 1: C
        Nc = subjects_range(c+1) - subjects_range(c);
        H1(end+1: end + Nc, end+1: end + Nc) = ones(Nc)/Nc;
    end
    N = subjects_range(C+1);
    H2 = ones(N)/N;
    H = 2 * eye(N)-2 * H1 + H2;
    L_D_hat = lambda3 * norm(H, 'fro');
    rho_1 = gamma_1 * L_D_hat;
    lambda4 = round(rho_1 * lambda4_coef , 2);

    [S, threshold_break] = deal(3, 1e-4); % number of datasets
    for i = 1:K
        ind(1+(i-1)*S : i*S) = i:K:(S-1)*K+i;
    end
    while true
        % read how many runs we have so far
        filename = "data/run_done.txt"; % define file name
        fid = fopen(filename, "r+");
        formatSpec = '%f';
        lines = fscanf(fid, formatSpec);
        seed = lines(end);
        if seed > numRuns
            break;
        end
        fprintf(fid, "\n%d", seed + 1);
        fclose(fid);
        fprintf('%s, lambda1 = %f, lambda2 = %f, lambda3 = %f, lambda4_coef = %f, alpha = %f, K = %d, seed = %d\n', E_id, lambda1, lambda2, lambda3, lambda4_coef, alpha, K, seed);
        rng(seed, 'twister'); 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Z_hat_ks = zeros(S, num_voxels, K);
        Z_check_ks = zeros(S, num_voxels, K);
        for k = 1: K
            Z_hat_ks(:, :, k) = normrnd(0, 1, [S, num_voxels]); 
            Z_check_ks(:, :, k) = normrnd(0, 1, [S, num_voxels]); 
        end
        Z_ks = Z_hat_ks + Z_check_ks;
        Z1 = squeeze(Z_ks(1, :, :))';
        Z2 = squeeze(Z_ks(2, :, :))';
        Z3 = squeeze(Z_ks(3, :, :))';
        
        D = normrnd(0, 1, [num_subjects, S * K]);
        D_hat = normrnd(0, 1, [num_subjects, S * K]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Iterations Begin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Delta_D = [];
        Delta_D_hat = [];
        Delta_Z_hat_ks = zeros(1, K);
        Delta_Z_check_ks = zeros(1, K);
        % [cost, cost_recon, cost_Z, cost_D_tilde, cost_D_hat] = deal(1e100, 1e100, 1e100, 1e100, 1e100);
        tic;
        [t, breakLoop] = deal(1, true);
        while breakLoop
            D_old = D;
            D_hat_old = D_hat;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update {Z_hat_k} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for k = 1: K
                D_k = blkdiag(D(:, k), D(:, K + k), D(:, 2 * K + k));
                L_hat_k = norm(D_k' * D_k, 'fro');
                rho_hat_k = gamma_k * L_hat_k;
                X_hat_k = X;  
                for k_prime = 1: k - 1
                    D_k_prime = blkdiag(D(:, k_prime), D(:, K + k_prime), D(:, 2 * K + k_prime));
                    X_hat_k = X_hat_k - D_k_prime * (Z_hat_ks(:, :, k_prime) + Z_check_ks(:, :, k_prime));
                end
                for k_prime = k + 1: K
                    D_k_prime = blkdiag(D(:, k_prime), D(:, K + k_prime), D(:, 2 * K + k_prime));
                    X_hat_k = X_hat_k - D_k_prime * (Z_hat_ks(:, :, k_prime) + Z_check_ks(:, :, k_prime));
                end
                X_hat_k = X_hat_k - D_k * Z_check_ks(:, :, k);
                U_hat_k = Z_hat_ks(:, :, k) - 1/rho_hat_k * D_k' * (D_k * Z_hat_ks(:, :, k) - X_hat_k); 
                for j = 1: num_voxels
                    Z_hat_ks(:, j, k) = max(0, 1 - lambda1/(rho_hat_k * norm(U_hat_k(:, j)))) * U_hat_k(:, j);
                end 
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update {Z_check_k} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for k = 1: K
                D_k = blkdiag(D(:, k), D(:, K + k), D(:, 2 * K + k));
                L_check_k = norm(D_k' * D_k, 'fro');
                rho_check_k = gamma_k * L_check_k;
                X_check_k = X;
                for k_prime = 1: k - 1
                    D_k_prime = blkdiag(D(:, k_prime), D(:, K + k_prime), D(:, 2 * K + k_prime));
                    X_check_k = X_check_k - D_k_prime * (Z_hat_ks(:, :, k_prime) + Z_check_ks(:, :, k_prime));
                end
                for k_prime = k + 1: K
                    D_k_prime = blkdiag(D(:, k_prime), D(:, K + k_prime), D(:, 2 * K + k_prime));
                    X_check_k = X_check_k - D_k_prime * (Z_hat_ks(:, :, k_prime) + Z_check_ks(:, :, k_prime));
                end
                X_check_k = X_check_k - D_k * Z_hat_ks(:, :, k);
                U_check_k = Z_check_ks(:, :, k) - 1/rho_check_k * D_k' * (D_k * Z_check_ks(:, :, k) - X_check_k);
                %%%%%%%%%%%%%%%%%%%%%%%% define tmp1 and tmp2 to reduce
                %%%%%%%%%%%%%%%%%%%%%%%% computation
                tmp1 = lambda2 * (1 - alpha) / rho_check_k;
                tmp2 = lambda2 * alpha / rho_check_k;
                parfor i = 1: S
                    U_check_prime_k = zeros(1, num_voxels);
                    for j = 1: num_voxels
                        U_check_prime_k(j) = max(0, 1 - tmp1 / abs(U_check_k(i, j))) * U_check_k(i, j);
                    end
                    Z_check_ks(i, :, k) = max(0, 1 - tmp2 / norm(U_check_prime_k)) * U_check_prime_k;
                end
           end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get {Z^s} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Z_ks = Z_hat_ks + Z_check_ks;
            for k = 1: K
                Z1(k, :) = Z_ks(1, :, k); 
                Z2(k, :) = Z_ks(2, :, k); 
                Z3(k, :) = Z_ks(3, :, k); 
            end 
 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update D_hat %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            W_hat = D_hat - lambda3/rho_1 * H * (D_hat - D);
            for j = 1: S * K
                D_hat(:, j) = max(0, 1 - lambda4/(rho_1 * norm(W_hat(:, j)))) * W_hat(:, j);
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Z = blkdiag(Z1, Z2, Z3);
            L_D = norm(Z * Z', 'fro') + lambda3 * norm(H, 'fro');
            rho_2 = gamma_2 * L_D;
            W = D - 1/rho_2 * ((D * Z - X_prime) * Z' + lambda3 * H * (D - D_hat));
            for j = 1: S * K
                D(:, j) = 1/max(1, norm(W(:, j))) * W(:, j); 
            end
           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute Changes And Costs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            change_D = sum(sum(abs(D - D_old)))/(num_subjects * S * K);
            Delta_D = [Delta_D, change_D];
            change_D_hat = sum(sum(abs(D_hat - D_hat_old)))/(num_subjects * S * K);
            Delta_D_hat = [Delta_D_hat, change_D_hat];
            if 0 == mod(t, 30)
                fprintf("Iteration: %d finished, time took: %f, change of D: %f, break threshold: %f, \n", t, toc/3600, change_D, threshold_break);
            end
            if change_D < threshold_break
                % [cost_recon, cost_Z, cost_D_tilde, cost_D_hat, cost] = cost_cal(X1, X2, X3, D, D_hat, Z1, Z2, Z3, K, Z_hat_ks, Z_check_ks, lambda1, lambda2, lambda3, lambda4, alpha, num_voxels, subjects_range);        
                fprintf("Training Took %f hours at iteration %d. \n", toc/3600, t);
                break; 
            end
            %{
            % save temp results
            if mod(t, 200) == 0
                path = E_id + "/" + string(seed) + "/" + string(t);
                mkdir("data", path);
                fn = fullfile("data", path, "D.mat");
                save(fn, "D", "num_subjects", "K", "subjects_range"); 
                fn = fullfile("data", path, "D_hat.mat");
                save(fn, "D_hat"); 
                fn = fullfile("data", path, "Z_hat_ks.mat");
                save(fn, "Z_hat_ks", "num_voxels");
                fn = fullfile("data", path, "Z_check_ks.mat");
                save(fn, "Z_check_ks"); 
                fn = fullfile("data", path, "Delta_D.mat");
                save(fn, "Delta_D"); 
                fn = fullfile("data", path, "Delta_D_hat.mat");
                save(fn, "Delta_D_hat"); 
                fn = fullfile("data", path, "Delta_Z_hat_ks.mat");
                save(fn, "Delta_Z_hat_ks"); 
                fn = fullfile("data", path, "Delta_Z_check_ks.mat");
                save(fn, "Delta_Z_check_ks"); 
                fn = fullfile("data", path, "costs.mat");
                save(fn, "cost", "cost_recon", "cost_Z", "cost_D_tilde", "cost_D_hat");
                fn = fullfile("data", path, "lambda4_" + string(lambda4) + "_.mat");
                save(fn, "lambda4");
                fn = fullfile("data", path, "idx_train.mat");
                save(fn, "idx_train");
                disp("Finished");
            end
            %}
            t = t + 1;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mkdir("data", E_id + "/" + string(seed));
        fn = fullfile("data", E_id, string(seed), "D.mat");
        save(fn, "D", "num_subjects", "K", "subjects_range"); 
        fn = fullfile("data", E_id, string(seed), "D_hat.mat");
        save(fn, "D_hat"); 
        fn = fullfile("data", E_id, string(seed), "Z_hat_ks.mat");
        save(fn, "Z_hat_ks", "num_voxels");
        fn = fullfile("data", E_id, string(seed), "Z_check_ks.mat");
        save(fn, "Z_check_ks"); 
        fn = fullfile("data", E_id, string(seed), "Delta_D.mat");
        save(fn, "Delta_D"); 
        fn = fullfile("data", E_id, string(seed), "Delta_D_hat.mat");
        save(fn, "Delta_D_hat"); 
        % fn = fullfile("data", E_id, string(seed), "costs.mat");
        % save(fn, "cost", "cost_recon", "cost_Z", "cost_D_tilde", "cost_D_hat");
        fn = fullfile("data", E_id, string(seed), "lambda4_" + string(lambda4) + "_.mat");
        save(fn, "lambda4");
        fn = fullfile("data", E_id, string(seed), "idx_train.mat");
        save(fn, "idx_train");
        disp("Finished");
    end
end 
 
function [cost_recon, cost_Z, cost_D_tilde, cost_D_hat, cost_t] = cost_cal(X1, X2, X3, D, D_hat, Z1, Z2, Z3, K, Z_hat_ks, Z_check_ks, lambda1, lambda2, lambda3, lambda4, alpha, num_voxels, subjects_range)
    D1 = D(:, 1: K);
    D2 = D(:, 1 + K: 2 * K);
    D3 = D(:, 1 + 2 * K: 3 * K); 
    cost_recon = 1/2 * norm(X1 - D1 * Z1, 'fro')^2 +  1/2 * norm(X2 - D2 * Z2, 'fro')^2 +  1/2 * norm(X3 - D3 * Z3, 'fro')^2;
    cost_Z = 0;
    parfor k = 1: K
        for j = 1: num_voxels
            cost_Z = cost_Z + lambda1 * norm(Z_hat_ks(:, j, k));
        end 
        for s = 1: size(Z_check_ks, 1)
            cost_Z = cost_Z + lambda2 * alpha * norm(Z_check_ks(s, :, k));
        end
        cost_Z = cost_Z + lambda2 * (1 - alpha) * sum(sum(abs(Z_check_ks(:, :, k)))); 
    end
    cost_D_tilde = 1/2 * lambda3 * FDDL_discriminative((D - D_hat)', subjects_range);
    cost_D_hat = 0;
    parfor i = 1: size(D, 2)
        cost_D_hat = cost_D_hat + lambda4 * norm(D_hat(:, i)); 
    end
    cost_t = cost_recon + cost_Z + cost_D_tilde + cost_D_hat;
end