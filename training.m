function training()
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
    if 0 == isSimulation
        datasets = ["New_AOD3", "New_SIRP", "New_SM"];
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
    else
        datasets = ["Synthesized_dataset1", "Synthesized_dataset2", "Synthesized_dataset3"];
        seed = 10;
        N_train = [110, 87]; % number of subjects to train per class
        N_vals = [28, 22]; % number of subjects to validate of [HC, SZ].
        [X1, ~, label_train, ~, ~] = pickTrainTest(datasets(1), N_train, N_vals, seed);
        subjects_range = label_to_range(label_train);
        [X2, ~, ~, ~, ~] = pickTrainTest(datasets(2), N_train, N_vals, seed);
        [X3, ~, ~, ~, ~] = pickTrainTest(datasets(3), N_train, N_vals, seed);
    end
   
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

    [S, threshold_break, MAXITER] = deal(3, 1e-5, 1200); % number of datasets
    for i = 1:K
        ind(1+(i-1)*S : i*S) = i:K:(S-1)*K+i;
    end
    indices = reshape(ind, S, K)';
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
        fprintf('lambda1 = %f, lambda2 = %f, lambda3 = %f, lambda4_coef = %f, alpha = %f, K = %d, seed = %d\n',lambda1, lambda2, lambda3, lambda4_coef, alpha, K, seed);
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
            X_check_k = X;
            Z = [Z1; Z2; Z3];
            D_tmp = blkdiag(D(:, 1: K),D(:, 1 + K: 2 * K),D(:, 1 + 2 * K: 3 * K));
            D_1 = D_tmp(:, indices(1, :));
            X_hat_k = X - D_tmp(:, ind) * Z(ind,:) + D_1 * Z_hat_ks(:, :, 1);
            for k = 1: K
                D_k = D_tmp(:, indices(k, :));
                L_hat_k = norm(D_k' * D_k, 'fro');
                rho_hat_k = gamma_k * L_hat_k;  
                if k > 1
                    D_k_minus_1 = D_tmp(:, indices(k - 1, :));
                    tmp = D_k_minus_1 * Z_hat_ks(:, :, k - 1) - X_hat_k;
                    X_hat_k = D_k * Z_hat_ks(:, :, k) - tmp;
                    U_hat_k = Z_hat_ks(:, :, k) - 1/rho_hat_k * D_k' * tmp;
                    Z_hat_ks(:, :, k) = max(0, 1 - (lambda1/rho_hat_k) ./ vecnorm(U_hat_k)) .* U_hat_k;
                else
                    U_hat_k = Z_hat_ks(:, :, k) - 1/rho_hat_k * D_k' * (D_k * Z_hat_ks(:, :, k) - X_hat_k);
                    Z_hat_ks(:, :, k) = max(0, 1 - (lambda1/rho_hat_k) ./ vecnorm(U_hat_k)) .* U_hat_k;
                end
                X_check_k = X_check_k - D_k * (Z_hat_ks(:, :, k) + Z_check_ks(:, :, k));
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update {Z_check_k} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            X_check_k = X_check_k + D_1 * Z_check_ks(:, :, 1);
            for k = 1: K
                D_k = D_tmp(:, indices(k, :));
                L_check_k = norm(D_k' * D_k, 'fro');
                rho_check_k = gamma_k * L_check_k;
                if k > 1
                    D_k_minus_1 = D_tmp(:, indices(k - 1, :));
                    tmp = D_k_minus_1 * Z_check_ks(:, :, k - 1) - X_check_k;
                    X_check_k = D_k * Z_check_ks(:, :, k) - tmp;
                        U_check_k = Z_check_ks(:, :, k) - 1/rho_check_k * D_k' * tmp;
                    tmp6 = lambda2 * (1 - alpha) / rho_check_k;
                    tmp7 = lambda2 * alpha / rho_check_k;
                    U_check_prime_k = max(0, 1- tmp6./abs(U_check_k)).* U_check_k;
                    Z_check_ks(:, :, k) = max(0, 1 - tmp7 ./ vecnorm(U_check_prime_k')').* U_check_prime_k;
                else
                    U_check_k = Z_check_ks(:, :, k) - 1/rho_check_k * D_k' * (D_k * Z_check_ks(:, :, k) - X_check_k);
                    tmp6 = lambda2 * (1 - alpha) / rho_check_k;
                    tmp7 = lambda2 * alpha / rho_check_k;
                    U_check_prime_k = max(0, 1- tmp6./abs(U_check_k)).* U_check_k;
                    Z_check_ks(:, :, k) = max(0, 1 - tmp7 ./ vecnorm(U_check_prime_k')').* U_check_prime_k;
                end
            end

            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get {Z^s} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Z_ks = Z_hat_ks + Z_check_ks;
            Z1 = squeeze(Z_ks(1, :, :))';
            Z2 = squeeze(Z_ks(2, :, :))';
            Z3 = squeeze(Z_ks(3, :, :))';
 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update D_hat %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            W_hat = D_hat - lambda3/rho_1 * H * (D_hat - D);
            D_hat = max(0, 1 - (lambda4/rho_1) ./ vecnorm(W_hat)).* W_hat;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Z = blkdiag(Z1, Z2, Z3);
            L_D = norm(Z * Z', 'fro') + lambda3 * norm(H, 'fro');
            rho_2 = gamma_2 * L_D;
            W = D - 1/rho_2 * ((D * Z - X_prime) * Z' + lambda3 * H * (D - D_hat));
            D = 1./max(1, vecnorm(W)).* W;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute Changes And Costs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            change_D = sum(sum(abs(D - D_old)))/(num_subjects * S * K);
            Delta_D = [Delta_D, change_D];
            
            if 0 == mod(t, 30)
                fprintf("Iteration: %d finished, time took: %f, change of D: %f, break threshold: %f, \n", t, toc/3600, change_D, threshold_break);
            end
            %%%%%%%%%%%%% Optimization, All empty, Early Stop %%%%%%%%%%%%%
            
            if t >= 500 && mod(t, 100) == 0 
                [a, b] = deal(0, 0);
                for k = 1: K
                    for s = 1: S
                         if 0 == sum(abs(Z_check_ks(s, :, k)))
                             b = b + 1; 
                         end
                         if 0 == sum(abs(D_hat(:, k + (s - 1) * S)))
                             a = a + 1; 
                         end
                    end
                end
                if b == 0 || b == K * S || a == 0 || a == K * S
                    fprintf("Break because all D_hat a == %f, Z_check b == %f\n", a, b);
                    [cost_recon, cost_Z, cost_D_tilde, cost_D_hat, cost] = cost_cal(X1, X2, X3, D, D_hat, Z1, Z2, Z3, K, Z_hat_ks, Z_check_ks, lambda1, lambda2, lambda3, lambda4, alpha, subjects_range);        
                    break;
                end
                if 0 == sum(sum(sum(abs(Z_hat_ks))))
                    fprintf("Break because sum(sum(abs(Z_hat_ks))) == %f\n", sum(sum(abs(Z_hat_ks))));
                    [cost_recon, cost_Z, cost_D_tilde, cost_D_hat, cost] = cost_cal(X1, X2, X3, D, D_hat, Z1, Z2, Z3, K, Z_hat_ks, Z_check_ks, lambda1, lambda2, lambda3, lambda4, alpha, subjects_range);        
                    break;
                end
            end
            
            if change_D < threshold_break || t == MAXITER
                fprintf("Training Took %f hours at iteration %d. change_D = %f \n", toc/3600, t, change_D);
                break; 
            end
            % save temp results
            %{
            if 0 == mod(t, 100)
                [cost_recon, cost_Z, cost_D_tilde, cost_D_hat, cost] = cost_cal(X1, X2, X3, D, D_hat, Z1, Z2, Z3, K, Z_hat_ks, Z_check_ks, lambda1, lambda2, lambda3, lambda4, alpha, subjects_range);
                mkdir("data", "results" + "/" + string(seed) + "/" + string(t));
                fn = fullfile("data", "results", string(seed), string(t), "D.mat");
                save(fn, "D", "num_subjects", "K", "subjects_range"); 
                fn = fullfile("data", "results", string(seed), string(t), "D_hat.mat");
                save(fn, "D_hat"); 
                fn = fullfile("data", "results", string(seed), string(t), "Z_hat_ks.mat");
                save(fn, "Z_hat_ks", "num_voxels");
                fn = fullfile("data", "results", string(seed), string(t), "Z_check_ks.mat");
                save(fn, "Z_check_ks"); 
                fn = fullfile("data", "results", string(seed), string(t), "Delta_D.mat");
                save(fn, "Delta_D"); 
                fn = fullfile("data", "results", string(seed), string(t), "costs.mat");
                save(fn, "cost", "cost_recon", "cost_Z", "cost_D_tilde", "cost_D_hat");
                fn = fullfile("data", "results", string(seed), string(t), "lambda4_" + string(lambda4) + "_.mat");
                save(fn, "lambda4"); 
            end
            %}
            t = t + 1;
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [cost_recon, cost_Z, cost_D_tilde, cost_D_hat, cost] = cost_cal(X1, X2, X3, D, D_hat, Z1, Z2, Z3, K, Z_hat_ks, Z_check_ks, lambda1, lambda2, lambda3, lambda4, alpha, subjects_range);
        mkdir("data", "results" + "/" + string(seed));
        fn = fullfile("data", "results", string(seed), "D.mat");
        save(fn, "D", "num_subjects", "K", "subjects_range"); 
        fn = fullfile("data", "results", string(seed), "D_hat.mat");
        save(fn, "D_hat"); 
        fn = fullfile("data", "results", string(seed), "Z_hat_ks.mat");
        save(fn, "Z_hat_ks", "num_voxels");
        fn = fullfile("data", "results", string(seed), "Z_check_ks.mat");
        save(fn, "Z_check_ks"); 
        fn = fullfile("data", "results", string(seed), "Delta_D.mat");
        save(fn, "Delta_D"); 
        fn = fullfile("data", "results", string(seed), "costs.mat");
        save(fn, "cost", "cost_recon", "cost_Z", "cost_D_tilde", "cost_D_hat");
        fn = fullfile("data", "results", string(seed), "lambda4_" + string(lambda4) + "_.mat");
        save(fn, "lambda4");
        %fn = fullfile("data", "results", string(seed), "idx_train.mat");
        %save(fn, "idx_train");
        disp("Finished");
    end
end 


function [cost_recon, cost_Z, cost_D_tilde, cost_D_hat, cost_t] = cost_cal(X1, X2, X3, D, D_hat, Z1, Z2, Z3, K, Z_hat_ks, Z_check_ks, lambda1, lambda2, lambda3, lambda4, alpha, subjects_range)
    D1 = D(:, 1: K);
    D2 = D(:, 1 + K: 2 * K);
    D3 = D(:, 1 + 2 * K: 3 * K); 
    cost_recon = 1/2 * norm(X1 - D1 * Z1, 'fro')^2 +  1/2 * norm(X2 - D2 * Z2, 'fro')^2 +  1/2 * norm(X3 - D3 * Z3, 'fro')^2;
    cost_Z = 0;
    for k = 1: K
        cost_Z = cost_Z + lambda1 * sum(vecnorm(Z_hat_ks(:, :, k)));
        cost_Z = cost_Z + lambda2 * alpha * sum(vecnorm(Z_check_ks(:, :, k)'));
        cost_Z = cost_Z + lambda2 * (1 - alpha) * sum(sum(abs(Z_check_ks(:, :, k)))); 
    end
    cost_D_tilde = 1/2 * lambda3 * FDDL_discriminative((D - D_hat)', subjects_range);
    cost_D_hat = lambda4 * sum(vecnorm(D_hat)); 
    cost_t = cost_recon + cost_Z + cost_D_tilde + cost_D_hat;
end