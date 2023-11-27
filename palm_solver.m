function [D, D_hat, Z1, Z2, Z3, Z_hat_ks, Z_check_ks] = palm_solver(lambda1, lambda2, lambda3, lambda4_coef, alpha, K, X1, X2, X3, subjects_range, seed, E_id)
    rng(seed, 'twister'); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [gamma_k, gamma_1, gamma_2] = deal(1, 1, 1); % Fix gamma_hat_k = gamma_check_k = gamma_k = gamma_1 = gamma_2 = 1;
    
    %%%%%%%%%%%%%%%%% Fetch Datasets Sizes. The Three Datasets Have The Same Sizes %%%%%%%%%%%%%%%%%
    [num_subjects, num_voxels] = size(X1);
    X = [X1', X2', X3']';
    X_prime = [X1, X2, X3];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S = 3; % number of datasets
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Iterations Begin %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic;
    [t, breakLoop, threshold_break] = deal(1, true, 1e-5);
    for i = 1:K
        ind(1+(i-1)*S : i*S) = i:K:(S-1)*K+i;
    end
    indices = reshape(ind, S, K)';
    while breakLoop
        D_old = D;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update {Z_hat_k} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        X_check_k = X;
        Z = [Z1; Z2; Z3];
        D_tmp = blkdiag(D(:, 1: K),D(:, 1 + K: 2 * K),D(:, 1 + 2 * K: 3 * K));
        D_1 = D_tmp(:, indices(1, :));
        X_hat_k = X - D_tmp(:,ind) * Z(ind,:) + D_1 * Z_hat_ks(:, :, 1);
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
        %{
        for k = 1: K
            Z1(k, :) = Z_ks(1, :, k); 
            Z2(k, :) = Z_ks(2, :, k); 
            Z3(k, :) = Z_ks(3, :, k); 
        end 
        %}
        %%% Optimization, replace for-loop %%% 
        Z1 = squeeze(Z_ks(1, :, :))';
        Z2 = squeeze(Z_ks(2, :, :))';
        Z3 = squeeze(Z_ks(3, :, :))';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update D_hat %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        W_hat = D_hat - lambda3/rho_1 * H * (D_hat - D);
        %{
        for j = 1: S * K
            D_hat(:, j) = max(0, 1 - lambda4/(rho_1 * norm(W_hat(:, j)))) * W_hat(:, j);
        end
        %}
        %%%%%%% Optimization: Replace for-loop with matrix opration. %%%%%%%
        D_hat = max(0, 1 - (lambda4/rho_1) ./ vecnorm(W_hat)).* W_hat;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Update D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Z = blkdiag(Z1, Z2, Z3);
        L_D = norm(Z * Z', 'fro') + lambda3 * norm(H, 'fro');
        rho_2 = gamma_2 * L_D;
        W = D - 1/rho_2 * ((D * Z - X_prime) * Z' + lambda3 * H * (D - D_hat));
        %{
        for j = 1: S * K
            D(:, j) = 1/max(1, norm(W(:, j))) * W(:, j); 
        end
        %}
        %%%%%%% Optimization: Replace for-loop with matrix opration. %%%%%%%
        D = 1./max(1, vecnorm(W)).* W;

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
                break;
            end
            if 0 == sum(sum(sum(abs(Z_hat_ks))))
                fprintf("Break because sum(sum(abs(Z_hat_ks))) == %f\n", sum(sum(abs(Z_hat_ks))));
                break;
            end
        end
        change_D = sum(sum(abs(D - D_old)))/(num_subjects * S * K);
        if change_D < threshold_break || 800 == t
            fprintf("Break at iteration %d because of convergence or t == 800.\n", t);
            break; 
        end
        %{  
        if mod(t, interval) == 0
            cost_t = cost(X1, X2, X3, D, D_hat, Z1, Z2, Z3, K, Z_hat_ks, Z_check_ks, lambda1, lambda2, lambda3, lambda4, alpha, num_voxels, subjects_range);
            % Don't use "abs(cost_t - costs_old) / costs_old" because
            % the converged values of 'cost_t' vary from thounsand to hundreds of
            % thousand, while the converged values of 'cost_t -
            % costs_old' are usually smaller than 5. If we use "abs(cost_t - costs_old) / costs_old",
            % sometimes it breaks quickly, while somethimes it never
            % breaks;
            if abs(cost_t - costs_old)/interval < threshold_break
                fprintf("Break at Iteration %d\n", t);
                breakLoop = false;
            end
            fprintf("Training Took %f hours at iteration %d. Change rate %f, threshold_break = %f\n", toc/3600, t, abs(cost_t - costs_old)/interval, threshold_break);
            costs_old = cost_t;
        end
        %}
        if 0 == mod(t, 30)
            disp(t);
        end 
           
        t = t + 1;
    end
    %mySave(lambda1, lambda2, lambda3, lambda4_coef, seed, E_id, Z_hat_ks, Z_check_ks, D_hat, D, num_subjects, num_voxels, K)
    fprintf("Training Took %f hours.\n", toc/3600);
end

function [cost_t] = cost(X1, X2, X3, D, D_hat, Z1, Z2, Z3, K, Z_hat_ks, Z_check_ks, lambda1, lambda2, lambda3, lambda4, alpha, num_voxels, subjects_range)
    D1 = D(:, 1: K);
    D2 = D(:, 1 + K: 2 * K);
    D3 = D(:, 1 + 2 * K: 3 * K);
    cost_recon = 1/2 * norm(X1 - D1 * Z1, 'fro')^2 +  1/2 * norm(X2 - D2 * Z2, 'fro')^2 +  1/2 * norm(X3 - D3 * Z3, 'fro')^2;
    cost_Z = 0;
    parfor k = 1: K
        for j = 1: num_voxels
            cost_Z = cost_Z + lambda1 * norm(Z_hat_ks(:, :, k));
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

function mySave(lambda1, lambda2, lambda3, lambda4_coef, seed, E_id, Z_hat_ks, Z_check_ks, D_hat, D, num_subjects, num_voxels, K)
    fileName = string(lambda1) + "_" + string(lambda2) + "_" + string(lambda3) + "_" + string(lambda4_coef) + "_" + string(seed);
    mkdir("data", E_id + "/" + fileName);
    fn = fullfile("data", E_id, fileName, "data.mat");
    save(fn, "Z_hat_ks", "Z_check_ks", "D_hat", "D", "num_subjects", "num_voxels", "K");
end