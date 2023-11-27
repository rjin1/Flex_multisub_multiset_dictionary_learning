function [] = tunning_worker_cv(E_id)
    addpath('utils');
    % datasets = ["Synthesized_dataset1", "Synthesized_dataset2", "Synthesized_dataset3"];
    datasets = ["New_AOD3", "New_SIRP", "New_SM"];
    [num_folds, num_inits] = deal(3, 1);
    N_train = [97, 76]; % number of subjects to train per class
    N_vals = [41, 33]; % number of subjects to validate of [HC, SZ].
    while true
        % find how many combinations of parameters are there
        fid = fopen("data/tasks_parameters.txt", "r"); 
        format = "%d %f %f %f %f %f %d";
        size_pars = [7 Inf];
        pars = fscanf(fid, format, size_pars); 
        num_par = size(pars, 2);
          
        %%%%%%%%%%%%%%%%%%%% read next task ID %%%%%%%%%%%%%%%%%%
        fid = fopen("data/tasks.txt", "r+");
        format = "%d";
        parIDs = fscanf(fid, format);
        parID = parIDs(end); 
        if parID > num_par
            break;
        end 
        fprintf(fid, "\n%d", parID + 1);
        fclose(fid);
        
        [task_id, lambda1, lambda2, lambda3, lambda4_coef, alpha, K] = ...
        deal(pars(1, parID), pars(2, parID), pars(3, parID), pars(4, parID), pars(5, parID), pars(6, parID), pars(7, parID));
        fprintf('%s, pars id = %d, task_id = %d, lambda1 = %f, lambda2 = %f, lambda3 = %f, lambda4_coef = %d, alpha = %d\n', E_id, parID, task_id, lambda1, lambda2, lambda3, lambda4_coef, alpha);
        tic;
        spmd(num_folds)
             [X1_train, X1_val, label_train, label_val, ~] = pickTrainTest(datasets(1), N_train, N_vals, labindex);
             [X2_train, X2_val, ~, ~, ~] = pickTrainTest(datasets(2), N_train, N_vals, labindex);
             [X3_train, X3_val, ~, ~, ~] = pickTrainTest(datasets(3), N_train, N_vals, labindex);
             subjects_range_train = label_to_range(label_train);
             [acc_NCs, acc_NNs, cost_recons] = deal(0, 0, 0);
             for id_init = 1: num_inits
                    seed = id_init + labindex;
                    [D_train, D_hat_train, Z1_train, Z2_train, Z3_train, ~, ~] = palm_solver(lambda1, lambda2, lambda3, lambda4_coef, alpha, K, X1_train, X2_train, X3_train, subjects_range_train, seed);
                    [D1_train, D2_train, D3_train] = deal(D_train(:, 1: K), D_train(:, K + 1: 2 * K), D_train(:, 2 * K + 1: 3 * K));
                    [D1_val, D2_val, D3_val] = deal(X1_val/Z1_train, X2_val/Z2_train, X3_val/Z3_train);
                    
                    [D_train_discriminative, D_val_discriminative] = deal([], []);
                    for k = 1: K
                        if 0 == norm(D_hat_train(:, k))
                             D_train_discriminative = [D_train_discriminative, D1_train(:, k)];
                             D_val_discriminative = [D_val_discriminative, D1_val(:, k)];
                        end
                        if 0 == norm(D_hat_train(:, k + K))
                             D_train_discriminative = [D_train_discriminative, D2_train(:, k)];
                             D_val_discriminative = [D_val_discriminative, D2_val(:, k)];
                        end
                        if 0 == norm(D_hat_train(:, k + 2*K))
                             D_train_discriminative = [D_train_discriminative, D3_train(:, k)];
                             D_val_discriminative = [D_val_discriminative, D3_val(:, k)];
                        end
                    end
                    [acc_NC, acc_NN] = pred(K, subjects_range_train, sum(N_vals),label_train, label_val, D_train_discriminative, D_val_discriminative);

                    cost_recons = cost_recons + 1/2 * norm(X1_val - D1_val * Z1_train, 'fro')^2 + 1/2 * norm(X2_val - D2_val * Z2_train, 'fro')^2 + 1/2 * norm(X3_val - D3_val * Z3_train, 'fro')^2;
                    %mySaveResults1(E_id, D_train, D_hat_train, Z_hat_ks_train, Z_check_ks_train, D1_val, D2_val, D3_val, labindex, id_init, task_id, "metric1");

                    [acc_NCs, acc_NNs] = deal(acc_NCs + acc_NC, acc_NNs + acc_NN);
             end
        end
        [acc_NC, acc_NN, cost_recon] = deal(0, 0, 0);
        for i = 1: num_folds
             [acc_NC, acc_NN, cost_recon] = deal(acc_NC + acc_NCs{i}, acc_NN + acc_NNs{i}, cost_recon + cost_recons{i});
        end
        [acc_NC, acc_NN, cost_recon] = deal(acc_NC/ num_folds/num_inits, acc_NN/ num_folds/num_inits, cost_recon / num_folds/num_inits);
        mkdir("data", E_id);
        fn = fullfile("data", E_id, string(lambda1) + "_" + string(lambda2) + "_" + string(lambda3) + "_" + string(lambda4_coef) + "_" + string(alpha) + "_" +  string(acc_NC) + "_" + string(acc_NN)+ "_" +  string(cost_recon) + "_" + string(task_id) + "_.mat");
        save(fn, "acc_NC"); 
        fprintf('%s, pars id = %d, task_id = %d Done Took %f hours!\n', E_id, parID, task_id, toc/3600);
    end
end


function [] = mySaveResults1(E_id, D, D_hat, Z_hat_ks, Z_check_ks, D1_val, D2_val, D3_val, labindex, id_init, task_id, metric)
    mkdir("data\" + E_id + "\" + string(task_id) + "\" +  metric + "\" + string(labindex)+ "\" + string(id_init));
    fn = fullfile("data", E_id, string(task_id), metric, string(labindex), string(id_init), "D_train.mat");
    save(fn, "D");
    fn = fullfile("data", E_id, string(task_id), metric, string(labindex), string(id_init), "D_hat_train.mat");
    save(fn, "D_hat");
    fn = fullfile("data", E_id, string(task_id), metric, string(labindex), string(id_init), "Z_hat_ks_train.mat");
    save(fn, "Z_hat_ks");
    fn = fullfile("data", E_id, string(task_id), metric, string(labindex), string(id_init), "Z_check_ks_train.mat");
    save(fn, "Z_check_ks");
    fn = fullfile("data", E_id, string(task_id), metric, string(labindex), string(id_init), "D1_val.mat");
    save(fn, "D1_val");
    fn = fullfile("data", E_id, string(task_id), metric, string(labindex), string(id_init), "D2_val.mat");
    save(fn, "D2_val");
    fn = fullfile("data", E_id, string(task_id), metric, string(labindex), string(id_init), "D3_val.mat");
    save(fn, "D3_val");
end
                   