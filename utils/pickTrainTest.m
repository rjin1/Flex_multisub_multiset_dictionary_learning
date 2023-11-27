function [Y_train,  Y_val, label_train, label_val, idx_train] = pickTrainTest(dataset, N_trains, N_vals, seed)
    rng(seed);
    data_fn = fullfile('data', strcat(dataset, '.mat'));
    load(data_fn);
    Y = Y';
    rows = size(Y, 1); % number of rows
    if ~exist('Y_range', 'var')
        Y_range = label_to_range(label); 
    end

    C = numel(Y_range) - 1; % numel: number of elements
    N_train = sum(N_trains);
    N_val = sum(N_vals);
    Y_train = zeros(rows, N_train);
    Y_val = zeros(rows, N_val);

    label_train = zeros(1, N_train);
    label_val = zeros(1, N_val);

    cur_train = 0; 
    cur_val = 0;
    idx_train = [];
    for c = 1: C 
        Yc = get_block_col(Y, c, Y_range);
        N_total_c = size(Yc, 2);
        
        label_train(:, cur_train + 1: cur_train + N_trains(c)) = c*ones(1, N_trains(c));
        label_val(:, cur_val + 1: cur_val + N_vals(c)) = c*ones(1, N_vals(c));
        idx = randperm(N_total_c);

        idx_train = [idx_train, Y_range(c) + idx(1: N_trains(c))];
        Y_train(:, cur_train + 1: cur_train + N_trains(c)) = Yc(:, idx(1: N_trains(c)));
        Y_val(:, cur_val + 1: cur_val + N_vals(c)) = Yc(:, idx(N_trains(c) + 1: N_trains(c) + N_vals(c)));
        
        cur_train = cur_train + N_trains(c);
        cur_val = cur_val + N_vals(c);
    end 
    %Y_train = normc(Y_train)'; 
    %Y_val = normc(Y_val)';
    Y_train = normc(Y_train)'; 
    Y_val = normc(Y_val)';
end 