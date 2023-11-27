 function [acc_NC, acc_NN] = pred(K, subjects_range_train, num_subjects_val,label_train, label_val, D_train, D_val)

    if isempty(D_train) || isempty(D_val) || size(D_train, 2) == 3*K
        [acc_NC, acc_NN] = deal(0, 0);
    else
        %% --------- classification by nearest centroid -------
        acc_NC = pred_NC(subjects_range_train, num_subjects_val, D_val, D_train, label_val);

        %% ----- alternative classificaiton by K-nearest neighbor ----------------
        acc_NN = pred_NN(num_subjects_val, D_val, D_train, label_train, label_val);
    end 
end 

function [acc_NC] = pred_NC(subjects_range_train, num_subjects_val, D_val, D_train, label_val)
    nClasses = 2;
    D_train_M = zeros(size(D_train', 1), nClasses);
    for c = 1:nClasses
        Xc = get_block_col(D_train', c, subjects_range_train);
        D_train_M(:, c) = mean(Xc, 2); % mean of each row(profiles)
    end
    
    E = zeros(nClasses, num_subjects_val);
    for c = 1: nClasses
        Mc = repmat(D_train_M(:, c), 1, num_subjects_val); 
        R1 = D_val' - Mc;
        E(c, :) = sum(R1.^2);
    end
    [~, pred] = min(E);
    acc_NC = double(sum(pred == label_val)) / num_subjects_val;
end

function [acc_NN] = pred_NN(num_subjects_val, D_val, D_train, label_train, label_val)
    Md = ExhaustiveSearcher(D_train);
    [idx,~] = knnsearch(Md, D_val,'K',1);
    pred = label_train(idx);
    acc_NN = double(sum(pred == label_val))/ num_subjects_val;
end