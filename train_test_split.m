seed = 10;
N_train = [110, 87]; % number of subjects to train per class
N_vals = [28, 22]; % number of subjects to validate of [HC, SZ].
datasets1 = ["Synthesized_dataset1", "Synthesized_dataset2", "Synthesized_dataset3"];
datasets2 = ["Synthesized_dataset_test_1", "Synthesized_dataset_test_2", "Synthesized_dataset_test_3"];
[X_train, X_val, ~, ~, ~] = pickTrainTest(datasets(1), N_train, N_vals, seed);

data_fn = fullfile('data', strcat(datasets(2), '.mat'));
load(data_fn);
X2 = normc(Y')';
data_fn = fullfile('data', strcat(datasets(3), '.mat'));
load(data_fn);
X3 = normc(Y')';