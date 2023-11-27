%%%%%%%%% Use a set of good estimated maps to check linear regression
seed = 1;
rng(seed);
K = 30;
V = 48546;
M = 200;
S = 3;

for i = 1:K
    ind(1+(i-1)*S : i*S) = i:K:(S-1)*K+i;
end

D = normrnd(0, 1, [M, S * K]);

Z_hat_ks = zeros(S, V, K);
Z_check_ks = zeros(S, V, K);
for k = 1: K
    Z_hat_ks(:, :, k) = normrnd(0, 1, [S, V]); 
    Z_check_ks(:, :, k) = normrnd(0, 1, [S, V]); 
end
Z_ks = Z_hat_ks + Z_check_ks;
Z1 = squeeze(Z_ks(1, :, :))';
Z2 = squeeze(Z_ks(2, :, :))';
Z3 = squeeze(Z_ks(3, :, :))';

D_tmp = blkdiag(D(:, 1: K),D(:, 1 + K: 2 * K),D(:, 1 + 2 * K: 3 * K));
Z = [Z1; Z2; Z3];
a = D_tmp(:,ind) * Z(ind,:);

b = 0;
for k = 1: K
    D_k = blkdiag(D(:, k), D(:, K + k), D(:, 2 * K + k));
    b = b + D_k * (Z_check_ks(:, :, k) + Z_hat_ks(:, :, k));
end

norm(a-b)
sum(sum(abs(a-b)))
