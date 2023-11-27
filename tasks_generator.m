fid = fopen("data/tasks_parameters.txt", "w"); 
idx = 1;
K = 30;
lambda1s = [0.008, 0.009, 0.01];
lambda2s = [0.34, 0.35, 0.36, 0.37, 0.4];
lambda3s = [0.35, 0.38, 0.4, 0.43, 0.45];
lambda4_coefs = [0.035, 0.04, 0.045, 0.05];
alphas = [0.99];
for i = 1:length(lambda1s)
    for j = 1:length(lambda2s)
        for k = 1:length(lambda3s)
            for l = 1: length(lambda4_coefs)
                for m = 1: length(alphas)
                    lambda1 = lambda1s(i);
                    lambda2 = lambda2s(j);
                    lambda3 = lambda3s(k);
                    lambda4_coef = lambda4_coefs(l);
                    alpha = alphas(m);
                    fprintf(fid, "%d %f %f %f %f %f %d\n", idx, lambda1, lambda2, lambda3, lambda4_coef, alpha, K);
                    idx = idx + 1;
                end
            end
        end
    end
end
fclose(fid);
disp("stop");