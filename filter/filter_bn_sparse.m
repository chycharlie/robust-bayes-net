% Experiments on robust learning of fixed-structure Bayes nets with random tree/graph structure.
% Compare the performance of MLE_noNoise, MLE, RANSAC, and Filtering (performance = error in total variation distance) for various d and m.
% Sparse matrices are used to allow better scalability (Fxq is N by m but each row is d-sparse).

% N = number of samples, eps = fraction of corruption, d = number of variables, m = number of conditional probabilities.

clear
rng(1, 'twister');  % Fix random seed for replication/debug.

% Whether the ground-truth Bayes net is a tree or a random graph.
bn_is_tree = 0;
fprintf('BN_is_tree = %d\n', bn_is_tree);

for target_m = 100:100:1000
	% For each value of m, repeat the experiments several times.
	for itr = 1:5
        eps = 0.1;
        if (bn_is_tree)
            d = target_m / 2;
        else
            d = 50;
        end
        % Construct the ground-truth Bayes Net (BN). 
        % Generate the in-degree of every node.
        parent = cell(d, 1);
        if (bn_is_tree)
            % Tree case: the degree is 0 for the root and 1 for other nodes.
            deg = [0; ones(d-1, 1)];
            m = d+d-1;
        else
            % Graph case: start with the empty graph; increase the degree of a random node until m reaches the targer m.
            deg = zeros(d, 1);
            m = d;
            while (m < target_m)
                i = randi(d);
                if (deg(i) < i-1)
                    m = m - 2^deg(i) + 2^(deg(i)+1);
                    deg(i) = deg(i) + 1;
                end
            end
        end
        % Generate the graph structure (edges go from nodes with smaller index to ones with larger index).
        % Generate the conditional probabilities p(i, a), drawn i.i.d from [0, 1/4] \cup [3/4, 1].
        p = zeros(m, 1);
        m = 0;
        for i = 1:d
            parent{i} = randsample(i-1, deg(i))';
            for j = 1:2^deg(i)
                m = m + 1;
                if (rand() > 0.5)
                    p(m) = rand(1)/4 + 3/4;
                else
                    p(m) = rand(1)/4;
                end
            end
        end
        fprintf('d = %d, m = %d, itr = %d\n', d, m, itr);

        % Take N samples from the ground-truth BN.
        % This implementation is much faster than looping over N.
        N = 10*floor(m / eps^2);
        X = zeros(round((1-eps)*N), d);
        k = 0;
        for i = 1:d
            for j = 1:2^deg(i)
                % For every node i and every parental configuration j, select the samples in which parent[i] = j, draw X_i randomly.
                parent_config = dec2bin(j-1, deg(i)) - '0';
                k = k + 1;
                matched_rows = all(bsxfun(@eq, X(:, parent{i}), parent_config), 2);
                X(matched_rows, i) = rand(sum(matched_rows), 1) < p(k);
            end
        end

        % Evaluate MLE without noise (gold standard).
        p_MLE_noNoise = empirical_bn(parent, X);
        fprintf('\td_TV(p, p_MLE_noNoise) = %f\n', dtv_bn(parent, p, p_MLE_noNoise));

        % eps-fraction of the samples are corrupted.
        if (bn_is_tree)
            % Draw eps*N corrupted samples from a product distribution when the ground-truth BN is a tree.
            % The mean of the product distribution is i.i.d in [0, 1].
            % eps*N samples can be arbitrarily corrupted. Feel free to use a different Y.
            p_noise = rand(1, d);
            Y = rand(round(eps*N), d);
            Y = bsxfun(@le, Y, p_noise);
        else
            % Draw corrupted samples from a tree BN when the ground-truth BN is a graph.
            % The conditional probabilities are drawn i.i.d from [0, 1/4] \cup [3/4, 1].
            Y = zeros(round(eps*N), d);
            k = 0;
            for i = 2:d
                parent_noise_i = randsample(i-1, min(i-1, 1))';
                for j = 1:2
                    k = k + 1;
                    if (rand() > 0.5)
                        p_noise_k = rand(1)/4 + 3/4;
                    else
                        p_noise_k = rand(1)/4;
                    end
                    % Pi(i, j) == 1 if (y(parent_noise_i) == j-1).
                    matched_rows = bsxfun(@eq, Y(:, parent_noise_i), j-1);
                    Y(matched_rows, i) = rand(sum(matched_rows), 1) < p_noise_k;
                end
            end
        end
        X = [X; Y];
        N = size(X, 1);

        % Evaluate MLE (i.e., empirical conditional mean) as baseline #1.
        p_MLE = empirical_bn(parent, X);
        fprintf('\td_TV(p, p_MLE) = %f\n', dtv_bn(parent, p, p_MLE));

        % Evaluate RANSAC as baseline #2.
        ransac_N = round(eps*N);
        nItr_ransac = 5;  % Full experiments set this to 100.
        min_dtv_ransac = 1;
        for ransac_itr = 1:nItr_ransac
            ransac_X = X(randsample(1:N, ransac_N), :);
            p_ransac = empirical_bn(parent, ransac_X);
            min_dtv_ransac = min(min_dtv_ransac, dtv_bn(parent, p, p_ransac));
        end
        fprintf('\td_TV(p, p_RANSAC) = %f\n', min_dtv_ransac);

        % Evaluate Filtering (our method).
        
        % Expand the sample dimension from d to m.  Fill in the unknown information with empirical conditional means.
        p_empirical_cond_prob = empirical_bn(parent, X);  % q = p_empirical_cond_prob.
        % Compute and store (F(X,q)-q) as a sparse matrix.
        % This may be slower, but an (N by m) matrix may not fit in memory, while its sparse version will fit.
        % Too slow to access a sparse matrix repeatly by index, create it all at once.
        % Fxq_minus_q = sparse(N, m);
        Fxq_q_i = zeros(N*d, 1);
        Fxq_q_j = zeros(N*d, 1);
        Fxq_q_v = zeros(N*d, 1);
        Fxq_q_nnz = 0;
        k = 0;
        for i = 1:d
            for j = 1:2^deg(i)
                parent_config = dec2bin(j-1, deg(i)) - '0';
                k = k + 1;
                if deg(i) == 0
                    Fxq_q_i(Fxq_q_nnz+1:Fxq_q_nnz+N) = (1:N)';
                    Fxq_q_j(Fxq_q_nnz+1:Fxq_q_nnz+N) = k;
                    Fxq_q_v(Fxq_q_nnz+1:Fxq_q_nnz+N) = X(:, i) - p_empirical_cond_prob(k);
                    Fxq_q_nnz = Fxq_q_nnz + N;
                else
                    matched_config = all(bsxfun(@eq, X(:, parent{i}), parent_config), 2);
                    sum_matched = sum(matched_config);
                    Fxq_q_i(Fxq_q_nnz+1:Fxq_q_nnz+sum_matched) = find(matched_config);
                    Fxq_q_j(Fxq_q_nnz+1:Fxq_q_nnz+sum_matched) = k;
                    Fxq_q_v(Fxq_q_nnz+1:Fxq_q_nnz+sum_matched) = X(matched_config, i) - p_empirical_cond_prob(k);
                    Fxq_q_nnz = Fxq_q_nnz + sum_matched;
                end
            end
        end
        Fxq_minus_q = sparse(Fxq_q_i, Fxq_q_j, Fxq_q_v, N, m);
        
        % The full version of the filtering algorithm should repeat until the top eigenvalue is small enough.
        % In our experiment, we only run one iteration, which is faster and gives reasonable results.
        
        % Compute the top eigenvalue and eigenvector of off-diag(cov(Fxq - q)).
        cov_Fxq_q = full(Fxq_minus_q' * Fxq_minus_q) / N;
        cov_Fxq_q = cov_Fxq_q - diag(diag(cov_Fxq_q));
        [v1, lambda1] = eigs(cov_Fxq_q, 1);
        % Project the (expanded) samples along the direction v1.
        projection_data_pair = [abs(Fxq_minus_q * v1) X];
        % Sort by the absolute value of the projection (first column).
        sorted_pair = sortrows(projection_data_pair);
        % Remove eps-fraction of the sample farthest from the projected mean.
        i = round((1-eps)*N);
        X = sorted_pair(1:i, 2:end);

        p_filter = empirical_bn(parent, X);
        fprintf('\td_TV(p, p_filter) = %f\n', dtv_bn(parent, p, p_filter))

        fprintf('\t\t%f fraction of the filtered samples are noise\n', 1 - size(intersect(X, Y, 'rows'), 1) / (eps*N));
	end
end