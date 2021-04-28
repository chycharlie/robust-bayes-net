% Estimate the total variation distance (via sampling) between two Bayes nets P and Q that share the same structure.
% Input: the common graph structure (parent) and two conditional probability tables (p and q).
% Output: the estimated total variation distance.

function [dtv] = dtv_bn(parent, p, q)
    % d_TV(p, q) = p(A) - q(A), where A = {x: p(x) > q(x)}.
	n = size(parent, 1);
    deg = zeros(n, 1);
    for i = 1:n
        deg(i) = numel(parent{i});
    end
    
    % Draw N samples from p to estimate p(A).
    % Use log-likelihood to decide if x is in A.
    N = 10^5;
    X = zeros(N, n);
    log_pX = zeros(N, 1);
    log_qX = zeros(N, 1);
    k = 0;
    for i = 1:n
        for j = 1:2^deg(i)
            parent_config = dec2bin(j-1, deg(i)) - '0';
            k = k + 1;
            matched_rows = all(bsxfun(@eq, X(:, parent{i}), parent_config), 2);
            X(matched_rows, i) = rand(sum(matched_rows), 1) < p(k);
            matched_rows_0 = (matched_rows & X(:, i)==0);
            matched_rows_1 = (matched_rows & X(:, i)==1);
            log_pX(matched_rows_1) = log_pX(matched_rows_1) + log(p(k));
            log_pX(matched_rows_0) = log_pX(matched_rows_0) + log(1-p(k));
            log_qX(matched_rows_1) = log_qX(matched_rows_1) + log(q(k));
            log_qX(matched_rows_0) = log_qX(matched_rows_0) + log(1-q(k));
        end
    end
    pA = sum(log_pX > log_qX) / N;
    
    % Estimate q(A) similarly.
    X = zeros(N, n);
    log_pX = zeros(N, 1);
    log_qX = zeros(N, 1);
    k = 0;
    for i = 1:n
        for j = 1:2^deg(i)
            parent_config = dec2bin(j-1, deg(i)) - '0';
            k = k + 1;
            matched_rows = all(bsxfun(@eq, X(:, parent{i}), parent_config), 2);
            X(matched_rows, i) = (rand(sum(matched_rows), 1) < q(k));
            matched_rows_0 = (matched_rows & X(:, i)==0);
            matched_rows_1 = (matched_rows & X(:, i)==1);
            log_pX(matched_rows_1) = log_pX(matched_rows_1) + log(p(k));
            log_pX(matched_rows_0) = log_pX(matched_rows_0) + log(1-p(k));
            log_qX(matched_rows_1) = log_qX(matched_rows_1) + log(q(k));
            log_qX(matched_rows_0) = log_qX(matched_rows_0) + log(1-q(k));
        end
    end
    qA = sum(log_pX > log_qX) / N;
    
    dtv = max(pA - qA, 0);
end