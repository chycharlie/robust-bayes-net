% Compute the empirical conditional probabilities of a known-structure Bayes net.
% Input: the graph structure (parent) and the samples (X).
% Output: the empirical conditional probability table (p).
% This function is randomized.

function [p] = empirical_bn(parent, X)
    % N = number of samples, d = number of variables, m = size of conditional probability table.
    N = size(X, 1);
    d = size(parent, 1);
    deg = zeros(d, 1);
    % Compute m.
    m = 0;
    for i = 1:d
        deg(i) = numel(parent{i});
        m = m + deg(i);
    end
    
    % Compute p.
    p = zeros(m, 1);
    k = 0;
    for i = 1:d
        for j = 1:2^deg(i)
            % parent_config = 
            parent_config = dec2bin(j-1, deg(i)) - '0';
            % cnt01 = the number of times that [parent(i) = j] occurred in X.
            if (deg(i) == 0)
                cnt01 = N;
                parent_config = [];
            else
                cnt01 = sum(all(bsxfun(@eq, X(:, parent{i}), parent_config), 2));
            end
            % cnt1 = the number of times that [parent(i) = j and X_i = 1] occurred in X.
            cnt1 = sum(all(bsxfun(@eq, X(:, [parent{i} i]), [parent_config 1]), 2));
            k = k + 1;
            if (cnt01 > 0)
                p(k) = cnt1 / cnt01;
            else
                p(k) = 0.5;
            end
        end
    end
end