% LogisticLearn - train a binary logistic classifier
%
% Arguments:
%   features - NxD matrix. Each row is a feature vector. There are N
%   training samples.
%   labels   - Nx1 vector of 0,1 ground truth labels
%   C        - the importance of sparseness of parameters. Higher value
%   give sparser parameters.
%   do_quadratic_kernel - First compute quadratic features and then learn
%   the classifier (false by default)
%
% Returns:
%   beta     - Dx1 vector of parameters of a logistic probability distribution
%   P(X=1) = 1 / (1 + exp(-beta' * feature))
%
function [beta, history] = LogisticLearn(features, labels, C, do_quadratic_kernel)
    if (nargin < 4)
        do_quadratic_kernel = false;
    end
    
    options = optimset('GradObj','on', 'LargeScale', 'off', 'MaxIter', 1000, ...
        'Display', 'iter', 'DerivativeCheck', 'off', 'outputfcn', @LogisticLogLikelihoodOutput);
    
    [num_samples, D] = size(features);
    
    if (do_quadratic_kernel)
        quad_ind = tril(ones(D,D)) > 0;
        quad_features = zeros(num_samples, sum(quad_ind(:)));
        
        for i = 1:num_samples
            cur_quad_features = features(i,:)' * features(i, :);
            quad_features(i, :) = cur_quad_features(quad_ind(:));
        end
        
        features = quad_features;
        [num_samples, D] = size(features);
    end
    
    init_beta = zeros(D, 1);
    percent_validation = 0.2;
    validation_idx = rand(num_samples, 1) < percent_validation;
%     validation_idx(1:floor(percent_validation*num_samples)) = true;
    train_idx = ~validation_idx;
    history.x = [];
    history.fval_train = [];
    history.fval_valid = [];
%     beta = fminunc(@LogisticLogLikelihood, init_beta, options, features, labels, C, train_idx);
    beta = minimize(init_beta, @LogisticLogLikelihood, 1500, features, labels, C, train_idx);

    function stop = LogisticLogLikelihoodOutput(x, optimValues, state, varargin)
        stop = false;
        history.x = [history.x, x];
        history.fval_train = [history.fval_train; optimValues.fval];
        fval_valid = LogisticLogLikelihood(x, features, labels, C, validation_idx);
        history.fval_valid = [history.fval_valid; fval_valid];
    
    end
end