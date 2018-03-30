function [negL, neg_grad_L] = LogisticLogLikelihood(beta, features, labels, C, sample_idx)

    features = features(sample_idx, :);
    labels = labels(sample_idx);
    
    logistic_arg = features * beta;
    logistic_arg(isinf(exp(-logistic_arg))) = -log(eps(realmax));
    
%     lc = max(la,lb)+log(1+exp(-|la-lb|))
    L = sum(logistic_arg .* (labels - 1) - log(1+exp(-logistic_arg))) - ...
        C*sum(abs(beta));   % L1 regularization
%         C*sum(beta.^2);   % L2 regularization
    negL = -L;
    
    if (isnan(L) || isinf(L))
        min(logistic_arg)
        max(logistic_arg)
        error('NAN or INF likelihood');
    end
    
    if (nargout > 1)
        logistic_prob = 1 ./ (1 + exp(-logistic_arg));
        D = size(features, 2);
        
        grad_L = sum(repmat(labels - logistic_prob, [1, D]) .* features, 1)';
        grad_L = grad_L - C*sign(beta);  % L1 regularization
%         grad_L = grad_L - 2*C*beta;  % L2 regularization
        neg_grad_L = -grad_L;
    end