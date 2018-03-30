function probabilities = LogisticPredict(features, beta, do_quadratic_kernel)

    if (nargin < 3)
        do_quadratic_kernel = false;
    end
    
    if (do_quadratic_kernel)
        [num_samples, D] = size(features);
    
        quad_ind = tril(ones(D,D)) > 0;
        quad_features = zeros(num_samples, sum(quad_ind(:)));
        
        for i = 1:num_samples
            cur_quad_features = features(i,:)' * features(i, :);
            quad_features(i, :) = cur_quad_features(quad_ind(:));
        end
        
        features = quad_features;
    end

    probabilities = 1 ./ (1 + exp(-features * beta));