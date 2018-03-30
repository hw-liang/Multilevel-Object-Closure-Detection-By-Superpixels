function [features, labels] = LogisticSampleData(num_samples)

    features = ones(num_samples, 3);
    labels = zeros(num_samples, 1);
    axis([0,1,0,1]);
    
    for i = 1:num_samples
        i
        [x,y,but] = ginput(1);
        
        if (but == 1)
            labels(i) = 1;
        end
        features(i, 1) = x;
        features(i, 2) = y;
        pos_ind = 1:i;
        pos_ind = pos_ind(labels(pos_ind) == 1);
        neg_ind = 1:i;
        neg_ind = neg_ind(labels(neg_ind) == 0);
        
        clf;
        scatter(features(pos_ind,1), features(pos_ind,2), 20, 'r', 'filled');
        hold on;
        scatter(features(neg_ind,1), features(neg_ind,2), 20, 'b', 'filled');
        axis([0,1,0,1]);
    end
    