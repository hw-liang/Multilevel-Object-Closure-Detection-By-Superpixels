function [E, P, frag_edge, frag_per] = SuperpixelAffinityEdgeContrastGivenFragments...
    (fragments, segment_fragments, fragment_segments, frag_angle, frag_sq_curv, image_data, edge_thresh, affinity)

    load edge_logistic.mat
    
    if (~exist('edge_thresh', 'var') || isempty(edge_thresh))
        edge_thresh = 0.05;
    end
    
    num_fragments = numel(fragments);
    num_sups = numel(segment_fragments);
    
    [M,N] = size(image_data.pb);
    
    E = zeros(num_sups, num_sups);
    P = zeros(num_sups, num_sups);
    frag_edge = zeros(num_fragments, 1);
    frag_per = zeros(num_fragments, 1);
    
    % Fill the edge and perimeter info for edges between segments
    for i = 1:num_fragments
        x = fragments{i}(:,1);
        y = fragments{i}(:,2);
        
        x = max(min(ceil(x), N), 1);
        y = max(min(ceil(y), M), 1);
        
        pix_index = sub2ind([M,N], y, x);
        num_pix = numel(pix_index);
        
        if (numel(fragment_segments{i}) == 2)
            [edge_dist_features, edge_strength_feature, angle_feature, curv_feature] = ... 
                ComputePixelEdgeFeatures(image_data, pix_index, frag_angle{i}, frag_sq_curv{i});
            feature = [edge_dist_features, edge_strength_feature, angle_feature, curv_feature];
%             feature = [edge_dist_features, angle_feature, curv_feature];
%             feature = [edge_dist_features, angle_feature, ];
%             feature = [edge_dist_features];

            edge = LogisticPredict([feature, ones(num_pix, 1)], edge_beta);
            
            if (exist('affinity', 'var') && ~isempty(affinity))
                a = affinity(fragment_segments{i}(1), fragment_segments{i}(2));
                weight = 0.6;
                edge = (1-weight)*edge+weight*a*edge;
            end
             
%             frag_edge(i) = sum(edge);
            frag_edge(i) = sum(min(double(edge > edge_thresh), 1));
        end
        
        frag_per(i) = num_pix;
        
        if (numel(fragment_segments{i}) == 2)
            sup1 = fragment_segments{i}(1);
            sup2 = fragment_segments{i}(2);
            
            E(sup1, sup2) = -frag_edge(i);
            E(sup2, sup1) = -frag_edge(i);
            P(sup1, sup2) = -frag_per(i);
            P(sup2, sup1) = -frag_per(i);
            
%             E(sup1, sup2) = -frag_edge(i)/(frag_per(i)+eps);
%             E(sup2, sup1) = -frag_edge(i)/(frag_per(i)+eps);
%             P(sup1, sup2) = -1;
%             P(sup2, sup1) = -1;
        end
    end
    
    % Fill the edge and perimeter info for superpixels
    for i = 1:num_sups
        E(i,i) = sum(frag_edge(segment_fragments{i}));
        P(i,i) = sum(frag_per(segment_fragments{i}));
%         E(i,i) = sum(frag_edge(segment_fragments{i})./(frag_per(segment_fragments{i})+eps));
%         P(i,i) = numel(segment_fragments{i});        
    end
    
    
    
