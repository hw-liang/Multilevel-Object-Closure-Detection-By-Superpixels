function [edge_dist_features, edge_strength_feature, angle_feature, curv_feature] = ... 
    ComputePixelEdgeFeatures(image_data, pix_index, pix_angle, pix_sq_curv)
                    
    edge_dist_features = image_data.pb_dist(pix_index);
    edge_strength_feature = image_data.pb(image_data.pb_closest(pix_index));
    angle_feature = abs(cos(pix_angle) .* cos(image_data.pb_theta(pix_index)) ...
        + sin(pix_angle) .* sin(image_data.pb_theta(pix_index)));
    curv_feature = pix_sq_curv;