function image_data = PrecomputeImageData(img, pb, pb_theta, tmap)
    lower_score = 0.02;
    higher_score = 0.05;
    thinned_score = hysthresh(pb, higher_score, lower_score);

    % Link the resulting branches into edges with a given minimum length
    min_edge_length = 10;
    edgelist = edgelink(thinned_score, min_edge_length);
    
    new_pb = zeros(size(pb));
    for e = 1:numel(edgelist)
        edge_idx = sub2ind(size(pb), edgelist{e}(:,1), edgelist{e}(:,2));
        mean_pb = mean(pb(edge_idx));
        new_pb(edge_idx) = mean_pb;
    end
    pb = new_pb;
    
    [pb_dist, pb_closest] = bwdist(pb > 0);
    
    dilated_pb = imdilate(pb, ones(5,5));
    
%     no = 6;
%     ss = 1;
%     ns = 2;
%     sc = sqrt(2);
%     el = 2;
%     k = 32;
%     fname = sprintf( ...
%         'unitex_%.2g_%.2g_%.2g_%.2g_%.2g_%d.mat',no,ss,ns,sc,el,k);
%     load(fname,'tex','fb','tsim'); % defines fb,tex,tsim
%     tmap = assignTextons(fbRun(fb,rgb2gray(img)),tex);
    
    image_data.pb = pb;
    image_data.dilated_pb = dilated_pb;
    image_data.img_r = img(:,:,1);
    image_data.img_g = img(:,:,2);
    image_data.img_b = img(:,:,3);
    hsv_img = rgb2hsv(img);
    image_data.img_h = hsv_img(:,:,1);
    image_data.img_s = hsv_img(:,:,2);
    image_data.img_v = hsv_img(:,:,3);
    image_data.cos_hs = cos(image_data.img_h * 2 * pi) .* image_data.img_s;
    image_data.sin_hs = sin(image_data.img_h * 2 * pi) .* image_data.img_s;

    image_data.pb_theta = pb_theta;
    image_data.tmap = tmap;
    image_data.pb_dist = pb_dist;
    image_data.pb_closest = pb_closest;
    