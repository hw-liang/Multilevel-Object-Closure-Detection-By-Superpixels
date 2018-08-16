% This function will find a single object with a certain mask for the
% range.
% Input: the image file name, original image file, superpixel file, image
% data, the original mask. The mask input at the very beginning is 0 for
% all elements. Its size is the same as the sup_image.
% Output: selected_labels, result_mask, all_selected_labels

% input: file name, superpixel labels, pb edge data, denominator for
% the cost function, the max number of solutions, the selected 
% superpixels (limit the range?) edge threshold for support.
% output: X is the first solution's selected superpixels, Xs is the
% selected superpixles for all possible solutions, cost is the result of the minimum
% cost function, the sup-image is the label image for superpixels.
% X is the new mask, and sup_image will also be updated.

function [selected_labels, result_mask, all_selected_labels] = find_an_object(img_filename, img, ...
    sup_image, image_data, mask, edge_thresh)  % Mask the area of insterest as 0 return extracted mask as 1.

    sup_image = sup_image .* (~mask); %%Interested target in mask is 0.
    if(numel(unique(sup_image(:)))<5)
       selected_labels = [];
       result_mask = [];
       all_selected_labels = [];
       return       
    end
    disp('Looking for closures using parametric maxflow');

    [X, Xs, cost, sup_image] = SuperpixelClosureGrouping(img_filename, sup_image, ...
        image_data, 'area', 4, [], edge_thresh);  % call the superpixel grouping on the ground level
    % In the returned sup_image, the masked area is marked as size(X,1)+1
    disp('Postprocessing the results');
    Xs = ChooseLargestSupComponent(sup_image, Xs);  % choose a group of large number of superpixels

    % Remove duplicates
    s = 1;
    while (s < size(Xs,2))
        duplicates = all(Xs(:,(s+1):end) == repmat(Xs(:,s), [1, size(Xs,2)-s]), 1);
        Xs(:, [false(1,s), duplicates]) = [];
        s = s + 1;
    end

    % Remove all-1 and all-0 solutions
    n = max(sup_image(:));
    all_1 = sum(Xs, 1) == n;
    Xs(:, all_1) = [];

    all_0 = sum(Xs, 1) == 0;
    Xs(:, all_0) = [];
    
    % !!Corner-case solution !!!result_mask = X;
    if (isempty(Xs))
        selected_labels = [];
        result_mask = [];
        return;
    end

    selected_labels = Xs(:, 1);  % default to select the first one
    all_selected_labels = Xs; % (new_num_sup,candidate)
    
    selected_sup = find(selected_labels);  % selected_sup is list of index of selected superpixels
    result_mask = zeros(size(mask));
    for i = 1:numel(selected_sup)
        % for each element index in the selected_sup
        result_mask(sup_image == selected_sup(i)) = 1;  % make selected superpixels to be 1.
    end

end

