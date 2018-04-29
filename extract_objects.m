% The function will extract all objects at different level from a single image.
% The input will be the original image, mask, and the tree node.
% The output will be all the masks saved to files.
% The function will also build a tree structure for the relationships.
function [] = extract_objects(img_filename, output_dir, img, num_sups, sup_image, image_data, ...
    mask, edge_thresh, counter)
    % Each element of the initial mask is 0, meaning none is selected.
    disp(['extract_objects is called with the counter', num2str(counter)]);
    s = sum(~mask(:));
    [m,n] = size(sup_image);
    div = max(sup_image(:));
    threshold = (m*n)/div;
    if (s < threshold)
        return;
    end
    
    % the masks_holder have all masks for the objects in this level.
    masks_holder = find_objects_same_level(img_filename, output_dir, img, threshold, sup_image, image_data, ...
    mask, edge_thresh, counter);
    
    counter = counter + 1;
    
    [a,b] = size(sup_image);
    while (~isempty(masks_holder))
        % when there is still an object to explore its next level objects
        object_mask = masks_holder(:, 1:b);
        masks_holder = masks_holder(:, b+1 : end);
        extract_objects(img_filename, output_dir, img, num_sups, sup_image, image_data, ~object_mask, edge_thresh, counter);
    end
end

