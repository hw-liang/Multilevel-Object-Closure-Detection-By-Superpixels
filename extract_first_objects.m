% The function will extract all objects at first level.
% The input will be the original image, mask(mask the area of interest as 0), and the tree node.
% The output will be all the masks saved to files.
% The function will also build a tree structure for the relationships.

function [] = extract_first_objects(img_filename, output_dir, img, num_sups, sup_image, image_data, ...
    mask, edge_thresh, counter,counter2, counter3)
    counter3 = 1;
    disp(['extract_objects is called with the counter', num2str(counter)]);
    
    area_interest = sum(~mask(:));
    [m,n] = size(sup_image);
    
    threshold = (m*n)/num_sups;
    if (area_interest < threshold)
        return;
    end
    
    % the masks_holder have all masks for the objects in this level.
    masks_holder = find_objects_same_level(img_filename, output_dir, img, threshold, sup_image, image_data, ...
    mask, edge_thresh, counter,counter2,counter3);
    
    counter = counter + 1;
    counter2 = 0;
    
    while (~isempty(masks_holder))
        % when there is still an object to explore its next level objects
        object_mask = masks_holder(:, 1:n); %% Target object is 1.
        masks_holder = masks_holder(:, n+1 : end);
        counter2 = counter2 +1;
        counter3 = counter3 +1;
        extract_first_objects(img_filename, output_dir, img, num_sups, sup_image, image_data, ~object_mask, edge_thresh, counter,counter2,counter3);    
    end
end

