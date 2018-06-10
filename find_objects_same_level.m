% This function will find all objects at the same level.
% The input will be:
% img_filename, original image, superpixel image, image data, the current mask, 
% and the current tree node.
% The output will the the masks for all the objects at the same level.
% An mask is an 2D array with the same size as the original image,
% but it will just have 0 and 1 as elements to indicate which area is
% valid, which area is not valid.
function [result_masks] = find_objects_same_level(img_filename, output_dir, img, threshold, sup_image, image_data, ...
    mask, edge_thresh, counter,counter2,counter3)   % Mask the area of insterest as 0
    disp(['find_objects_same_level is called with the counter', num2str(counter)]);
    objects_same_level = [];
    [selected_labels, a_mask] = find_an_object(img_filename, img, sup_image, image_data, mask, edge_thresh);
                                % Mask the area of insterest as 0
    objects_same_level = [objects_same_level, a_mask];
    object_size = sum(~mask(:)) - sum(a_mask(:));

    if((object_size ==0 || sum(a_mask(:))==0))
        result_masks = [];
        return    
    end
    
    while (object_size > threshold)  % the condition can be proportional to the size of the image (e.g. (total_size/100)) or a fixed threshold
        mask = mask | a_mask;  % may not be true
        [selected_labels, a_mask] = find_an_object(img_filename, img, sup_image, image_data, mask, edge_thresh);
        object_size = sum(a_mask(:));
        objects_same_level = [objects_same_level, a_mask];
    end   
    result_masks = objects_same_level;
    
    % warp things here
    disp('Saving solutions');
    core_name = img_filename(1:end-4);
    % Save the figure images into a files
    [pathstr, name, ext] = fileparts(img_filename);
    [m,n] = size(result_masks);
    [a,b] = size(sup_image);
    s = n/b  % The number of found objects.
    Xs = zeros(max(sup_image(:)), s); % to hold all objects' labels for the original sup_image.
    results_img_file2 = [output_dir,'/',core_name,'/',name,'_', num2str(counter),'_',num2str(counter3),'_',num2str(counter2), '_multiplesolutions.jpg'];  % the file name to hold multiple solutions

    for sol = 1:s
        % For each solution, write the foreground (white and black) to an
        % image.
        results_img_file = [output_dir,'/',core_name,'/',name, '_', num2str(counter), '_solution_',num2strPad(sol,3),'.jpg'];
        fg = result_masks(:, (sol-1)*b + 1 : sol*b);
        imwrite(fg, results_img_file, 'jpg');  % the white and black images
        object_sup = sup_image .* fg;
        temp_holder = unique(object_sup);
        temp_h = Xs(:, sol);
        temp_holder = temp_holder(temp_holder ~= 0);
        temp_h(temp_holder) = 1;
        Xs(:, sol) = temp_h;
    end
    Xs_size = size(Xs);
    if (Xs_size(2) ~= 0)
        results_img = DrawSuperpixelsAreaIterationsSingleFigure(img, sup_image, Xs(:,1:s));  % get multiple solutions into an 3D array
        imwrite(results_img, results_img_file2, 'jpg');  % write the image to the file
    else
        disp('No result.');
    end
    
end

