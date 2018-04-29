% This function will find a single object with a certain mask for the
% range.
% Input: the image file name, original image file, superpixel file, image
% data, the original mask. The mask input at the very beginning is 0 for
% all elements. Its size is the same as the sup_image.
% Output: selected_labels, result_mask, all_selected_labels
function [selected_labels, result_mask, all_selected_labels] = find_an_object(img_filename, img, ...
    sup_image, image_data, mask, edge_thresh)

    sup_image = sup_image .* (~mask);
    
    disp('Looking for closures using parametric maxflow');
    % input: file name, superpixel labels, pb edge data, denominator for
    % the cost function, the max number of solutions, the selected 
    % superpixels (limit the range?) edge threshold for support.
    % output: X is the first solution's selected pixels, Xs might be the
    % selected pixles for all possible solutions, cost is the result of the minimum
    % cost function, the sup-image is the label image for superpixels.
    % X is the new mask, and sup_image will also be updated.
    [X, Xs, cost, sup_image] = SuperpixelClosureGrouping(img_filename, sup_image, ...
        image_data, 'area', 10, [], edge_thresh);  % call the superpixel grouping on the ground level
    
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
    
%     result_mask = X;
    if (isempty(Xs))
        selected_labels = [];
        result_mask = [];
        return;
    end
    selected_labels = Xs(:, 1);  % default to select the first one
    all_selected_labels = Xs;
    
    selected_sup = find(selected_labels);  % find the selected superpixels
    result_mask = zeros(size(mask));
    for i = 1:numel(selected_sup)
        % for each element index in the selected_sup
        result_mask(sup_image == selected_sup(i)) = 1;  % make them 1.
    end
    
    
    
    
%     %% Save the solutions out to files
%     disp('Saving solutions');
%     % Save the figure images into a files
%     [pathstr, name, ext] = fileparts(img_filename);
%     results_img_file = [output_dir,'/',core_name,'/',name,'_multiplesolutions.jpg'];  % the file name to hold multiple solutions
% 
%     s = min([size(Xs,2), num_solutions]);  % the number of solutions might be smaller than the num_solutions
%     results_img = DrawSuperpixelsAreaIterationsSingleFigure(img, sup_image, Xs(:,1:s));  % get multiple solutions into an 3D array
%     imwrite(results_img, results_img_file, 'jpg');  % write the image to the file
% 
%     for sol = 1:s
%         % For each solution, write the foreground (white and black) to an
%         % image.
%         results_img_file = [output_dir,'/',core_name,'/',name,'_solution_',num2strPad(sol,3),'.jpg'];
%         fg = SupValueImage_MEX(sup_image, double(Xs(:,sol)));  % get the foreground image
%         imwrite(fg, results_img_file, 'jpg');  % the white and black images
%     end

end

