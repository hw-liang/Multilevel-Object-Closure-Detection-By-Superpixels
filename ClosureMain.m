% Main function for finding closure in images:
%
% Arguments:
%   img_filename    - filename of the image to process (any format readable
%                     by Matlab
%   output_dir      - the output directory for the solutions
%   num_sups        - number of superpixels to extract
%   num_solutions   - number of closure solutions
%   edge_threshold  - Pb edge threshold (smaller values result in more
%                     solutions)
%   sup_method      - The method used to generate superpixels.
%                     Valid values: 'ncuts', 'turbo' (default)
%   Example: ClosureMain(img_filename, '.', 100, 10, 0.05, 'turbo');
%            This would open the image in img_file, extract 100 superpixels
%            using the Turbopixels method, use an edge threshold of 0.05,
%            generate at most 10 solutions and will store them in the
%            current directory.
%
function ClosureMain(img_filename, output_dir, num_sups, num_solutions, ...
    edge_thresh, sup_algorithm)
    
    %% setup files
    core_name = img_filename(1:end-4);
    if (~exist(core_name,'dir'))
        mkdir(core_name)
    end
    
    use_gpb = false;  % if use the global pb
    if (ispc)  % if it is a PC
        % Global Pb is not supported on windows
        use_gpb = false;
    end
    
    img = im2double(imread(img_filename));  % read a image and convert it to double precision
    img_copy = img;  % get a copy of the original image
    % Add base case here
    [m,n] = size(img);
    if (m*n) < 1000
        % If the total size of the image is less than 1000 pixels, return.
        disp('The input image is too small. The function does not continue.');
       return; 
    end
    
    %% Compute Pb
    disp('Computing Pb');
    pb_file = [core_name,'/',img_filename(1:end-4),'_pb.mat'];  % make a new name for the data file for the image
    if (~exist(pb_file, 'file'))  % if the pb file does not exist
        [pb, theta, tmap] = pbCGTG(img);  % get the boundary edge for the image
        save(pb_file, 'pb', 'theta', 'tmap');  % save the data to the pb_file, not need to calculate again next time
    else
        load(pb_file);  % If the pb_file exists, just load it.
    end

    if (use_gpb)  % if use the global pb edge detector
        disp('Computing gPb');
        cur_dir = pwd;
        cd globalPb;
        gpb_file = [core_name,'/',img_filename(1:end-4),'_gpb.mat'];
        rsz = 0.5;
        if (true || ~exist(gpb_file, 'file'))
            [gPb_thin, gPb, maxo] = globalPb(img_filename, 'temp.bmp', rsz);
            save(gpb_file, 'gPb_thin', 'gPb', 'maxo');
        end
        cd(cur_dir);
        if (exist(gpb_file, 'file'))
            load(gpb_file);
            pb = gPb_thin;
        end
    end
    
    disp('Postprocessing Pb results');
    image_data_file = [core_name,'/',img_filename(1:end-4),'_image_data.mat'];
    if (~exist(image_data_file, 'file'))
        % If the data does not exist, save the data to the file.
        image_data = PrecomputeImageData(img, pb, theta, tmap);
        save(image_data_file, 'image_data');
    else
        load(image_data_file);
    end
    
    
    %% Compute superpixels
    disp('Computing Superpixels');
    superpixels_file = [core_name,'/',img_filename(1:end-4),'_num_sups_',num2str(num_sups),'.seg'];
    if (~exist(superpixels_file, 'file'))
        % If the superpixel data has not once been saved to the
        % superpixels_file, it will be computed.
        switch (sup_algorithm)
            case 'ncuts'
                sup_image = Superpixels_Ncuts_Pb(img, num_sups, image_data_file);
            case 'turbo'                
                sup_image = superpixels_pb(img, num_sups, [], [], [], image_data.pb);
            case 'slic'
                sup_image = vl_slic(single(img), num_solutions, edge_thresh) + 1;
                % The returning value is a UINT32 array containing the
                % superpixel identifier (label) for each image pixel.
                % It is a matrix for labelling each pixel in the original
                % image.
            otherwise
                error('Unsupported superpixel method');
        end
        writeSeg(sup_image, superpixels_file);
    else
        % Otherwise, it will just read the data.
        sup_image = readSeg(superpixels_file);
    end
    sup_image = CleanSupImage(sup_image);  % fill the empty pixel locations

    
    %% Parametric Maxflow
    disp('Looking for closures using parametric maxflow');
    % input: file name, superpixel labels, pb edge data, denominator for
    % the cost function, the max number of solutions, the selected 
    % superpixels (limit the range?) edge threshold for support.
    % output: X is the first solution's selected pixels, Xs might be the
    % selected pixles for all possible solutions, cost is the result of the minimum
    % cost function, the sup-image is the label image for superpixels.
    
    [X2, Xs2, cost, sup_image] = SuperpixelClosureGrouping(img_filename, sup_image, ...
        image_data, 'area', num_solutions, [], edge_thresh);  % call the superpixel grouping on the ground level
    
    % test code:
    unselected_sup = find(X2);  % find the selected superpixels
    for i = 1:numel(unselected_sup)
        % for each element index in the selected_sup
        sup_image(sup_image == unselected_sup(i)) = 0;  % make them 0.
    end
    
    [X, Xs, cost, sup_image] = SuperpixelClosureGrouping(img_filename, sup_image, ...
        image_data, 'area', num_solutions, [], edge_thresh);  % call the superpixel grouping on the ground level
    
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

    %% Save the solutions out to files
    disp('Saving solutions');
    % Save the figure images into a files
    [pathstr, name, ext] = fileparts(img_filename);
    results_img_file = [output_dir,'/',core_name,'/',name,'_multiplesolutions.jpg'];  % the file name to hold multiple solutions

    s = min([size(Xs,2), num_solutions]);  % the number of solutions might be smaller than the num_solutions
    results_img = DrawSuperpixelsAreaIterationsSingleFigure(img, sup_image, Xs(:,1:s));  % get multiple solutions into an 3D array
    imwrite(results_img, results_img_file, 'jpg');  % write the image to the file

    % Crop the image here, use the first solution as the target
    % Use the original image, superpixel image, and the binary labels
%     [m,n,c] = size(img);
%     h_min = intmax; h_max = intmin; v_min = intmax; v_max = intmin;
%     X_labels = Xs(:,1);
%     X_labels = X;
%     for i = 1:m
%         for j = 1:n
%             label = sup_image(i,j);
%             if (label > 0)
%                if X_labels(label) == 1
%                     % in the group of selected superpixels
%                     if (j < h_min)
%                         h_min = j;
%                     end
%                     if (j > h_max)
%                         h_max = j;
%                     end
%                     if (i < v_min)
%                         v_min = i;
%                     end
%                     if (i > v_max)
%                         v_max = i;
%                     end
%                end
%             end
%         end
%     end
%     v_min, h_min, v_max, h_max
%     test_crop = imcrop(img, [h_min v_min (h_max-h_min) (v_max-v_min)]); 
%     imshow(test_crop);
%     imwrite(test_crop, 'test.jpg');
    for sol = 1:s
        % For each solution, write the foreground (white and black) to an
        % image.
        results_img_file = [output_dir,'/',core_name,'/',name,'_solution_',num2strPad(sol,3),'.jpg'];
        fg = SupValueImage_MEX(sup_image, double(Xs(:,sol)));  % get the foreground image
        imwrite(fg, results_img_file, 'jpg');  % the white and black images
    end
%     ClosureMain('test.jpg', '.', 100, 10, 0.05, 'slic');
