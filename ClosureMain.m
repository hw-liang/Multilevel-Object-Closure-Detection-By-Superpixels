% Main function for finding closure in images:
%   Arguments:
%   img_filename    - filename of the image to process (any format readable
%                     by Matlab
%   output_dir      - the output directory for the solutions
%   num_sups        - number of superpixels to extract
%   edge_threshold  - Pb edge threshold (smaller values result in more
%                     solutions)
%   sup_method      - The method used to generate superpixels.
%                     Valid values: 'ncuts', 'turbo' (default)
%   Example: ClosureMain(img_filename, '.', 100, 0.05, 'ncuts');
%            This would open the image in img_file, extract 100 superpixels
%            using the Turbopixels method, use an edge threshold of 0.05,
%            generate at most 10 solutions and will store them in the
%            current directory.

function ClosureMain(img_filename, output_dir, num_sups, edge_thresh, sup_algorithm)
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
    % Add base case here
    [m,n] = size(img);
    if ((m*n) < 1000)
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
            otherwise
                error('Unsupported superpixel method');
        end
        writeSeg(sup_image, superpixels_file);
    else
        % Otherwise, it will just read the data.
        sup_image = readSeg(superpixels_file);
    end
    sup_image = CleanSupImage(sup_image);  % fill the empty pixel locations
    % Add mask
    original_mask = zeros(size(sup_image));  % The original mask is all 0s, and the dimension is the same as the sub_image.
    counter = 0;
    counter2 = 1;
    counter3 = 1;
    extract_first_objects(img_filename, output_dir, img, num_sups, sup_image, image_data, original_mask, edge_thresh, counter,counter2,counter3);
    
