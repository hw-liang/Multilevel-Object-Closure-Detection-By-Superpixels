% Detect closures for all the images in a given directory
%
function GenerateResultsDirectory(img_dir, num_sups, num_solutions, edge_thresh, sup_method, output_dir)
    if (nargin < 2 || isempty(num_sups))
        num_sups = 100;
    end
    
    if (nargin < 3 || isempty(num_solutions))
        num_solutions = 5;
    end
    
    if (nargin < 4 || isempty(edge_thresh))
        edge_thresh = 0.05;
    end
    
    if (nargin < 5 || isempty(sup_method))
        sup_method = 'turbo';
    end
    
    if (nargin < 6 || isempty(output_dir))
        output_dir = img_dir;
    end
    
    files = dir([img_dir,'/*.jpg']);
    
    img_idx = 1:numel(files);
    
    num_images = numel(img_idx);
    
    for i = 1:num_images
        disp(['Processing image ',num2str(i),' out of ',num2str(num_images)]);

        img_filename = [img_dir,files(img_idx(i)).name];

        ClosureMain(img_filename, output_dir, num_sups, num_solutions, edge_thresh, sup_method);
    end