% Main function for computing closure
% Most of the function contains experimental stuff.
% Tested setting for parameters are:
% recursive_on = false; do_padding = true; alg_type = 'area';

function [X, Xs, cost, sup_image, costs] = SuperpixelClosureGrouping(img_filename, ...
    sup_image, image_data, alg_type, num_solutions, selected_sups, edge_thresh)
       
    display_on = false;  % default false
    recursive_on = false;  % default false. This function has not been completed yet.
    
    try
        img = im2double(imread(img_filename));
    catch
        img = img_filename;
    end


    if (~exist('edge_thresh', 'var'))
        edge_thresh = [];
    end
   
    num_sups = max(sup_image(:));  % number of superpixels
    % computer the adjacent matrix that shows the adjaceny relationships
    % among superpixels.
    adjmat = SuperpixelAdjacencyGraphNcuts(sup_image); 
        
    if (exist('selected_sups', 'var') && ~isempty(selected_sups))
        % If there are selected superpixels
        num_sups = sum(selected_sups);  % total number of selected superpixels
        adjmat = adjmat(selected_sups, selected_sups);  % make a adjacent matrix only for the selected superpixels
        % SupValueImage_MEX will return the foreground (white) and
        % background image (black)
        sup_image(SupValueImage_MEX(sup_image, double(selected_sups)) == 0) = 0;  % make the previous selected pixels to be 0
        % [b,m,new_sups] = unique(sup_image(:));
        % new_sups = unique(sup_image(:)); % These two lines may not be needed.
        % sup_image = reshape(new_sups, size(sup_image)) - 1;
    end    

    % Compute all the edge fragments between superpixels
    
    % Include (or not) the fragments around the image border
    do_padding = true;  % default as true
    if (do_padding)
        s = padarray(sup_image, [1,1]); % add 0s to both horizontal and vertical boundaries
        s(s == 0) = num_sups+1;  % use another label for all pixels with label 0
        % input: padded new label image, original image padded, 
        [fragments, junctions, neighbor_data, sup_image, avgcolor, ...
            polyfragments, poly_params, T] = seg2fragments(s, ...
            padarray(img,[2,2,0]), 10, 3);
        % The first two outputs are the angles and tangents for the edges
        % of superpixels.
        [poly_angle, poly_tangent, poly_sq_curv] = ...
            polyfragments_angle_and_curvature(T, poly_params);

        fragments = cellfun(@(x) (x - 1), fragments, 'UniformOutput', false);
        sup_image = sup_image(2:end-1, 2:end-1);
        sup_image(sup_image == (num_sups+1)) = 0;
        segment_fragments = neighbor_data.segment_fragments(1:end-1);
        num_sups = numel(segment_fragments);

        fragment_segments = cellfun(@(x) (x(x <= num_sups)), ...
            neighbor_data.fragment_segments, 'UniformOutput', false);
    else
        [fragments, junctions, neighbor_data, sup_image, avgcolor, polyfragments, poly_params, T] = seg2fragments(sup_image, img, 10, 3);
        [poly_angle, poly_tangent, poly_sq_curv] = polyfragments_angle_and_curvature(T, poly_params);

        segment_fragments = neighbor_data.segment_fragments;
        num_sups = numel(segment_fragments);

        fragment_segments = neighbor_data.fragment_segments;
    end
    
    STATS = regionprops(sup_image, 'Area');  % calculate the areas for each superpixel 
    A = diag([STATS(:).Area]);  % create a matrix with the diagonal as the areas of superpixels
    
    gap_multiplier = ones(size(A));  % for G
    edge_affinity = ones(size(A));  % for E
    
    % Choose the closure cost type (currently only 'area' is well tested)
    try
        if (exist('alg_type', 'var') && ~isempty(alg_type))
            switch alg_type
                case 'area'
                    affinity = A;
                case 'appearance_affinity'
                    fg = double(rand(size(sup_image)) > 0.5);
                    [fg_prob, bg_prob, rgbhist] = ForegroundBackgroundProbColorHist(img, fg, linspace(0,1,5), linspace(0,1,5), linspace(0,1,5));
                    hist_affinities = SuperpixelAffinityHistogram(sup_image, adjmat, rgbhist, 125);
                    affinity = hist_affinities;
                    % gap_multiplier(eye(size(adjmat)) == 0) = affinity(eye(size(adjmat)) == 0);
                case 'junction_affinity'
                    [affinity_junc, affinity_junc_on] = Junc3AffinityOne(img, sup_image, image_data);
                    affinity = (affinity_junc .* adjmat);
                    % contour_affinity = ContourAffinity(affinity, neighbor_data, fragment_segments);
                    % edge_affinity = contour_affinity;
                    % [D, P, L, maxEdge] = ShortestPaths(affinity);
                    % maxD = max(D(~isinf(D)));
                    % D(isinf(D)) = maxD;
                    % D = D / maxD;
                    % affinity = exp(-D*10);
                    % affinity(eye(size(affinity))==1) = 0;

                case 'junction_and_area'
                    affinity = affinity_junc .* (diag(A)*diag(A)');
                case 'junction_area_weighted_gap'
                    affinity = (affinity_junc .* adjmat);
                    gap_multiplier = (diag(A)*diag(A)');
                otherwise
                    error('Unknown affinity');
            end
        else
            affinity = A;
        end
    catch
        error('Unsupported affinity type');
    end
    
    [E, P, frag_edge, frag_per] = SuperpixelAffinityEdgeContrastGivenFragments(...
        fragments, segment_fragments, fragment_segments, ...
        poly_angle, poly_sq_curv, image_data, edge_thresh, edge_affinity);
    G = P - E;
    if (size(G) ~= size(gap_multiplier))
        gap_multiplier(end,:) = [];
        gap_multiplier(:, end) = [];
        affinity(end,:) = [];
        affinity(:, end) = [];
    end
    G = G .* gap_multiplier; % G is the edge gap matrix
    D = sum(affinity,2);  % a column with sums of rows; each element is the area of a superpixel
    % X = MinimizeRatioBinary(double(rand(num_sups,1)>0.5),G, affinity);
    
    % Can call closure extraction recursively on each solution and its
    % complement (used for recursively extract multiple objects)
    if (recursive_on)
        if (num_solutions == 1)
            Xs = X;
        else
            Xs1 = [];
            if (sum(X > 0.5) > 1)  % if there is grouping of superpixels
                [X1, Xs1, cost1, sup_image1] = SuperpixelClosureGrouping(img_filename, ...
                    sup_image, image_data, num_solutions-1, X > 0.5);  % recursively call the self function
                repetitions = sum((Xs1 == repmat(X, 1, size(Xs1, 2))), 1) == numel(X);   % find the repetitions
                Xs1 = Xs1(:, ~repetitions);  % get the groupings not repeated
            end
            Xs2 = [];
            if (sum(X < 0.5) > 1)
                [X2, Xs2, cost2, sup_image2] = PartGrouping2(img_filename, sup_image, image_data, num_solutions-1, X < 0.5);
                repetitions = sum((Xs2 == repmat(X, 1, size(Xs2, 2))), 1) == numel(X);
                Xs2 = Xs2(:, ~repetitions);
            end
            Xs = [X, Xs1, Xs2];
            e = EvaluateMultipleSolutions(Xs, @CostGapOverAffinity, G, diag(D));
            e(isnan(e) | isinf(e)) = inf;
            [costs,idx] = sort(e);
            Xs = Xs(:,idx(1:min(size(Xs,2), num_solutions)));
            costs = costs(1:min(size(Xs,2), num_solutions));
        end
    else
        % if recursion is not on; default option
        % use parametric maxflow to find the grouping the the superpixels
        % Xs is the matrix with each column as a grouping; lambda is a
        % vector with all possible lambda values.
        [Xs,lambda] = ParametricMaxflow_MEX([D';diag(G)'], 2*G, -1, 0, 100);
        % The second argument is a function pointer.
        e = EvaluateMultipleSolutions(Xs, @CostGapOverAffinity, G, diag(D));  % Computes the costs of multiple solutions given a cost function
        e(isnan(e) | isinf(e)) = inf;
        [costs,idx] = sort(e);  % sort the cost results
        Xs = Xs(:,idx(1:min(size(Xs,2), num_solutions)));  % ordered the groupings by the ascending cost order
        costs = costs(1:min(size(Xs,2), num_solutions));  % also order the costs
        X = Xs(:,1);  % pick the first column as the best grouping
    end
    
    if (display_on)  % If the display is on, it will show the results.
        figure(1); clf;
        DrawSuperpixelsFigure(img, sup_image, X);
        figure(2); clf;
        DrawEdgeFragments(img, fragments, [0,1,0], 5*(frag_edge + eps)./(frag_per + 10*eps)+eps);
        figure(3); clf;
        SuperpixelAdjacencyDisplay(img, sup_image, adjmat, affinity);
        figure(4); clf;
        DrawSuperpixelsAreaIterationsSingleFigure(img, sup_image, Xs);
        pause;
    end
    
    cost = costs(1);  % pick the first cost
    
    if (exist('selected_sups', 'var') && ~isempty(selected_sups))
        X_new = zeros(numel(selected_sups), 1);
        Xs_new = zeros(numel(selected_sups), size(Xs,2));
        X_new(selected_sups) = X;
        Xs_new(selected_sups, :) = Xs;
        X = X_new;
        Xs = Xs_new;
    end
