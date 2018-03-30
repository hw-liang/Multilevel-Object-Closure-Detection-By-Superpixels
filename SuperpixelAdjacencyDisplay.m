% SuperpixelAdjacencyDisplay - Display the adjacency graph for superpixels
% overlaid on the superpixel segmentation.
%
% Arguments:
%   img - MxNx3 or MxN image.
%   sup_image - superpixel segmentation with 0 elements indicating boundary
%   or background pixels. MxN matrix.
%   adjmat - num_sups x num_sups adjacency matrix. adjmat(i,j) =
%   adjmat(j,1) = 1 iff superpixels i and j are adjacent.
%   affinity - num_sups x num_sups matrix of affinities between
%   superpixels. Pairs with higher affinities will be represented by
%   thicker lines. Default: all affinities are 1.
%
function SuperpixelAdjacencyDisplay(img, sup_image, adjmat, affinity)

    if (nargin < 4)
        affinity = ones(size(adjmat));
    end
    max_affinity = max(affinity(:));
    max_width = 10;
    min_width = 2;
    affinity = affinity ./ max_affinity;
    stats = regionprops(sup_image, 'Centroid');
    num_sups = numel(stats);
    
    center = reshape([stats(:).Centroid], [2, num_sups]);
    
    boundary = sup_image==0;
    if (sum(boundary(:)) == 0);
        boundary = seg2bmap(sup_image);
    end
    
    boundary = bwmorph(boundary, 'dilate', 0);
    
    imagesc(display_logical(img, boundary, [0,1,0]));
    axis off;
    axis image;
%     title(['Affinity between superpixels (max\_affinity = ',num2str(max_affinity), ')']);
%     hold on;
    
    [sup1, sup2] = find(tril(adjmat, -1));
    num_edges = numel(sup1);
    
    for i = 1:num_edges
        width = min_width + (max_width - min_width) * abs(affinity(sup1(i), sup2(i)));
        width = min(max(width, min_width), max_width);
             
        % Display in Green if the affinity is undefined, Blue if defined
        if (affinity(sup1(i), sup2(i)) < 0)
            line(center(1, [sup1(i),sup2(i)]), center(2, [sup1(i),sup2(i)]), 'LineWidth', width, 'Color', 'g');
        else
            line(center(1, [sup1(i),sup2(i)]), center(2, [sup1(i),sup2(i)]), 'LineWidth', width, 'Color', 'y')
        end
    end