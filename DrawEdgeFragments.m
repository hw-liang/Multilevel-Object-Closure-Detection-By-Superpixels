% Draw edge fragments in a given color with a given weight
%
% fragments is a cell array of fragments with every cell storing an Mx2
% matrix of (x,y) fragment coordinates
% color is a 3x1 vector representing the color of the fragments
% weight is a vector os  weights for each fragment (will affect the color)
function DrawEdgeFragments(img, fragments, color, weight)

    num_fragments = numel(fragments);
    
    if (nargin < 4)
        weight = ones(size(fragments));
    end
    
    if (isscalar(weight))
        weight = repmat(weight, size(fragments));
    end
    
    imshow(img);
    hold on;
    
    for i = 1:num_fragments
        plot(fragments{i}(:,1), fragments{i}(:,2), 'Color', color, 'LineWidth', weight(i));
    end
    
    hold off;