% SuperpixelAdjacencyGraphNcuts - Compute the adjacency matrix for
% superpixels for a seg_image computer with Ncuts (every pixel is in range
% 1..K where K is the number of labels). If a pixel is 0, it is unlabeled
%
% Arguments:
%   seg_image - superpixel segmentation. MxN matrix.
%
% Returns:
%   adjmat - num_sups x num_sups adjacency matrix. adjmat(i,j) =
%   adjmat(j,1) = 1 iff superpixels i and j are adjacent.
%
function adjmat = SuperpixelAdjacencyGraphNcuts(seg_image)

    num_sups = max(seg_image(:));
    adjmat = zeros(num_sups, num_sups);
    [M,N] = size(seg_image);
    
    for x = 1:M-1
        for y = 1:N-1
            s = seg_image(x,y);
            s1 = seg_image(x+1, y);
            s2 = seg_image(x, y+1);
            s3 = seg_image(x+1, y+1);
            
            if (s > 0 && s1 > 0)
                adjmat(s, s1) = 1;                
                adjmat(s1, s) = 1;
            end
            
            if (s > 0 && s2 > 0)
                adjmat(s, s2) = 1;
                adjmat(s2, s) = 1;
            end
            
            if (s > 0 && s3 > 0)
                adjmat(s, s3) = 1;
                adjmat(s3, s) = 1;
            end
        end
    end