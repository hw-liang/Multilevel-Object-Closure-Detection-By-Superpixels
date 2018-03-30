% Ensures that 4-connectivity of all superpixels
%
function sup_image = CleanSupImage(sup_image)

    [M,N] = size(sup_image);  % return the numbers of rows and columns
    new_sup_image = zeros(M,N);

    num_sups = max(sup_image(:));
    
    for sup = 1:num_sups
        cur_sup_BW = (sup_image == sup);
        labels = bwlabel(cur_sup_BW,4);
        if (max(labels(:)) > 1)
            stats = regionprops(labels, 'Area', 'PixelIdxList');
            [s, ind] = sort([stats.Area]);
            new_sup_image(stats(ind(end)).PixelIdxList) = sup;
        else
            new_sup_image(cur_sup_BW) = sup;
        end
    end
    
    % Find unlabeled pixels and label them with their surrounding
    % superpixel labels
    [r,c] = find(new_sup_image == 0);  % find all zero pixels in the 
    for p = 1:numel(r)
        neigh_sups = zeros(4,1);
        if (r(p)-1 > 0)
            neigh_sups(1) = sup_image(r(p)-1, c(p));
        end

        if (r(p)+1 <= M)
            neigh_sups(2) = sup_image(r(p)+1, c(p));
        end

        if (c(p)-1 > 0)
            neigh_sups(3) = sup_image(r(p), c(p)-1);
        end

        if (c(p)+1 <= N)
            neigh_sups(4) = sup_image(r(p), c(p)+1);
        end

        neigh_sups = neigh_sups(neigh_sups ~= 0);

        new_sup_image(r(p),c(p)) = mode(neigh_sups);  % use the most frequent number to fill the empty pixel
    end
    
    sup_image = new_sup_image;