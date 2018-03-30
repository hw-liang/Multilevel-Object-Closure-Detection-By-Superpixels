function seg_image = Superpixels_Ncuts_Pb(img, num_sups, image_data_file)

    if (nargin >= 3)
        [Sp,NcutDiscrete,NcutEigenvectors] = NcutImage(img, num_sups, image_data_file);
    else
        [Sp,NcutDiscrete,NcutEigenvectors] = NcutImage(img, num_sups);
    end
    
    % Connected components analysis, cut apart disconnected clusters
    % Clean up Sp.  Merge small segments, break disconnected clusters.
    MIN_SIZE = 30;
    min_c = min(Sp(:)); max_c = max(Sp(:));
    new_num = max_c+1;
    for c_i=min_c:max_c
        [L,num] = bwlabel(Sp==c_i,4);
        if num > 1 || sum(L(:)) < MIN_SIZE
            for n_i=1:num
                the_inds = find(L==n_i);
                if length(the_inds) < MIN_SIZE
                    % Merge each small segment with a nearby one.
                    the_perim = bwperim(imdilate(L==n_i,strel('diamond',1)));
                    the_perim(the_inds)=0;  % bwperim has errors on boundary.
                    nhbs = find(the_perim);
                    the_dists = dist2(NcutEigenvectors(the_inds,:),NcutEigenvectors(nhbs,:));
                    the_dists = min(the_dists,[],1);
                    [mind,mind] = min(the_dists);
                    Sp(the_inds) = Sp(nhbs(mind));
                else
                    % Renumber large segments.
                    Sp(the_inds) = new_num;
                    new_num = new_num+1;
                end
            end
        end
    end

    % Renumber segments.
    Spv = Sp';
    Spv = Spv(:);
    uu = unique(Spv);
    N_seg = length(uu);
    the_map(uu) = 1:N_seg;
    Spv = the_map(Spv);
    seg_image = reshape(Spv,[size(img,2) size(img,1)])';

