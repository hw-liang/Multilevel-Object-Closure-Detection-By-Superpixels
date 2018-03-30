function Xs = ChooseLargestSupComponent(sup_image, Xs)

    num_solutions = size(Xs, 2);
    for s = 1:num_solutions
        fg = SupValueImage_MEX(sup_image, double(Xs(:,s)));
        label_img = bwlabel(fg);
        if (max(label_img(:))<1)
            fg = SupValueImage_MEX(sup_image, double(Xs(:,s-1)));
            label_img = bwlabel(fg);        
        end
        region_stats = regionprops(label_img, 'Area', 'PixelIdxList');

        % Find the region with the largest area
        [area, max_idx] = max([region_stats(:).Area]);
        new_fg_img = zeros(size(fg));
        
        new_fg_img(region_stats(max_idx).PixelIdxList) = 1;
        Xs(:,s) = Mean_Superpixels_Image_MEX(sup_image, new_fg_img) > 0.5;
        
    end