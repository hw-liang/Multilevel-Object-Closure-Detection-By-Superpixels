% Given an image, a superpixel image and which sups are selected in the
% foreground, draw the foreground region.

function disp_img = DrawSuperpixelsFigure(img, sup_image, X, filename)

    fg = SupValueImage_MEX(sup_image, double(X));  % get the foreground

    img(repmat(fg, [1,1,size(img,3)]) == 0) = ...
        img(repmat(fg, [1,1,size(img,3)]) == 0).^0.3;
    
    
    disp_img = display_logical(im2double(display_logical(img, ...
        bwmorph(seg2bmap(fg),'dilate', 5), [1,0,0])), bwmorph(seg2bmap(sup_image), 'dilate', 1), [0,1,0]);
    
    if (nargout == 0)
        imshow(disp_img);
    end
    
    if (exist('filename', 'var') && ~isempty(filename))
        imwrite(disp_img, filename, 'jpg');
    end