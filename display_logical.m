function disp_img = display_logical(img, contour, color, brightness)
    
    if (size(img,3) == 1)
        img = repmat(img,[1,1,3]);
    end
    disp_img = img;
    
    if (nargin < 4)
        brightness = 1;
    end
    
    brightness = min(1, max(0, brightness));
    gamma = 1-brightness;
    img = img .^ gamma;

    ind = logical(zeros(size(img,1),size(img,2),3));
    ind(:,:,1) = contour;
    disp_img(ind) = color(1);
    ind = logical(zeros(size(img,1),size(img,2),3));
    ind(:,:,2) = contour;
    disp_img(ind) = color(2);
    ind = logical(zeros(size(img,1),size(img,2),3));
    ind(:,:,3) = contour;
    disp_img(ind) = color(3);
    disp_img = disp_img ./ max(disp_img(:));
    disp_img = uint8(disp_img * 255);
