function disp_img_total = DrawSuperpixelsAreaIterationsSingleFigure(img, ...
    sup_image, Xs, num_rows, num_cols, number_on)

    img = im2double(img);
    
    num_solutions = size(Xs, 2);
    
    if (~exist('num_rows', 'var') || isempty(num_rows))
        num_rows = floor(sqrt(num_solutions));  % the number of rows for multiple solutions in one image
        num_cols = ceil(num_solutions / num_rows);  % the number of columns for multiple solutions in one image
    end
    
    if (~exist('number_on', 'var') || isempty(number_on))
        number_on = true;
    end
    
    [M,N,C] = size(img);
    disp_img_total = zeros(M*num_rows, N*num_cols, C);  % the size of the single image, which holds multiple solutions.
    
    for i = 1:num_solutions
        % for each solution
        row = floor((i-1) / num_cols);
        col = mod((i-1), num_cols);
        
        if (nargout == 0)  % if the number of output is 0
            subplot('Position',[col/num_cols (num_rows-row-1)/num_rows 1/num_cols 1/num_rows]);
            iptsetpref('ImshowBorder', 'tight');
            DrawSuperpixelsFigure(img, sup_image, Xs(:,i));
        else
            disp_img = DrawSuperpixelsFigure(img, sup_image, Xs(:,i));  % draw the superpixel groups
            start_row = M*row+1;
            start_col = N*col+1;
            % assign one solution image to the total image array, which is
            % the output of this function
            disp_img_total(start_row:(start_row+M-1),start_col:(start_col+N-1),:) = im2double(disp_img);
        end
        
        
        if (number_on && nargout == 0)
            title(num2str(i));
        end
    end