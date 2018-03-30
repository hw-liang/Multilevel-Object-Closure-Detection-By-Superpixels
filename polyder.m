% Given a cell array of polynomial fragments parameters, compute the
% polynomial parameters of the derivative of a given order
function poly_params_der = polyder(poly_params, der_order) 

    order = cellfun(@(x) ((size(x,1):-1:1)'-1), poly_params, 'UniformOutput', false);
    
    assert(der_order == 1);
    poly_params_der = cellfun(@(x,y) (repmat(x, [1,2]) .* y), order, poly_params, 'UniformOutput', false);
%     poly_params_der = cellfun(@(x,y) (repmat(factorial(x)./factorial(max(x-der_order,0)), [1,2]) .* y), order, poly_params, 'UniformOutput', false);
    
    poly_params_der = cellfun(@(x) (x(1:end-der_order, :)), poly_params_der, 'UniformOutput', false);