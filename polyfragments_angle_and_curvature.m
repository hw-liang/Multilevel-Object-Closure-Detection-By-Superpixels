% Computes the angles, tangent angles, and squared curvature for all the
% points on all the given fragments
%
% Arguments:
%   T  -  a cell array of arclength parameters for each edge fragments
%   poly_params  - a cell array of polynomial parameters for each edge
%                  fragment
%   Thus the i-th edge fragment can be approximated with
%   T{i}*poly_params{i}
%
%  Returns:
%  poly_angle - A cell array of tangent angles (mod pi) for each point on
%  each edge fragment
%  poly_tangent - A cell array of tangent angles for each point on
%  each edge fragment
%  poly_sq_curv - A cell array of squared curvature at each point on each
%  edge fragment
%
function [poly_angle, poly_tangent, poly_sq_curv] = polyfragments_angle_and_curvature(T, poly_params)

    pder = polyder(poly_params, 1);
    pder2 = polyder(pder, 1);
    
    xy_der = cellfun(@(t,p) (t(:,2:end) * p), T, pder, 'UniformOutput', false);
    xy_der2 = cellfun(@(t,p) (t(:,3:end) * p), T, pder2, 'UniformOutput', false);
    
    poly_tangent = cellfun(@(xy) (angle(complex(xy(:,1), xy(:,2)))), xy_der, 'UniformOutput', false);
    poly_angle = cellfun(@(a) (mod(a, pi)), poly_tangent, 'UniformOutput', false);
    
    poly_sq_curv = cellfun(@(x,xx) ...
        ((xx(:,2).*x(:,1) - xx(:,1).*x(:,2)) ./ (x(:,1).^2 + x(:,2).^2).^1.5), ...
        xy_der, xy_der2, 'UniformOutput', false);
    
    poly_sq_curv = cellfun(@(x) (x.^2), poly_sq_curv, 'UniformOutput', false);
    