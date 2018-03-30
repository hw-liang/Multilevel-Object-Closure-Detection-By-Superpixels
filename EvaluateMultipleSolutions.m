% Computes the costs of multiple solution given a function that can compute
% a cost for a single solution
function costs = EvaluateMultipleSolutions(Xs, eval_func, varargin)

    num_evals = size(Xs,2);
    costs = zeros(num_evals, 1);
    
    for i=1:num_evals
        try
            costs(i) = eval_func(Xs(:,i), varargin{:});
        catch
            costs(i) = inf;
        end
    end