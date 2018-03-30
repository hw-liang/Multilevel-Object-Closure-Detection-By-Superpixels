function cost = CostCutOverAffinity(X, W)
    D = diag(sum(W,2));
    L = D - W;
    
    A = double(W > 0);
    DA = diag(sum(A,2));
    LA = DA - A;
    cost = ((X'*L*X) / (X'*LA*X) ) / ((X'*W*X) / (X'*A*X) );
%    cost = (X'*L*X) / (X'*W*X);