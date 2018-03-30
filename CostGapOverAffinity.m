function cost = CostGapOverAffinity(X, GAP, W)
    cost = (X'*GAP*X) / (X'*W*X);