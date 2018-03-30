#include <stdio.h>
#include <math.h>
#include "parametric.h"
#include "mex.h"

// print cut C_k (0<=k<=N where N is the number returned by Solve()).
void PrintCut(Parametric* p, int numNodes, int k)
{
	int i;

	printf("Cut (");
	for (i=0; i<numNodes; i++)
	{
		int x = (p->GetRegionCount(p->GetRegion(i)) >= k) ? 1 : 0; // label of node i in C_k
		printf(" %d ", x);
	}
	printf(")");
}

//Assigns the cut C_k to the matrix of returned cuts
void AssignCutToSolution(Parametric* p, int numNodes, int k, mxLogical* selectedNodes, int num_solutions, int num_breakpoints)
{
	int i;

    if (k >= num_breakpoints-num_solutions) {
        for (i=0; i<numNodes; i++)
        {
            *(selectedNodes + (k-(num_breakpoints-num_solutions))*numNodes + i) = (p->GetRegionCount(p->GetRegion(i)) >= k); // label of node i in C_k
        }
    }
}
// Minimize energy of the type 
// sum_i((nodeCost(1,i)*lambda+nodeCost(2,i))*X_i) + sum_ij(X_i*X_j*edgeCost(i,j))
//
// Parameters:
//    nodeCost  - 2xN unary node costs including lambda's as described above
//    edgeCost  - pairwise costs (can be either an NxN matrix or 3xE matrix
//                were each row corresponds to an edge (i,j, and cost)
//                i and j are 0-based
//    lambdaMin - minimum lambda for the search
//    lambdaMax - maximum lambda for the search
//    num_sol   - number of requested solutions
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
		const mxArray *prhs[])
{
	double *nodeCost;
	mxLogical *selectedNodes;
	int numNodes;
	int num_solutions = 100000;
    double *edgeCost;
    double *lambdas;
    
	/* The input must be a noncomplex scalar double.*/
	numNodes = (int)mxGetN(prhs[0]);  // get the number of columns in an array, which is the number of nodes.
    // Each node represents a superpixel.
    printf("Number of nodes: %d \n", numNodes);
    
	/* Assign pointers to each input and output. */
	nodeCost = mxGetPr(prhs[0]);
	edgeCost = mxGetPr(prhs[1]);
    
    if (nrhs > 4) {
        num_solutions = (int)mxGetScalar(prhs[4]);
    }
    
    // Compute the number of edges
    int numEdges = 0;
    if (mxGetM(prhs[1]) == 3) {
        numEdges = mxGetN(prhs[1]);
        printf("Sparse edge representation \n");
    }
    else {
        printf("Full edge representation \n");
        for (int i = 0; i < numNodes; i++)
        {
            for (int j = 0; j < i; j++)
            {

                if (fabs(*(edgeCost + j*numNodes + i)) > 0)
                {
                    numEdges++;
                }
            }
        }
    }
    
    printf("Number of edges: %d \n", numEdges);
    
    

    Parametric::RegionId r;
	Parametric* p = new Parametric(/*maximum # of nodes*/ numNodes, /*maximum # of edges*/ numEdges); 
	Parametric::REAL lambdaMin = (*mxGetPr(prhs[2]));
	Parametric::REAL lambdaMax = (*mxGetPr(prhs[3]));

    printf("lambdaMin: %f \n", lambdaMin);
    printf("lambdaMax: %f \n", lambdaMax);
    
    p -> AddNode(numNodes); 

    
	for (int i = 0; i < numNodes; i++)
	{
        if (*(nodeCost + i*2 + 0) < 0) {
//            printf("Adding term : (%f * lambda + %f) * X_%d \n", *(nodeCost + i*2 + 0), *(nodeCost + i*2 + 1), i+1);
            delete p;
            mexErrMsgTxt("Negative lambda terms detected");
        }
            
        p -> AddUnaryTerm(i,  *(nodeCost + i*2 + 0), *(nodeCost + i*2 + 1)); // (nodeCost(1,i)*lambda+nodeCost(2,i))*X_i
//        printf("Adding term : (%f * lambda + %f) * X_%d \n", *(nodeCost + i*2 + 0), *(nodeCost + i*2 + 1), i+1);
	}

    
    if (mxGetM(prhs[1]) == 3) {
        for (int e = 0; e < numEdges; e++) {
            int i = (int)(*(edgeCost + e*3 + 0));
            int j = (int)(*(edgeCost + e*3 + 1));
            double cost = *(edgeCost + e*3 + 2);
            if (fabs(cost) > 0)
            {
                if (cost > 0) {
                    printf("Adding term: %f * X_%d * X_%d \n", cost, j+1, i+1);
                    delete p;
                    mexErrMsgTxt("Irregular terms detected");
                }
                
                p -> AddPairwiseTerm(/* NodeId */ i, /* NodeId */ j,     0, cost*(-0.5), cost*(-0.5), 0);  // X_i*X_j*edgeCost(i,j)
                p -> AddUnaryTerm(/* NodeId */ i,    0,  cost*0.5);
                p -> AddUnaryTerm(/* NodeId */ j,    0, cost*0.5);
                //                printf("Adding term: %f * X_%d * X_%d \n", *(edgeCost + j*numNodes + i), j+1, i+1);
            }
        }
    }
    else
    {
        for (int i = 0; i < numNodes; i++)
        {
            for (int j = 0; j < i; j++)
            {

                if (fabs(*(edgeCost + j*numNodes + i)) > 0)
                {
                    if (*(edgeCost + j*numNodes + i) > 0) {                    
                        delete p;
                        mexErrMsgTxt("Irregular terms detected");
                    }

                    p -> AddPairwiseTerm(/* NodeId */ i, /* NodeId */ j,     0, *(edgeCost + j*numNodes + i)*(-0.5), *(edgeCost + j*numNodes + i)*(-0.5), 0);  // X_i*X_j*edgeCost(i,j)
                    p -> AddUnaryTerm(/* NodeId */ i,    0, *(edgeCost + j*numNodes + i)*0.5); 
                    p -> AddUnaryTerm(/* NodeId */ j,    0, *(edgeCost + j*numNodes + i)*0.5); 
    //                printf("Adding term: %f * X_%d * X_%d \n", *(edgeCost + j*numNodes + i), j+1, i+1);
                }

            }
        }
    }

    printf("Solving parametric maxflow... \n");
    
    int N = p -> Solve(lambdaMin, lambdaMax);

    num_solutions = N < num_solutions? N : num_solutions;
    
    printf("Solved (%d breakpoints) ! \n", N);
    
	plhs[0] = mxCreateLogicalMatrix(numNodes, num_solutions+1);
    plhs[1] = mxCreateDoubleMatrix(N, 1, mxREAL);

	selectedNodes = mxGetLogicals(plhs[0]);
    lambdas = mxGetPr(plhs[1]);
	
	r = p->GetFirstRegion();
	if (p->GetRegionLambda(r) > lambdaMin)
	{
//		PrintCut(p, numNodes, 0);
		AssignCutToSolution(p, numNodes, 0, selectedNodes, num_solutions, N);
//		printf(" : optimal for lambda's in (-infinity; %f]\n", p->GetRegionLambda(r));
	}

	// Consider cuts C_0,...,C_{N-1}
	for (int k=1; k<N; k++, r=p->GetNextRegion(r))
	{
		lambdas[k-1] = p->GetRegionLambda(r);
//		PrintCut(p, numNodes, k);
		AssignCutToSolution(p, numNodes, k, selectedNodes, num_solutions, N);
//		printf(" : optimal for lambda's in [%f; %f]\n", p->GetRegionLambda(r), p->GetRegionLambda(p->GetNextRegion(r)));
	}

	// consider cut C_N
	r = p->GetLastRegion();
	lambdas[N-1] = p->GetRegionLambda(r);
	if (p->GetRegionLambda(r) < lambdaMax)
	{
//		PrintCut(p, numNodes, N);
		AssignCutToSolution(p, numNodes, N, selectedNodes, num_solutions, N);
//		printf(" : optimal for lambda's in [%f; +infinity)\n", p->GetRegionLambda(r));
	}
	delete p;

}
