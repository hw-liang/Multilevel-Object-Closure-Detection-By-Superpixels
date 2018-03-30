#include <stdio.h>
#include "parametric.h"



// Example2: Minimize energies
// E(x,y ; lambda) = (lambda - 1) x + (lambda + 2) y - xy
// for all values of lambda in the range [-10;10]


// print cut C_k (0<=k<=N where N is the number returned by Solve()).
void PrintCut(Parametric* p, int nodeNum, int k)
{
	int i;

	printf("Cut (");
	for (i=0; i<nodeNum; i++)
	{
		int x = (p->GetRegionCount(p->GetRegion(i)) >= k) ? 1 : 0; // label of node i in C_k
		printf(" %d ", x);
	}
	printf(")");
}



int main()
{
	int k, N, nodeNum = 2;
	Parametric::RegionId r;
	Parametric* p = new Parametric(/*maximum # of nodes*/ 2, /*maximum # of edges*/ 1); 
	Parametric::REAL lambdaMin = -10, lambdaMax = 10;

	p -> AddNode(2); // add 3 nodes; their id's are 0,1,2

	p -> AddUnaryTerm(/* NodeId */ 0,    1, -1); // add term   (lambda - 1) x
	p -> AddUnaryTerm(/* NodeId */ 1,    1, 2);  // add term   (lambda + 2) y

	//p -> AddPairwiseTerm(/* NodeId */ 0, /* NodeId */ 1,     0, 0, 0, -0.5);  // add term   - x y
	p -> AddPairwiseTerm(/* NodeId */ 1, /* NodeId */ 0,     0, 0, 0, -1);  // add term   - x y

	N = p -> Solve(lambdaMin, lambdaMax);



	// consider cut C_0
	r = p->GetFirstRegion();
	if (p->GetRegionLambda(r) > lambdaMin)
	{
		PrintCut(p, nodeNum, 0);
		printf(" : optimal for lambda's in (-infinity; %f]\n", p->GetRegionLambda(r));
	}

	// Consider cuts C_0,...,C_{N-1}
	for (k=1; k<N; k++, r=p->GetNextRegion(r))
	{
		PrintCut(p, nodeNum, k);
		printf(" : optimal for lambda's in [%f; %f]\n", p->GetRegionLambda(r), p->GetRegionLambda(p->GetNextRegion(r)));
	}

	// consider cut C_N
	r = p->GetLastRegion();
	if (p->GetRegionLambda(r) < lambdaMax)
	{
		PrintCut(p, nodeNum, N);
		printf(" : optimal for lambda's in [%f; +infinity)\n", p->GetRegionLambda(r));
	}



	delete p;

	return 0;
}
