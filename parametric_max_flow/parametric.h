// Copyright (C) 2007 Vladimir Kolmogorov <v.kolmogorov@cs.ucl.ac.uk>

#ifndef __PARAMETRIC_H__
#define __PARAMETRIC_H__

#include "maxflow//graph.h"
#include "maxflow//block.h"
#include "assert.h"

// Software for minimizing functions of the form
//
//   E(x_1,...,x_n ; lambda) = \sum_i E^{lambda}_i(x_i) + \sum_{(i,j)} E_{ij}(x_i,x_j)
// 
// for all values of lambda in the range [lambdaMin,LambdaMax]
//
// The unary terms E^{lambda}_i(x_i) depend on the global parameter lambda as follows:
//   E^{lambda}_i(x_i) = (A_i*lambda + B_i) * x_i
// 
// All coefficients A_i must be non-negative. 
// Pairwise terms must be submodular, i.e. E_{ij}(0,0)+E_{ij}(1,1)<=E_{ij}(0,1)+E_{ij}(1,0)
//
// For example usage, see example.cpp

class Parametric
{
private:
	struct Region;
public:
	typedef double REAL; // can be changed to double

	typedef Graph<REAL,REAL,REAL> GraphType;
	typedef int NodeId;
	typedef Region* RegionId;

	Parametric(int nodeNumMax, int edgeNumMax);
	~Parametric();

	// Adds node(s) to the graph. By default, one node is added (num=1); then first call returns 0, second call returns 1, and so on. 
	// If num>1, then several nodes are added, and node_id of the first one is returned.
	// Can be called at most nodeNumMax times.
	NodeId AddNode(int num = 1);

	// Adds term E^{lambda}_i(x_i) = (A*lambda + B) * x_i (where lambda is the global parameter).
	// A must be non-negative.
	// Can be called only after AddNode(). Can be called several times per node.
	void AddUnaryTerm(NodeId i, REAL A, REAL B);
	// Can be called at most edgeNumMax times.
	// The term must be submodular, i.e. E00+E11<=E01+E10.
	void AddPairwiseTerm(NodeId i, NodeId j, REAL E00, REAL E01, REAL E10, REAL E11);

	// Solves the problem for all lambda's in the range [lambdaMin,lambdaMax].
	// After the call the nodes are partitioned into N non-emptry ordered regions R_0, ..., R_{N-1}
	// (number N is returned by Solve()).
	// For each region k there is corresponding threshold lambda_k, which
	// is the value at which the region's label switch from 1 to 0 as lambda increases.
	// (Exception: if lambda_0==lambdaMin, then the lambda of region 0 is not
	// fixed yet, the interval [lambdaMin,lambdaMax] would have to be enlarged.
	// Similarly if lambda_{N-1}==lambdaMax).
	//
	// The sequence (lambdaMin, lambda_0, ..., lambda_{N-1}, lambdaMax) is non-decreasing.
	// 
	// There may also be two special regions R0 and R1 which never switch their label:
	// R0 has label 0 and R1 has label 1 for all lambda's in (-infinity,+infinity).
	// These regions are not counted by N. 
	// (Note that R0 and R1 may only contain nodes for which A>0).
	//
	// There holds 0<=N<nodeNum.
	//
	// (Partitions and thresholds can be accessed via functions GetRegion(),...,GetRegionCount() below).
	// -------------------------------------------------------------------------------------------------
	// Partitions define N+1 cuts C_0, ..., C_N as follows: In cut C_k
	//    Regions R1,  R_0, ..., R_{k-1} are labeled as 1
	//    Regions R0,  R_k, ..., R_{N-1} are labeled as 0
	// Explicit formula for C_k:
	//    x_i = 1   if   GetRegionCount(GetRegion(i)) >= k
	//        = 0   otherwise
	// -------------------------------------------------------------------------------------------------
	// Properties of cuts:
	//   0<k<N:  Cut C_k is optimal for the values of lambda in the range [lambda_{k-1},lambda_k].
	//   k=0:    If lambda_0     > lambdaMin then cut C_0 is optimal for (-\infty, lambda_0].
	//   k=N:    If lambda_{N-1} < lambdaMax then cut C_N is optimal for [lambda_{N-1},+\infty).
	// -------------------------------------------------------------------------------------------------
	// There must hold lambdaMin <= lambdaMax.
	// 
	// Solve() can be called multiple times with the following restriction:
	// the interval [lambdaMin,lambdaMax] must be the same or larger than that for the previous call.
	// Suppose, for example, that you want to compute the cut for a specific value of lambda,
	// check whether it's "good" or not, and if not compute the entire sequence of cuts.
	// Then you can start with interval [lambda,lambda] and then, if necessary,
	// iteratively increase the interval until there holds lambdaMin < lambda_0 <= lambda_{N-1} < lambdaMax.
	//
	// If you want modify the energy after calling Solve(), see ResetRegions().
	int Solve(REAL lambdaMin, REAL lambdaMax);


	/////////////////////////////////////////////////////////////////////////////////////////////////
	// Functions GetRegion(),...,GetRegionCount() can be called only after Solve()

	// Returns region to which node i belongs.
	RegionId GetRegion(NodeId i);

	// The list of regions can also be accessed sequentially using the following four functions.
	// Note: R0 and R1 are not included in this list.
	RegionId GetFirstRegion(); // returns region whose count is 0
	RegionId GetLastRegion(); // returns region whose count is N-1
	RegionId GetPrevRegion(RegionId r); // for region with count k, returns region k-1. Can be called for k=0, but the result is not specified.
	RegionId GetNextRegion(RegionId r); // similar

	// returns lambda corresponding to region r. 
	// Cannot be called for R0 or R1.
	REAL GetRegionLambda(RegionId r);
	// returns integer k such that -1 <= k <= N.
	// (Regions -1 and N are special: they correspond to R0 and R1, respectively).
	int GetRegionCount(RegionId r);

	/////////////////////////////////////////////////////////////////////////////////////////////////
	// If you want to modify energy after calling Solve(), you must call ResetRegions() first.
	// After that functions above (GetRegion(), ..., GetRegionCount()) cannot be called, until the next call to Solve().
	// Calling ResetRegions() removes the restriction on the interval for Solve().
	void ResetRegions();









	/////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////


private:
	struct Node;
	struct Arc;

	struct Node
	{
		Arc*		m_first;
		REAL		m_A, m_B;
		Node*		m_next; // next node from the same region
		Region*		m_region;
		GraphType::node_id	m_node;
	};

	struct Arc
	{
		Node*		m_to;
		Arc*		m_sister;
		Arc*		m_next;
		REAL		m_rcap;
	};

	struct Region
	{
		Node*		m_firstNode;
		Region*		m_prev;
		Region*		m_next;
		REAL		m_lambda;
		bool		m_isFixed;
		int			m_regionCount;
	};

	int		m_nodeNum;
	int		m_edgeNum;
	int		m_nodeNumMax;
	int		m_edgeNumMax;
	Node*	m_nodes;
	Arc*	m_arcs;
	REAL	m_lambdaMin, m_lambdaMax;

	Block<Region>*	m_regionBlocks;
	Region*	m_regionFirst;
	Region*	m_regionLast;
	int		m_regionNum;
	Region	m_region0;
	Region	m_region1;
	GraphType*	m_graph;

	void ProcessRegion(Region* r);
	void ReassignToR0();
	void ReassignToR1();
};



inline Parametric::NodeId Parametric::AddNode(int num)
{
	assert(!m_regionBlocks);
	assert(m_nodeNum + num <= m_nodeNumMax);
	int i = m_nodeNum;
	m_nodeNum += num;
	return i;
}

inline void Parametric::AddUnaryTerm(NodeId i, REAL A, REAL B)
{
	assert(!m_regionBlocks);
	assert(i>=0 && i<m_nodeNum);
	assert(A >= 0);
	m_nodes[i].m_A += A;
	m_nodes[i].m_B += B;
}

inline void Parametric::AddPairwiseTerm(NodeId from, NodeId to, REAL A, REAL B, REAL C, REAL D)
{
	assert(!m_regionBlocks);
	assert(from>=0 && from<m_nodeNum);
	assert(to>=0 && to<m_nodeNum);
	assert(m_edgeNum < m_edgeNumMax);

	Arc* a = m_arcs + 2*m_edgeNum;
	Arc* arev = a + 1;
	m_edgeNum ++;

	a->m_to = m_nodes+to;
	arev->m_to = m_nodes+from;

	a->m_sister = arev;
	arev->m_sister = a;

	a->m_next = m_nodes[from].m_first;
	m_nodes[from].m_first = a;

	arev->m_next = m_nodes[to].m_first;
	m_nodes[to].m_first = arev;

	// set capacities
	/* 
	   E = A A  +  0   B-A
	       D D     C-D 0
	   Add edges for the first term
	*/
	m_nodes[from].m_B += A - D;
	B -= A; C -= D;

	/* now need to represent
	   0 B
	   C 0
	*/

	assert(B + C >= 0); /* check submodularity */
	if (B < 0)
	{
		/* Write it as
		   B B  +  -B 0  +  0   0
		   0 0     -B 0     B+C 0
		*/
		m_nodes[from].m_B += B; /* first term */
		m_nodes[to].m_B -= B; /* second term */
		a->m_rcap = 0;
		arev->m_rcap = B+C;
	}
	else if (C < 0)
	{
		/* Write it as
		   -C -C  +  C 0  +  0 B+C
		    0  0     C 0     0 0
		*/
		m_nodes[from].m_B -= C; /* first term */
		m_nodes[to].m_B += C; /* second term */
		a->m_rcap = B+C;
		arev->m_rcap = 0;
	}
	else /* B >= 0, C >= 0 */
	{
		a->m_rcap = B;
		arev->m_rcap = C;
	}
}

inline Parametric::RegionId Parametric::GetRegion(NodeId i)
{
	assert(m_regionBlocks);
	return m_nodes[i].m_region;
}

inline Parametric::RegionId Parametric::GetFirstRegion()
{
	assert(m_regionBlocks);
	return m_regionFirst;
}

inline Parametric::RegionId Parametric::GetLastRegion()
{
	assert(m_regionBlocks);
	return m_regionLast;
}

inline Parametric::RegionId Parametric::GetPrevRegion(RegionId r)
{
	assert(m_regionBlocks);
	return r->m_prev;
}

inline Parametric::RegionId Parametric::GetNextRegion(RegionId r)
{
	assert(m_regionBlocks);
	return r->m_next;
}

inline Parametric::REAL Parametric::GetRegionLambda(RegionId r)
{
	assert(m_regionBlocks);
	return r->m_lambda;
}

inline int Parametric::GetRegionCount(RegionId r)
{
	assert(m_regionBlocks);
	return r->m_regionCount;
}

#endif
