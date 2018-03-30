// Copyright (C) 2007 Vladimir Kolmogorov <v.kolmogorov@cs.ucl.ac.uk>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parametric.h"

Parametric::Parametric(int nodeNumMax, int edgeNumMax)
	: m_nodeNumMax(nodeNumMax),
	  m_nodeNum(0),
	  m_edgeNumMax(edgeNumMax),
	  m_edgeNum(0),
	  m_regionFirst(NULL),
	  m_regionBlocks(NULL),
	  m_regionNum(0),
	  m_graph(NULL)
{
	m_nodes = new Node[m_nodeNumMax];
	memset(m_nodes, 0, m_nodeNumMax*sizeof(Node));
	m_arcs = new Arc[2*m_edgeNumMax];

	m_region0.m_prev = m_region0.m_next = NULL;
	m_region0.m_isFixed = true;
	m_region1.m_prev = m_region1.m_next = NULL;
	m_region1.m_isFixed = true;
}

Parametric::~Parametric()
{
	delete [] m_nodes;
	delete [] m_arcs;
	if (m_regionBlocks) delete m_regionBlocks;
	if (m_graph) delete m_graph;
}

void Parametric::ResetRegions()
{
	if (m_regionBlocks)
	{
		delete m_regionBlocks;
		m_regionBlocks = NULL;
		delete m_graph;
		m_graph = NULL;
	}
}

void Parametric::ReassignToR0()
{
	Node* i;
	Region* r;

	while ( (r = m_regionFirst) )
	{
		assert(r->m_isFixed);
		if (r -> m_lambda > m_lambdaMin) break;

		for (i=r->m_firstNode; i; i=i->m_next)
		{
			if (i->m_A > 0) break;
		}
		if (i) break;

		// reassign
		m_regionFirst = r->m_next;
		if (m_regionFirst) m_regionFirst->m_prev = NULL;
		else m_regionLast = NULL;
		m_regionNum --;
		for (i=r->m_firstNode; i; i=i->m_next)
		{
			i->m_region = &m_region0;
		}
	}
}

void Parametric::ReassignToR1()
{
	Node* i;
	Region* r;

	while ( (r = m_regionLast) )
	{
		assert(r->m_isFixed);
		if (r -> m_lambda < m_lambdaMax) break;

		for (i=r->m_firstNode; i; i=i->m_next)
		{
			if (i->m_A > 0) break;
		}
		if (i) break;

		// reassign
		m_regionLast = r->m_prev;
		if (m_regionLast) m_regionLast->m_next = NULL;
		else m_regionFirst = NULL;
		m_regionNum --;
		for (i=r->m_firstNode; i; i=i->m_next)
		{
			i->m_region = &m_region1;
		}
	}
}

void Parametric::ProcessRegion(Region* r)
{
	Node* i;
	Node* j;
	Arc* a;
	Region* r_new[2];
	REAL lambda, lambdaMin, lambdaMax;
	REAL A, B;
	int s;

	// choose threshold
	A = B = 0;
	for (i=r->m_firstNode; i; i=i->m_next)
	{
		assert(i->m_region == r);
		A += i->m_A;
		B += i->m_B;
	}
	assert(A >= 0);
	lambdaMin = r->m_lambda;
	lambdaMax = (r->m_next) ? r->m_next->m_lambda : m_lambdaMax;
	if (A <= 0)
	{
		lambda = (B >= 0) ? lambdaMin : lambdaMax;
	}
	else
	{
		lambda = - B / A;
		if (lambda < lambdaMin) lambda = lambdaMin;
		if (lambda > lambdaMax) lambda = lambdaMax;
	}

	// build m_graph
	for (i=r->m_firstNode; i; i=i->m_next)
	{
		assert(i->m_region == r);
		i->m_node = m_graph->add_node();
		m_graph->add_tweights(i->m_node, i->m_A*lambda + i->m_B, 0);

		for (a=i->m_first; a; a=a->m_next)
		{
			j = a->m_to;
			if (j > i || j->m_region != r) continue;

			m_graph->add_edge(i->m_node, j->m_node, a->m_rcap, a->m_sister->m_rcap);
		}
	}

	// compute maxflow
	m_graph->maxflow();

	// update flow
	GraphType::arc_id ga = m_graph->get_first_arc();
	for (i=r->m_firstNode; i; i=i->m_next)
	{
		// old slope: i->m_A*lambda + i->m_B
		// new slope: ((Graph::node*)(i->m_node))->tr_cap = i->m_A*lambda + Bnew
		i->m_B = m_graph->get_trcap(i->m_node) - i->m_A*lambda;

		for (a=i->m_first; a; a=a->m_next)
		{
			j = a->m_to;
			if (j > i || j->m_region != r) continue;

			a->m_rcap = m_graph->get_rcap(ga);
			ga = m_graph->get_next_arc(ga);
			a->m_sister->m_rcap = m_graph->get_rcap(ga);
			ga = m_graph->get_next_arc(ga);
		}
	}

	// find out whether the region should be split or not
	s = m_graph->what_segment(r->m_firstNode->m_node);
	for (i=r->m_firstNode; i; i=i->m_next)
	{
		if (s != m_graph->what_segment(i->m_node)) break;
	}

	if (i == NULL)
	{
		r->m_isFixed = true;
		r->m_lambda = lambda;
	}
	else
	{
		// split r
		r_new[0] = r;
		r_new[1] = m_regionBlocks->New();
		r_new[1]->m_prev = r_new[0];
		r_new[1]->m_next = r_new[0]->m_next;
		if (r_new[0]->m_next) r_new[0]->m_next->m_prev = r_new[1];
		r_new[0]->m_next = r_new[1];
		r_new[1]->m_lambda = lambda;
		r_new[1]->m_isFixed = false;
		m_regionNum ++;
		if (r == m_regionLast) m_regionLast = r_new[1];

		Node* i_prev[2] = { NULL, NULL };
		for (i=r->m_firstNode; i; i=i->m_next)
		{
			s = m_graph->what_segment(i->m_node);
			i->m_region = r_new[s];
			if (i_prev[s]) i_prev[s]->m_next     = i;
			else           r_new[s]->m_firstNode = i;
			i_prev[s] = i;
		}
		assert(i_prev[0] && i_prev[1]);
		i_prev[0]->m_next = i_prev[1]->m_next = NULL;
	}

	m_graph->reset();
}

int Parametric::Solve(REAL lambdaMin, REAL lambdaMax)
{
	Node* i;
	Node* currentNode;
	Region* r;

	assert(lambdaMin <= lambdaMax);

	if (m_regionBlocks == NULL)
	{
		m_regionBlocks = new Block<Region>(1024);
		m_regionFirst = m_regionLast = m_regionBlocks->New();
		m_regionFirst->m_prev = m_regionFirst->m_next = NULL;
		m_regionFirst->m_lambda = lambdaMin;
		m_regionFirst->m_isFixed = false;
		m_regionFirst->m_firstNode = m_nodes;
		m_regionNum = 1;

		for (i=m_nodes; i<m_nodes+m_nodeNum; i++)
		{
			i->m_region = m_regionFirst;
			i->m_next = i+1;
		}
		(i-1)->m_next = NULL;

		m_region0.m_firstNode = m_region1.m_firstNode = NULL;

		m_graph = new GraphType(m_nodeNum, m_edgeNum);
	}
	else
	{
		assert(lambdaMin <= m_lambdaMin && lambdaMax >= m_lambdaMax);
		if (!m_regionFirst)
		{
			// only R0 and R1 left
			m_lambdaMin = lambdaMin;
			m_lambdaMax = lambdaMax;
			return 0; 
		}
		if (lambdaMin < m_lambdaMin && m_regionFirst->m_lambda == m_lambdaMin)
		{
			m_regionFirst->m_lambda = lambdaMin;
			m_regionFirst->m_isFixed = false;
		}
		if (lambdaMax > m_lambdaMax && m_regionLast->m_lambda == m_lambdaMax)
		{
			m_regionLast->m_isFixed = false;
		}
	}
	m_lambdaMin = lambdaMin;
	m_lambdaMax = lambdaMax;


	// main loop
	currentNode = m_nodes;
	while (currentNode < m_nodes + m_nodeNum)
	{
		// process region starting at currentNode
		r = currentNode->m_region;
		if (r->m_isFixed) { currentNode ++; continue; }
		assert(currentNode == r->m_firstNode);

		ProcessRegion(r);
	}

	ReassignToR0();
	ReassignToR1();

	// set correct counts
	int k = 0;
	for (r=m_regionFirst; r->m_next; r=r->m_next)
	{
		r->m_regionCount = k ++;
	}
	r->m_regionCount = k++;
	assert(k == m_regionNum);
	m_region0.m_regionCount = -1;
	m_region1.m_regionCount = m_regionNum;

	return m_regionNum;
}
