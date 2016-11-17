#include "stdafx.h"
#include "mesh.h"

//
// Constructor
//
Mesh::Mesh(void)
{
}

//
// Destructor
//
Mesh::~Mesh(void)
{
}


//
// Defines the mesh.
// x      - vertex coordinates (order ascending)
// degree - polynomial degrees
// "degree.size() == x.size() - 1"
//
void Mesh::Set(const std::vector<double>& x, const std::vector<size_t>& degree)
{
	assert(degree.size() == x.size() - 1);
	assert(x.size() >= 2);

	// copies the vector
	m_x = x;

	const size_t N = degree.size();
	m_elt.resize(N);

	for(size_t n = 0; n < N; n++)
		m_elt[n].Set(x[n], x[n + 1], degree[n]);
}

//
// Genarates the mesh on the interval [a, b].
// Mesh has $nodeNo$ nodes (it means $nodeNo-1$ elements).
// Each element has $degree$
//
void Mesh::GenLin(double a, double b, size_t nodeNo, size_t degree)
{
const double dx = (b - a) / (nodeNo - 1);
std::vector<double> x(nodeNo);
std::vector<size_t> deg(nodeNo - 1, degree);

	assert(b > a);
	assert(nodeNo >= 2);

	for(size_t i = 0; i < nodeNo - 1; i++)
		x[i] = a + i * dx;

	// To avoid the rounding errors
	x.back() = b; 

	Set(x, deg);
}

//
// Generates the mesh on the interval [a, b]
//
void Mesh::GenLinExp(double a, double b, int m, size_t degree)
{
const size_t linNodeNo = (1 << m);
const double h = 1. / (double)linNodeNo;
std::vector<double> x;
double r;

	assert(b > a);

	r = a;
	for(size_t i = 0; i < linNodeNo; i++)
	{
		x.push_back(r);
		r += h;
	}

	while(r < b)
	{
		x.push_back(r);
		r *= (1 + h); 
	}

        // The last node must be on the end of interval
	x.push_back(b); 

	std::vector<size_t> deg(x.size() - 1, degree);

	Set(x, deg);

}

/*
//!
//! Create connectivity array
//!
void Mesh::CreateCnnt()
{
int count = 0;
size_t m;

	assert(!m_elt.empty());

	// Visiting vertex basis functions on the first element
	if(m_elt[0].m_vertBnd[0] == BndType_Dir)
		m_elt[0].m_vertDof[0] = -1;
	else
	{
		m_elt[0].m_vertDof[0] = count;
		count++;
	}
	m_elt[0].m_vertDof[1] = count;

	// Visiting vertex basis functions on interior elements 
	for(m = 1; m < m_elt.size() - 1; m++)
	{
		m_elt[m].m_vertDof[0] = count;
		count++;
		m_elt[m].m_vertDof[1] = count;
	}

	// Visiting vertex basis functions on the last element
	m = m_elt.size() - 1;

	m_elt[m].m_vertDof[0] = count;
	count++;
	if(m_elt[m].m_vertBnd[1] == BndType_Dir)
		m_elt[m].m_vertDof[1] = -2;
	else
	{
		m_elt[m].m_vertDof[1] = count;
		count++;
	}

	// Visiting buble basis functions on all elements
	for(m = 0; m < m_elt.size(); m++)
	{
		for(size_t j = 0; j < m_elt[m].m_p - 1; j++)
		{
			m_elt[m].m_bubbDof[j] = count;
			count++;
		}
	}	
}
*/


//
// Create connectivity array
//
void Mesh::CreateCnnt(BndrType left, BndrType right)
{
const size_t N = m_elt.size();
int idx;
size_t n, j;

	assert(!m_elt.empty());

	// Left end
	if(left  == BndrType_Dir)
		idx = -1;
	else
		idx = 0;

	for(n = 0; n < N; n++)
	{
		Element& e = m_elt[n];

		for(j = 0; j < e.m_dof.size(); j++)
			e.m_dof[j] = idx++;

		// The last basis function of the last element must be the first basis function of the next element.
		idx--;
	}	

	// Right end
	if(right == BndrType_Dir)
		m_elt.back().m_dof.back() = -2;
}

//
// Calculate the dimension of the finite element space
//
size_t Mesh::Dim(BndrType left, BndrType right) const
{
size_t M = 0;

	assert(m_elt.size() > 0);
	for(size_t n = 0; n < m_elt.size(); n++)
		M += m_elt[n].P();

	if(left == BndrType_Neu)
		M++;

	if(right == BndrType_Neu)
		M++;

	return (M - 1);
}

//
// Returns the bandwith of band matrix
//
size_t Mesh::GetBand(void) const
{
size_t pMax = 1;
	
	for(size_t n = 0; n < m_elt.size(); n++)
	{
		if(pMax < m_elt[n].P())
			pMax = m_elt[n].P();
	}
	return pMax;
}

//
// Returns "true" if "x" belongs to interval
//
bool Mesh::IsInRange(double x) const
{
	const bool b1 = (x >= m_x.front());
	const bool b2 = (x <= m_x.back());

	return b1 && b2;
}

//
// Searching the interval containg the value $x$: $x_{m} <= x <= x_{m+1}$
// returns index $m$
//
size_t Mesh::FindElt(double x) const
{
size_t n;
	for(n = 0; n < m_x.size() - 1; n++)
	{
		if(m_x[n] <= x && x <= m_x[n + 1])
			break;
	}
	return n;
}

//
// Adds element to the mesh
//
void Mesh::AddToMesh(const std::vector<size_t>& eltToSplit)
{
size_t i, n;
std::vector<double> newX(m_x);
std::vector<size_t> newDegree(m_x.size() - 1 + eltToSplit.size(), m_elt[0].P()); // All elements has the same degree

	for(i = 0; i < eltToSplit.size(); i++)
	{
		n = eltToSplit[i];
		double tmp = (X(n) + X(n + 1)) / 2;
		newX.push_back(tmp);
	}
	std::sort(newX.begin(), newX.end());

	Set(newX, newDegree);
}
