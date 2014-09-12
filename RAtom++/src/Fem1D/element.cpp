#include "stdafx.h"
#include "element.h"


//
// Constructor
//
Element::Element(void)
{
	m_c1 = 0;
	m_c2 = 0;
}

//
// Destructor
//
Element::~Element(void)
{
}

//
// Write the element to file
//
void Element::Write(FILE* out) const
{
//	fprintf(out, "p = %ld, vertDir={%ld, %ld}, vertDof={%ld, %ld}, bubbDof={", 
//		m_p, m_vertBnd[0], m_vertBnd[1], m_vertDof[0], m_vertDof[1]);
//	for(size_t i = 0; i < m_p - 1; i++)
//		fprintf(out, "%ld, ", m_bubbDof[i]);
//	fprintf(out, "}\n\n");
}

/*
//!
//! Returns "true", if the one of element vercies has the Neumann boundary condition.
//!
bool Element1D::IsNeumann() const
{
const bool b0 = (m_vertBnd[0] == BndType_Neu);
const bool b1 = (m_vertBnd[1] == BndType_Neu);

	return (b0 || b1);
}
*/

//
// Returns ID of referenced Lobatto basis function.
// Returns function cooperates with function Mesh::CreateCnnt()
//
size_t Element::PsiId(size_t i) const
{
	assert(i < m_dof.size());

	if(i == 0)
		return 0;

	if(i == m_dof.size() - 1)
		return 1;

	return i + 1;
}

//
// Sets (defines) the element woth ends [x0, x1] and degree "p".
//
void Element::Set(double x0, double x1, size_t p)
{
	assert(x1 > x0);

	// Calculates the transformation coefficients form interval $[x_m, x_{m+1}]$ to reference interval $[-1, 1]$.
	m_c1 = (x1 + x0) / 2;
	m_c2 = (x1 - x0) / 2;

	m_dof.resize(p + 1);
}
