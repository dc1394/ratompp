#include "stdafx.h"
#include "prob.h"

//
// Constructor
//
Prob::Prob() : m_gamma(0)
{
}

//
// Constructor
//
Prob::Prob(Bndr left, Bndr right, double gamma, std::shared_ptr<const util::Fun1D> const & g, std::shared_ptr<const util::Fun1D> const & f)
{
	DefineProb(left, right, gamma, g, f);
}

//
// Defines the problem
//
void Prob::DefineProb(Bndr left, Bndr right, double gamma, std::shared_ptr<const util::Fun1D> const & g, std::shared_ptr<const util::Fun1D> const & f)
{
	m_left = left;
	assert(left.m_type != BndrType_Emp);

	m_right = right;
	assert(right.m_type != BndrType_Emp);

	m_gamma = gamma;
	m_g = g;
	m_f = f;
}



//
// Returns the element $o[i][j]$ overlap matrix element
// The elements are read from precomputed array.
//
double Prob::CalcO(const Element& e, size_t i, size_t j) const
{
	return e.Jac() * Memi(i, j);
}

//
// Returns the element $b[i]$ of load matrix element.
// Gauss quadrature applied.
//
double Prob::CalcB(const Element& e, size_t i) const
{
double w, s, b = 0;

	assert(m_f);
	assert(m_xGauss.size() == m_wGauss.size());

	for(size_t n = 0; n < m_xGauss.size(); n++)
	{
		s = m_xGauss[n];
		w = m_wGauss[n];
		b += w * Basis(i, s) * m_f->Get(e.X(s));
	}
	return e.Jac() * b;
}

//
// Returns the element $e[i][j]$ stiffness matrix element.
// The elements are read from precomputed array.
//
double Prob::CalcS(const Element& e, size_t i, size_t j) const
{
const double v1 = m_gamma * Mesi(i, j);
double v0 = 0, s, w;

        if(m_g) // If function "g" is defined
	{
		for(size_t n = 0; n < m_xGauss.size(); n++)
		{
			s = m_xGauss[n];
			w = m_wGauss[n];
			v0 += w * Basis(i, s) * Basis(j, s) * m_g->Get(e.X(s));
		}
	}

	const double jac = e.Jac();
	return v1 / jac + v0 * jac;
}

//
// Creates linear mesh
//
void Prob::GenMeshLin(double a, double b, size_t nodeNo, size_t degree)
{
	Mesh::GenLin(a, b, nodeNo, degree);
	CreateCnnt(m_left.m_type, m_right.m_type);
}

//
// Creates linear-exponential mesh
//
void Prob::GenLinExp(double a, double b, int m, size_t degree)
{
	Mesh::GenLinExp(a, b, m, degree);
	CreateCnnt(m_left.m_type, m_right.m_type);
}

//
// Defines the mesh
//
void Prob::SetMesh(const std::vector<double>& x, const std::vector<size_t>& degree)
{
	Mesh::Set(x, degree);
	CreateCnnt(m_left.m_type, m_right.m_type);
}

//
// Add elements to the mesh
//
void Prob::AddToMesh(const std::vector<size_t>& eltToSplit)
{
	Mesh::AddToMesh(eltToSplit);
	CreateCnnt(m_left.m_type, m_right.m_type);
}

//
// Returns dimensinality of the problem
//
size_t Prob::Dim() const
{
	return Mesh::Dim(m_left.m_type, m_right.m_type);
}


