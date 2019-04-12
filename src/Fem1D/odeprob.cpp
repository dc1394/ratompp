#include "stdafx.h"
#include "odeprob.h"


//
// Constructor
//
OdeProb::OdeProb(Bndr left, Bndr right, double gamma, std::shared_ptr<const util::Fun1D> const & g, std::shared_ptr<const util::Fun1D> const & f)
	: Prob(left, right, gamma, g, f)
{
	assert(g);
	assert(f);
}

//
// Defines the ODE problem
//
void OdeProb::Define(Bndr left, Bndr right, double gamma, std::shared_ptr<const util::Fun1D> const & g, std::shared_ptr<const util::Fun1D> const & f)
{
	DefineProb(left, right, gamma, g, f);
}

//
// Solves the problem (WITHOUT adaptation)
//
void OdeProb::Solve()
{
//	// Writes mesh on the screan
//	for(size_t m = 0; m < m_elt.size(); m++)
//		m_elt[m].Write(stdout);

	Malloc();
	Assemble();
	m_s->SolveSymPos(m_b, m_y);
	// m_y->Write("Ysol.dat");
}

//
// Adaptive solution of the problem
//
void OdeProb::SolveAdapt(double absMaxCoef)
{
EltInfo eltInfo;
std::vector<size_t> eltToSplit(1);

	while(true)
	{
//		printf("\n\nSTEP %d\n", step++);
//		printf("DIM = %d\n", Dim());

		Solve();
		MaxMinCoef(eltInfo);

		if(eltInfo.m_maxMinCoef < absMaxCoef)
			break;

		eltToSplit[0] = eltInfo.m_eltId;
		AddToMesh(eltToSplit);
	}
}

//
// Finds the element with the largest coefficient contained in the minimal coefficients
// BUBBLE functins are considered only.
//
void OdeProb::MaxMinCoef(EltInfo& eltInfo) const
{
double coef, minCoef;

	eltInfo.m_maxMinCoef = 0; // Initialization
	for(size_t n = 0; n < EltNo(); n++) // For each element
	{
		const Element& e = Elt(n);
		minCoef = DBL_MAX;

		for(size_t j = 1; j < e.DofNo() - 1; j++) // For each BUBBLE DOF in element
		{
			const int dof = e.m_dof[j];
			if(dof < 0) // Skip Dirichlet boundary conditions
				continue;

			// Find the smallest coefficient for element "e"
			coef = fabs(m_y->Get(dof));
			if(coef < minCoef)
				minCoef = coef;
		}

		// Set the largest coefficient with minimal coefficients
		if(minCoef > eltInfo.m_maxMinCoef)
		{
			eltInfo.m_eltId = n;
			eltInfo.m_maxMinCoef = minCoef;
			eltInfo.m_eigVal = 0;
		}
	}
}



//
// Allocates the required memory
//
void OdeProb::Malloc(void)
{
const size_t M = Dim();
const size_t band = GetBand();

	m_s.reset();
	m_s = std::make_unique<ClpMtxBand>(M, band, 0);
	m_s->Zero();

	m_b.reset();
	m_b = std::make_unique<Vec>(M);
	m_b->Zero();

	m_y.reset();
	m_y = std::make_unique<Vec>(M);
	m_y->Zero();
}


//
// Assembling algorithm for equation solving
//
void OdeProb::Assemble()
{
size_t i, j, psiI, psiJ;
int ni; // row position in matrix S
int nj; // column position in matrix S
// const size_t M = Dim();
const size_t N = EltNo(); // Number of elements

// const double bndr[3] = {0, m_left.m_val, m_right.m_val};

	// Element loop
	for(size_t n = 0; n < N; n++)
	{
		const Element& e = Elt(n);
		const size_t DofNo = e.DofNo();

		// Loop over basis functions
		for(i = 0; i < DofNo; i++)
		{
			ni = e.m_dof[i];
			if(ni < 0)
				continue;

			psiI = e.PsiId(i);

			// Loop over basis functions
			for(j = i; j < DofNo; j++)
			{
				psiJ = e.PsiId(j);

				nj = e.m_dof[j];
				if(nj > -1)
					m_s->Set(ni, nj) += CalcS(e, psiI, psiJ);
				//else // Dirichlet boundary conditions are ZERO, hence it can be skiped
				//	m_b->Set(ni) -= bndr[-nj] * CalcS(e, psiI, psiJ);
			}

			// Contribution of the vertex basis function $v_{m_1}$ to the right hand side $b$
			m_b->Set(ni) += CalcB(e, psiI);
		}
	}
}


/*
//!
//! Assembling algorithm for equation solving
//!
void OdeProb::Assemble()
{
size_t i, j;
int m1; // m1 - row position in matrix S
int m2; // m2 - column position in matrix S
size_t m;
const size_t N = Dim(m_left.m_type, m_right.m_type);
const size_t M = EltNo(); // Number of elements


	// Element loop
	for(m = 0; m < M; m++)
	{
		const Element& e = Elt(m);

		// Loop over vertex basis functions
		for(i = 0; i < 2; i++)
		{
			if(e.m_vertBnd[i] == BndType_Dir)
				continue;

			m1 = e.m_vertDof[i];
			// Loop over vertex basis functions
			for(j = 0; j < 2; j++)
			{
				m2 = e.m_vertDof[j];
				if(m2 > -1)
					m_s->Set(m1, m2) += CalcS(e, i, j);				
				else
					m_b->Set(m1) -= m_bndDir[-m2] * CalcS(e, i, j);
			}

			// Loop over bubble basis functions
			for(j = 0; j < e.m_p - 1; j++)
			{	
				m2 = e.m_bubbDof[j];
				assert(m2 > -1);
				m_s->Set(m1, m2) += CalcS(e, i, j + 2);
			}

			// Contribution of the vertex basis function $v_{m_1}$ to the right hand side $b$
			m_b->Set(m1) += CalcB(e, i);
		}


		// Loop over bubble basis functions
		for(i = 0; i < e.m_p - 1; i++)
		{
			m1 = e.m_bubbDof[i];
			assert(m1 > -1);

			// Loop over vertex basis functions
			for(j = 0; j < 2; j++)
			{
				m2 = e.m_vertDof[j];
				if(m2 > -1)
					m_s->Set(m1, m2) += CalcS(e, i + 2, j);
				else
					m_b->Set(m1) -= m_bndDir[-m2] * CalcS(e, i + 2, j);
			}

			// Loop over bubble basis functions
			for(j = 0; j < e.m_p - 1; j++)
			{
				m2 = e.m_bubbDof[j];
				assert(m2 > -1);
				m_s->Set(m1, m2) += CalcS(e, i + 2, j + 2);
			}

			// Contribution of the bubble basis function $v_{m_1}$ to the right hand side $b$
			m_b->Set(m1) += CalcB(e, i + 2);
		}
	}

	// Applying the Neumanna boundary conditions
	if(EltFront().m_vertBnd[0] == BndType_Neu)
		m_b->Set(0) -= m_neumann[0];

	if(EltBack().m_vertBnd[1] == BndType_Neu)
	{
		m = EltNo() - 1;
		m_b->Set(m) += m_neumann[1];
	}


	// m_s->Write("Smtx.dat");
	// m_b->Write("Bvec.dat");
}
*/


//
// Returns the value of the solution at $x$
//
double OdeProb::GetSol(double x) const
{
int m1;

	// The equation must be solved!
	assert(m_y);

	assert(IsInRange(x));
	const size_t n = FindElt(x);
	const Element& e = Elt(n);

	// s - Locat coordiante for element "e"
	const double s = std::min(std::max(e.Xinv(x), -1.0), 1.0); // MIN, MAX - To avoid the rounding errors

	// Sum over all basis function with support on the element $e$
	double val = 0;

/*
	// Left vertex basis function
	m1 = e.m_dof.front();
	if(m1 > -1)
		val += m_y->Get(m1) * Basis(0, s);
	else
	{
		// Apply the Dirichlet boundary conditions
		if(m_left.m_type == BndrType_Dir)
			val += m_left.m_val * Basis(0, s);
	}

	// Right vertex basis function
	m1 = e.m_dof.back();
	if(m1 > -1)
		val += m_y->Get(m1) * Basis(1, s);
	else
	{
		// Apply the Dirichlet boundary conditions
		if(m_right.m_type == BndrType_Dir)
			val += m_right.m_val * Basis(1, s);
	}

	// Bubble basis functions
	for(size_t j = 1; j < e.m_dof.size() - 1; j++)
	{
		m1 = e.m_dof[j];
		val += m_y->Get(m1) * Basis(j, s);
	}
*/

	// It works only with zero Dirichlet bpundary conditions
	for(size_t i = 0; i < e.m_dof.size(); i++)
	{
		const int m = e.m_dof[i];
		if(m < 0)
			continue;
		const size_t psiI = e.PsiId(i);

		val += m_y->Get(m) * Basis(psiI, s);
	}

	// Non-zero Dirichlet bpundary conditions are applied
	{
		// Left vertex basis function
		assert(m_left.m_type == BndrType_Dir);
		m1 = e.m_dof.front();
		if(m1 < 0)
		{
			const size_t psiI = e.PsiId(0);
			val += m_left.m_val * Basis(psiI, s);
		}

		// Right vertex basis function
		assert(m_right.m_type == BndrType_Dir);
		m1 = e.m_dof.back();
		if(m1 < 0)
		{
			const size_t psiI = e.PsiId(1);
			val += m_right.m_val * Basis(psiI, s);
		}
	}

	return val;
}

//
// Writes the solution the the file $path$.
// The solution is written for $pointNo$ points.
//
void OdeProb::WriteSol(const char* path, size_t pointNo) const
{
    auto const dx = (XBack() - XFront()) / (pointNo -1);

	auto out = std::unique_ptr<FILE, decltype(&std::fclose)>(std::fopen(path, "wt"), std::fclose);

	auto x = XFront();
	for (std::size_t i = 0UL; i < pointNo - 1; i++)
	{
		std::fprintf(out.get(), "%20lf\t%20lf\n", x, GetSol(x));
		x += dx;
	}

	// The last point must be written (to avoid the rounding errors)
	x = XBack();
	std::fprintf(out.get(), "%20lf\t%20lf\n", x, GetSol(x));
}
