#include "stdafx.h"
#include "eigprob.h"
#include <memory>           // for std::unique_ptr    
#include <boost/format.hpp> // for boost::format

//
// Defines the problem
//
void EigProb::Define(double gamma, std::shared_ptr<util::Fun1D> const & g)
{
	DefineProb(Bndr(BndrType_Dir, 0), Bndr(BndrType_Dir, 0), gamma, g, nullptr);
}


//
// Solves the eigenproblem, WITHOUT adaptive procedure
//
void EigProb::Solve(size_t eigNo, double abstol)
{
	Malloc();
	Assemble();
	m_s->EigenGen(eigNo, abstol, *m_w, *m_z, *m_o);
//	for(size_t i = 0; i < eigNo; i++)
//	{
//		double v = GetEigVal(i);
//		printf("%ld %.16lf\n", i, v);
//	}
}

//
// Solve the eigenproblem adatively
//
void EigProb::SolveAdapt(size_t eigNo, double abstol, double absMaxCoef)
{
std::vector<EltInfo> eltInfo(eigNo);
std::vector<EltInfo>::iterator ii, newEnd;
std::vector<size_t> eltToSplit;
// size_t n,  step = 0;
size_t i;
double maxCoef;


	while(true)
	{
//		printf("\n\nSTEP %d\n", step++);

		Solve(eigNo, abstol);
		MaxMinCoef(eltInfo);
		std::sort(eltInfo.begin(), eltInfo.end());
		newEnd = std::unique(eltInfo.begin(), eltInfo.end());

		maxCoef = 0;
		for(i = 0; i < eltInfo.size(); ++i)
		{
			if(eltInfo[i].m_maxMinCoef > maxCoef)
				maxCoef = eltInfo[i].m_maxMinCoef;
		}
		if(maxCoef < absMaxCoef)
			break;
		
# ifdef _DEBUG
		printf("Elt splitted. DIM = %d\n", Dim());
		printf("MAX COEF = %15.6E\n", maxCoef);
# endif
		
		eltToSplit.clear();
		for(ii = eltInfo.begin(); ii != newEnd; ++ii) 
		{
# ifdef _DEBUG
			printf("Eig = %d, EltId = %d, maxMinCoef = %15.6E\n", (*ii).m_eigVal, (*ii).m_eltId, (*ii).m_maxMinCoef);
# endif
			eltToSplit.push_back((*ii).m_eltId);
		}
		AddToMesh(eltToSplit);
	}
	

/*
	eltInfo.erase(newEnd, eltInfo.end());
	for(i = 0; i < eltInfo.size(); ++i) 
		printf("Eig = %d, EltId = %d, maxMinCoef = %15.6E\n", eltInfo[i].m_eigVal, eltInfo[i].m_eltId, eltInfo[i].m_maxMinCoef);

	for(i = 0; i < eltInfo.size(); i++)
	{
		n = eltInfo[i].m_eltId;
		double tmp = (X(n - 1) + X(n)) / 2;
		m_x.push_back(tmp);
	}
	std::sort(m_x.begin(), m_x.end());
*/
}


//
// Allocates the required memory
//
void EigProb::Malloc(void)
{
const size_t M = Dim();
const size_t band = GetBand();

	m_s.reset();
	m_s = std::make_shared<ClpMtxBand>(M, band, 0);
	m_s->Zero();

	m_w.reset();
	m_w = std::make_shared<Vec>(M);

	m_z.reset();
	m_z = std::make_shared<ClpMtx>(M, M);

	m_o.reset();
	m_o = std::make_shared<ClpMtxBand>(M, band, 0);
}

//
// Assembling algorithm for eigenvalue problem
//
void EigProb::Assemble()
{
size_t i, j;
int ni; // row position in matrix S and O
int nj; // column position in matrix S and O
// const size_t M = Dim();
const size_t N = EltNo(); // Number of elements

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

			const size_t psiI = e.PsiId(i);

			// Loop over basis functions
			for(j = i; j < DofNo; j++)
			{
				const size_t psiJ = e.PsiId(j);

				nj = e.m_dof[j];
				if(nj > -1)
				{
					m_s->Set(ni, nj) += CalcS(e, psiI, psiJ);
					m_o->Set(ni, nj) += CalcO(e, psiI, psiJ);
				}
			}
		}
	}

//	m_s->Write("Smtx.dat");
//	m_o->Write("Omtx.dat");
}

/*
//!
//! Assembling algorithm for eigenvalue problem
//!
void EigProb::Assemble()
{
size_t i, j;
int m1; // m1 - row position in matrix S and O
int m2; // m2 - column position in matrix S and O
size_t m;
const size_t N = Dim();
const size_t M = EltNo(); // Number of elements

	// Element loop
	for(m = 0; m < M; m++)
	{
		const Element1D& e = Elt(m);

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
				{
					m_s->Set(m1, m2) += CalcS(e, i, j);
					m_o->Set(m1, m2) += CalcO(e, i, j);
				}
			}

			// Loop over bubble basis functions
			for(j = 0; j < e.m_p - 1; j++)
			{	
				m2 = e.m_bubbDof[j];
				assert(m2 > -1);
				m_s->Set(m1, m2) += CalcS(e, i, j + 2);
				m_o->Set(m1, m2) += CalcO(e, i, j + 2);
			}
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
				{
					m_s->Set(m1, m2) += CalcS(e, i + 2, j);
					m_o->Set(m1, m2) += CalcO(e, i + 2, j);
				}
			}

			// Loop over bubble basis functions
			for(j = 0; j < e.m_p - 1; j++)
			{
				m2 = e.m_bubbDof[j];
				assert(m2 > -1);
				m_s->Set(m1, m2) += CalcS(e, i + 2, j + 2);
				m_o->Set(m1, m2) += CalcO(e, i + 2, j + 2);
			}
		}
	}

	// Tylko waunki brzegowe Dirichleta sa zaimplementowane
	assert(EltFront().m_vertBnd[0] == BndType_Dir);
	assert(EltBack().m_vertBnd[1] == BndType_Dir);


	// m_s->Write("Smtx.dat");
	// m_o->Write("Omtx.dat");
}
*/

//
// Returns the value of $eig$ eigenfunction at point $x$
//
double EigProb::GetEigFun(size_t eig, double x) const
{
// int m1;

	assert(eig < m_w->size());
	assert(IsInRange(x));

	const size_t n = FindElt(x);
	const Element& e = Elt(n);

	// s - local variable for element "e"
	const double s = std::min(std::max(e.Xinv(x), -1.0), 1.0); // MIN, MAX - To avoid the rounding errors

	// Sum over all basis function with support on the element $e$
	double val = 0;

	for(size_t i = 0; i < e.m_dof.size(); i++)
	{
		const int m = e.m_dof[i];
		if(m < 0)
			continue;
		const size_t psiI = e.PsiId(i);

		val += m_z->Get(m, eig) * Basis(psiI, s);
	}
	return val;
}

//
// Returns the value of $eig$ first derivative of eigenfunction at point $x$
//
double EigProb::GetEigFunD1(size_t eig, double x) const
{
// int m1;

	assert(eig < m_w->size());
	assert(IsInRange(x));

	const size_t n = FindElt(x);
	const Element& e = Elt(n);

	// s -  local variable for element "e"
	const double s = std::min(std::max(e.Xinv(x), -1.0), 1.0); // MIN, MAX - To avoid the rounding errors

	// Sum over all basis function with support on the element $e$
	double val = 0;

	for(size_t i = 0; i < e.m_dof.size(); i++)
	{
		const int m = e.m_dof[i];
		if(m < 0)
			continue;
		const size_t psiI = e.PsiId(i);

		val += m_z->Get(m, eig) * BasisD1(psiI, s);
	}
	return val;
}


//
// Returns the value of $eig$ second derivative of eigenfunction at point $x$
//
double EigProb::GetEigFunD2(size_t eig, double x) const
{
// int m1;

	assert(eig < m_w->size());
	assert(IsInRange(x));

	const size_t n = FindElt(x);
	const Element& e = Elt(n);

	// s -  local variable for element "e"
	const double s = std::min(std::max(e.Xinv(x), -1.0), 1.0); // MIN, MAX - To avoid the rounding errors

	// Sum over all basis function with support on the element $e$
	double val = 0;

	for(size_t i = 0; i < e.m_dof.size(); i++)
	{
		const int m = e.m_dof[i];
		if(m < 0)
			continue;
		const size_t psiI = e.PsiId(i);

		val += m_z->Get(m, eig) * BasisD2(psiI, s);
	}
	return val;
}



//
// Returns the value of the $eig$ eigenvalue
//
double EigProb::GetEigVal(size_t eig) const
{
	if(eig > m_w->size())
		return 0;

	return m_w->Get(eig);
}



//
// Writes eigen-function $eig$ to file
// If "pointNo == 0", then eigenfunction is stored in mesh nodes only.
// Argument "pointNo" determines number of addtional points netween mesh nodes
// where the eigenfunction is stored.
//
void EigProb::WriteEigFun(const char* path, size_t eig, size_t pointNo) const
{
    double x;

	auto out = std::unique_ptr<FILE, decltype(&std::fclose)>(std::fopen(path, "wt"), std::fclose);
	if (!out)
	{
        throw std::invalid_argument((boost::format("Cannot open file '%s' for write.") % path).str());
	}


	for (std::size_t n = 0UL; n < XNo() - 1; n++)
	{
		x = X(n);
		auto const dx = (X(n + 1) - X(n)) / (pointNo + 1);
		for (std::size_t i = 0UL; i < pointNo + 1; i++)
		{
			std::fprintf(out.get(), "%20lf\t%20lf\n", x, GetEigFun(eig, x));
			x += dx;
		}
	}

	// The last point must be written (to avoid the rounding errors)
	x = XBack();
	std::fprintf(out.get(), "%20lf\t%20lf\n", x, GetEigFun(eig, x));
}


//
// Writes eigen-function $eig$ to file
// If "pointNo == 0", then eigenfunction is stored in mesh nodes only.
// Argument "pointNo" determines number of addtional points netween mesh nodes
// where the eigenfunction is stored.
// ll - angular quantum number, used for extralopation at r==0 only
//
void EigProb::WriteEigFunEx(const char* path, size_t ll,  size_t eig, size_t pointNo) const
{
	auto out = std::unique_ptr<FILE, decltype(&std::fclose)>(std::fopen(path, "wt"), std::fclose);
	if(!out)
	{
        throw std::invalid_argument((boost::format("Cannot open file '%s' for write.") % path).str());
	}

	std::fprintf(out.get(), "# %14s %18s %18s %18s %18s\n", "r = radius", "R(r)", "R(r)/r", "R'(r)", "R''(r)");
    
    double radr;
	for (std::size_t n = 0; n < XNo() - 1; n++)
	{
		auto x = X(n);
		auto const dx = (X(n + 1) - X(n)) / (pointNo + 1);
		for (std::size_t i = 0; i < pointNo + 1; i++)
		{
			auto const rad = GetEigFun(eig, x);
			auto const radD1 = GetEigFunD1(eig, x);
			auto const radD2 = GetEigFunD2(eig, x);

			if(x == 0) // Extrapolacja liniowa w x == 0
			{
				const double h = 1E-4;
				if(ll == 0)
					radr = 2 * GetEigFun(eig, h) / h - GetEigFun(eig, 2 * h) / (2 * h);
				else  // if ll>0, then eigenfunctions are 0 for r=0
					radr = 0;
			}
			else
			{
				radr = rad / x;
			}

			std::fprintf(out.get(), "%18.8E %18.8E %18.8E %18.8E %18.8E\n", x, rad, radr, radD1, radD2);
			x += dx;
		}
	}

	// The last point must be written (to avoid the rounding errors)
	auto const x = XBack();
	std::fprintf(out.get(), "%18.8E %18.8E %18.8E %18.8E %18.8E\n", x, GetEigFun(eig, x), GetEigFun(eig, x) / x, GetEigFunD1(eig, x), GetEigFunD2(eig, x));
}


//
// Writes memmbers of eginevectors into file
//
void EigProb::WriteEigCoef(const char* path, size_t eigNo) const
{
    auto out = std::unique_ptr<FILE, decltype(&std::fclose)>(std::fopen(path, "wt"), std::fclose);

    std::fprintf(out.get(), "EIGENVALUES:\n");
    for (std::size_t i = 0UL; i < eigNo; i++)
    {
       std::fprintf(out.get(), "Lambda[%2lu] = %20.11E\n", static_cast<unsigned long>(i), (*m_w)[i]);
    }
	//////////////////////////////////////////////////////////////////////////////////////////
	std::fprintf(out.get(), "\n\nEIGENVECTOR COEFFICIENTS\n");	
	std::fprintf(out.get(), "   DOF"); 
	for (std::size_t i = 0UL; i < eigNo; i++)
	{
	    std::fprintf(out.get(), "     Eig[%2lu]      ", static_cast<unsigned long>(i));
	}
    
    std::fprintf(out.get(), "\n");
	std::fprintf(out.get(), "======");
	
    for (std::size_t i = 0UL; i < eigNo; i++)
    {
        std::fprintf(out.get(), "==================");
    }
    std::fprintf(out.get(), "\n");


	std::size_t row = 0;
	for (std::size_t n = 0UL; n < EltNo(); n++) // For each element
	{
		auto const & e = Elt(n);

		std::fprintf(out.get(), "Elt:%06lu = [ %14.6E, %14.6E ]\n", static_cast<unsigned long>(n), e.X(-1), e.X(1));
		
		for(size_t j = 0; j < e.DofNo(); j++) // For each DOF in elemt
		{
			auto const dof = e.m_dof[j];
			if(dof < 0)
				continue;
			std::fprintf(out.get(), "%6d", dof);
            for (std::size_t i = 0UL; i < eigNo; i++)
            {
                std::fprintf(out.get(), "%18.6E", m_z->Get(dof, i));
            }
			
		    std::fprintf(out.get(), "\n");
			row++;
		}
		std::fprintf(out.get(), "\n");

	}
/*
	for(row = 0; row < m_z->RowNo(); row++)
	{
		fprintf(out, "%6d", row);
		for(i = 0; i < eigNo; i++)
			fprintf(out, "%18.6E", m_z->Get(row, i));
		fprintf(out, "\n");
	}
*/
	//////////////////////////////////////////////////////////////////////////////////////////
	size_t coefNo = EltBack().P();
	std::fprintf(out.get(), "\n\nABSOLUTE SUM OF LAST '%lu' COEFFICIENTS\n", static_cast<unsigned long>(coefNo));
	
    for(std::size_t i = 0UL; i < eigNo; i++)
    {
        std::fprintf(out.get(), "     Eig[%2lu]      ", static_cast<unsigned long>(i));
    }
    std::fprintf(out.get(), "\n");
	
    for (std::size_t i = 0; i < eigNo; i++)
	{
	    std::fprintf(out.get(), "==================");
	}
    std::fprintf(out.get(), "\n");

	std::vector<double> sum(eigNo);
	for (std::size_t i = 0UL; i < eigNo; i++)
	{
        for (row = m_z->RowNo() - coefNo; row < m_z->RowNo(); row++)
        {
            sum[i] += fabs(m_z->Get(row, i));
        }
	}
    
    for (std::size_t i = 0UL; i < eigNo; i++)
    {
        std::fprintf(out.get(), "%18.6E", sum[i]);
    }
	std::fprintf(out.get(), "\n");
}

//
// Finds the element with the largest coefficient contained in the minimal coefficients
// BUBBLE functins are considered only.
//
void EigProb::MaxMinCoef(std::vector<EltInfo>& eltInfo) const
{
const size_t eigNo = eltInfo.size();
double coef, minCoef;
int dof;

        for(size_t i = 0; i < eigNo; i++) // For each eigenfunction
	{
		eltInfo[i].m_maxMinCoef = 0; // Inicjalizacja
		for(size_t n = 0; n < EltNo(); n++) // For each element
		{
			const Element& e = Elt(n);
			minCoef = DBL_MAX;

			for(size_t j = 1; j < e.DofNo() - 1; j++) // For each BUBBLE DOF at element
			{
				dof = e.m_dof[j];
				if(dof < 0) // Skip Dirichlet boundary conditions
					continue;

				// Find the minimal coefficient for element "e"
				coef = fabs(m_z->Get(dof, i));
				if(coef < minCoef)
					minCoef = coef;
			}

			// Set the largest coefficient for the smallest coefficients
			if(minCoef > eltInfo[i].m_maxMinCoef)
			{
				eltInfo[i].m_eltId = n;
				eltInfo[i].m_maxMinCoef = minCoef;
				eltInfo[i].m_eigVal = i;
			}
		}
	}
}

