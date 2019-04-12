#include "stdafx.h"
#include "clpmtx.h"
// #include "Except.h"
#include <memory>       // for std::unique_ptr

#ifndef USE_MKL
extern "C"
{
    std::int32_t dgesv_(int *n, int *nrhs, double *a, int *lda, 
         std::int32_t *ipiv, double *b, int *ldb, int *info);
	
    std::int32_t dsysv_(char *uplo, int *n, int *nrhs, double *a, 
	std::int32_t *lda, int *ipiv, double *b, int *ldb, 
		double *work, std::int32_t *lwork, int *info);
}
#endif

//!
//! Konstruktor 
//!
ClpMtx::ClpMtx()
{
	m_rowNo = 0;
	m_colNo = 0;
	// comment out by dc1394 - Jan/14/2014
	//m_array = nullptr;
}


//!
//! Konstruktor
//! rowNo - liczba wierszy
//! colNo - liczba kolumn
//!
ClpMtx::ClpMtx(size_t rowNo, size_t colNo)
{
    //delete[] m_array;
	m_rowNo = rowNo;
	m_colNo = colNo;
	//m_array = new double[rowNo * colNo];
	m_array.clear();
    m_array.resize(rowNo * colNo);

	Zero();
}

//!
//! Definiuje nowy rozmiar macierzy. Poprzednia zawartosc jest niszczona
//! rowNo - liczba wierszy
//! colNo - liczba kolumn
//!
void ClpMtx::SetSize(size_t rowNo, size_t colNo)
{
    //delete[] m_array;
	m_rowNo = rowNo;
	m_colNo = colNo;
	//m_array = new double[rowNo * colNo];
	m_array.clear();
    m_array.resize(rowNo * colNo);

	Zero();
}

//!
//! Zwraca element (row, col)
//!
double ClpMtx::Get(size_t row, size_t col) const
{
	return m_array[Elt(row, col)];
}

//!
//! Ustala wartosc elementu (row, col)
//!
double& ClpMtx::Set(size_t row, size_t col)
{
	return m_array[Elt(row, col)];
}

//!
//! Zwraca indeks elementu w tablicy "m_array"
//!
size_t ClpMtx::Elt(size_t row, size_t col) const
{
	assert(row < m_rowNo);
	assert(col < m_colNo);

	return col * m_rowNo + row;
}

//!
//! Rozwiazuje uklad rowan
//!	(*)	A x = b
//! Macierz A musi byc kwadratowa
//!
//! b - [IN/OUT] prawa strona rowania (*). Rozmiar rego wektora musi byc rowny rozmiarowi macierzy A.
//! On exit zawiera rozwiazanie.
//! 
void ClpMtx::Dgesv(const Vec& b, Vec& x)
{
    auto n = static_cast<std::int32_t>(m_colNo);
    std::int32_t nrhs = 1;
    std::int32_t lda = n;
    std::int32_t ldb = n;
    std::int32_t info;

	assert(m_colNo == m_rowNo);

	// Kopiowanie. Rozwiazanie zwracane jest na wekotrze "x"
	x = b;

#ifdef USE_MKL
    std::vector<std::int32_t, util::mkl_allocator<std::int32_t> > ipiv(n);
#else
    std::vector<std::int32_t> ipiv(n);
#endif
    dgesv_(&n, &nrhs, m_array.data(), &lda, ipiv.data(), x.data(), &ldb, &info);

	if (info != 0)
    {
		throw std::invalid_argument("Error in 'ClpMtx::Solve'");
    }

//  assert(ret == 0);
//	assert(info == 0);
}

//  DSYSV computes the solution to a real system of linear equations
//      A * X = B,
//  where A is an N-by-N symmetric matrix and X and B are vectors.
//  Matrix A must be upper triangular matrix
//  
void ClpMtx::Dsysv(const Vec& b, Vec& x)
{
    char uplo = 'U';  // Upper triangle of A is stored

    auto n = static_cast<std::int32_t>(m_colNo); // The number of linear equations, i.e., the order of the matrix A.
    std::int32_t nrhs = 1; 
    std::int32_t lda = n;
    std::int32_t * ipiv = nullptr;
    std::int32_t ldb = n;
    double tmp[2];
    std::int32_t lwork;
    std::int32_t info;

	// Kopiowanie. Rozwiazanie zwracane jest na wekotrze "x"
	x = b;

	// Zapytanie o wymagana pamiec
	lwork = -1;
	dsysv_(&uplo, &n, &nrhs, m_array.data(), &lda, ipiv, x.data(), &ldb, tmp, &lwork, &info);
	lwork = static_cast<std::int32_t>(tmp[0]);

#ifdef USE_MKL
    std::vector<double, util::mkl_allocator<double> > work(lwork);
    std::vector<std::int32_t, util::mkl_allocator<std::int32_t> > ipiv2(n);
#else
    std::vector<double> work(lwork);
    std::vector<std::int32_t> ipiv2(n);
#endif

    dsysv_(&uplo, &n, &nrhs, m_array.data(), &lda, ipiv2.data(), x.data(), &ldb, work.data(), &lwork, &info);

	if (info != 0)
    {
		throw std::invalid_argument("Error in 'ClpMtx::Dsysv'");
    }
}


//!
//! Zerowanie macierzy 
//!
void ClpMtx::Zero()
{
	// Wszystkie elemnty sa rowne zero
	for(size_t i = 0; i < m_rowNo * m_colNo; i++)
    {
		m_array[i] = 0;
    }
}

//!
//! Zapisuje macierz do pliku "path"
//!
void ClpMtx::Write(const char* path, bool rowId) const
{
    auto out = std::unique_ptr<FILE, decltype(&std::fclose)>(std::fopen(path, "wt"), std::fclose);
	for (auto row = 0U; row < RowNo(); row++)
	{
		if(rowId)
			std::fprintf(out.get(), "ROW=%2lu ", static_cast<unsigned long>(row));

		for (auto col = 0U; col < ColNo(); col++)
			std::fprintf(out.get(), "%15.6E", Get(row, col));
        std::fprintf(out.get(), "\n");
	}
}

