#include "stdafx.h"
#include "clpmtx.h"
// #include "Except.h"


//extern "C"
//{
//    int dgesv_(int *n, int *nrhs, double *a, int *lda, 
//         int *ipiv, double *b, int *ldb, int *info);
//	
//    int dsysv_(char *uplo, int *n, int *nrhs, double *a, 
//	int *lda, int *ipiv, double *b, int *ldb, 
//		double *work, int *lwork, int *info);
//}

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
int n = static_cast<int>(m_colNo);
int nrhs = 1;
int lda = n;
int ldb = n;
int info;

	assert(m_colNo == m_rowNo);

	// Kopiowanie. Rozwiazanie zwracane jest na wekotrze "x"
	x = b;

	auto ipiv = reinterpret_cast<int *>(mkl_malloc(n * sizeof(int), 64));
	dgesv_(&n, &nrhs, m_array.data(), &lda, ipiv, x.data(), &ldb, &info);
	mkl_free(ipiv);
	
	if (info != 0)
    {
		throw std::invalid_argument("Error in 'ClpMtx::Solve'");
    }
//	assert(ret == 0);
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

int n = static_cast<int>(m_colNo); // The number of linear equations, i.e., the order of the matrix A.
int nrhs = 1; 
int lda = n;
int * ipiv = nullptr;
int ldb = n;
double tmp[2];
int lwork;
int info;

	// Kopiowanie. Rozwiazanie zwracane jest na wekotrze "x"
	x = b;

	// Zapytanie o wymagana pamiec
	lwork = -1;
	dsysv_(&uplo, &n, &nrhs, m_array.data(), &lda, ipiv, x.data(), &ldb, tmp, &lwork, &info);
	lwork = static_cast<int>(tmp[0]);

	auto work = reinterpret_cast<double *>(mkl_malloc(lwork * sizeof(double), 64));
	ipiv = reinterpret_cast<int *>(mkl_malloc(n * sizeof(int), 64));

	dsysv_(&uplo, &n, &nrhs, m_array.data(), &lda, ipiv, x.data(), &ldb, work, &lwork, &info);

	mkl_free(ipiv);
	mkl_free(work);

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
		m_array[i] = 0;
}

//!
//! Zapisuje macierz do pliku "path"
//!
void ClpMtx::Write(const char* path, bool rowId) const
{
    auto out = std::unique_ptr<FILE, decltype(&std::fclose)>(fopen(path, "wt"), std::fclose);
	for (auto row = 0U; row < RowNo(); row++)
	{
		if(rowId)
			fprintf(out.get(), "ROW=%2lu ", static_cast<unsigned long>(row));

		for (auto col = 0U; col < ColNo(); col++)
			fprintf(out.get(), "%15.6E", Get(row, col));
		fprintf(out.get(), "\n");
	}
}


