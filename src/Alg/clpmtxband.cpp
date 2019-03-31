#include "stdafx.h"
#include "clpmtxband.h"

// March 22nd, 2014 Added by dc1394
#include <mkl.h>

// March 22nd, 2014 Modified by dc1394
//extern "C"
//{
//void dsbevx_(char *jobz, char *range, char *uplo, int *n, 
//	int *kd, double *ab, int *ldab, double *q, int *
//	ldq, double *vl, double *vu, int *il, int *iu, 
//	double *abstol, int *m, double *w, double *z__, 
//	int *ldz, double *work, int *iwork, int *ifail, 
//	int *info);
//	
//void dsbgvx_(char *jobz, char *range, char *uplo, int *n, 
//	int *ka, int *kb, double *ab, int *ldab, double *
//	bb, int *ldbb, double *q, int *ldq, double *vl, 
//	double *vu, int *il, int *iu, double *abstol, int 
//	*m, double *w, double *z__, int *ldz, double *work, 
//	int *iwork, int *ifail, int *info);	
//	
//void dpbsvx_(char *fact, char *uplo, int *n, int *kd, 
//	int *nrhs, double *ab, int *ldab, double *afb, 
//	int *ldafb, char *equed, double *s, double *b, int *
//	ldb, double *x, int *ldx, double *rcond, double *ferr,
//	 double *berr, double *work, int *iwork, int *info);	
//}


// Constructor
// n - number of columns
// kl - number of subdiagonals
// ku - number of superdiagonals
//
ClpMtxBand::ClpMtxBand(size_t n, size_t ku, size_t kl) : m_mtx(kl + ku + 1, n), m_zero(0)
{
	m_kl = kl;
	m_ku = ku;
}

//
// Destructor
//
ClpMtxBand::~ClpMtxBand(void)
{
}

//
// Returns element (row, col)
//
double ClpMtxBand::Get(size_t row, size_t col) const
{
	if(InBand(row, col))
		return m_mtx.Get(RowEx(row, col), col);

	return 0;
}

//
// Sets value of element (row, col)
//
double& ClpMtxBand::Set(size_t row, size_t col)
{
	if(InBand(row, col))
		return m_mtx.Set(RowEx(row, col), col);

	return m_zero;
}


//
// Returns row id in rectangular matrix "m_mtx"
// for element (row, col) in band matrix "this"
//
size_t ClpMtxBand::RowEx(size_t row, size_t col) const
{
	return m_ku + row - col;
}

//
// Returns "true", if element (row, col) belongs to band matrix "this"
//
bool ClpMtxBand::InBand(size_t row, size_t col) const
{
const bool b1 = (col <= row + m_ku);
const bool b2 = (row <= col + m_kl);

	return (b1 && b2);
}

//
// returns number of columns in matrix
//
size_t ClpMtxBand::ColNo() const
{
	return m_mtx.m_colNo;
}


//
// WRAPPER for "dsbevx" procedure from LAPACK
//
// Oblicza kilka najmniejszych wartosci wlasnych i wektorow wlasnych zagadnienia wlasnego
//
//		A x = \lambda x
//
//	gdzie macierz A pasmowa, symetryczna. Macierze A podawane jest jako "trojkatna gorna" orozmiarze $N$
//
// eigNo  - [IN] liczba obliczanych najmniejszych wartosci wlasnych
// abstol - [IN] The absolute error tolerance for the eigenvalues.
// w  - [OUT] wektor z obliczonymi wartosciami wlasnymi. MUSI BYS rozmiaru $N$.
// z -  [OUT] wektor z obliczonymi wektorami wlasnymi. MUSI BYS rozmiaru $N x N$.
//
void ClpMtxBand::Eigen(size_t eigNo, double abstol, Vec& w, ClpMtx& z)
{
int m;

int n = static_cast<int>(m_mtx.m_colNo);
int ldz = n; // The leading dimension of the array Z.  LDZ >= 1, and if JOBZ = 'V', LDZ >= max(1,N).
int ldq = n; // The leading dimension of the array Q.  If JOBZ = 'N', LDQ >= 1. If JOBZ = 'V', LDQ >= max(1,N).
int ku = static_cast<int>(m_ku);
int ldab = ku + 1; // The leading dimension of the array AB.  LDAB >= KA+1.
int il = 1; // If RANGE='I', the indice of the smallest eigenvalues to be returned.
int iu = static_cast<int>(eigNo);

double vl = 0, vu = 0; // Not referenced if RANGE = 'A' or 'I'.

char jobz = 'V';  // Compute eigenvalues and eigenvectors
char range = 'I'; // the IL-th through IU-th eigenvalues will be found
char uplo = 'U';  // Upper triangles of A and B are stored;

double *q, *work;
int *iwork, *ifail, info;


	// Tylko "trojkatna gorna" jest zdefiniowana
	assert(m_kl == 0);
	
	// March 22nd, 2014 Modified by dc1394
	//q = (double*)malloc(n * n * sizeof(double));
	//work = (double*)malloc(7 * n * sizeof(double));
	//iwork = (int*)malloc(5 * n * sizeof(int));
	//ifail = (int*)malloc(n * sizeof(int));
	q = reinterpret_cast<double*>(mkl_malloc(n * n * sizeof(double), 64));
	work = reinterpret_cast<double*>(mkl_malloc(7 * n * sizeof(double), 64));
	iwork = reinterpret_cast<int*>(mkl_malloc(5 * n * sizeof(int), 64));
	ifail = reinterpret_cast<int*>(mkl_malloc(n * sizeof(int), 64));

	dsbevx_(&jobz, &range, &uplo, &n, &ku, m_mtx.m_array.get(), &ldab, q, &ldq, 
		&vl, &vu, &il, &iu, &abstol, &m, 
		&w.front(), 
		z.m_array.get(), &ldz, work, iwork, ifail, &info);

	// March 22nd, 2014 Modified by dc1394
	//free(work);
	//free(iwork);
	//free(ifail);
	//free(q);
	mkl_free(work);
	mkl_free(iwork);
	mkl_free(ifail);
	mkl_free(q);

	if(info != 0)
		throw std::invalid_argument("Error in 'ClpMtxBand::EigenGen'");
//	assert(ret == 0);
//	assert(info == 0);

}

//
// WRAPPER for "dsbgvx" procedure from LAPACK
//
// Oblicza kilka najmniejszych wartosci wlasnych i wektorow wlasnych uogolnionego zagadnienia wlasnego
//
//		A x = \lambda B x
//
//	gdzie A, B pasmowe, symetryczne oraz B jest dodatnio okreslona.
// Macierze A i B podawane sa jako trojkatne gorne.
//
// eigNo  - [IN] liczba obliczanych najmniejszych wartosci wlasnych
// abstol - [IN] The absolute error tolerance for the eigenvalues.
// w - [OUT] wektor z obliczonymi wartosciami wlasnymi. MUSI BYS rozmiaru $N$.
// z - [OUT] wektor z obliczonymi wektorami wlasnymi. MUSI BYS rozmiaru $N x N$.
// b - [IN]  macierz B, pasmowa, symetryczna, dodatnio okreslona
//
void ClpMtxBand::EigenGen(size_t eigNo, double abstol, Vec& w, ClpMtx& z, ClpMtxBand& b)
{
int m;

int n = static_cast<int>(m_mtx.m_colNo);
int ldz = n; // The leading dimension of the array Z.  LDZ >= 1, and if JOBZ = 'V', LDZ >= max(1,N).
int ldq = n; // The leading dimension of the array Q.  If JOBZ = 'N', LDQ >= 1. If JOBZ = 'V', LDQ >= max(1,N).

int ka = static_cast<int>(m_ku);
int ldab = ka + 1; // The leading dimension of the array AB.  LDAB >= KA+1.

int kb = static_cast<int>(b.m_ku);
int ldbb = kb + 1; // The leading dimension of the array BB.  LDBB >= KB+1. 

int il = 1; // If RANGE='I', the indice of the smallest eigenvalues to be returned.

double vl = 0, vu = 0; // Not referenced if RANGE = 'A' or 'I'.
int iu = static_cast<int>(eigNo);

char jobz = 'V';  // Compute eigenvalues and eigenvectors
char range = 'I'; // the IL-th through IU-th eigenvalues will be found
char uplo = 'U';  // Upper triangles of A and B are stored;

double *q, *work;
int *iwork, *ifail, info;

// Tylko "trojkatna gorna" jest zdefiniowana
	assert(m_kl == 0);

	// March 22nd, 2014 Modified by dc1394
	//q = (double*)malloc(n * n * sizeof(double));
	//work = (double*)malloc(7 * n * sizeof(double));
	//iwork = (int*)malloc(5 * n * sizeof(int));
	//ifail = (int*)malloc(iu * sizeof(int));
	q = reinterpret_cast<double*>(mkl_malloc(n * n * sizeof(double), 64));
	work = reinterpret_cast<double*>(mkl_malloc(7 * n * sizeof(double), 64));
	iwork = reinterpret_cast<int*>(mkl_malloc(5 * n * sizeof(int), 64));
	ifail = reinterpret_cast<int*>(mkl_malloc(iu * sizeof(int), 64));

	dsbgvx_(&jobz, &range, &uplo, &n, &ka, &kb, m_mtx.m_array.get(), &ldab, 
		b.m_mtx.m_array.get(), &ldbb, q, &ldq, &vl, 
		&vu, &il, &iu, &abstol, &m, 
		&w.front(), 
		z.m_array.get(), &ldz, work, iwork, ifail, &info);

	// March 22nd, 2014 Modified by dc1394
	//free(work);
	//free(iwork);
	//free(ifail);
	//free(q);
	mkl_free(work);
	mkl_free(iwork);
	mkl_free(ifail);
	mkl_free(q);

	if(info != 0)
		throw std::invalid_argument("Error in 'ClpMtxBand::EigenGen'");
//	assert(ret == 0);
//	assert(info == 0);
}

//
// WRAPPER for "dpbsvx" procedure from LAPACK
//
// Rozwiazuje uklad rownan:
//		A x = b
// $A$ jest macierza symetryczna, dodatnia okreslona.
// $A$ jest przechowywana jako trojkatna gorna.
//
void ClpMtxBand::SolveSymPos(std::unique_ptr<Vec> const & b, std::unique_ptr<Vec> & x)
{
char fact = 'N';  // The matrix A will be copied to AFB and factored.
char equed = 'N'; // Specifies the form of equilibration that was done.   
char uplo = 'U';  // Upper triangles of A and B are stored;

int nrhs = 1; // The number of right-hand sides
int info; 

int n = static_cast<int>(m_mtx.m_colNo);
int ldb = n; // The leading dimension of the array B
int ldx = n; // The leading dimension of the array X

int kd = static_cast<int>(m_ku);
int ldab = kd + 1; // The leading dimension of the array AB.  LDAB >= KA+1.
int ldafb = kd + 1; // The leading dimension of the array AFB.

double rcond; // The estimate of the reciprocal condition number
double *afb, *work;
int* iwork;
double ferr[1]; // The estimated forward error bound
double berr[1]; // The componentwise relative backward error

	// Macierz musi byc trojkatna gorna
	assert(m_kl == 0);

	// March 22nd, 2014 Modified by dc1394
	//afb = (double*)malloc(ldafb * n * sizeof(double));
	//work = (double*)malloc(3 * n * sizeof(double));
	//iwork = (int*)malloc(n * sizeof(int));
	afb = reinterpret_cast<double*>(mkl_malloc(ldafb * n * sizeof(double), 64));
	work = reinterpret_cast<double*>(mkl_malloc(3 * n * sizeof(double), 64));
	iwork = reinterpret_cast<int*>(mkl_malloc(n * sizeof(int), 64));

	dpbsvx_(&fact, &uplo, &n, &kd, 
		&nrhs, m_mtx.m_array.get(), &ldab, afb, &ldafb, 
		&equed, NULL, 
		(double*)(&b->front()), 
		&ldb, 
		&x->front(), 
		&ldx, 
		&rcond, ferr, berr, work, iwork, 
		&info);

	// March 22nd, 2014 Modified by dc1394
	//free(afb);
	//free(work);
	//free(iwork);
	mkl_free(afb);
	mkl_free(work);
	mkl_free(iwork);

	if(info != 0)
		throw std::invalid_argument("Error in 'ClpMtxBand::Solve'");
//	assert(ret == 0);
//	assert(info == 0);
}

//
// Zeroeing all elements in band matrix
//
void ClpMtxBand::Zero()
{
	m_mtx.Zero();
}

//
// Writes matrix into file "path"
//
void ClpMtxBand::Write(const char* path) const
{
FILE* out;
size_t kl, row, ku;

	out = fopen(path, "wt");
	fprintf(out, "DIAGONAL\n");
	for(row = 0; row < ColNo(); row++)
		fprintf(out, "%4lu %lf\n", static_cast<unsigned long>(row), Get(row, row));

	// Subdiagonals
	for(kl = 1; kl <= m_kl; kl++)
	{
		fprintf(out, "SUB-DIAGONAL %lu\n", static_cast<unsigned long>(kl));
		for(row = kl; row < ColNo(); row++)
			fprintf(out, "%4lu %lf\n", static_cast<unsigned long>(row), Get(row, row - kl));
	}

	// Superdiagonals
	for(ku = 1; ku <= m_ku; ku++)
	{
		fprintf(out, "SUPER-DIAGONAL %lu\n", static_cast<unsigned long>(ku));
		for(row = 0; row < ColNo(); row++)
		{
			if(row + ku < ColNo())
				fprintf(out, "%4lu %lf\n", static_cast<unsigned long>(row), Get(row, row + ku));
		}
	}

	fclose(out);
}



