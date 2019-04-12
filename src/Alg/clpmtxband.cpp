#include "stdafx.h"
#include "clpmtxband.h"

// March 22nd, 2014 Modified by dc1394
#ifndef USE_MKL
extern "C"
{
    void dsbevx_(char *jobz, char *range, char *uplo, std::int32_t *n, 
	    std::int32_t *kd, double *ab, std::int32_t *ldab, double *q, std::int32_t *
	    ldq, double *vl, double *vu, std::int32_t *il, std::int32_t *iu, 
	    double *abstol, std::int32_t *m, double *w, double *z__, 
	    std::int32_t *ldz, double *work, std::int32_t *iwork, std::int32_t *ifail, 
	    std::int32_t *info);
	
    void dsbgvx_(char *jobz, char *range, char *uplo, std::int32_t *n, 
	    std::int32_t *ka, std::int32_t *kb, double *ab, std::int32_t *ldab, double *
	    bb, std::int32_t *ldbb, double *q, std::int32_t *ldq, double *vl, 
	    double *vu, std::int32_t *il, std::int32_t *iu, double *abstol, std::int32_t 
	    *m, double *w, double *z__, std::int32_t *ldz, double *work, 
	    std::int32_t *iwork, std::int32_t *ifail, std::int32_t *info);	
	
    void dpbsvx_(char *fact, char *uplo, std::int32_t *n, std::int32_t *kd, 
	    std::int32_t *nrhs, double *ab, std::int32_t *ldab, double *afb, 
	    std::int32_t *ldafb, char *equed, double *s, double *b, std::int32_t *
	    ldb, double *x, std::int32_t *ldx, double *rcond, double *ferr,
	    double *berr, double *work, std::int32_t *iwork, std::int32_t *info);	
}
#endif

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
// Returns element (row, col)
//
double ClpMtxBand::Get(size_t row, size_t col) const
{
    if (InBand(row, col))
        return m_mtx.Get(RowEx(row, col), col);

    return 0;
}

//
// Sets value of element (row, col)
//
double& ClpMtxBand::Set(size_t row, size_t col)
{
    if (InBand(row, col))
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
    std::int32_t m;

    auto n = static_cast<std::int32_t>(m_mtx.m_colNo);
    std::int32_t ldz = n; // The leading dimension of the array Z.  LDZ >= 1, and if JOBZ = 'V', LDZ >= max(1,N).
    std::int32_t ldq = n; // The leading dimension of the array Q.  If JOBZ = 'N', LDQ >= 1. If JOBZ = 'V', LDQ >= max(1,N).
    auto ku = static_cast<std::int32_t>(m_ku);
    std::int32_t ldab = ku + 1; // The leading dimension of the array AB.  LDAB >= KA+1.
    std::int32_t il = 1; // If RANGE='I', the indice of the smallest eigenvalues to be returned.
    auto iu = static_cast<std::int32_t>(eigNo);

    double vl = 0, vu = 0; // Not referenced if RANGE = 'A' or 'I'.

    char jobz = 'V';  // Compute eigenvalues and eigenvectors
    char range = 'I'; // the IL-th through IU-th eigenvalues will be found
    char uplo = 'U';  // Upper triangles of A and B are stored;

    std::int32_t info;

    // Tylko "trojkatna gorna" jest zdefiniowana
    assert(m_kl == 0);

    // March 22nd, 2014 Modified by dc1394
#ifdef USE_MKL
    std::vector<double, util::mkl_allocator<double> > q(n * n);
    std::vector<double, util::mkl_allocator<double> > work(7 * n);
    std::vector<std::int32_t, util::mkl_allocator<std::int32_t> > iwork(5 * n);
    std::vector<std::int32_t, util::mkl_allocator<std::int32_t> > ifail(n);
#else
    std::vector<double> q(n * n);
    std::vector<double> work(7 * n);
    std::vector<std::int32_t> iwork(5 * n);
    std::vector<std::int32_t> ifail(n);
#endif   
    dsbevx_(&jobz, &range, &uplo, &n, &ku, m_mtx.m_array.data(), &ldab, q.data(), &ldq,
        &vl, &vu, &il, &iu, &abstol, &m,
        w.data(),
        z.m_array.data(), &ldz, work.data(), iwork.data(), ifail.data(), &info);

    // March 22nd, 2014 Modified by dc1394

    if (info != 0)
        throw std::invalid_argument("Error in 'ClpMtxBand::EigenGen'");
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
    std::int32_t m;

    auto n = static_cast<std::int32_t>(m_mtx.m_colNo);
    std::int32_t ldz = n; // The leading dimension of the array Z.  LDZ >= 1, and if JOBZ = 'V', LDZ >= max(1,N).
    std::int32_t ldq = n; // The leading dimension of the array Q.  If JOBZ = 'N', LDQ >= 1. If JOBZ = 'V', LDQ >= max(1,N).

    auto ka = static_cast<std::int32_t>(m_ku);
    std::int32_t ldab = ka + 1; // The leading dimension of the array AB.  LDAB >= KA+1.

    auto kb = static_cast<std::int32_t>(b.m_ku);
    std::int32_t ldbb = kb + 1; // The leading dimension of the array BB.  LDBB >= KB+1. 

    std::int32_t il = 1; // If RANGE='I', the indice of the smallest eigenvalues to be returned.

    double vl = 0, vu = 0; // Not referenced if RANGE = 'A' or 'I'.
    auto iu = static_cast<std::int32_t>(eigNo);

    char jobz = 'V';  // Compute eigenvalues and eigenvectors
    char range = 'I'; // the IL-th through IU-th eigenvalues will be found
    char uplo = 'U';  // Upper triangles of A and B are stored;

    std::int32_t info;

    // Tylko "trojkatna gorna" jest zdefiniowana
    assert(m_kl == 0);

#ifdef USE_MKL
    std::vector<double, util::mkl_allocator<double> > q(2 * n * n);
    std::vector<double, util::mkl_allocator<double> > work(7 * n);
    std::vector<std::int32_t, util::mkl_allocator<std::int32_t> > iwork(5 * n);
    std::vector<std::int32_t, util::mkl_allocator<std::int32_t> > ifail(iu);
#else
    std::vector<double> q(2 * n * n);
    std::vector<double> work(7 * n);
    std::vector<std::int32_t> iwork(5 * n);
    std::vector<std::int32_t> ifail(iu);
#endif

    dsbgvx_(&jobz, &range, &uplo, &n, &ka, &kb, m_mtx.m_array.data(), &ldab,
        b.m_mtx.m_array.data(), &ldbb, q.data(), &ldq, &vl,
        &vu, &il, &iu, &abstol, &m,
        w.data(),
        z.m_array.data(), &ldz, work.data(), iwork.data(), ifail.data(), &info);
    
    if (info != 0)
        throw std::invalid_argument("Error in 'ClpMtxBand::EigenGen'");
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

    std::int32_t nrhs = 1; // The number of right-hand sides
    std::int32_t info;

    auto n = static_cast<std::int32_t>(m_mtx.m_colNo);
    std::int32_t ldb = n; // The leading dimension of the array B
    std::int32_t ldx = n; // The leading dimension of the array X

    auto kd = static_cast<std::int32_t>(m_ku);
    std::int32_t ldab = kd + 1; // The leading dimension of the array AB.  LDAB >= KA+1.
    std::int32_t ldafb = kd + 1; // The leading dimension of the array AFB.

    double rcond; // The estimate of the reciprocal condition number
    double ferr[1]; // The estimated forward error bound
    double berr[1]; // The componentwise relative backward error

        // Macierz musi byc trojkatna gorna
    assert(m_kl == 0);

    // March 22nd, 2014 Modified by dc1394

#ifdef USE_MKL
    std::vector<double, util::mkl_allocator<double> > afb(ldafb * n);
    std::vector<double, util::mkl_allocator<double> > work(3 * n);
    std::vector<std::int32_t, util::mkl_allocator<std::int32_t> > iwork(n);
#else
    std::vector<double> afb(ldafb * n);
    std::vector<double> work(3 * n);
    std::vector<std::int32_t> iwork(n);
#endif

    dpbsvx_(
        &fact,
        &uplo,
        &n,
        &kd,
        &nrhs,
        m_mtx.m_array.data(),
        &ldab,
        afb.data(),
        &ldafb,
        &equed,
        nullptr,
        b->data(),
        &ldb,
        x->data(),
        &ldx,
        &rcond,
        ferr,
        berr,
        work.data(),
        iwork.data(),
        &info);

    if (info != 0)
        throw std::invalid_argument("Error in 'ClpMtxBand::Solve'");
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
    auto out = std::unique_ptr<FILE, decltype(&std::fclose)>(std::fopen(path, "wt"), std::fclose);
    fprintf(out.get(), "DIAGONAL\n");
    for (std::size_t row = 0UL; row < ColNo(); row++)
        fprintf(out.get(), "%4lu %lf\n", static_cast<unsigned long>(row), Get(row, row));

    // Subdiagonals
    for (std::size_t kl = 1UL; kl <= m_kl; kl++)
    {
        fprintf(out.get(), "SUB-DIAGONAL %lu\n", static_cast<unsigned long>(kl));
        for (std::size_t row = kl; row < ColNo(); row++)
            fprintf(out.get(), "%4lu %lf\n", static_cast<unsigned long>(row), Get(row, row - kl));
    }

    // Superdiagonals
    for (std::size_t ku = 1UL; ku <= m_ku; ku++)
    {
        fprintf(out.get(), "SUPER-DIAGONAL %lu\n", static_cast<unsigned long>(ku));
        for (std::size_t row = 0; row < ColNo(); row++)
        {
            if (row + ku < ColNo())
                fprintf(out.get(), "%4lu %lf\n", static_cast<unsigned long>(row), Get(row, row + ku));
        }
    }
}

