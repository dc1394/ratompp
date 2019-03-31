#ifndef __RATOM_CLPMTXBAND_H__
#define __RATOM_CLPMTXBAND_H__


/** \brief Wrapper for band matrix from LAPACK library
*
* \author Zbigniew Romanowski [ROMZ]
*
*/

#include "clpmtx.h"


class ClpMtxBand
{
public:
	ClpMtxBand(size_t n, size_t ku, size_t kl);
	~ClpMtxBand(void);

	double Get(size_t row, size_t col) const;
	double& Set(size_t row, size_t col);

	void Eigen(size_t eigNo, double abstol, Vec& w, ClpMtx& z);
	void EigenGen(size_t eigNo, double abstol, Vec& w, ClpMtx& z, ClpMtxBand& b);

	void SolveSymPos(std::unique_ptr<Vec> const & b, std::unique_ptr<Vec> & x);

	size_t ColNo() const;


	void Zero();

	void Write(const char* path) const;

private:
	size_t RowEx(size_t row, size_t col) const;
	bool InBand(size_t row, size_t col) const;

private:
        // Reactangular matrix used for strong band matrix
	ClpMtx m_mtx;

        // kl - number of subdiagonals
	size_t m_kl;

        // ku - number of superdiagonals
	size_t m_ku;

        // Auxiliary mememer, works as "zero".
	double m_zero;
};

#endif

