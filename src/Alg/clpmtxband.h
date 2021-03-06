#ifndef __RATOM_CLPMTXBAND_H__
#define __RATOM_CLPMTXBAND_H__

/** \brief Wrapper for band matrix from LAPACK library
*
* \author Zbigniew Romanowski [ROMZ]
*
*/

#include "clpmtx.h"
#include <memory>   // for std::unique_ptr

class ClpMtxBand final
{
public:
	ClpMtxBand(std::size_t n, std::size_t ku, std::size_t kl);
	~ClpMtxBand() = default;

	double Get(std::size_t row, std::size_t col) const;
	double& Set(std::size_t row, std::size_t col);

	void Eigen(std::size_t eigNo, double abstol, Vec& w, ClpMtx& z);
	void EigenGen(std::size_t eigNo, double abstol, Vec& w, ClpMtx& z, ClpMtxBand& b);

	void SolveSymPos(std::unique_ptr<Vec> const & b, std::unique_ptr<Vec> & x);

	std::size_t ColNo() const;


	void Zero();

	void Write(const char* path) const;

private:
	std::size_t RowEx(std::size_t row, std::size_t col) const;
	bool InBand(std::size_t row, std::size_t col) const;

private:
        // Reactangular matrix used for strong band matrix
	ClpMtx m_mtx;

        // kl - number of subdiagonals
	std::size_t m_kl;

        // ku - number of superdiagonals
	std::size_t m_ku;

        // Auxiliary mememer, works as "zero".
	double m_zero;
};

#endif

