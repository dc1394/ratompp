#ifndef __RATOM_CLPMTX_H__
#define __RATOM_CLPMTX_H__

#include "../Util/vec.h"

/** \brief Wrapper for rectangular matrix from LAPACk library
*
* \author Zbigniew Romanowski [ROMZ]
*
*/


class ClpMtx final
{
	friend class ClpMtxBand;
public:
	ClpMtx();
	ClpMtx(std::size_t rowNo, std::size_t colNo);
	~ClpMtx() = default;

	void SetSize(std::size_t rowNo, std::size_t colNo);

	double Get(std::size_t row, std::size_t col) const;
	double& Set(std::size_t row, std::size_t col);

	std::size_t ColNo() const { return m_colNo; }
	std::size_t RowNo() const { return m_rowNo; }

	void Dgesv(const Vec& b, Vec& x);
	void Dsysv(const Vec& b, Vec& x);

	void Zero();

	void Write(const char* path, bool rowId = true) const;

private:
	std::size_t Elt(std::size_t row, std::size_t col) const;

private:
        // Number of columns
	std::size_t m_colNo;

        // Number of rows
	std::size_t m_rowNo;

        // Array with data
#ifdef USE_MKL
    std::vector<double, util::mkl_allocator<double> > m_array;
#else
    std::vector<double> m_array;
#endif
};

#endif

