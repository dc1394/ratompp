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
	ClpMtx(size_t rowNo, size_t colNo);
	~ClpMtx() = default;

	void SetSize(size_t rowNo, size_t colNo);

	double Get(size_t row, size_t col) const;
	double& Set(size_t row, size_t col);

	size_t ColNo() const { return m_colNo; }
	size_t RowNo() const { return m_rowNo; }

	void Dgesv(const Vec& b, Vec& x);
	void Dsysv(const Vec& b, Vec& x);

	void Zero();

	void Write(const char* path, bool rowId = true) const;

private:
	size_t Elt(size_t row, size_t col) const;

private:
        // Number of columns
	size_t m_colNo;

        // Number of rows
	size_t m_rowNo;

        // Array with data
	// Arranged by dc1394 - Jan/14/2014
	//double* m_array;
    std::vector<double> m_array;
};

#endif

