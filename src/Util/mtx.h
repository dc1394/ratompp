#ifndef __RATOM_MTX_H__
#define __RATOM_MTX_H__


/** \brief Rectangular matrix
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/

class Mtx
{
public:
	Mtx(size_t rowNo, size_t colNo);
	virtual ~Mtx(void);

	double Get(size_t row, size_t col) const;
	double& Set(size_t row, size_t col);

	size_t RowNo() const { return m_rowNo; }
	size_t ColNo() const { return m_colNo; }

private:
	size_t EltIdx(size_t row, size_t col) const;

private:
	size_t m_colNo;

	size_t m_rowNo;

	double* m_array;
};


#endif

