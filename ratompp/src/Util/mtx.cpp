#include "stdafx.h"
#include "mtx.h"

//
// Constructor
//
Mtx::Mtx(size_t rowNo, size_t colNo) : m_colNo(colNo), m_rowNo(rowNo)
{
	m_array = new double[m_rowNo * m_colNo];
}

//
// Destructor
//
Mtx::~Mtx(void)
{
	delete [] m_array;
}

//
// Returns element (row, col)
//
double Mtx::Get(size_t row, size_t col) const
{
	return m_array[EltIdx(row, col)];
}

//
// Sets element (row, col)
//
double& Mtx::Set(size_t row, size_t col)
{
	return m_array[EltIdx(row, col)];
}

//
// Returns index of element (row, col)
//
size_t Mtx::EltIdx(size_t row, size_t col) const
{
	return row + col * m_rowNo;
}


