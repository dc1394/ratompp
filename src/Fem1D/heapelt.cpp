#include "stdafx.h"
#include "heapelt.h"

//
// Constructor
//
HeapElt::HeapElt(void)  : m_left(0), m_right(0), m_delta(0), m_c1(0), m_c2(0)
{
}

//
// Constructor
//
HeapElt::HeapElt(double left, double right, double delta, const std::vector<double>& c) : 
		m_left(left), m_right(right), m_delta(delta), m_coef(c)
{
	assert(right > left); 

	m_c1 = 0.5 * (m_right + m_left);
	m_c2 = 0.5 * (m_right - m_left);
}

//
// Destructor
//
HeapElt::~HeapElt(void)
{
}

//
// Comparison operator required for heap operations
//
bool HeapElt::operator<(const HeapElt& e) const
{
	return (m_delta < e.m_delta);
}

//
// Returns local variable "s" corresponding to global variable "x"
//
double HeapElt::Xinv(double x) const
{
	return (x - m_c1) / m_c2; 
}

// March 16th, 2014 Added by dc1394
//
// Returns c2
//
double HeapElt::Getc2() const
{
	return m_c2;
}
