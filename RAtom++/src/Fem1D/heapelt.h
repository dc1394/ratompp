#ifndef __RATOM_HEAPELT_H__
#define __RATOM_HEAPELT_H__


/** \brief Element of heap. Used for adaptive algorithm.
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/
// March 16th, 2014 Modified by dc1394

class HeapElt
{
public:
	HeapElt(void);
	HeapElt(double left, double right, double delta, const std::vector<double>& c);
	~HeapElt(void);

	double Xinv(double x) const;
	// March 16th, 2014 Added by dc1394
	double Getc2() const;

	bool operator<(const HeapElt& e) const;

public:
        // Left end of interval
	double m_left;

        // Right end of inteval
	double m_right;

        // Approximation error for interval
	double m_delta;

        // Coefficientt of linear transformation (mapping) X
	double m_c1, m_c2;


        // Vector of interpolation coefficients
	std::vector<double> m_coef;
};

#endif

