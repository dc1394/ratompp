#ifndef __RATOM_FUN2D_H__
#define __RATOM_FUN2D_H__


/** \brief Represents function from R^2 into R.
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/

class Fun2D
{
public:
	Fun2D(void) { }
	virtual ~Fun2D(void) { }

	virtual double Get(double x, double y) const = 0;
};

#endif

