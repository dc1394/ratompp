#ifndef __RATOM_FUN3D_H__
#define __RATOM_FUN3D_H__


/** \brief Represents function from R^3 into R.
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/

class Fun3D
{
public:
	Fun3D(void) { }
	virtual ~Fun3D(void) { }

	virtual double Get(double x, double y, double z) const = 0;
};


#endif

