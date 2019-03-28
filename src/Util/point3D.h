#ifndef __RATOM_POINT3D_H__
#define __RATOM_POINT3D_H__


/** \brief Point in R^3 space
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/

class Point3D
{
public:
	Point3D(void) : m_x(0), m_y(0), m_z(0) { }
	Point3D(double x, double y, double z) : m_x(x), m_y(y), m_z(z) { }
	virtual ~Point3D(void) { }

	// Zwraca odleglosc miedzy punktem "a" i "b".
	static double Dist(const Point3D* a, const Point3D* b)
	{
	double x, y, z;

		x = a->m_x - b->m_x;
		y = a->m_y - b->m_y;
		z = a->m_z - b->m_z;

		return sqrt(x * x + y * y + z * z);
	}

public:
	double m_x;
	double m_y;
	double m_z;

};


#endif

