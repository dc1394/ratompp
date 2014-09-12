#ifndef __RATOM_INT1DDISCRET_H__
#define __RATOM_INT1DDISCRET_H__


/** \brief Algorytm obliczania calki 1D i 3D z funkcji jednej zmiennej
*
* \author Zbigniew Romanowski [ROMZ]
*
* \version 17-Sty-2007 (romz) - first release
*
*/

class Int1DDiscret
{
public:
	Int1DDiscret(void);
	~Int1DDiscret(void);

//	static double Spher(const Vec& f, const Mesh* mesh);
//	static double Spher(const Vec& f, const Vec& g, const Mesh* mesh);

	static double Order4a(const Vec& f, double h);
	static double Order3(const Vec& f, double h);
	static double Order4(const Vec& f, double h);

	static double NewtonCotes3(const Vec& f, double h);
	static double NewtonCotes4(const Vec& f, double h);
	static double NewtonCotes5(const Vec& f, double h);
	static double NewtonCotes6(const Vec& f, double h);

	static double NewtonCotes6Ex(const Vec& f, double h);
};


#endif

