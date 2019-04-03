#ifndef __RATOM_GAUSS2D_H__
#define __RATOM_GAUSS2D_H__


/**
*	\brief Oblicza calke w dwoch wymiarach po prostkoacie. 
*		Prostokat ma boki rownolegle do osi ukladu wspolrzednych. 
*		Algorytm oparty na iloczynie kwadratur Gaussa w jednym wymiarze
*
*		Algorytm adaptacyjny dzieli prostokaty na CZTERY mniejsze prostokaty.
*
*
*	\author Zbigniew Romanowski, ROMZ
*
*	\version 19-Cze-2007 [romz]
*
*/



class Gauss2D
{
public:
	Gauss2D(double epsAbs, size_t order);
	~Gauss2D(void);

	double ErrAbs() const;
	double ErrRel() const;
	size_t RecNo() const;
	size_t FunCall() const;


	double CalcAdapt(const Fun2D& f, double x0, double x1, double y0, double y1);
	double Calc(const Fun2D& f, double x1, double x2, double y1, double y2);


private:
	double Gauss(const Fun2D& f) const;
	

	double Run(const Fun2D& f, double x0, double x1, double y0, double y1);
	static void GauLeg(double x1, double x2, std::vector<double>& x, std::vector<double>& w);

private:
	// Wagi wezlow
	std::vector<double> m_w;

	// Wspolrzedne wezlow
	std::vector<double> m_x;

	// Wspolrzedne wezlow
	std::vector<double> m_y;

	// Liczba obliczen funkcji
	size_t m_funCall;

	// Koncowa liczba prostokatow
	size_t m_recNo;

	// Dokladnosc z jaka ma byc obliczona calka
	double m_epsAbs;
};

#endif

