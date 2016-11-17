#ifndef __RATOM_BRICK_H__
#define __RATOM_BRICK_H__



/**
*	\brief Reprezentuje prostopadloscian. Sciany prostopadloscianu 
*		rownolegle do pwoerzchni wyznzconych przez osie ukladu wspolrzednych.
*
*	\author Zbigniew Romanowski, ROMZ
*
*	\version 07-Lip-2007 [romz]
*
*/

class Brick
{
public:
	Brick(void);
	~Brick(void);

	void HalfXup(void);
	void HalfXdown(void);

	void HalfYup(void);
	void HalfYdown(void);

	void HalfZup(void);
	void HalfZdown(void);

	bool operator<(const Brick& r) const;

	// Wykorzystywane przez kolejke priorytetowa
	double Pri() const { return m_errAbs; }

public:
	// Wspolrzedne prostokata
	double m_x0;

	// Wspolrzedne prostokata
	double m_x1;

	// Wspolrzedne prostokata
	double m_y0;

	// Wspolrzedne prostokata
	double m_y1;

	// Wspolrzedne prostokata
	double m_z0;

	// Wspolrzedne prostokata
	double m_z1;

	// Absolutny Blad kwadratury
	double m_errAbs;

//	// WZGLEDNY blad kwadratury
//	double m_errRel;

	// Wartosc kwadratury
	double m_quad;

	// Numer wspolrzednej o nazwiekszej roznicy dzielonej. 
	// $0$ oznacza os $X$. $1$ oznacza os $Y$. $2$ oznacza os $Z$.
	int m_divAx;

};

#endif

