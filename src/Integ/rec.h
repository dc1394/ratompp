#ifndef __RATOM_REC_H__
#define __RATOM_REC_H__



class Rec
{
public:
	Rec(void);
	~Rec(void);

	void HalfLeft();
	void HalfRight();

	void HalfTop();
	void HalfBottom();

//	void Quad(int i);

	bool operator<(const Rec& r) const;

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

	// BEZWZGLEDNY (Absolutny) blad kwadratury
	double m_errAbs;

//	// WZGLEDNY blad kwadratury
//	double m_errRel;

	// Wartosc kwadratury
	double m_quad;

	// Numer wspolrzednej o nazwiekszej roznicy dzielonej.
	// $0$ oznacza os $X$. $1$ oznacza os $Y$.
	int m_divAx;

};

#endif

