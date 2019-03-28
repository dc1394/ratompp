#ifndef __RATOM_SHREAL_H__
#define __RATOM_SHREAL_H__


/** \brief Real spherical harmonic
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/


class ShReal
{
public:
	ShReal(void);
	~ShReal(void);
	double Get(int l, int m, double x, double y, double z, double r) const;
	double Get(int l, int m, double theta, double phi) const;

	double Clm(int l, int m) const;


	double Plm(int l, int m, double theta) const;
	double Qm(int m, double phi) const;

	void CheckNorm(void) const;

private:
	int Idx(int l, int m) const;
	double Prod(int la, int ma, int lb, int mb) const;

private:
	// Wspolczynniki normalizacji w kartezjanskim ukladzie wspolrzednych
	static double m_clm[16];

	// Wspolczynniki normalizacji w sferycznym ukladzie wspolrzednych
	static double m_slm[16];
};




 #endif
 
