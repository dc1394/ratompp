#ifndef __RATOM_EXCHS_H__
#define __RATOM_EXCHS_H__


/** \brief Slater approximation for echange and correlation.
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/



class ExchS //: public Xc
{
public:
	ExchS(void);
	virtual ~ExchS(void);

	virtual double E(double rhoa, double rhob, double gaa, double gab, double gbb) const;

	virtual double Vrhoa(double rhoa, double rhob, double gaa, double gab, double gbb) const;
	virtual double Vrhob(double rhoa, double rhob, double gaa, double gab, double gbb) const;

	virtual double Vgaa(double rhoa, double rhob, double gaa, double gab, double gbb) const;
	virtual double Vgbb(double rhoa, double rhob, double gaa, double gab, double gbb) const;
	virtual double Vgab(double rhoa, double rhob, double gaa, double gab, double gbb) const;


private:
	static const double m_alpha;
	static const double m_c;
};


#endif

