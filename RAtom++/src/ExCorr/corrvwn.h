#ifndef __RATOM_CORRVWN_H__
#define __RATOM_CORRVWN_H__


/** \brief VWN approximation.
*
*	\note 1) S. J. Vosko, L. Wilk, M. Nusair "Accurate spin dependent electron liquid correlation
*		energies for local spin dnsity calculations: A critical analysis",
*		Can. J. Phys. vol. 58, 1200-1211 (1980)
*		2) Vosko, Wilk, Phys. Rev. B, vol. 22, 3812 (1980)
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/

// March 18th, 2014 Modified by dc1394
#include "xc.h"



class CorrVwn : public Xc
{
public:
	//CorrVwn(void);
	CorrVwn(std::function<double(double)> rhoTilde);
	virtual ~CorrVwn(void);

	virtual double V(double rho, double gRho) const;
	virtual double E(double rho, double gRho) const;
	virtual double EdiffV(double rho, double gRho) const;
	
	// March 18th, 2014 Added by dc1394
	virtual double V(double r) const;
	virtual double E(double r) const;

	virtual const char* Name() const
	{
		// April 4th, 2014 Modified by dc1394
		//return "vwn";
		return pxcfunc_->info->name;
	}

private:
	virtual void Help(double rho, double* ec, double* vc) const;
};


#endif

