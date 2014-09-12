#ifndef __RATOM_CORRPBE_H__
#define __RATOM_CORRPBE_H__


/** PBE correlation.
*
* \author @dc1394
*
* March 8th, 2014
*/

#include "xc.h"

class CorrPbe : public Xc
{
public:
	CorrPbe(std::function<double(double)> rhoTilde,
		    std::function<double(double)> rhoTildeDeriv,
			std::function<double(double)> rhoTildeLapl);
	virtual ~CorrPbe(void);

	virtual double V(double r) const;
	virtual double E(double r) const;

	virtual const char* Name() const
	{
		// April 4th, 2014 Modified by dc1394
		//return "PBE Correlation";
		return pxcfunc_->info->name;
	}
};

#endif	// __RATOM_CORRPBE_H__
