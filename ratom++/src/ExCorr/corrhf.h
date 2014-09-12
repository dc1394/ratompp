#ifndef __RATOM_CORRHF_H__
#define __RATOM_CORRHF_H__


/** \brief HF Correlation
*
* \author @dc1394
* March 28th, 2014
*/


#include "xc.h"



class CorrHf : public Xc
{
public:
	CorrHf(void);
	virtual ~CorrHf(void);

	virtual double V(double r) const;
	virtual double E(double r) const;

	virtual const char* Name() const
	{
		return "hartree-fock";
	}
};

#endif	// __RATOM_CORRHF_H__

