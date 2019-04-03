#ifndef __RATOM_FUN1D_H__
#define __RATOM_FUN1D_H__


/** \brief Represents function from R into R.
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/
// March 16th, 2014 Modified by dc1394

namespace util {
    class Fun1D
    {
    public:
        Fun1D() = default;
        virtual ~Fun1D() = default;

        virtual double Get(double x) const = 0;

        // April 3rd, 2014 Added by dc1394
        virtual double GetDeriv(double) const
        {
            return 0.0;
        }
        virtual double Get2ndDeriv(double) const
        {
            return 0.0;
        }
    };
}

#endif

