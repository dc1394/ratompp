#ifndef __RATOM_XC_H__
#define __RATOM_XC_H__


/** \brief Represents any echang-coreelation potential.
*
* \author Zbigniew Romanowski [ROMZ@wp.pl]
*
*/
// March 7th, 2014	Modified by dc1394

// March 7th, 2014	Added by dc1394
#include "../Util/spin.h"
#include <functional>

// May 24th, 2014 Modified by dc1394
namespace excorr {
    template <util::Spin Spin>
    class Xc final {
        const std::function<double(double)> V_;
        const std::function<double(double)> E_;
        const std::function<char const *()> Name_;

        //std::function<double(double)> rhoTilde_;
        //std::function<double(double)> rhoTildeDeriv_
        //std::function<double(double)> rhoTildeLapl_;

        Xc(Xc const &) = delete;
        Xc & operator=(const Xc &) = delete;

    public:
        // Constructor
        template <typename T>
        Xc(const T & obj)
            : V_([&obj](double r) { return obj.xc_vxc<Spin>(r); }),
              E_([&obj](double r) { return obj.xc_exc(r); }),
              Name_([&obj]() { return obj.Name(); })
        {}
        // Destructor
        ~Xc(void) {}

        double V(double r) const;
        double E(double r) const;
        double EdiffV(double r) const;

        char const * Name() const;

        // March 8th, 2014	added by dc1394
        //double my_xc_gga_exc(double r) const;
        //double my_xc_gga_vxc(double r) const;
    };
}

#endif

