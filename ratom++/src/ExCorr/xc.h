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
#include <functional>       // for std::function
#include <string>           // 

// May 24th, 2014 Modified by dc1394
namespace excorr {
    template <util::Spin S>
    class Xc final {
        std::function<double(double)> const V_;
        std::function<double(double)> const E_;
        std::function<std::string()> const Name_;

        //std::function<std::pair<double, double>(double)> const rhoTilde_;
        //std::function<double(double)> rhoTildeDeriv_
        //std::function<double(double)> rhoTildeLapl_;

        Xc(Xc const &) = delete;
        Xc & operator=(const Xc &) = delete;

    public:
        // Constructor
        template <typename T>
        Xc(const T & obj)
            :   V_([this, obj](double r) { return obj.xc_vxc<S>(r); }),
                E_([&obj](double r) { return obj.xc_exc<S>(r); }),
                Name_([&obj]() { return obj.Name(); })
        {}
        // Destructor
        ~Xc() = default;

        double V(double r) const
        {
            return V_(r);
        }
        double E(double r) const
        {
            return E_(r);
        }
        double EdiffV(double r) const
        {
            return E_(r) - V_(r);
        }

        std::string Name() const
        {
            return Name_();
        }

        // March 8th, 2014	added by dc1394
        //double my_xc_gga_exc(double r) const;
        //double my_xc_gga_vxc(double r) const;
    };
}

#endif

