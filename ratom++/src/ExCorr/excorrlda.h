#ifndef __RATOM_EXCORR_LDA_H__
#define __RATOM_EXCORR_LDA_H__

/** \brief Represents any LDA Exchange-Correlation potential.
*
* \author @dc1394
* \date May 24th, 2014
*
*/

#include "../Util/isspin.h"
#include <array>
#include <memory>
#include <utility>
#include <functional>
#include <type_traits>
#include <cstdint>
#include <xc.h>

extern void* enabler;

namespace excorr {
    template <util::IsSpin T, typename X = void>
    using alpha_enabler = typename std::enable_if<T == util::IsSpin::Alpha, X>::type*&;
    template <util::IsSpin T, typename X = void>
    using beta_enabler = typename std::enable_if<T == util::IsSpin::Beta, X>::type*&;

    class ExCorrLDA final {
        ExCorrLDA(ExCorrLDA const &) = delete;
        ExCorrLDA & operator=(ExCorrLDA const &) = delete;

        const std::function<std::pair<double, double> (double)> rhoTilde_;
        const std::function<void(xc_func_type * xcfunc)> deleter_;
        const std::unique_ptr<xc_func_type, decltype(deleter_)> pxcfunc_;

        std::array<double, 2> xc_vxc_impl(double r) const;

    public:
        ExCorrLDA(std::function<std::pair<double, double>(double)> rhoTilde,
                  std::uint32_t xc_type);
        ~ExCorrLDA();

        double xc_exc(double r) const;
        template <util::IsSpin Spin, alpha_enabler<Spin> T = enabler>
        double xc_vxc(double r) const;
        template <util::IsSpin Spin, beta_enabler<Spin> T = enabler>
        double xc_vxc(double r) const;

        char const * Name() const;
    };
}

#endif  // __RATOM_EXCORR_LDA_H__
