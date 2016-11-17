#include "exchhf.h"
#include <cmath>

namespace excorr {
    ExchHf::ExchHf(std::function<double(double)> const & Vh, double Z)
        :   Vh_(Vh),
            Z_(Z)
    {
    }
    
    double ExchHf::xc_exc(double r) const
    {
        if (std::fabs(Z_ - 1.0) < ZERO) {
            return -0.5 * Vh_(r);
        }
        else if (std::fabs(Z_ - 2.0) < ZERO) {
            return -0.25 * Vh_(r);
        }
        else {
            return NAN;
        }
    }

    //
    // Potencial
    //
    /*double ExchHf::V(double r) const
    {
        if (std::fabs(Z_ - 1.0) < ZERO)
            return -Vh_(r);
        else if (std::fabs(Z_ - 2.0) < ZERO)
            return -0.5 * Vh_(r);
        else
            throw std::runtime_error("–¢ŽÀ‘•");
    }*/
        
}