#include "stdafx.h"
#include "exchslater.h"

// May 25th, 2014 Added by dc1394
namespace excorr {
    //
    // Constructor
    //
    // March 8th, 2014	Modified by dc1394
    //ExchSlater::ExchSlater(void) : m_c(-pow(1.5 / M_PI, 2. / 3.))
    ExchSlater::ExchSlater(std::function<double(double)> rhoTilde)
        : ExCorrLDA(std::move(rhoTilde))
    {
        // March 8th, 2014	Added by dc1394
        xc_func_init(pxcfunc_.get(), XC_LDA_X, XC_POLARIZED);
    }


    // March 8th, 2014	Modified by dc1394
    //
    // Destructor
    //
    ExchSlater::~ExchSlater(void)
    {
    }

    //
    // Potencial
    //
    //double ExchSlater::V(double rho, double /* gRho */) const
    //{
    //    if (rho == 0)
    //        return 0;

    //    return m_c / Rs(rho);
    //}


    ////
    //// Energy density
    ////
    //double ExchSlater::E(double rho, double /* gRho */) const
    //{
    //    if (rho == 0)
    //        return 0;

    //    return 0.75 * m_c / Rs(rho);
    //}


    ////
    //// Difference. Optimized version :-)
    ////
    //double ExchSlater::EdiffV(double rho, double /* gRho */) const
    //{
    //    if (rho == 0)
    //        return 0;

    //    return -0.25 * m_c / Rs(rho);
    //}
}

