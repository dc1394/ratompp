//#ifndef __RATOM_EXCHSLATER_H__
//#define __RATOM_EXCHSLATER_H__
//
//
///** \brief Slater approximation.
//*
//* \author Zbigniew Romanowski [ROMZ@wp.pl]
//*
//*/
//
//// May 25th, 2014 Modified by dc1394
////#include "xc.h"
//#include "excorrlda.h"
//
//
//// May 25th, 2014 Added by dc1394
//namespace excorr {
//    // May 25th, 2014 Modified by dc1394
//    class ExchSlater : public ExCorrLDA
//    {
//        ExchSlater(ExchSlater const &) = delete;
//        ExchSlater & operator=(const ExchSlater &) = delete;
//
//    public:
//        // March 7th, 2014	Modified by dc1394
//        //ExchSlater(void);
//        ExchSlater(std::function<double(double)> rhoTilde);
//        virtual ~ExchSlater(void);
//
//        //virtual double V(double rho, double gRho) const;
//        // March 7th, 2014	Modified by dc1394
//        virtual double V(double r) const;
//        //virtual double E(double rho, double gRho) const;
//        // March 7th, 2014	Modified by dc1394
//        virtual double E(double r) const;
//        //virtual double EdiffV(double rho, double gRho) const;
//
//        virtual const char* Name() const
//        {
//            // April 4th, 2014 Modified by dc1394
//            //return "slater";
//            return pxcfunc_->info->name;
//        }
//
//
//    private:
//        double m_c;
//    };
//}
//
//#endif
//
