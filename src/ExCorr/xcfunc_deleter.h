/*! \file xcfunc_deleter.h
    \brief xcfuncを解放する

    Copyright ©  2016 @dc1394 All Rights Reserved.
    This software is released under the GNU GPL v3.
*/

#ifndef __RATOM_XCFUNC_DELETER_H__
#define __RATOM_XCFUNC_DELETER_H__

#include <xc.h>                         // for xc_func_end
#include <boost/checked_delete.hpp>     // for boost::checked_delete

namespace excorr {
    //! A lambda expression.
    /*!
        xc_func_typeへのポインタを解放するラムダ式
        \param xcfunc xc_func_type へのポインタ
    */
    static auto const xcfunc_deleter = [](xc_func_type * xcfunc) {
        xc_func_end(xcfunc);
        boost::checked_delete(xcfunc);
    };
}

#endif  // __RATOM_XCFUNC_DELETER_H__
