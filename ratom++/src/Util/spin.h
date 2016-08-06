#ifndef __RATOM_SPIN_H__
#define __RATOM_SPIN_H__

#include <cstdint>  // for std::int32_t

namespace util {
    enum class Spin : std::int32_t {
        Alpha = 1,
        Beta = 2
    };
}

#endif  // __RATOM_SPIN_H__
