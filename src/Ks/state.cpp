#include "stdafx.h"
#include "state.h"
#include <cstdint>      // for std::int32_t
#include <sstream>      // for std::istringstream
#include <stdexcept>    // for std::invalid_argument

namespace ks {
    //
    // Constructor
    //
    State::State(std::string const & name) : m_name(name)
    {
        std::int32_t ee, nn;
        char ll;

        assert(!name.empty());
        std::istringstream iss(name);
        iss >> nn >> ll >> ee;

        //if (std::sscanf(name.c_str(), "%d%c%d", &nn, &ll, &ee) == EOF)
        if (iss.fail() ||   // すべて正常に読み込めたかチェック
            !iss.eof())     // 文字列の最後まで処理したかチェックを行う。
        {
            throw std::invalid_argument("name (str) is incorrect!");
        }
        
        // printf("nn = %d, ll = %c, ee = %d\n", nn, ll, ee);

        m_l = GetL(ll);
        m_n = static_cast<std::size_t> (nn - m_l - 1);
        m_eigVal = 0.0;
        //auto const occ = static_cast<std::std::size_t>(ee);
        m_occ = static_cast<std::size_t>(ee);

        // printf("name = %s, n = %lu, l = %lu, occ = %lu\n", name, (unsigned long)m_n, (unsigned long)m_l, (unsigned long)m_occ);
        // fflush(stdout);
    }

    std::size_t State::GetL(char ll) const
    {
        std::size_t ret = 0;
        switch (ll)
        {
        case 's': ret = 0; break;
        case 'p': ret = 1; break;
        case 'd': ret = 2; break;
        case 'f': ret = 3; break;
        default: assert(0); break;
        }
        return ret;
    }

    //
    // Comparision operator
    // States are sorder according to its eigenvalues
    //
    bool State::operator< (const State& s) const
    {
        return (m_eigVal < s.m_eigVal);
    }
}

