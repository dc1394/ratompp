/*! \file mkl_allocator.h
    \brief MKL用のアロケータクラスの宣言と実装

    Copyright ©  2015 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#ifndef _MKL_ALLOCATOR_H_
#define _MKL_ALLOCATOR_H_

#pragma once

#include <cstddef>  // for std::ptrdiff_t
#include <new>      // for placement new,std::bad_alloc
#include <limits>   // for std::numeric_limits
#include <mkl.h>    // for mkl_allocator

namespace util {
    template <typename T>
    //! A class.
    /*!
        MKL用のアロケータクラス
    */
    struct mkl_allocator {
        // #region 型エイリアス

        using size_type = std::size_t;
        using difference_type = std::ptrdiff_t;
        using pointer = T *;
        using const_pointer = T const *;
        using reference = T &;
        using const_reference = T const &;
        using value_type = T;
            
        // #endregion 型エイリアス

        // #region 構造体内構造体

        template <class U>
        //! A templater struct.
        /*!
            convert an allocator<T> to allocator<U>
        */
        struct rebind {
            typedef mkl_allocator<U> other;
        };
        
        // #endregion 構造体内構造体

        // #region コンストラクタ・デストラクタ
        
        //! A constructor.
        /*!
            デフォルトコンストラクタ
        */
        mkl_allocator() noexcept = default;
        
        //! A constructor.
        /*!
            デフォルトコピーコンストラクタ
        */
        mkl_allocator(mkl_allocator const &) noexcept = default;

        template <typename U>
        //! A constructor (template).
        /*!
            コピーコンストラクタ
        */
        mkl_allocator(mkl_allocator<U> const &) noexcept
        {
        }

        //! A destructor.
        /*!
            デフォルトデストラクタ
        */
        ~mkl_allocator() = default;

        // #endregion コンストラクタ・デストラクタ

        // #region メンバ関数

        //! A public member function.
        /*!
            メモリを割り当てる
            \param size 割り当てるメモリのサイズ
            \param 未使用
            \return 割り当てたメモリの先頭アドレス
        */
        pointer allocate(size_type size, const_pointer hint = 0)
        {
            auto const p = mkl_malloc(2 * size * sizeof(T), 64);
            
            if (!p) {
                throw std::bad_alloc();
            }

            return reinterpret_cast<pointer>(p);
        }

        //! A public member function.
        /*!
            割当て済みの領域を初期化する
            \param p 割り当て済みのメモリの先頭アドレス
            \param val 初期化に使う値
        */
        void construct(pointer p, T const & val)
        {
            new (reinterpret_cast<void *>(p)) T(val);
        }

        //! A public member function.
        /*!
            割り当て済みのメモリを解放する
            \param p 割り当て済みのメモリの先頭アドレス
            \param 未使用 
        */
        void deallocate(pointer p, size_type)
        {
            mkl_free(reinterpret_cast<void *>(p));
        }

        //! A public member function.
        /*!
            初期化済みの領域を削除する
            \param p 初期化済みの領域の先頭アドレス
        */
        void destroy(pointer p)
        {
            p->~T();
        }

        //! A public member function (const).
        /*!
            アドレスを返す
            \param value アドレスを返す対象の変数
            \return アドレス
        */
        pointer address(reference value) const
        {
            return &value;
        }

        //! A public member function (const).
        /*!
            アドレスを返す（constポインタバージョン）
            \param value アドレスを返す対象の変数
            \return アドレス
        */
        const_pointer address(const_reference value) const
        {
            return &value;
        }
        
        //! A public member function (const).
        /*!
            割当てることができる最大の要素数を返す
            \return 割当てることができる最大の要素数
        */
        size_type max_size() const noexcept
        {
            return (std::numeric_limits<std::size_t>::max()) / sizeof(T);
        }
    };

    template <typename T, typename U>
    //! A global function (template function).
    /*!
        operator==の宣言と実装
        \param 未使用
        \param 未使用
        \return 左辺と右辺が等値かどうか
    */
    inline bool operator ==(mkl_allocator<T> const &, const mkl_allocator<U>)
    {
        return true;
    }

    template <typename T, typename U>
    //! A global function (template function).
    /*!
        operator!=の宣言と実装
        \param 未使用
        \param 未使用
        \return 左辺と右辺が等値かどうか
    */
    inline bool operator !=(mkl_allocator<T> const &, const mkl_allocator<U>)
    {
        return false;
    }
}

#endif  // _MKL_ALLOCATOR_H_
