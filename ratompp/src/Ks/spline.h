#ifndef _SPLINE_H_
#define _SPLINE_H_

#ifdef _MSC_VER
	#pragma once
#endif

#include <vector>

#if !defined(__INTEL_COMPILER) || !defined(__GXX_EXPERIMENTAL_CXX0X__) || (_MSC_VER < 1800)
	#include <boost/noncopyable.hpp>
#endif

namespace Thomas_Fermi {
	class Spline
#if !defined(__INTEL_COMPILER) || !defined(__GXX_EXPERIMENTAL_CXX0X__) || (_MSC_VER < 1800)
		: private boost::noncopyable
#endif
	{
#if defined(__INTEL_COMPILER) || defined(__GXX_EXPERIMENTAL_CXX0X__) || (_MSC_VER >= 1800)
		Spline(const Spline &) = delete;
		Spline & operator=(const Spline &) = delete;
		Spline() = delete;
#endif

#if defined(__INTEL_COMPILER) || defined(__GXX_EXPERIMENTAL_CXX0X__)
		static constexpr double TINY = 1.0E-30;
#else
		static const double TINY;
#endif

		typedef std::vector<double> dvector;

		const std::size_t size;

		const dvector x_;
		const dvector y_;
		dvector y2_;
		
	public:
		Spline(const dvector & x, const dvector & y);
		double operator()(double x) const;
		double df_dx(double x) const;
	};
}

#endif	// _SPLINE_H_
