#include "Spline.h"
#include <cstdint>
#include <boost/cast.hpp>

namespace Thomas_Fermi {
#if !defined(__INTEL_COMPILER) && !defined(__GXX_EXPERIMENTAL_CXX0X__)
	const double Spline::TINY = 1.0E-30;
#endif	
	Spline::Spline(const dvector & x, const dvector & y)
		:	size(x.size()), x_(x), y_(y), y2_(size, 0.0)
	{
		if (x_.size() != y_.size())
			BOOST_ASSERT(!"xとyのサイズが異なります！");

		dvector u(size - 1, 0.0);

		y2_[0] = 0.0;

		// 3重対角のアルゴリズムの分解ループ
		// y2, uに分解済みの因子を仮に保管する
		for (std::size_t i = 1; i < size - 1; i++) {
			const double sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
			const double p = sig * y2_[i - 1] + 2.0;
			y2_[i] = (sig - 1.0) / p;
			u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) -
				   (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
			u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) -
					sig * u[i - 1]) / p;
		}

		double qn_1 = 0.0, un_1 = 0.0;

		y2_[size - 1] = (un_1 - qn_1 * u[size - 2]) /
						(qn_1 * y2_[size - 2] + 1.0);

		for (std::int32_t k = boost::numeric_cast<std::int32_t>(size - 2); k >= 0; k--)	// 3重対角アルゴリズムの後退代入ループ
			y2_[k] = y2_[k] * y2_[k + 1] + u[k];
	}

	double Spline::operator()(double x) const
	{
		std::size_t klo = 0;												// 表の中の正しい位置を二分探索で求める
		std::size_t khi = size - 1;

		while (khi - klo > 1) {
			const std::size_t k = (khi + klo) >> 1;

			if (x_[k] > x)
				khi = k;
			else
				klo = k;
		}

		const double h = x_[khi] - x_[klo];
		if (std::fabs(h) < TINY)
			BOOST_ASSERT(!"Bad xa input to routine splint()");				// xaの各値は異ならなければならない

		const double a = (x_[khi] - x) / h;
		const double b = (x - x_[klo]) / h;									// 3次スプラインの値を評価

		const double y = a * y_[klo] + b * y_[khi] +
						 ((a * a * a - a) * y2_[klo] +
						  (b * b * b - b) * y2_[khi]) *
						 (h * h) / 6.0;
		return y;
	}

	double Spline::df_dx(double x) const
	{
		std::uint32_t klo = 0;												// 表の中の正しい位置を二分探索で求める
		std::uint32_t khi = static_cast<std::uint32_t>(size - 1);

		while (khi - klo > 1) {
			const std::uint32_t k = static_cast<std::uint32_t>((khi + klo) >> 1);

			if (x_[k] > x)
				khi = k;
			else
				klo = k;
		}

		const double h = x_[khi] - x_[klo];
		if (std::fabs(h) < TINY)
			BOOST_ASSERT(!"Bad xa input to routine splint()");				// xaの各値は異ならなければならない

		const double a = (x_[khi] - x) / h;
		const double b = (x - x_[klo]) / h;									// 3次スプラインの値を評価

		const double dy_dx = (y_[khi] - y_[klo]) / (x_[khi] - x_[klo]) +
							 (- (3.0 * a * a - 1.0) * y2_[klo] + 
							  (3.0 * b * b - 1.0) * y2_[khi]) *
							 h / 6.0;
		return dy_dx;
	}
}
