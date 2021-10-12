#pragma once
#include <iostream>
#include "Matrix.h"
#include "MatrixExceptions.h"

enum class SIDE {LEFT, RIGHT};

template<typename TYPE>
class HouseHolderReflectors
{
protected:
	Matrix<TYPE> const& _mat;
	Vector<TYPE> _tau;
public:
	HouseHolderReflectors(Matrix<TYPE> const& mat) : _mat(mat), _tau(std::min(mat.shape().i - 1, mat.shape().j)) {}
	void emplaceScalarFactor(TYPE tau, size_t index) { _tau(index) = tau; }
	Shape shape() const { return { _mat.shape().i, _mat.shape().i }; }

	template<typename TYPE, typename E1>
	friend Matrix<TYPE> operator*(HouseHolderReflectors<TYPE> const& Q, MatrixExpr<TYPE, E1> const& expr);

	template<typename TYPE, typename E1>
	friend Matrix<TYPE> operator*(MatrixExpr<TYPE, E1> const& expr, HouseHolderReflectors<TYPE> const& Q);
};

template<typename TYPE, typename E1>
Matrix<TYPE> operator*(HouseHolderReflectors<TYPE> const& Q, MatrixExpr<TYPE, E1> const& expr)
{
	Matrix<TYPE> copy = expr;//Force evaluation and copy

	size_t m = Q._mat.shape().i;
	if (m != copy.shape().i)
		throw IncompatibleShapesException(Q._mat.shape(), copy.shape(), OPERATION::MULT);

	size_t s = std::min(Q._mat.shape().i - 1, Q._mat.shape().j);

	for (size_t k = s; k > 0; --k)
	{
		size_t km = k - 1;
		auto col_v = column(Q._mat, km);
		auto v = subVector(col_v, k, m - 1);

		TYPE tau = std::abs(Q._tau(km));
		TYPE sign = Q._tau(km) / tau;

		for (size_t kk = 0; kk < copy.shape().j; ++kk)
		{
			auto col_A = column(copy, kk);
			auto vec_A = subVector(col_A, m - 1 - v.size(), m - 1);

			auto sub_vec_A = subVector(vec_A, 1, vec_A.size() - 1);

			TYPE proj = dot(sub_vec_A, v) + vec_A(0); //(""" v(0) = 1 """)

			vec_A(0) = sign * (proj / tau - vec_A(0));
			for (size_t kkk = 0; kkk < sub_vec_A.size(); ++kkk)
			{
				sub_vec_A(kkk) = sign * (proj / tau * v(kkk) - sub_vec_A(kkk));
			}
		}
	}

	return copy;
}

template<typename TYPE, typename E1>
Matrix<TYPE> operator*(MatrixExpr<TYPE, E1> const& expr, HouseHolderReflectors<TYPE> const& Q)
{
	Matrix<TYPE> copy = expr;//Force evaluation and copy

	size_t m = Q._mat.shape().i;
	if (m != copy.shape().i)
		throw IncompatibleShapesException(copy.shape(), Q._mat.shape(), OPERATION::MULT);

	size_t s = std::min(Q._mat.shape().i - 1, Q._mat.shape().j);

	for (size_t k = 0; k < s; ++k)
	{
		auto col_v = column(Q._mat, k);
		auto v = subVector(col_v, k + 1, m - 1);

		TYPE tau = std::abs(Q._tau(k));
		TYPE sign = Q._tau(k) / tau;

		for (size_t kk = 0; kk < copy.shape().i; ++kk)
		{
			auto row_A = row(copy, kk);
			auto vec_A = subVector(row_A, k, copy.shape().i - 1);

			auto sub_vec_A = subVector(vec_A, 1, vec_A.size() - 1);

			TYPE proj = dot(v, sub_vec_A) + vec_A(0); //(""" v(0) = 1 """)

			vec_A(0) = sign * (proj / tau - vec_A(0));
			for (size_t kkk = 0; kkk < sub_vec_A.size(); ++kkk)
			{
				sub_vec_A(kkk) = sign * (proj / tau * v(kkk) - sub_vec_A(kkk));
			}
		}
	}

	return copy;
}

template<typename TYPE>
class HouseHolderReflectors<std::complex<TYPE>>
{
protected:
	Matrix<std::complex<TYPE>> const& _mat;
	Vector<std::complex<TYPE>> _tau;
public:
	HouseHolderReflectors(Matrix<std::complex<TYPE>> const& mat) : _mat(mat), _tau(std::min(mat.shape().i - 1, mat.shape().j)) {}
	void emplaceScalarFactor(std::complex<TYPE> tau, size_t index) { _tau(index) = tau; }
	Shape shape() const { return { _mat.shape().i, _mat.shape().i }; }

	template<typename TYPE, typename E1>
	friend Matrix<std::complex<TYPE>> operator*(HouseHolderReflectors<std::complex<TYPE>> const& Q, MatrixExpr<std::complex<TYPE>, E1> const& expr);

	template<typename TYPE, typename E1>
	friend Matrix<std::complex<TYPE>> operator*(MatrixExpr<std::complex<TYPE>, E1> const& expr, HouseHolderReflectors<std::complex<TYPE>> const& Q);
};

template<typename TYPE, typename E1>
Matrix<std::complex<TYPE>> operator*(HouseHolderReflectors<std::complex<TYPE>> const& Q, MatrixExpr<std::complex<TYPE>, E1> const& expr)
{
	Matrix<std::complex<TYPE>> copy = expr;//Force evaluation and copy

	size_t m = Q._mat.shape().i;
	if (m != copy.shape().i)
		throw IncompatibleShapesException(Q._mat.shape(), copy.shape(), OPERATION::MULT);

	size_t s = std::min(Q._mat.shape().i - 1, Q._mat.shape().j);

	for (size_t k = s; k > 0; --k)
	{
		size_t km = k - 1;
		auto col_v = column(Q._mat, km);
		auto v = subVector(col_v, k, m - 1);

		TYPE tau = std::abs(Q._tau(km));
		std::complex<TYPE> sign = Q._tau(km) / tau;

		for (size_t kk = 0; kk < copy.shape().j; ++kk)
		{
			auto col_A = column(copy, kk);
			auto vec_A = subVector(col_A, m - 1 - v.size(), m - 1);

			auto sub_vec_A = subVector(vec_A, 1, vec_A.size() - 1);

			std::complex<TYPE> proj = dot(sub_vec_A, v) + vec_A(0); //(""" v(0) = 1 """)

			vec_A(0) = sign * (proj / tau - vec_A(0));
			for (size_t kkk = 0; kkk < sub_vec_A.size(); ++kkk)
			{
				sub_vec_A(kkk) = sign * (proj / tau * v(kkk) - sub_vec_A(kkk));
			}
		}
	}

	return copy;
}

template<typename TYPE, typename E1>
Matrix<std::complex<TYPE>> operator*(MatrixExpr<std::complex<TYPE>, E1> const& expr, HouseHolderReflectors<std::complex<TYPE>> const& Q)
{
	Matrix<std::complex<TYPE>> copy = expr;//Force evaluation and copy

	size_t m = Q._mat.shape().i;
	if (m != copy.shape().i)
		throw IncompatibleShapesException(copy.shape(), Q._mat.shape(), OPERATION::MULT);

	size_t s = std::min(Q._mat.shape().i - 1, Q._mat.shape().j);

	for (size_t k = 0; k < s; ++k)
	{
		auto col_v = column(Q._mat, k);
		auto v = subVector(col_v, k + 1, m - 1);

		TYPE tau = std::abs(Q._tau(k));
		std::complex<TYPE> sign = Q._tau(k) / tau;

		for (size_t kk = 0; kk < copy.shape().i; ++kk)
		{
			auto row_A = row(copy, kk);
			auto vec_A = subVector(row_A, k, copy.shape().i - 1);

			auto sub_vec_A = subVector(vec_A, 1, vec_A.size() - 1);

			std::complex<TYPE> proj = dot(v, conj(sub_vec_A)) + vec_A(0); //(""" v(0) = 1 """)

			vec_A(0) = sign * (proj / tau - vec_A(0));
			for (size_t kkk = 0; kkk < sub_vec_A.size(); ++kkk)
			{
				sub_vec_A(kkk) = sign * (proj / tau * std::conj(v(kkk)) - sub_vec_A(kkk));
			}
		}
	}

	return copy;
}

template<typename TYPE>
class HouseHolderQR
{
public:
	HouseHolderQR(Matrix<TYPE>& mat) : 
		_mat(mat), 
		_s(std::min(_mat.shape().i - 1, _mat.shape().j)), 
		_Q(_mat),
		_R(_mat)
	{
		size_t m = _mat.shape().i;
		size_t n = _mat.shape().j;

		if (_s >= 1)
		{
			for (size_t k = 0; k < _s; ++k)
			{
				auto first_column = column(_mat, k);
				auto a1 = subVector(first_column, k, m-1);
				TYPE normA1sq = dot(a1, a1);
				TYPE normA1 = std::sqrt(normA1sq);
				TYPE signA1 = a1(0) / std::abs(a1(0));
				Vector<TYPE> v = a1;
				v(0) += signA1 * normA1;

				TYPE tempV0 = v(0);
				v /= tempV0;//Store "normalized" vector such that v(0)=1

				TYPE nVsq = dot(v, v);

				_Q.emplaceScalarFactor(signA1 * nVsq / TYPE(2), k);

				for (size_t kk = k+1; kk < n; ++kk)
				{
					auto kk_column = column(_mat, kk);
					auto kk_vec = subVector(kk_column, k, m-1);
					
					TYPE proj = dot(kk_vec, v) / nVsq;

					for (size_t kkk = 0; kkk < kk_vec.size(); ++kkk)
					{
						kk_vec(kkk) = signA1 * (TYPE(2) * proj * v(kkk)-kk_vec(kkk));
					}
				}

				a1(0) = normA1;

				for (size_t kkk = 1; kkk < a1.size(); ++kkk)
				{
					a1(kkk) = v(kkk);
				}
			}
		}		
	}
protected:
	Matrix<TYPE>& _mat;
	size_t _s;
public:
	HouseHolderReflectors<TYPE> _Q;
	UpperTriangular<TYPE, Matrix<TYPE>> _R;
	
};

template<typename TYPE>
class HouseHolderQR<std::complex<TYPE>>
{
public:
	HouseHolderQR(Matrix<std::complex<TYPE>>& mat) : 
		_mat(mat), 
		_s(std::min(_mat.shape().i - 1, _mat.shape().j)), 
		_Q(_mat), 
		_R(_mat)
	{
		size_t m = _mat.shape().i;
		size_t n = _mat.shape().j;

		if (_s >= 1)
		{
			for (size_t k = 0; k < _s; ++k)
			{
				auto first_column = column(_mat, k);
				auto a1 = subVector(first_column, k, m - 1);
				std::complex<TYPE> normA1sq = dot(a1, a1);
				std::complex<TYPE> normA1 = std::sqrt(normA1sq);
				std::complex<TYPE> signA1 = a1(0) / std::abs(a1(0));
				Vector<std::complex<TYPE>> v = a1;
				v(0) += signA1 * normA1;

				std::complex<TYPE> tempV0 = v(0);
				v /= tempV0;//Store "normalized" vector such that v(0)=1

				std::complex<TYPE> nVsq = dot(v, v);

				_Q.emplaceScalarFactor(signA1 * nVsq / TYPE(2), k);

				for (size_t kk = k+1; kk < n; ++kk)
				{
					auto kk_column = column(_mat, kk);
					auto kk_vec = subVector(kk_column, k, m - 1);

					std::complex<TYPE> proj = dot(kk_vec, v) / nVsq;

					for (size_t kkk = 0; kkk < kk_vec.size(); ++kkk)
					{
						kk_vec(kkk) = std::conj(signA1) * (TYPE(2) * proj * v(kkk) - kk_vec(kkk));
					}
				}

				a1(0) = normA1;

				for (size_t kkk = 1; kkk < a1.size(); ++kkk)
				{
					a1(kkk) = v(kkk);
				}
			}
		}
	}
protected:
	Matrix<std::complex<TYPE>>& _mat;
	size_t _s;
public:
	HouseHolderReflectors<std::complex<TYPE>> _Q;
	UpperTriangular<std::complex<TYPE>, Matrix<std::complex<TYPE>>> _R;
};

