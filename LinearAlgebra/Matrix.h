#pragma once
#include <iostream>
#include <vector>
#include <complex>
#include <type_traits>
#include <functional>
#include "Traits.h"
#include "Vector.h"
#include "Shape.h"
#include "MatrixExceptions.h"
#include "VectorExceptions.h"


template<typename TYPE, typename E>
class MatrixExpr
{	
public:
	TYPE operator()(size_t i, size_t j) const { return static_cast<E const&>(*this)(i,j); }
	Shape shape() const { return static_cast<E const&>(*this).shape(); }
}; 

//Stored by rows
template<typename TYPE>
class Matrix : public MatrixExpr<TYPE, Matrix<TYPE>>
{
public:

	Matrix(size_t m , size_t n) : _shape(m, n), _mat(_shape.size()) {}

	template<typename E>
	Matrix(MatrixExpr<TYPE, E> const& expr) : _shape(expr.shape()), _mat(expr.shape().i* expr.shape().j)
	{
		size_t m = expr.shape().i;
		size_t n = expr.shape().j;

		for (size_t i = 0; i < m; ++i)
		{
			for (size_t j = 0; j < n ; ++j)
			{			
				_mat[i*n + j] = expr(i, j);
			}
		}
	}

	//Possible amélioration sur l'affectation de MatrixExpr en laissant la gestion de l'aliasing à l'utilisateur (B=B*A nécessite seulement de copier la première ligne de B en temporaire, même chose avec B=A*B sauf qu'il faut copier les colonnes)
	//Ici quoiqu'il arrive, l'expression B=B*A ou B=A*B se déroule comme suit : (Matrix) M=B*A ou M=A*B puis, B=M (operator=) donc avec une copie globale
	Matrix<TYPE>& operator=(Matrix<TYPE> const& mat)
	{ 
		size_t m = mat.shape().i;
		size_t n = mat.shape().j;

		_mat.resize(m*n);

		//if expr depends on mat, then we change expr when changing mat
		for (size_t i = 0; i < m; ++i)
		{
			for (size_t j = 0; j < n; ++j)
			{
				_mat[i * n + j] = mat(i, j);
			}
		}		
		return *this; 
	}

	TYPE operator()(size_t i, size_t j) const { return _mat[i * _shape.j + j]; }
	TYPE& operator()(size_t i, size_t j) { return _mat[i * _shape.j + j]; }

	Shape shape() const { return _shape; }
protected:
	Shape _shape;
	std::vector<TYPE> _mat;	

	template<typename TYPE>
	friend std::ostream& operator<<(std::ostream& os, Matrix<TYPE> const& mat);
};

template<typename TYPE>
std::ostream& operator<<(std::ostream& os, Matrix<TYPE> const& mat)
{
	if (mat._shape.notNull())
	{
		for (size_t i = 0; i < mat._shape.i - 1; ++i)
		{
			for (size_t j = 0; j < mat._shape.j; ++j)
			{
				os << mat(i, j) << ',';
			}
			os << std::endl;
		}
		for (size_t j = 0; j < mat._shape.j-1; ++j)
		{
			os << mat(mat._shape.i - 1, j) << ',';
		}
		os << mat(mat._shape.i - 1, mat._shape.j - 1);
	}	
	return os;
}

template<typename TYPE, typename E1, typename E2>
class MatrixSum : public MatrixExpr<TYPE, MatrixSum<TYPE, E1, E2>>
{
public:
	MatrixSum(E1 const& e1, E2 const& e2) : _e1(e1), _e2(e2)
	{
		if (_e1.shape() != _e2.shape())
			throw IncompatibleShapesException(_e1.shape(), _e2.shape(), OPERATION::ADD);
	}

	TYPE operator()(size_t i, size_t j) const { return _e1(i,j) + _e2(i,j); }
	Shape shape() const { return _e1.shape(); }
protected:
	E1 const& _e1;
	E2 const& _e2;
};

template<typename TYPE, typename E1, typename E2>
MatrixSum<TYPE, E1, E2> operator+(MatrixExpr<TYPE, E1> const& e1, MatrixExpr<TYPE, E2> const& e2)
{
	return MatrixSum<TYPE, E1, E2>(*static_cast<const E1*>(&e1), *static_cast<const E2*>(&e2));
}

template<typename TYPE, typename E1, typename E2>
class MatrixDiff : public MatrixExpr<TYPE, MatrixDiff<TYPE, E1, E2>>
{
public:
	MatrixDiff(E1 const& e1, E2 const& e2) : _e1(e1), _e2(e2)
	{
		if (_e1.shape() != _e2.shape())
			throw IncompatibleShapesException(_e1.shape(), _e2.shape(), OPERATION::SUB);
	}

	TYPE operator()(size_t i, size_t j) const { return _e1(i, j) - _e2(i, j); }
	Shape shape() const { return _e1.shape(); }
protected:
	E1 const& _e1;
	E2 const& _e2;
};

template<typename TYPE, typename E1, typename E2>
MatrixDiff<TYPE, E1, E2> operator-(MatrixExpr<TYPE, E1> const& e1, MatrixExpr<TYPE, E2> const& e2)
{
	return MatrixDiff<TYPE, E1, E2>(*static_cast<const E1*>(&e1), *static_cast<const E2*>(&e2));
}

template<typename TYPE, typename E1, typename E2>
class MatrixProd : public MatrixExpr<TYPE, MatrixProd<TYPE, E1, E2>>
{
public:
	MatrixProd(E1 const& e1, E2 const& e2) : _e1(e1), _e2(e2)
	{
		if (_e1.shape().j != _e2.shape().i)
			throw IncompatibleShapesException(_e1.shape(), _e2.shape(), OPERATION::MULT);
	}

	TYPE operator()(size_t i, size_t j) const {
		TYPE s{};
		for (size_t k = 0; k < _e1.shape().j; ++k)
		{
			s += _e1(i, k) * _e2(k, j);
		}
		return s;
	}
	Shape shape() const { return { _e1.shape().i, _e2.shape().j }; }
protected:
	E1 const& _e1;
	E2 const& _e2;
};

template<typename TYPE, typename E1, typename E2>
MatrixProd<TYPE, E1, E2> operator*(MatrixExpr<TYPE, E1> const& e1, MatrixExpr<TYPE, E2> const& e2)
{
	return MatrixProd<TYPE, E1, E2>(*static_cast<const E1*>(&e1), *static_cast<const E2*>(&e2));
}

template<typename TYPE, typename E1>
class MatrixMult : public MatrixExpr<TYPE, MatrixMult<TYPE, E1>>
{
public:
	MatrixMult(E1 const& e1, TYPE const& t) : _e1(e1), _t(t)
	{
	}

	TYPE operator()(size_t i, size_t j) const {
		return _t*_e1(i,j);
	}
	Shape shape() const { return { _e1.shape().i, _e1.shape().j }; }
protected:
	E1 const& _e1;
	TYPE const& _t;
};

template<typename TYPE, typename E1>
MatrixMult<TYPE, E1> operator*(MatrixExpr<TYPE, E1> const& e1, TYPE const& t)
{
	return MatrixMult<TYPE, E1>(*static_cast<const E1*>(&e1), t);
}

template<typename TYPE, typename E1>
MatrixMult<TYPE, E1> operator*(TYPE const& t, MatrixExpr<TYPE, E1> const& e1)
{
	return MatrixMult<TYPE, E1>(*static_cast<const E1*>(&e1), t);
}

template<typename TYPE, typename E1>
	class MatrixConj : public MatrixExpr<TYPE, MatrixConj<TYPE, E1>>
{
public:
	MatrixConj(E1 const& e1) : _e1(e1)
	{
	}

	TYPE operator()(size_t i, size_t j) const {
		return std::conj(_e1(i, j));
	}
	Shape shape() const { return { _e1.shape().i, _e1.shape().j }; }
protected:
	E1 const& _e1;
};

template<typename TYPE, typename E1, typename std::enable_if < is_complex<TYPE>{}>::type>
MatrixConj<TYPE, E1> conj(MatrixExpr<TYPE, E1> const& e1)
{
	return MatrixConj<TYPE, E1>(*static_cast<const E1*>(&e1));
}

template<typename TYPE, typename E1>
class MatrixTranspose : public MatrixExpr<TYPE, MatrixTranspose<TYPE, E1>>
{
public:
	MatrixTranspose(E1 const& e1) : _e1(e1)
	{
	}

	TYPE operator()(size_t i, size_t j) const {
		return _e1(j, i);
	}
	Shape shape() const { return { _e1.shape().j, _e1.shape().i }; }
protected:
	E1 const& _e1;
};

template<typename TYPE, typename E1>
MatrixTranspose<TYPE, E1> transpose(MatrixExpr<TYPE, E1> const& e1)
{
	return MatrixTranspose<TYPE, E1>(*static_cast<const E1*>(&e1));
}

template<typename TYPE, typename E1>
class MatrixAdjoint : public MatrixExpr<TYPE , MatrixAdjoint<TYPE, E1>>
{
public:
	MatrixAdjoint(E1 const& e1) : _e1(e1)
	{
	}

	TYPE operator()(size_t i, size_t j) const {
		return std::conj(_e1(j, i));
	}
	Shape shape() const { return { _e1.shape().j, _e1.shape().i }; }
protected:
	E1 const& _e1;
};

template<typename TYPE, typename E1>
MatrixTranspose<TYPE, E1> adjoint(MatrixExpr<TYPE, E1> const& e1)
{
	return MatrixTranspose<TYPE, E1>(*static_cast<const E1*>(&e1));
}

template<typename TYPE, typename E1>
MatrixAdjoint<std::complex<TYPE>, E1> adjoint(MatrixExpr<std::complex<TYPE>, E1> const& e1)
{
	return MatrixAdjoint<std::complex<TYPE>, E1>(*static_cast<const E1*>(&e1));
}

template<typename TYPE, typename E1>
class SubMatrix : public MatrixExpr<TYPE, SubMatrix<TYPE, E1>>
{
public:
	SubMatrix(E1 & e1, size_t i, size_t j, size_t m, size_t n) : _e1(e1), _i(i), _j(j), _m(m), _n(n) {}

	TYPE operator()(size_t i, size_t j) const {
		return _e1(_i + i, _j + j);
	}

	template<typename = std::enable_if_t<has_ref_access_operator<E1, TYPE>::value>>
	TYPE& operator()(size_t i, size_t j) {
		return _e1(_i + i, _j + j);
	}

	Shape shape() const { return { _m, _n }; }

	template<typename E2>
	void operator+=(MatrixExpr<TYPE, E2> const& mat)
	{
		if (_m != mat.shape().i || _n != mat.shape().j)
			throw IncompatibleShapesException(Shape(_m, _n), mat.shape(), OPERATION::ADD);
		for (size_t i = 0; i < _m; ++i)
		{
			for (size_t j = 0; j < _n; ++j)
			{
				_e1(_i + i, _j + j) += mat(i, j);
			}
		}
	}

	template<typename E2>
	void operator-=(MatrixExpr<TYPE, E2> const& mat)
	{
		if(_m != mat.shape().i ||_n != mat.shape().j)
			throw IncompatibleShapesException(Shape(_m, _n), mat.shape(), OPERATION::SUB);
		for (size_t i = 0; i < _m; ++i)
		{
			for (size_t j = 0; j < _n; ++j)
			{
				_e1(_i + i, _j + j) -= mat(i, j);
			}
		}
	}
protected:
	E1 & _e1;
	size_t _i;
	size_t _j;
	size_t _m;
	size_t _n;
};

template<typename TYPE, typename E1>
SubMatrix<TYPE, E1> subMatrix(MatrixExpr<TYPE, E1> & mat, size_t i, size_t j, size_t m, size_t n)
{
	return SubMatrix<TYPE, E1>(static_cast<E1&>(mat), i, j, m, n);
}

template<typename TYPE, typename E1>
class ConstSubMatrix : public MatrixExpr<TYPE, ConstSubMatrix<TYPE, E1>>
{
public:
	ConstSubMatrix(E1 const& e1, size_t i, size_t j, size_t m, size_t n) : _e1(e1), _i(i), _j(j), _m(m), _n(n) {}

	TYPE operator()(size_t i, size_t j) const {
		return _e1(_i + i, _j + j);
	}

	Shape shape() const { return { _m, _n }; }

protected:
	E1 const& _e1;
	size_t _i;
	size_t _j;
	size_t _m;
	size_t _n;
};

template<typename TYPE, typename E1>
ConstSubMatrix<TYPE, E1> subMatrix(MatrixExpr<TYPE, E1> const& mat, size_t i, size_t j, size_t m, size_t n)
{
	return SubMatrix<TYPE, E1>(static_cast<E1 const&>(mat), i, j, m, n);
}

template<typename TYPE, typename E1>
class MatrixColumn : public VectorExpr<TYPE, MatrixColumn<TYPE, E1>>
{
public:
	MatrixColumn(E1 & e1, size_t j) : _e1(e1), _j(j) {
	}

	TYPE operator()(size_t i) const { return _e1(i, _j); }

	template<typename = std::enable_if_t<has_ref_access_operator<E1, TYPE>::value, bool>>
	TYPE& operator()(size_t i) { return _e1(i, _j); }

	size_t size() const { return _e1.shape().i; };

	template<typename E2>
	MatrixColumn<TYPE, E1>& operator+=(VectorExpr<TYPE, E2> const& e2) {
		if (this->size() != e2.size())
			throw IncompatibleVectorSizeException(this->size(), e2.size());

		for (size_t i = 0; i < this->size(); ++i)
		{
			_e1(i, _j) += e2(i, _j);
		}
		return *this;
	}

	template<typename E2>
	MatrixColumn<TYPE, E1>& operator-=(VectorExpr<TYPE, E2> const& e2) {
		if (this->size() != e2.size())
			throw IncompatibleVectorSizeException(this->size(), e2.size());

		for (size_t i = 0; i < this->size(); ++i)
		{
			_e1(i, _j) -= e2(i, _j);
		}
		return *this;
	}

	MatrixColumn<TYPE, E1>& operator*=(TYPE t) {
		for (size_t i = 0; i < this->size(); ++i)
		{
			_e1(i, _j) *= t;
		}
		return *this;
	}

	MatrixColumn<TYPE, E1>& operator/=(TYPE t) {
		for (size_t i = 0; i < this->size(); ++i)
		{
			_e1(i, _j) /= t;
		}
		return *this;
	}

protected:
	E1 & _e1;
	size_t _j;
};

template<typename TYPE, typename E1>
MatrixColumn<TYPE, E1> column(MatrixExpr<TYPE, E1>& mat, size_t j)
{
	return MatrixColumn<TYPE, E1>(static_cast<E1&>(mat), j);
}


template<typename TYPE, typename E1>
class ConstMatrixColumn : public VectorExpr<TYPE, ConstMatrixColumn<TYPE, E1>>
{
public:
	ConstMatrixColumn(E1 const& e1, size_t j) : _e1(e1), _j(j) {
	}

	TYPE operator()(size_t i) const { return _e1(i, _j); }

	size_t size() const { return _e1.shape().i; };

protected:
	E1 const& _e1;
	size_t _j;
};

template<typename TYPE, typename E1>
ConstMatrixColumn<TYPE, E1> column(MatrixExpr<TYPE, E1> const& mat, size_t j)
{
	return ConstMatrixColumn<TYPE, E1>(static_cast<E1 const&>(mat), j);
}

template<typename TYPE, typename E1>
class MatrixRow : public VectorExpr<TYPE, MatrixRow<TYPE, E1>>
{
public:
	MatrixRow(E1& e1, size_t i) : _e1(e1), _i(i) {
	}

	TYPE operator()(size_t j) const { return _e1(_i, j); }

	template<typename = std::enable_if_t<has_ref_access_operator<E1, TYPE>::value, bool>>
	TYPE& operator()(size_t j) { return _e1(_i, j); }

	size_t size() const { return _e1.shape().j; };

	template<typename E2>
	MatrixRow<TYPE, E1>& operator+=(VectorExpr<TYPE, E2> const& e2) {
		if (this->size() != e2.size())
			throw IncompatibleVectorSizeException(this->size(), e2.size());

		for (size_t j = 0; j < this->size(); ++j)
		{
			_e1(_i, j) += e2(j);
		}
		return *this;
	}

	template<typename E2>
	MatrixRow<TYPE, E1>& operator-=(VectorExpr<TYPE, E2> const& e2) {
		if (this->size() != e2.size())
			throw IncompatibleVectorSizeException(this->size(), e2.size());

		for (size_t j = 0; j < this->size(); ++j)
		{
			_e1(_i, j) -= e2(j);
		}
		return *this;
	}

	MatrixRow<TYPE, E1>& operator*=(TYPE t) {
		for (size_t j = 0; j < this->size(); ++j)
		{
			_e1(_i, j) *= t;
		}
		return *this;
	}

	MatrixRow<TYPE, E1>& operator/=(TYPE t) {
		for (size_t j = 0; j < this->size(); ++j)
		{
			_e1(_i, j) /= t;
		}
		return *this;
	}

protected:
	E1& _e1;
	size_t _i;
};

template<typename TYPE, typename E1>
MatrixRow<TYPE, E1> row(MatrixExpr<TYPE, E1>& mat, size_t i)
{
	return MatrixRow<TYPE, E1>(static_cast<E1&>(mat), i);
}

template<typename TYPE, typename E1>
class ConstMatrixRow : public VectorExpr<TYPE, ConstMatrixRow<TYPE, E1>>
{
public:
	ConstMatrixRow(E1 const& e1, size_t i) : _e1(e1), _i(i) {
	}

	TYPE operator()(size_t j) const { return _e1(_i, j); }
protected:
	E1 const& _e1;
	size_t _i;
};

template<typename TYPE, typename E1>
ConstMatrixRow<TYPE, E1> row(MatrixExpr<TYPE, E1> const& mat, size_t j)
{
	return ConstMatrixRow<TYPE, E1>(static_cast<E1 const&>(mat), j);
}

template<typename TYPE, typename E1>
class UpperTriangular : public MatrixExpr<TYPE, UpperTriangular<TYPE, E1>>
{
public:
	UpperTriangular(E1& e1) : _e1(e1) {}

	TYPE operator()(size_t i, size_t j) const {
		if (i > j)
			return TYPE(0);
		else
			return _e1(i, j);
	}

	TYPE& operator()(size_t i, size_t j) {
		if (i > j)
			throw BadAccesException(i, j);
		else
			return _e1(i, j);
	}

	Shape shape() const { return _e1.shape(); }
protected:
	E1& _e1;
};

template<typename TYPE, typename E1>
UpperTriangular<TYPE, E1> upper(MatrixExpr<TYPE, E1>& e1)
{
	return UpperTriangular<TYPE, E1>(static_cast<E1&>(e1));
}

template<typename TYPE, typename E1>
class ConstUpperTriangular : public MatrixExpr<TYPE, ConstUpperTriangular<TYPE, E1>>
{
public:
	ConstUpperTriangular(E1 const& e1) : _e1(e1) {}

	TYPE operator()(size_t i, size_t j) const {
		if (i > j)
			return TYPE(0);
		else
			return _e1(i, j);
	}

	Shape shape() const { return _e1.shape(); }
protected:
	E1 const& _e1;
};

template<typename TYPE, typename E1>
ConstUpperTriangular<TYPE, E1> upper(MatrixExpr<TYPE, E1> const& e1)
{
	return ConstUpperTriangular<TYPE, E1>(static_cast<E1 const&>(e1));
}

template<typename TYPE, typename E1, typename TYPE_EXP>
class PowMatrix : public MatrixExpr<TYPE, PowMatrix<TYPE, E1, TYPE_EXP>>
{
public:
	PowMatrix(E1 const& e1, TYPE_EXP exp) : _e1(e1), _exp(exp) {}

	TYPE operator()(size_t i, size_t j) const {
		return std::pow(_e1(i, j), _exp);
	}

	Shape shape() const { return _e1.shape(); }
protected:
	E1 const& _e1;
	TYPE_EXP _exp;
};

template<typename TYPE, typename E1, typename TYPE_EXP>
PowMatrix<TYPE, E1, TYPE_EXP> pow(MatrixExpr<TYPE, E1> const& e1, TYPE_EXP exp)
{
	return PowMatrix<TYPE, E1, TYPE_EXP>(static_cast<E1 const&>(e1), exp);
}

template<typename TYPE, typename E1, typename RETURN_TYPE>
class ByElementFunctionMatrix : public MatrixExpr<TYPE, ByElementFunctionMatrix<TYPE, E1, RETURN_TYPE>>
{
public:
	ByElementFunctionMatrix(E1 const& e1, std::function<RETURN_TYPE(TYPE)> byElementFunction) : _e1(e1), _byElementFunction(byElementFunction) {}

	RETURN_TYPE operator()(size_t i, size_t j) const {
		return _byElementFunction(_e1(i, j));
	}

	Shape shape() const { return _e1.shape(); }
protected:
	E1 const& _e1;
	std::function<RETURN_TYPE(TYPE)> _byElementFunction;
};

template<typename TYPE, typename E1, typename RETURN_TYPE>
ByElementFunctionMatrix<TYPE, E1, RETURN_TYPE> vectorize(MatrixExpr<TYPE, E1> const& e1, std::function<RETURN_TYPE(TYPE)> byElementFunction)
{
	return ByElementFunctionMatrix<TYPE, E1, RETURN_TYPE>(static_cast<E1 const&>(e1), byElementFunction);
}

template<typename TYPE, typename E1>
TYPE norm(MatrixExpr<TYPE, E1> const& e1)//Frobenius norm
{
	auto sq = pow(e1, 2);
	TYPE s{};

	for (size_t i = 0; i < e1.shape().i; ++i)//Amélioration si implémentation d'itérateurs, notamment pour les matrices triangulaires supérieures où il y a une moitié de 0
	{
		for (size_t j = 0; j < e1.shape().j; ++j)
		{
			s += sq(i, j);
		}
	}
	return std::sqrt(s);
}

template<typename TYPE, typename E1>
TYPE norm(MatrixExpr<std::complex<TYPE>, E1> const& e1)//Frobenius norm
{
	auto f = [](std::complex<TYPE> elem) { return std::real(elem) * std::real(elem) + std::imag(elem) * std::imag(elem); };
	auto sq = vectorize<std::complex<TYPE>, E1, TYPE>(e1, f);	
	TYPE s{};

	for (size_t i = 0; i < e1.shape().i; ++i)//Amélioration si implémentation d'itérateurs, notamment pour les matrices triangulaires supérieures où il y a une moitié de 0
	{
		for (size_t j = 0; j < e1.shape().j; ++j)
		{
			s += sq(i, j);
		}
	}
	return std::sqrt(s);
}

typedef Matrix<float> MatrixF;
typedef Matrix<double> MatrixD;
typedef Matrix<std::complex<float>> MatrixCF;
typedef Matrix<std::complex<double>> MatrixCD;