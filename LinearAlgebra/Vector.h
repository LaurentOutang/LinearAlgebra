#pragma once
#include <vector>
#include <iostream>
#include <complex>
#include <type_traits>
#include "Traits.h"
#include "VectorExceptions.h"

template<typename TYPE, typename E>
class VectorExpr
{
public:
	TYPE operator()(size_t i) const
	{
		return static_cast<E const&>(*this)(i);
	}

	size_t size() const
	{
		return static_cast<E const&>(*this).size();
	}
};

template<typename TYPE>
class Vector : public VectorExpr<TYPE, Vector<TYPE>>
{
public:

	Vector(std::initializer_list<TYPE> list) : _vec(list) {}

	Vector(size_t i) : _vec(i) {}

	template<typename E>
	Vector(VectorExpr<TYPE, E> const& expr) : _vec(expr.size()) {
		for (size_t i = 0; i != expr.size(); ++i) {
			_vec[i] = expr(i);
		}
	}

	Vector(Vector<TYPE> const& v) : _vec(v.size()) {
		for (size_t i = 0; i != v.size(); ++i) {
			_vec[i] = v(i);
		}
	}

	TYPE operator()(size_t i) const { return _vec[i]; }

	TYPE& operator()(size_t i) { return _vec[i]; }

	size_t size() const { return _vec.size(); }

	template<typename E1>
	Vector& operator+=(VectorExpr<TYPE, E1> const& expr)
	{
		for (size_t k = 0; k < this->size(); ++k)
		{
			_vec[k] += expr(k);
		}
		return *this;
	}

	template<typename E1>
	Vector& operator-=(VectorExpr<TYPE, E1> const& expr)
	{
		for (size_t k = 0; k < this->size(); ++k)
		{
			_vec[k] -= expr(k);
		}
		return *this;
	}

	Vector& operator*=(TYPE t)
	{
		for (size_t k = 0; k < this->size(); ++k)
		{
			_vec[k] *= t;
		}
		return *this;
	}

	Vector& operator/=(TYPE t)
	{
		for (size_t k = 0; k < this->size(); ++k)
		{
			_vec[k] /= t;
		}
		return *this;
	}
protected:
	std::vector<TYPE> _vec;
};

template<typename TYPE, typename E1, typename E2>
TYPE dot(VectorExpr<TYPE, E1> const& u, VectorExpr<TYPE, E2> const& v)
{
	TYPE s{};
	for (size_t i = 0; i < u.size(); ++i)
	{
		s += u(i) * v(i);
	}
	return s;
}

template<typename TYPE, typename E1, typename E2>
std::complex<TYPE> dot(VectorExpr<std::complex<TYPE>, E1> const& u, VectorExpr<std::complex<TYPE>, E2> const& v)
{
	std::complex<TYPE> s{};
	for (size_t i = 0; i < u.size(); ++i)
	{
		s += u(i) * std::conj(v(i));
	}
	return s;
}

template<typename TYPE>
std::ostream& operator<<(std::ostream& os, Vector<TYPE> const& v)
{
	size_t vsize = v.size();
	if (vsize > 0)
	{
		for (size_t i = 0; i < vsize-1; ++i)
		{
			os << v(i) << ',';
		}
		os << v(vsize - 1);
	}
	return os;
}

template<typename TYPE, typename E1, typename E2>
class VectorSum : public VectorExpr<TYPE, VectorSum<TYPE, E1, E2>>
{
protected:
	static TYPE test;
	E1 const& _e1;
	E2 const& _e2;
public:
	VectorSum(E1 const& e1, E2 const& e2) : _e1(e1), _e2(e2) {
		if (_e1.size() != _e1.size())
			throw IncompatibleVectorSizeException(_e1.size(), _e2.size())
	}

	TYPE operator()(size_t i) const { return _e1(i) + _e2(i); }
	size_t size() const { return _e1.size(); }
};

template<typename TYPE, typename E1, typename E2>
TYPE VectorSum<TYPE, E1, E2>::test{};

template<typename TYPE, typename E1, typename E2>
VectorSum<TYPE, E1, E2> operator+(VectorExpr<TYPE, E1> const& e1, VectorExpr<TYPE, E2> const& e2)
{
	return VectorSum<TYPE, E1, E2>(*static_cast<const E1*>(&e1), *static_cast<const E2*>(&e2));
}

template<typename TYPE, typename E1, typename E2>
class VectorDiff : public VectorExpr<TYPE, VectorDiff<TYPE, E1, E2>>
{
protected:
	E1 const& _e1;
	E2 const& _e2;
public:
	VectorDiff(E1 const& e1, E2 const& e2) : _e1(e1), _e2(e2) {
		if (_e1.size() != _e2.size())
			throw IncompatibleVectorSizeException(_e1.size(), _e2.size());
	}

	TYPE operator()(size_t i) const { return _e1(i) - _e2(i); }
	size_t size() const { return _e1.size(); }
};

template<typename TYPE, typename E1, typename E2>
VectorDiff<TYPE, E1, E2> operator-(VectorExpr<TYPE, E1> const& e1, VectorExpr<TYPE, E2> const& e2)
{
	return VectorDiff<TYPE, E1, E2>(*static_cast<const E1*>(&e1), *static_cast<const E2*>(&e2));
}

template<typename TYPE, typename E1>
class VectorMult : public VectorExpr<TYPE, VectorMult<TYPE, E1>>
{
protected:
	E1 const& _e1;
	TYPE const& _t;
public:
	VectorMult(E1 const& e1, TYPE const& t) : _e1(e1), _t(t) {
	}

	TYPE operator()(size_t i) const { return _t*_e1(i); }
	size_t size() const { return _e1.size(); }
};

template<typename TYPE, typename E1>
VectorMult<TYPE, E1> operator*(VectorExpr<TYPE, E1> const& e1, TYPE const& t)
{
	return VectorMult<TYPE, E1>(*static_cast<const E1*>(&e1), t);
}

template<typename TYPE, typename E1>
VectorMult<TYPE, E1> operator*(TYPE const& t, VectorExpr<TYPE, E1> const& e1)
{
	return VectorMult<TYPE, E1>(*static_cast<const E1*>(&e1), t);
}

template<typename TYPE, typename E1, typename = std::enable_if_t < is_complex<TYPE>{}>>
class VectorConj : public VectorExpr<TYPE, VectorConj<TYPE, E1>>
{
protected:
	E1 const& _e1;
public:
	VectorConj(E1 const& e1) : _e1(e1) {}

	TYPE operator()(size_t i) const { return std::conj(_e1(i)); }
	size_t size() const { return _e1.size(); }
};

template<typename TYPE, typename E1, typename = std::enable_if_t<is_complex<TYPE>{}>>
VectorConj<TYPE, E1> conj(VectorExpr<TYPE, E1> const& e1)
{
	return VectorConj<TYPE, E1>(*static_cast<const E1*>(&e1));
}

template<typename TYPE, typename E1>
class SubVector : public VectorExpr<TYPE, SubVector<TYPE, E1>>
{
public:
	SubVector(E1& e1, size_t start, size_t end) : _e1(e1), _start(start), _end(end) {
		if (_end < _start)
			throw SubVectorSizeException(_start, _end);
	}

	TYPE operator()(size_t i) const { return _e1(i + _start); }

	template<typename = std::enable_if_t<has_ref_access_operator<E1, TYPE>::value, bool>>
	TYPE& operator()(size_t i) { return _e1(i + _start); }

	size_t size() const { return _end - _start + 1; };

	template<typename E2>
	SubVector& operator+=(VectorExpr<TYPE, E2> const& vec)
	{
		if (this->size() != vec.size())
			throw IncompatibleVectorSizeException(this->size(), vec.size());

		for (size_t k = 0; k < this->size(); ++k)
		{
			this->operator()(k) += vec(k);
		}
		return *this;
	}

	template<typename E2>
	SubVector& operator-=(VectorExpr<TYPE, E2> const& vec)
	{
		if (this->size() != vec.size())
			throw IncompatibleVectorSizeException(this->size(), vec.size());

		for (size_t k = 0; k < this->size(); ++k)
		{
			this->operator()(k) -= vec(k);
		}
		return *this;
	}

	SubVector& operator*=(TYPE t)
	{
		for (size_t k = 0; k < this->size(); ++k)
		{
			this->operator()(k) *= t;
		}
		return *this;
	}

	SubVector& operator/=(TYPE t)
	{
		for (size_t k = 0; k < this->size(); ++k)
		{
			this->operator()(k) /= t;
		}
		return *this;
	}

protected:
	E1& _e1;
	size_t _start;
	size_t _end;
};

template<typename TYPE, typename E1>
SubVector<TYPE, E1> subVector(VectorExpr<TYPE, E1>& e1, size_t start, size_t end)
{
	return SubVector<TYPE, E1>(static_cast<E1&>(e1), start, end);
}

template<typename TYPE, typename E1>
class ConstSubVector : public VectorExpr<TYPE, ConstSubVector<TYPE, E1>>
{
public:
	ConstSubVector(E1 const& e1, size_t start, size_t end) : _e1(e1), _start(start), _end(end) {
		if (_end < _start)
			throw SubVectorSizeException(_start, _end);
	}

	TYPE operator()(size_t i) const { return _e1(i + _start); }

	size_t size() const { return _end - _start + 1; };
protected:
	E1 const& _e1;
	size_t _start;
	size_t _end;
};

template<typename TYPE, typename E1>
ConstSubVector<TYPE, E1> subVector(VectorExpr<TYPE, E1> const& e1, size_t start, size_t end)
{
	return ConstSubVector<TYPE, const E1>(static_cast<E1 const&>(e1), start, end);
}

typedef Vector<double> VectorD;
typedef Vector<float> VectorF;
typedef Vector<std::complex<double>> VectorCD;
typedef Vector<std::complex<float>>  VectorCF;