#pragma once
#include <complex>

template<class T> struct is_complex : std::false_type {};
template<class T> struct is_complex<std::complex<T>> : std::true_type {};

template<class T, class U> 
class has_ref_access_operator
{
private:
	template<typename V> 
	static std::true_type test(U& (V::*)(size_t)) {
		return std::true_type{};
	};

	template<typename V>
	static std::true_type test(U& (V::*)(size_t, size_t)) {
		return std::true_type{};
	};

	//Functions like <enable_if> U& operator()(...)
	template<typename V>
	static auto test(void*, void*) -> decltype(test(&V::template operator()<void>))
	{
		return test(&V::template operator()<void>);
	}

	template<typename V>
	static auto test2(void*, void*) -> decltype(test(&V::operator()))
	{
		return test(&V::operator());
	}


	//Sink hole for template function operator() & non-template function if return type & arguments do not match
	template<typename V> 
	static std::false_type test(...)
	{
		return std::false_type{};
	}

	//Sink hole for non-template function operator()
	template<typename V>
	static std::false_type test2(...)
	{
		return std::false_type{};
	}

public:
	static constexpr bool value = std::is_same<decltype(test<T>(0, 0)), std::true_type>::value || std::is_same<decltype(test2<T>(0, 0)), std::true_type>::value;
};