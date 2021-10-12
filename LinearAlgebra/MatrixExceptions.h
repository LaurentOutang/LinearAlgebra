#pragma once
#include <exception>
#include "Shape.h"

enum class OPERATION { MULT, ADD, SUB };

class IncompatibleShapesException : public std::exception
{
public:
	IncompatibleShapesException(Shape shape1, Shape shape2, OPERATION operation)
	{
		std::stringstream ss;
		ss << "Shape " << shape1 << " is incompatible with " << shape2 << "for ";
		if (operation == OPERATION::MULT)
		{
			ss << "multiplication.";
		}
		else if (operation == OPERATION::ADD)
		{
			ss << "addition.";
		}
		else if (operation == OPERATION::SUB)
		{
			ss << "substraction.";
		}
		msg = ss.str();
	}

	const char* what() const noexcept override 
	{ 
		return msg.c_str();
	}
protected:
	std::string msg;
};

class BadAccesException : public std::exception
{
public:
	BadAccesException(size_t i, size_t j)
	{
		std::stringstream ss;
		ss << "You are trying to access element " << Shape(i, j) << '.';
		msg = ss.str();
	}

	const char* what() const noexcept override
	{
		return msg.c_str();
	}
protected:
	std::string msg;
};

