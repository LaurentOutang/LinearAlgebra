#pragma once
#include <exception>
#include <sstream>

class IncompatibleVectorSizeException : public std::exception
{
public:
	IncompatibleVectorSizeException(size_t l1, size_t l2)
	{
		std::stringstream ss;
		ss << "Vector of size " << l1 << " is incompatible with vector of size " << l2 << "for this operation.";
		msg = ss.str();
	}
	const char* what() const noexcept override
	{
		return msg.c_str();
	}
protected:
	std::string msg;
};

class SubVectorSizeException : public std::exception
{
public:
	SubVectorSizeException(size_t start, size_t end)
	{
		std::stringstream ss;
		ss << "Trying to extract vector of size " << end - start + 1 << " from index " << start << " to " << end << '.';
		msg = ss.str();
	}
	const char* what() const noexcept override
	{
		return msg.c_str();
	}
protected:
	std::string msg;
};