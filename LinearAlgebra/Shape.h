#pragma once

struct Shape
{
	Shape(size_t i, size_t j) : i(i), j(j) {}
	bool operator==(Shape const& s) const { return this->i == s.i && this->j == s.j; }
	bool operator!=(Shape const& s) const { return this->i != s.i || this->j != s.j; }
	bool notNull() const { return i != 0 && j != 0; }
	size_t size() const { return i * j; }

	size_t i;
	size_t j;
};

std::ostream& operator<<(std::ostream& os, Shape const& shape)
{
	os << '(' << shape.i << ',' << shape.j << ')';
	return os;
}

