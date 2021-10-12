#include <iostream>
#include "Vector.h"
#include "Matrix.h"
#include "HouseHolderQR.h"

VectorD give()
{
	return { 1.0, 2.0, 3.0, 4.0 };
}

int main()
{
	MatrixCD Id5(5, 5);
	Id5(0, 0) = 1.0;
	Id5(1, 1) = 1.0;
	Id5(2, 2) = 1.0;
	Id5(3, 3) = 1.0;
	Id5(4, 4) = 1.0;

	MatrixCD C(5, 4);
	C(0, 0) = std::complex<double>(1.0, 1.0);
	C(0, 1) = std::complex<double>(2.0, 2.0);
	C(0, 2) = std::complex<double>(3.0, 3.0);
	C(0, 3) = std::complex<double>(4.0, 4.0);

	C(1, 0) = std::complex<double>(5.0, 5.0);
	C(1, 1) = std::complex<double>(6.0, 6.0);
	C(1, 2) = std::complex<double>(7.0, 7.0);
	C(1, 3) = std::complex<double>(8.0, 8.0);

	C(2, 0) = std::complex<double>(9.0, 9.0);
	C(2, 1) = std::complex<double>(10.0, 10.0);
	C(2, 2) = std::complex<double>(11.0, 11.0);
	C(2, 3) = std::complex<double>(12.0, 12.0);

	C(3, 0) = std::complex<double>(13.0, 13.0);
	C(3, 1) = std::complex<double>(14.0, 14.0);
	C(3, 2) = std::complex<double>(15.0, 15.0);
	C(3, 3) = std::complex<double>(16.0, 16.0);

	C(4, 0) = std::complex<double>(17.0, 17.0);
	C(4, 1) = std::complex<double>(18.0, 18.0);
	C(4, 2) = std::complex<double>(19.0, 19.0);
	C(4, 3) = std::complex<double>(20.0, 20.0);

	MatrixCD C_copy = C;

	HouseHolderQR<std::complex<double>> qr2(C);

	std::cout << C << std::endl << std::endl;


	MatrixCD Q = Id5 * qr2._Q;

	std::cout << Q << std::endl << std::endl;

	MatrixCD QQh = qr2._Q * adjoint(Q);

	MatrixCD QhQ = adjoint(Q) * qr2._Q;

	std::cout << norm(QQh - Id5) << std::endl << std::endl;

	std::cout << norm(QhQ - Id5) << std::endl << std::endl;

	MatrixCD C_recombined = qr2._Q * qr2._R;

	std::cout << norm(C_recombined - C_copy) << std::endl << std::endl;

	/*MatrixD Id(5, 5);
	Id(0, 0) = 1.0;
	Id(1, 1) = 1.0;
	Id(2, 2) = 1.0;
	Id(3, 3) = 1.0;
	Id(4, 4) = 1.0;

	MatrixD C(5, 4);
	C(0, 0) = 1.0;
	C(0, 1) = 2.0;
	C(0, 2) = 3.0;
	C(0, 3) = 4.0;

	C(1, 0) = 5.0;
	C(1, 1) = 6.0;
	C(1, 2) = 7.0;
	C(1, 3) = 8.0;

	C(2, 0) = 9.0;
	C(2, 1) = 10.0;
	C(2, 2) = 11.0;
	C(2, 3) = 12.0;

	C(3, 0) = 13.0;
	C(3, 1) = 14.0;
	C(3, 2) = 15.0;
	C(3, 3) = 16.0;

	C(4, 0) = 17.0;
	C(4, 1) = 18.0;
	C(4, 2) = 19.0;
	C(4, 3) = 20.0;

	MatrixD E = adjoint(C);
	std::cout << E << std::endl << std::endl;

	std::cout << C << std::endl << std::endl;

	HouseHolderQR<double> qr2(C);

	std::cout << C << std::endl << std::endl;


	MatrixD Q = qr2._Q * Id;
	
	std::cout << Q << std::endl << std::endl;

	MatrixD QQt = qr2._Q * transpose(Q);

	MatrixD QtQ = transpose(Q) * qr2._Q;

	std::cout << QQt << std::endl << std::endl;*/
	return 0;
}