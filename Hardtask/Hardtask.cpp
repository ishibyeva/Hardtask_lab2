// Hardtask.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#include "Hardtask.h"

using namespace std;

int main()
{
	int maxnum;
	double h, xn=100;
    std::cout << "Enter maximum n"<< std::endl; 
	std::cin >> maxnum;
	std::cout << "Enter h" << std::endl;
	std::cin >> h;
	//std::cout << "Enter xn" << std::endl;
	//std::cin >> xn;
	std::cout << "**********Implicit Runge-Kutta method with h = " << h << ", border = "<<xn<<", max step count = "<<maxnum<<"**********" << std::endl;
	std::cout << "num\t" << "x\t\t" << "u point\t\t\t" << "v point\t\t\t" << "error\t\t" << std::endl;
	RezOsn2(maxnum,h,xn);
	system("pause");
	return 0;
}




