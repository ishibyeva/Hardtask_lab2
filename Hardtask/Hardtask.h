#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <iomanip>


struct Point
{
	int num;
	double x, y1, y2;
};
struct twopoint
{
	double x1, x2;
};
double osnov2y1(double x, double u1, double u2)//f1'
{
	double y1 = (-1)*500.005*u1 + 499.995*u2;
	return y1;
}
double osnov2y2(double x, double u1, double u2)//f2'
{
	double y2 = 499.995*u1 + (-1)*500.005*u2;
	return y2;
}
Point metodRKS_N(int num, double h, double x, double u1, double u2)
{
	double v1 = u1;
	double v2 = u2;
	double a1 = -500.005, a2 = 499.995;
	double det = 1.0 / abs(pow((1 - h * 0.5*a1), 2) - pow(h*0.5*a2, 2));
	double k1 = 0, k2 = 0;
	//
	k1 = osnov2y1(x, u1, u2);
	k2 = det*((1-pow(h*0.5*a1,2)+pow(h*0.5*a2,2))*k1+h*a2*osnov2y2(x, u1, u2));
	v1 = u1 + (h*0.5) * (k1 + k2);

	k1 = 0, k2 = 0;
	k1 = osnov2y2(x, u1, u2);
	k2 = det*(h*a2*osnov2y1(x, u1, u2)+(1+pow(h*0.5*a2,2)-pow(h*0.5*a1,2))*k1);
	v2 = u2 + (h*0.5) * (k1 + k2);

	x += h;

	Point st;
	st.num = num;
	st.x = x;
	st.y1 = v1;
	st.y2 = v2;

	return st;
}
Point Toch_R(int num, double h, double x)
{
	Point t;
	t.x = x + h;
	t.y1 = 10 * exp(-0.01*(x+h)) - 3 * exp(-1000 * (x+h));
	t.y2 = 10 * exp(-0.01*(x+h)) + 3 * exp(-1000*(x+h));
	t.num = num;
	return t;
}

void RezOsn2(int maxnum, double h00, double xn)
{
	double E = 0.005;
	Point *mas;
	mas = new Point[maxnum];
	Point *obh;
	obh = new Point[maxnum];
	
	double *e = new double[maxnum]; 
	double *e1 = new double[maxnum]; 
	double *e2 = new double[maxnum];
	 e1[0] = 0;
	 e2[0] = 0;

	double  x0 = 0.0, h0 = h00, u01 = 7.0, u02 = 13.0;
	//xn - граница отрезка интегрирования
	int n = maxnum;


	Point t;
	t.num = 0;
	t.x = x0;
	t.y1 = u01;
	t.y2 = u02;
	mas[0] = t;
	obh[0] = t;
	double hrez = 0;

	int i = 1;
	while (i < n)
	{
		Point t1, t12, t2;

		//(x(n+1),v(n+1))
		x0 = mas[i - 1].x;
		u01 = mas[i - 1].y1;
		u02 = mas[i - 1].y2;
		//
		if (x0 < 0.035)
		{
			t1 = metodRKS_N(i, h0, x0, u01, u02);
			t12 = metodRKS_N(i, h0*0.5, x0, u01, u02);
			double xx = t12.x, uu1 = t12.y1, uu2 = t12.y2;
			t2 = metodRKS_N(i, h0*0.5, xx, uu1, uu2);

			e1[i] = abs(t2.y1 - t1.y1);
			e2[i] = abs(t2.y2 - t1.y2);
			double S = abs((e1[i] > e2[i]) ? e1[i] : e2[i]);

			int p = 2;
			e[i] = S * pow(2, p);

			if (S < E / (pow(2, p + 1)))
			{
				h0 = 2.0 * h0;
				mas[i] = t1;
				obh[i] = Toch_R(i, h0, x0);

				if (mas[i].x > xn)
					break;
				i++;

			}
			if (S > E)
			{
				h0 = h0 * (0.5);

			}
			if ((S > E / (pow(2, p + 1))) && (S < E))
			{
				mas[i] = t1;
				obh[i] = Toch_R(i, h0, x0);

				if (mas[i].x > 0.035)
					break;
				i++;


			}
		}
		else
		{
			h0 = h00;
			t1 = metodRKS_N(i, h0, x0, u01, u02);
			t2= Toch_R(i, h0, x0);
			e1[i] = abs(t2.y1 - t1.y1);
			e2[i] = abs(t2.y2 - t1.y2);
			mas[i] = t1;
			obh[i] = t2;
			i++;
		}
	}

	for (int d = 0; (d < i) && (mas[d].x < xn); d++)
	{

		std::cout << std::setw(7) << std::left<< d  <<std::setw(7) << mas[d].x << std::setw(20) << std::right
			<< "(" << std::setprecision(7) << obh[d].y1 << " ; " << std::setprecision(7) << obh[d].y2 << ")"
			<< std::setw(20) << std::right<< "(" << std::setprecision(7) << mas[d].y1 << " ; " << std::setprecision(7) << mas[d].y2 << ")"
			<< std::setw(20) << std::right
			<< "("<< std::setprecision(7) << e1[d]<<"__"<< std::setprecision(7) << e2[d]<<")"<< std::endl;
		std::cout << std::endl;
	}

}