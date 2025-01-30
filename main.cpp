/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include"opt_alg.h"
#include <fstream>

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();

int main()
{
	try
	{
		lab1();
	}
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	system("pause");
	return 0;
}

void lab0()
{
	//Funkcja testowa
	double epsilon = 1e-2;
	int Nmax = 0;
	matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
	solution opt;
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Wahadlo
	Nmax = 10000;
	epsilon = 1e-2;
	lb = 0;
	ub = 4;
	double teta_opt = 1.5;
	opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Zapis symulacji do pliku csv
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(opt.x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
	ofstream Sout("symulacja_lab0.csv");
	Sout << hcat(Y[0], Y[1]);
	Sout.close();
	Y[0].~matrix();
	Y[1].~matrix();
}

void lab1()
{
	// Ekspansja testowa
	int Nmax {1000};
	double step{1}, x0;
	double epsilon {1e-5}, gamma{1e-20};
	srand(time(NULL));
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> x0r(-100, 100);

	// Test 1zm

	//x0 = x0r(gen);
	//auto xp = expansion(ff1T, x0, step, alpha, Nmax);
	//printf("\nX_start: %f Boundaries: %f %f iters %i\n", x0 ,xp[0], xp[1], solution::f_calls);
	//solution::clear_calls();

	//auto fb = fib(ff1T, -100, 100, epsilon);
	//printf("\n X: %f Y: %f Fcalls: %i", fb.x(), fb.y(), solution::f_calls);
	//solution::clear_calls();

	//auto lg = lag(ff1T, -100, 100, epsilon, gamma, Nmax);
	//printf("\n X: %f Y: %f Fcalls: %i", lg.x(), lg.y(), solution::f_calls);

	// Loop

	ofstream Sout("symulacja_lab1.csv");
	double a[] = {1.25, 2.8, 6.5};
	for (double alpha : a) {
		for(int i = 0; i < 100; i++) {
			x0 = x0r(gen);
			auto p = expansion(ff1T, x0, step, alpha, Nmax);
			//printf("\n%f %f %f %i ", x0, p[0], p[1], solution::f_calls);
			Sout << x0 << "," << p[0] << "," << p[1] << "," << solution::f_calls << ",";
			solution::clear_calls();
			solution fb = fib(ff1T, p[0], p[1], epsilon);
			Sout << fb.x() << "," << fb.y() << "," << solution::f_calls << ",";
			solution::clear_calls();
			solution lg = lag(ff1T, p[0], p[1], epsilon, gamma, Nmax);
			Sout << lg.x() << "," << lg.y() << "," << solution::f_calls << endl;
			solution::clear_calls();
		}
	}
	Sout.close();
	// FF1R
	auto fibR = fib(ff1R, 1e-5, 1e-2, epsilon);
	std::cout << "FIB" << fibR << std::endl;
	auto lagR = lag(ff1R, 1e-5, 1e-2, epsilon, gamma, Nmax);
	std::cout << "LAG" << lagR << std::endl;
}

void lab2()
{

}

void lab3()
{

}

void lab4()
{

}

