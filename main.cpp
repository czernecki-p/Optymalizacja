/*********************************************
Kod stanowi uzupe�nienie materia��w do �wicze�
w ramach przedmiotu metody optymalizacji.
Kod udost�pniony na licencji CC BY-SA 3.0
Autor: dr in�. �ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia G�rniczo-Hutnicza
Data ostatniej modyfikacji: 19.09.2023
*********************************************/

#include"opt_alg.h"

void lab0();
void lab1();
void lab2();
void lab3();
void lab4();

int main()
{
	try
	{
		lab0();
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
	int Nmax = 10000;
	matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
	solution opt;
	a(0) = -1;
	a(1) = 2;
	opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
	cout << opt << endl << endl;
	solution::clear_calls();

	//Wahadlo
	Nmax = 1000;
	epsilon = 1e-2;
	lb = 0;
	ub = 5;
	double teta_opt = 1;
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
	// Parametry
    double a = -100, b = 100;
    double epsilon = 1e-3, gamma = 1e-3;
    int Nmax = 1000;
    double alpha = 2.0, d = 1.0;

    // Wartości optymalizacji
    matrix ud1, ud2;

    // Funkcja testowa
    double* bounds = expansion(ff1T, 0.0, d, alpha, Nmax, ud1, ud2);
    solution fib_opt = fib(ff1T, bounds[0], bounds[1], epsilon, ud1, ud2);
    solution lag_opt = lag(ff1T, bounds[0], bounds[1], epsilon, gamma, Nmax, ud1, ud2);

    cout << "Wyniki optymalizacji funkcji testowej:" << endl;
    cout << "Metoda Fibonacciego: " << fib_opt << endl;
    cout << "Metoda Lagrange'a: " << lag_opt << endl;

    delete[] bounds;
}

void lab2()
{
// Parametry
double a = 1, b = 100;
double epsilon = 1e-3, gamma = 1e-3;
int Nmax = 2000;

// Problem rzeczywisty
matrix ud1, ud2;
double* bounds = expansion(ff0R, 50.0, 1.0, 2.0, Nmax, ud1, ud2);

solution fib_opt = fib(ff0R, bounds[0], bounds[1], epsilon, ud1, ud2);
solution lag_opt = lag(ff0R, bounds[0], bounds[1], epsilon, gamma, Nmax, ud1, ud2);

cout << "Wyniki optymalizacji problemu rzeczywistego:" << endl;
cout << "Metoda Fibonacciego: " << fib_opt << endl;
cout << "Metoda Lagrange'a: " << lag_opt << endl;

delete[] bounds;
}

void lab3()
{

}

void lab4()
{

}

