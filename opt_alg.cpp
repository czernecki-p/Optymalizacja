#include"opt_alg.h"

using std::vector;

// lb, ub - dolna i górna macierz do ograniczania
// epison - dok?adno?? oblicze?

solution MC(matrix(*ff)(matrix, matrix, matrix), int N, matrix lb, matrix ub, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		while (true)
		{
			Xopt = rand_mat(N);
			for (int i = 0; i < N; ++i)
				Xopt.x(i) = (ub(i) - lb(i)) * Xopt.x(i) + lb(i);
			Xopt.fit_fun(ff, ud1, ud2);
			if (Xopt.y < epsilon)
			{
				Xopt.flag = 1;
				break;
			}
			if (solution::f_calls > Nmax)
			{
				Xopt.flag = 0;
				break;
			}
		}
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution MC(...):\n" + ex_info);
	}
}

double* expansion(matrix(*ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2){
	try{
		double* p = new double[2]{0,0};
		int i =0;
		solution X0(x0), X1(x0+d);
		X0.fit_fun(ff, ud1, ud2);
		X1.fit_fun(ff, ud1, ud2);

		bool cond{false};
		if (X1.y == X0.y) cond = true;
		else if (X1.y > X0.y) {
			d *= -1;
			X1.x = x0 +d;
			X1.fit_fun(ff, ud1, ud2);
			if(X1.y >= X0.y) cond = true;
		}

		if(cond) {
			p[0] = m2d((d>0) ? X0.x : X1.x);
			p[1] = m2d((d>0) ? X1.x : X0.x -d);
		}

		solution X2;
		while (solution::f_calls <= Nmax) {
			++i;
			X2.x = x0 + pow(alpha,i) * d;
			X2.fit_fun(ff, ud1, ud2);
			if(X2.y >= X1.y) break;
			X0 = X1;
			X1 = X2;
		}

		p[0] = m2d((d>0) ? X0.x : X2.x);
		p[1] = m2d((d>0) ? X2.x : X0.x);
		return p;
	}
	catch (string ex_info)
	{
		throw ("double* expansion(...):\n" + ex_info);
	}
}

solution fib(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		double iterations = (b - a) / epsilon;
		std::vector <int> fib = {1,1};
		while (fib.back() <= iterations) {
			fib.push_back(fib[fib.size() - 1] + fib[fib.size() - 2]);
		}
		int minFib = fib.size();

		double a_i = a, b_i = b;
		double c_i = b_i - ((double)fib[minFib-2] / fib[minFib-1]) * (b_i - a_i);
		double d_i = a_i + ((double)fib[minFib - 2] / fib[minFib - 1]) * (b_i - a_i);

		matrix c_mat(1,1,c_i), d_mat(1,1,d_i);
		matrix fc = ff(c_mat, ud1, ud2);
		matrix fd = ff(d_mat, ud1, ud2);

		for (int i = minFib - 2; i > 1;--i) {
			if (fc<fd) {
				b_i = d_i;
				d_i = c_i;
				fd = fc;
				c_i = b_i - ((double)fib[i - 1] / fib[i]) * (b_i - a_i);
				c_mat = matrix(1,1,c_i);
				fc = ff(c_mat, ud1, ud2); Xopt.f_calls++;
			}
			else {
				a_i = c_i;
				c_i = d_i;
				fc = fd;
				d_i = a_i + ((double)fib[i -1] / fib[i]) * (b_i - a_i);
				d_mat = matrix(1,1, d_i);
				fd = ff(d_mat, ud1, ud2); Xopt.f_calls++;
			}
			//printf("%f\n", b_i - a_i);
		}
		double x_opt = (a_i + b_i) / 2;
		Xopt.x = matrix(1,1,x_opt);
		Xopt.y = ff(Xopt.x, ud1, ud2);
		Xopt.flag = 1;

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution fib(...):\n" + ex_info);
	}

}

solution lag(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, double gamma, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt{}, A{a}, B{b}, C{(a+b)/2}, D{a}, D_prev{};
		double l{}, m{};
		A.fit_fun(ff, ud1, ud2);
		B.fit_fun(ff, ud1, ud2);
		C.fit_fun(ff, ud1, ud2);
		do {
			if(solution::f_calls > Nmax) {
				Xopt = D;
				Xopt.flag = 2;
				return Xopt;
			}
			l = m2d(A.y * (pow(B.x) - pow(C.x)) + B.y() * (pow(C.x()) - pow(A.x())) + C.y() * (pow(A.x()) - pow(B.x())));
			m = m2d(A.y * (B.x - C.x) + B.y*(C.x - A.x) + C.y*(A.x-B.x));
			if(m<=0) {
				Xopt.flag = 2;
				return D;
			}
			D_prev = D;
			D.x=0.5*l/m;
			D.fit_fun(ff, ud1, ud2);
			//printf("%f\n", B.x()-A.x());
			if(A.x <= D.x && D.x <= C.x) {
				if(D.y < C.y) {
					B=C;
					C=D;
				} else {
					A=D;
				}
			}else if (C.x <= D.x && D.x <= B.x) {
				if(D.y < C.y) {
					A=C;
					C=D;
				} else {
					B=D;
				}
			} else {
				//printf("\nAlgorytm nie jest zbiezny.");
				Xopt = D;
				Xopt.flag = 2;
				return Xopt;
			}
		}while(B.x()-A.x() >= epsilon || abs(D.x() - D_prev.x()) >= gamma);
		Xopt = D_prev;
		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution lag(...):\n" + ex_info);
	}
}

solution HJ(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution HJ(...):\n" + ex_info);
	}
}

solution HJ_trial(matrix(*ff)(matrix, matrix, matrix), solution XB, double s, matrix ud1, matrix ud2)
{
	try
	{
		//Tu wpisz kod funkcji

		return XB;
	}
	catch (string ex_info)
	{
		throw ("solution HJ_trial(...):\n" + ex_info);
	}
}

solution Rosen(matrix(*ff)(matrix, matrix, matrix), matrix x0, matrix s0, double alpha, double beta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Rosen(...):\n" + ex_info);
	}
}

solution pen(matrix(*ff)(matrix, matrix, matrix), matrix x0, double c, double dc, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try {
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution pen(...):\n" + ex_info);
	}
}

solution sym_NM(matrix(*ff)(matrix, matrix, matrix), matrix x0, double s, double alpha, double beta, double gamma, double delta, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution sym_NM(...):\n" + ex_info);
	}
}

solution SD(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution SD(...):\n" + ex_info);
	}
}

solution CG(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution CG(...):\n" + ex_info);
	}
}

solution Newton(matrix(*ff)(matrix, matrix, matrix), matrix(*gf)(matrix, matrix, matrix),
	matrix(*Hf)(matrix, matrix, matrix), matrix x0, double h0, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution Newton(...):\n" + ex_info);
	}
}

solution golden(matrix(*ff)(matrix, matrix, matrix), double a, double b, double epsilon, int Nmax, matrix ud1, matrix ud2)
{
	try
	{
		solution Xopt;
		//Tu wpisz kod funkcji

		return Xopt;
	}
	catch (string ex_info)
	{
		throw ("solution golden(...):\n" + ex_info);
	}
}

