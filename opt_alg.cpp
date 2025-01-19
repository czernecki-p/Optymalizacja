#include"opt_alg.h"

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


double* expansion(matrix(ff)(matrix, matrix, matrix), double x0, double d, double alpha, int Nmax, matrix ud1, matrix ud2){
    try{
        double *p = new double[2] {0,0};
        solution X0{x0}, X1 {x0 + d}, X2{};

        X0.fit_fun(ff,ud1, ud2);
        X1.fit_fun(ff,ud1, ud2);

        if (X0.y() == X1.y()) {
            p[0] = x0, p[1] = x0+d;
            return p;
        }
        if (X0.y() < X1.y()) {
            d = -d;
            X1.x=x0+d;
            X1.fit_fun(ff,ud1, ud2);
            if (X0.y() <= X1.y) {
                p[0] = X1.x(), p[1] = x0-d;
                return p;
            }
        }
        int i = 0;
        bool loop = true;
        while(loop){
            if(solution::f_calls > Nmax) {
                X1.flag = 1;
                throw std::range_error(" iter > Nmax ");
            } else {
                ++i;
                X2.x = x0 + pow(alpha,i) * d;
                X2.fit_fun(ff, ud1, ud2);
                if (X1.y() <= X2.y()) {
                    loop = false;
                } else {
                    X0 = X1, X1 = X2;
                }
            }
        }
        if ( d > 0) {
            p[0] = X0.x(), p[1] = X2.x();
        } else {
            p[0] = X2.x(), p[1] = X0.x();
        }
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
                fc = ff(c_mat, ud1, ud2);
            }
            else {
                a_i = c_i;
                c_i = d_i;
                fc = fd;
                d_i = a_i + ((double)fib[i -1] / fib[i]) * (b_i - a_i);
                d_mat = matrix(1,1, d_i);
                fd = ff(d_mat, ud1, ud2);
            }
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
		solution Xopt;
		//Tu wpisz kod funkcji

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

