#include"user_funs.h"


matrix ff0R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{ m2d(x),0.5 });
	matrix* Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);
	int n = get_len(Y[0]);
	double teta_max = Y[1](0, 0);
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));
	Y[0].~matrix();
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);
	double m = 1, l = 0.5, b = 0.5, g = 9.81;
	double I = m*pow(l, 2);
	dY(0) = Y(1);
	dY(1) = ((t <= ud2(1))*ud2(0) - m*g*l*sin(Y(0)) - b*Y(1)) / I;
	return dY;
}
matrix ff1_real(double t, matrix Y, matrix ud1, matrix ud2) {
    double a = 0.98;
    double b = 0.63;
    double g = 9.81;
    double P_A = 0.5, P_B = 1.0;
    double V0_A = 5.0, V0_B = 1.0;
    double T0_A = 90.0, T0_B = 20.0;
    double Tin_B = 20.0, Fin_B = 10.0 / 1000.0;
    double D_B = 36.5665 / 10000.0;
    
    double V_A = Y(0), V_B = Y(1), T_B = Y(2);
    double D_A = m2d(ud2);
    D_A /= 10000.0;

    double Fout_A = V_A > 0 ? a * b * D_A * sqrt(2 * g * V_A / P_A) : 0;
    double Fout_B = V_B > 0 ? a * b * D_B * sqrt(2 * g * V_B / P_B) : 0;

    double dV_A_dt = -Fout_A;
    double dV_B_dt = Fout_A + Fin_B - Fout_B;

    double dT_B_dt = (Fin_B / V_B) * (Tin_B - T_B) + (Fout_A / V_B) * (Y(3) - T_B);

    matrix dY(3, 1);
    dY(0) = dV_A_dt; dY(1) = dV_B_dt; dY(2) = dT_B_dt;

    return dY;
}