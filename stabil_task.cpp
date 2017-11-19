// (c) 2017 Alexander Pchelintsev, pchelintsev.an@yandex.ru

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <cmath>
#include <sys/time.h>

#include "rk4.hpp"

using namespace std;

typedef boost::numeric::ublas::matrix<double> matrix;
typedef boost::numeric::ublas::identity_matrix<double> identity;
typedef boost::numeric::ublas::vector<double> vector_col;

#define sys_alpha  0.2
#define sys_beta   6
#define sys_gamma  1
#define sys_b      0.5
#define sys_g      9.8

#define q_1        3
#define q_2        3
#define q_3        3
#define r          0.05
#define d_z0       -1
#define k_z0       1
#define k_P        10

#define FORMAT     "png"

inline void create_matrix(const double *arr, int rows, int cols, matrix &m)
{
	m.resize(rows, cols);
	for(int i = 0; i < rows; i++)
	{
		int ind = i * cols;
		for(int j = 0; j < cols; j++)
			m(i, j) = arr[ind + j];
	}
}

matrix a, a_tr, b, prod_rb, prod_brb, Q;

class riccati: public rk4 <matrix, matrix>
{
public:
	matrix f(const double &t, const matrix &x, const matrix *x_prev,
	         const matrix *u_prev, const matrix *K, const matrix *H)
	{
		matrix tmp = prod(x, prod_brb);
		tmp = prod(tmp, x);
		return -prod(x, a) - prod(a_tr, x) + tmp - Q;
	}
};

inline vector_col calc_z_vect(const double &t)
{
	vector_col z_vect(3);
	z_vect(0) = 1.15;
	z_vect(1) = 0;
	z_vect(2) = 11;
	return z_vect;
}

inline double f1()
{
	return 0;
}

inline double f2(const double &phi, const double &omega)
{
	return sys_alpha * omega * omega * sin(phi) * cos(phi) - sys_g * sin(phi);
}

inline double f3(const double &phi)
{
	return sys_gamma * cos(phi);
}

inline vector_col calc_f_vect(const vector_col &x)
{
	vector_col f_vect(3);
	f_vect(0) = f1();
	f_vect(1) = f2(x(0), x(2));
	f_vect(2) = f3(x(0));
	return f_vect;
}

class eq_H: public rk4 <vector_col, matrix>
{
public:
	vector_col f(const double &t, const vector_col &x, const vector_col *x_prev,
	             const vector_col *u_prev, const matrix *K, const vector_col *H)
	{
		return -prod(trans(a - prod(prod_brb, *K)), x) - prod(Q, calc_z_vect(t))+
		        prod(*K, calc_f_vect(*x_prev));
	}
};

class eq_X: public rk4 <vector_col, matrix>
{
public:
	vector_col f(const double &t, const vector_col &x, const vector_col *x_prev,
	             const vector_col *u_prev, const matrix *K, const vector_col *H)
	{
		vector_col u;
		calc_u(*H, *K, x, u);
		return prod(a, x) + prod(b, u) + calc_f_vect(*x_prev);
	}
	void calc_u(const vector_col &H, const matrix &K, const vector_col &x, vector_col &u)
	{
		u = prod(prod_rb, H - prod(K, x));
	}
};

inline double norm(const vector<vector_col> &v1, const vector<vector_col> &v2)
{
	double res = 0;
	for(int i = 0; i < v1.size(); i++)
	{
		int sz = v1[i].size();
		double sum = 0;
		for(int j = 0; j < sz; j++)
			sum += pow(v1[i](j)-v2[i](j), 2);
		sum = sqrt(sum);
		if(sum > res)
			res = sum;
	}
	return res;
}

void open_file(const string &coord, int num, ofstream &f)
{
	char fc[2], ch[2];
	fc[0] = coord[0];
	ch[0] = num + 49;
	fc[1] = ch[1] = '\0';
	string name = coord + ch;

	f.open((name + ".txt").c_str());
	f << "set term " << FORMAT << " size 1100, 700 font 20" << "\nset output \""<<
	     name << "." << FORMAT << "\"\n";
	f << "set xlabel \"t\"\n";
	f << "set ylabel \"" << name << "\"\n";

	f << "plot \"-\" with line lc 1 title \"" << fc << ch << "\"";
	if(coord.length() == 2)
		f << ", \"-\" with line lc rgb '#0' lw 2 title \"z" << ch << "\"";
	f << "\n";
}

inline void write_to_stream(ostringstream *ss, const double &t, const vector_col &coord)
{
	for(int i = 0, s = coord.size(); i < s; i++)
		ss[i] << t << " " << coord[i] << "\n";
}

inline void write_to_file(ofstream &f, const ostringstream &ss1, const ostringstream *ss2 = NULL)
{
	f << ss1.str();
	if(ss2 != NULL)
		f << "e\n" << (*ss2).str();
}

unsigned long GetTickCount()
{
	timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec * 1000 + tv.tv_usec / 1000;
}

int main()
{
	cout << "\nT > ";
	double T;
	cin >> T;
	cout << "n > ";
	int n;
	cin >> n;
	double h = T / n;
	cout << "eps > ";
	double eps;
	cin >> eps;

	double A[] = { 0,         1, 0,
	               0, -sys_beta, 0,
	               0,         0, 0 };
	create_matrix(A, 3, 3, a);
	a_tr = trans(a);

	double B[] = {      0,
	                    0,
	               -sys_b };
	create_matrix(B, 3, 1, b);
	matrix b_tr = trans(b);

	matrix inv_r(1, 1);
	inv_r(0, 0) = 1.0 / r;

	prod_rb  = prod(inv_r, b_tr);
	prod_brb = prod(b, prod_rb);

	identity I(3);
	Q = I;
	Q(0, 0) = q_1;
	Q(1, 1) = q_2;
	Q(2, 2) = q_3;

	matrix P = k_P * I;
	vector_col delta_z0(3);
	delta_z0(0) = delta_z0(1) = delta_z0(2) = d_z0;

	riccati ob_eq_ric;
	vector<matrix> K_t;
	cout << "\nSolving the Riccati equation\n";
	unsigned long t1 = GetTickCount();
	ob_eq_ric.solve(T, n, P, K_t, h, -1, false);
	cout << "K(0) = " << K_t[0] << "\n\n";

	try
	{
		vector_col z0 = calc_z_vect(0);
		vector_col zT = calc_z_vect(T);
		vector_col c  = k_z0 * z0 + delta_z0;

		vector<vector_col> x_prev_t, h_prev_t, u_prev_t, Z_t;
		vector<matrix>::iterator pK = K_t.begin();
		for(int num_t = 0; pK != K_t.end(); num_t++, pK++)
		{
			x_prev_t.push_back(c);
			h_prev_t.push_back(prod(P, zT));
			Z_t.push_back(calc_z_vect(num_t * h));
			u_prev_t.push_back(prod(prod_rb, prod(P, zT) - prod(*pK, c)));
		}

		eq_H ob_eq_H;
		vector<vector_col> H_t;
		eq_X ob_eq_X;
		vector<vector_col> X_t, u_t;
		int iter = 0;
		while(true)
		{
			iter++;
			cout << "Iteration " << iter << "\n";
			ob_eq_H.solve(T, n, prod(P, zT), H_t, h, -1, false, &x_prev_t, &u_prev_t, &K_t);
			ob_eq_X.solve(0, n, c, X_t, h, 1, false, &x_prev_t, &u_prev_t, &K_t, &H_t, &u_t);

			double error_x = norm(X_t, x_prev_t), error_h = norm(H_t, h_prev_t),
			       error_u = norm(u_t, u_prev_t);
			cout << "error_x = " << error_x << "\nerror_h = " << error_h <<
			        "\nerror_u = " << error_u << "\n\n";
			if(error_x < eps && error_h < eps && error_u < eps)
				break;

			copy(X_t.begin(), X_t.end(), x_prev_t.begin());
			copy(H_t.begin(), H_t.end(), h_prev_t.begin());
			copy(u_t.begin(), u_t.end(), u_prev_t.begin());
			X_t.clear();
			H_t.clear();
			u_t.clear();
		}
		cout << "The number of iterations = " << iter << "\n\n";
		double Time = (GetTickCount() - t1) / 60000.0;

		ostringstream sx[3], sz[3], su[1];
		double t;
		for(int i = 0; i <= n; i++)
		{
			t = i * h;
			cout << "X(" << t << ") = " << X_t[i];
			cout << "\nZ(" << t << ") = " << Z_t[i];
			cout << "\nu(" << t << ") = " << u_t[i] << "\n\n";
			write_to_stream(sx, t, X_t[i]);
			write_to_stream(sz, t, Z_t[i]);
			write_to_stream(su, t, u_t[i]);
		}
		cout << "Computing time is " << Time << " min.\n\n";

		ofstream xz[3], u[1];
		for(int i = 0; i < 3; i++)
		{
			open_file("xz", i, xz[i]);
			write_to_file(xz[i], sx[i], sz + i);
			xz[i].close();
		}
		for(int i = 0; i < 1; i++)
		{
			open_file("u", i, u[i]);
			write_to_file(u[i], su[i]);
			u[i].close();
		}
	}
	catch(const char *mes)
	{
		cout << mes << "\n\n";
	}
	return 0;
}
