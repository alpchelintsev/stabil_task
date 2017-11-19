// (c) 2017 Alexander Pchelintsev, pchelintsev.an@yandex.ru

#ifndef _RK4_H_INCLUDED_
#define _RK4_H_INCLUDED_

#include <iostream>
#include <vector>

template <class TypeSol, class TypeK>
class rk4
{
public:
	void solve(const double &t0, int n, const TypeSol &x0,
	           std::vector<TypeSol> &x, const double &h, int d = 1, bool info = true,
	           std::vector<TypeSol> *x_prev_t = NULL, std::vector<TypeSol> *u_prev_t = NULL,
	           std::vector<TypeK> *K_t = NULL, std::vector<TypeSol> *H_t = NULL,
	           std::vector<TypeSol> *u_t = NULL)
	{
		x.push_back(x0);
		TypeSol k1, k2, k3, k4, x_c, x_new;
		double t = t0, s = d*h, s2 = s/2;
		TypeSol *x_prev, *u_prev, *H = NULL;
		TypeK *K;
		typename std::vector<TypeSol>::iterator px, pu, pH;
		typename std::vector<TypeK>::iterator pK;
		bool prev = x_prev_t != NULL && u_prev_t != NULL && K_t != NULL;
		bool isH = d == 1 && H_t != NULL && u_t != NULL;
		if(prev)
		{
			if(d == 1)
			{
				px = (*x_prev_t).begin();
				pu = (*u_prev_t).begin();
				pK = (*K_t).begin();
				if(isH)
				{
					pH = (*H_t).begin();
					H = new TypeSol(*pH);
				}
			}
			else
			{
				px = (*x_prev_t).end() - 1;
				pu = (*u_prev_t).end() - 1;
				pK = (*K_t).end() - 1;
			}
			x_prev = new TypeSol(*px);
			u_prev = new TypeSol(*pu);
			K = new TypeK(*pK);
		}
		if(isH)
		{
			TypeSol u;
			calc_u(*H, *K, x0, u);
			(*u_t).push_back(u);
		}
		for(int i = 1; i <= n; i++)
		{
			if(d == 1)
				x_c = *(x.end() - 1);
			else
				x_c = *x.begin();
			if(prev)
			{
				k1 = f(t, x_c, x_prev, u_prev, K, H);
				if(d == 1)
				{
					px++; pu++; pK++;
					if(isH)
						pH++;
					bool not_end = px != (*x_prev_t).end();
					if(not_end)
					{
						*x_prev = (*x_prev + *px) / 2;
						*u_prev = (*u_prev + *pu) / 2;
						*K = (*K + *pK) / 2;
						if(isH)
							*H = (*H + *pH) / 2;
					}
					k2 = f(t + s2, x_c + s2 * k1, x_prev, u_prev, K, H);
					k3 = f(t + s2, x_c + s2 * k2, x_prev, u_prev, K, H);
					if(not_end)
					{
						*x_prev = *px;
						*u_prev = *pu;
						*K = *pK;
						if(isH)
							*H = *pH;
					}
					k4 = f(t + s, x_c + s * k3, x_prev, u_prev, K, H);
				}
				else
				{
					if(px != (*x_prev_t).begin())
					{
						px--; pu--; pK--;
						*x_prev = (*x_prev + *px) / 2;
						*u_prev = (*u_prev + *pu) / 2;
						*K = (*K + *pK) / 2;
					}
					k2 = f(t + s2, x_c + s2 * k1, x_prev, u_prev, K);
					k3 = f(t + s2, x_c + s2 * k2, x_prev, u_prev, K);
					*x_prev = *px;
					*u_prev = *pu;
					*K = *pK;
					k4 = f(t + s, x_c + s * k3, x_prev, u_prev, K);
				}
			}
			else
			{
				k1 = f(t, x_c);
				k2 = f(t + s2, x_c + s2 * k1);
				k3 = f(t + s2, x_c + s2 * k2);
				k4 = f(t + s, x_c + s * k3);
			}
			x_new = x_c + s/6 * (k1 + 2*k2 + 2*k3 + k4);
			if(d == 1)
				x.push_back(x_new);
			else
				x.insert(x.begin(), x_new);
			if(isH)
			{
				TypeSol u;
				calc_u(*H, *K, x_new, u);
				(*u_t).push_back(u);
			}
			t += s;
			if(info)
				std::cout << "Calculated " << t << "\n";
		}
		if(prev)
		{
			delete x_prev;
			delete u_prev;
			delete K;
			if(isH)
				delete H;
		}
	}
	virtual TypeSol f(const double &t, const TypeSol &x,
	                  const TypeSol *x_prev = NULL,
	                  const TypeSol *u_prev = NULL,
	                  const TypeK *K = NULL,
	                  const TypeSol *H = NULL) = 0;
	virtual void calc_u(const TypeSol &H, const TypeK &K, const TypeSol &x, TypeSol &u)
	{
	}
};

#endif //_RK4_H_INCLUDED_
