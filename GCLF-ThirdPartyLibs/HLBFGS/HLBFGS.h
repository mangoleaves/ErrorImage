///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// HLBFGS                                                                    //
// http://www.loria.fr/~liuyang/software/HLBFGS/							 //
//                                                                           //
// HLBFGS is a hybrid L-BFGS optimization framework which unifies L-BFGS     //
// method, Preconditioned L-BFGS method and                                  //
// Preconditioned Conjugate Gradient method.                                 //
//                                                                           //
// Version 1.2                                                               //
// March 09, 2010                                                            //
//                                                                           //
// Copyright (C) 2009--2010                                                  //
// Yang Liu                                                                  //
//																			 //
// xueyuhanlang@gmail.com                                                    //
//                                                                           //
// HLBFGS is HLBFGS is freely available for non-commercial purposes.		 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef HLBFGS_H
#define HLBFGS_H

#include <vector>

//////////////////////////////////////////////////////////////////////////

//! ICFS_INFO stores ICFS's working arrays
class ICFS_INFO
{
public:
	ICFS_INFO()
	{
		p = 15;
	}
	~ICFS_INFO()
	{
	}
	void allocate_mem(int N)
	{
		if (N > 0)
		{
			lcol_ptr.resize(N + 1);
			ldiag.resize(N);
			iwa.resize(3 * N);
			wa1.resize(N);
			wa2.resize(N);
			r.resize(N);
			p = 15;
			CGP.resize(N);
			CGR.resize(N);
			CGQ.resize(N);
			CGZ.resize(N);
		}
	}
	int * get_lcol_ptr()
	{
		return &lcol_ptr[0];
	}
	int * get_lrow_ind()
	{
		return &lrow_ind[0];
	}
	double * get_ldiag()
	{
		return &ldiag[0];
	}
	double * get_l()
	{
		return &l[0];
	}
	int * get_iwa()
	{
		return &iwa[0];
	}
	double * get_wa1()
	{
		return &wa1[0];
	}
	double * get_wa2()
	{
		return &wa2[0];
	}
	int & get_p()
	{
		return p;
	}
	double * get_r()
	{
		return &r[0];
	}
	double * get_CGP()
	{
		return &CGP[0];
	}
	double * get_CGQ()
	{
		return &CGQ[0];
	}
	double * get_CGR()
	{
		return &CGR[0];
	}
	double * get_CGZ()
	{
		return &CGZ[0];
	}
	double & get_icfs_alpha()
	{
		return icfs_alpha;
	}
	void set_lrow_ind_size(int size)
	{
		lrow_ind.resize(size);
	}
	void set_l_size(int size)
	{
		l.resize(size);
	}
private:
	std::vector<int> lcol_ptr;
	std::vector<int> lrow_ind;
	std::vector<double> ldiag;
	std::vector<double> l;
	std::vector<int> iwa;
	std::vector<double> wa1;
	std::vector<double> wa2;
	int p;
	std::vector<double> r;
	std::vector<double> CGP;
	std::vector<double> CGQ;
	std::vector<double> CGR;
	std::vector<double> CGZ;
	double icfs_alpha;
};
//////////////////////////////////////////////////////////////////////////
//! Stores the pointers of hessian matrix
class HESSIAN_MATRIX
{
public:
	HESSIAN_MATRIX(int N)
	{
		n = N;
		nnz = 0;
		values = 0;
		rowind = 0;
		colptr = 0;
		diag = 0;
	}
	~HESSIAN_MATRIX()
	{

	}
	void set_dimension(int dim)
	{
		n = dim;
	}
	void set_nonzeros(int nz)
	{
		nnz = nz;
	}
	void set_values(double *array_values)
	{
		values = array_values;
	}
	void set_rowind(int *array_rowind)
	{
		rowind = array_rowind;
	}
	void set_colptr(int *array_colptr)
	{
		colptr = array_colptr;
	}
	void set_diag(double *array_diag)
	{
		diag = array_diag;
	}
	int get_dimension()
	{
		return n;
	}
	int get_nonzeros()
	{
		return nnz;
	}
	double * get_values()
	{
		return values;
	}
	int * get_rowind()
	{
		return rowind;
	}
	int * get_colptr()
	{
		return colptr;
	}
	double * get_diag()
	{
		return diag;
	}
	ICFS_INFO& get_icfs_info()
	{
		return l_info;
	}
private:
	int n;
	int nnz;
	double *values;
	int *rowind;
	int *colptr;
	double *diag;
	ICFS_INFO l_info;
};
//////////////////////////////////////////////////////////////////////////
//! HLBFGS initialization
//Dimension of arrays: 20, 20
void INIT_HLBFGS(double PARAMETERS[], int INFO[]);

void HLBFGS_MESSAGE(bool print, int id, const double PARAMETERS[]);
//////////////////////////////////////////////////////////////////////////
void HLBFGS_UPDATE_First_Step(int N, int M, double *q, double *s, double *y,
							  double *rho, double *alpha, int bound, int cur_pos, int iter);

void HLBFGS_UPDATE_Hessian(int N, int M, double *q, double *s, double *y,
						   int cur_pos, double *diag, int INFO[]);

void HLBFGS_UPDATE_Second_Step(int N, int M, double *q, double *s, double *y,
							   double *rho, double *alpha, int bound, int cur_pos, int iter);

void CONJUGATE_GRADIENT_UPDATE(int N, double *q, double *prev_q_update,
							   double *prev_q_first_stage, int INFO[]);
//////////////////////////////////////////////////////////////////////////
void HLBFGS_BUILD_HESSIAN_INFO(HESSIAN_MATRIX& m_hessian, int INFO[]);
//////////////////////////////////////////////////////////////////////////
//! HLBFGS functions
void HLBFGS(int N, int M, double *x, void EVALFUNC(int, double*, double*,
			double*, double*), void EVALFUNC_H(int, double*, double*, double*,
			double*, HESSIAN_MATRIX&), void USER_DEFINED_HLBFGS_UPDATE_H(int, int,
			double*, double*, double*, int, double*, int[]), void NEWITERATION(int,
			int, double*, double*, double*, double*), double PARAMETERS[],
			int INFO[]);
//////////////////////////////////////////////////////////////////////////

#include "LineSearch.h"
#include "ICFS.h"

template<typename EVALFUNC, typename NEWITERATION>
void HLBFGS_template(
	int N, int M, double* x,
	EVALFUNC& evalfunc,
	void EVALFUNC_H(int, double*, double*, double*, double*, HESSIAN_MATRIX&),
	void USER_DEFINED_HLBFGS_UPDATE_H(int, int, double*, double*, double*, int, double*, int[]),
	NEWITERATION& newiteration, double PARAMETERS[], int INFO[])
{
	int T = INFO[6];
	if (N < 1 || M < 0 || T < 0 || INFO[4] < 1)
	{
		HLBFGS_MESSAGE(INFO[5] != 0, 0, PARAMETERS);
		return;
	}
	//allocate mem
	std::vector<double> q_vec(N), g_vec(N), alpha_vec(M <= 0 ? 0 : M), rho_vec(
		M <= 0 ? 0 : M), s_vec(M <= 0 ? 0 : M * N), y_vec(M <= 0 ? 0 : M
		* N), prev_x_vec(N), prev_g_vec(N), diag_vec, wa_vec(N);
	double *q = &q_vec[0];
	double *g = &g_vec[0];
	double *alpha = M <= 0 ? 0 : &alpha_vec[0];
	double *rho = M <= 0 ? 0 : &rho_vec[0];
	double *s = M <= 0 ? 0 : &s_vec[0];
	double *y = M <= 0 ? 0 : &y_vec[0];
	double *prev_x = &prev_x_vec[0];
	double *prev_g = &prev_g_vec[0];
	double *diag = 0;
	double *wa = &wa_vec[0];
	double update_alpha = 1;
	HESSIAN_MATRIX m_hessian(N);
	if (INFO[3] == 1)
	{
		diag_vec.resize(N);
		diag = &diag_vec[0];
		std::fill(diag_vec.begin(), diag_vec.end(), 1.0);
	}
	std::vector<double> prev_q_first_stage_vec, prev_q_update_vec;
	double *prev_q_first_stage = 0;
	double *prev_q_update = 0;
	double scale = 0.0;
	double cg_dginit = 0;
	if (INFO[10] == 1)
	{
		if (INFO[11] == 1)
		{
			prev_q_first_stage_vec.resize(N);
			prev_q_first_stage = &prev_q_first_stage_vec[0];
		}
		prev_q_update_vec.resize(N);
		prev_q_update = &prev_q_update_vec[0];
	}

	//initialize
	INFO[1] = 0;
	INFO[2] = 0;
	double f = 0;
	int maxfev = INFO[0], bound = 0, nfev = 0, cur_pos = 0, start = 0;
	//line search parameters
	double stp, ftol = PARAMETERS[0], xtol = PARAMETERS[1], gtol =
		PARAMETERS[2], stpmin = PARAMETERS[3], stpmax = PARAMETERS[4];
	int info, keep[20];
	double gnorm, rkeep[40];
	std::fill(&rkeep[0], &rkeep[40], 0.0);
	std::fill(&keep[0], &keep[20], 0);

	m_hessian.get_icfs_info().allocate_mem(N);
	char task1 = 'N';
	char task2 = 'T';
	double prev_f;
	int i;

	//////////////////////////////////////////////////////////////////////////
	do
	{
		if (INFO[7] == 1 && ((T == 0) || (INFO[2] % T == 0)))
		{
			EVALFUNC_H(N, x, INFO[2] == 0 ? 0 : prev_x, &f, g, m_hessian);
			HLBFGS_BUILD_HESSIAN_INFO(m_hessian, INFO);
		}
		else if (INFO[2] == 0)
		{
			evalfunc(N, x, 0, &f, g);
			INFO[1]++;
		}

		if (INFO[2] > 0 && M > 0)
		{
			//compute s and y
			start = cur_pos * N;
#ifdef USE_OPENMP
#pragma omp parallel for private(i)
#endif
			for (i = 0; i < N; i++)
			{
				s[start + i] = x[i] - prev_x[i];
				y[start + i] = g[i] - prev_g[i];
			}
			rho[cur_pos] = 1.0 / HLBFGS_DDOT(N, &y[start], &s[start]);
			if (INFO[13] == 1)
			{
				update_alpha = 1.0 / (rho[cur_pos] * 6 * (prev_f - f
					+ HLBFGS_DDOT(N, g, &s[start])) - 2.0);
			}
			else if (INFO[13] == 2)
			{
				update_alpha = 1.0 / (rho[cur_pos] * 2 * (prev_f - f
					+ HLBFGS_DDOT(N, g, &s[start])));
			}
			else if (INFO[13] == 3)
			{
				update_alpha = 1.0 / (1 + rho[cur_pos] * (6 * (prev_f - f) + 3
					* (HLBFGS_DDOT(N, g, &s[start]) + HLBFGS_DDOT(N,
					prev_g, &s[start]))));
			}
			if (INFO[13] != 0)
			{
				if (update_alpha < 0.01)
				{
					update_alpha = 0.01;
				}
				else if (update_alpha > 100)
				{
					update_alpha = 100;
				}
				rho[cur_pos] *= update_alpha;
			}
		}
#ifdef USE_OPENMP
#pragma omp parallel for private(i)
#endif
		for (i = 0; i < N; i++)
		{
			q[i] = -g[i];
		}

		if (INFO[2] > 0 && M > 0)
		{
			bound = INFO[2] > M ? M - 1 : INFO[2] - 1;
			HLBFGS_UPDATE_First_Step(N, M, q, s, y, rho, alpha, bound, cur_pos,
				INFO[2]);
		}

		if (INFO[10] == 0)
		{
			if (INFO[7] == 1)
			{
				ICFS_INFO& l_info = m_hessian.get_icfs_info();
				dstrsol_(&N, l_info.get_l(), l_info.get_ldiag(),
					l_info.get_lcol_ptr(), l_info.get_lrow_ind(), q, &task1);
				dstrsol_(&N, l_info.get_l(), l_info.get_ldiag(),
					l_info.get_lcol_ptr(), l_info.get_lrow_ind(), q, &task2);
			}
			else
			{
				USER_DEFINED_HLBFGS_UPDATE_H(N, M, q, s, y, cur_pos, diag, INFO);
			}
		}
		else
		{
			if (INFO[7] == 1)
			{
				ICFS_INFO& l_info = m_hessian.get_icfs_info();
				dstrsol_(&N, l_info.get_l(), l_info.get_ldiag(),
					l_info.get_lcol_ptr(), l_info.get_lrow_ind(), q, &task1);
				CONJUGATE_GRADIENT_UPDATE(N, q, prev_q_update,
					prev_q_first_stage, INFO);
				cg_dginit = -HLBFGS_DDOT(N, q, q);
				dstrsol_(&N, l_info.get_l(), l_info.get_ldiag(),
					l_info.get_lcol_ptr(), l_info.get_lrow_ind(), q, &task2);
			}
			else
			{
				INFO[12] = 0;
				USER_DEFINED_HLBFGS_UPDATE_H(N, M, q, s, y, cur_pos, INFO[3]
				== 0 ? (&scale) : diag, INFO);
				if (INFO[3] == 0)
				{
					if (M > 0 && INFO[2] > 0 && scale != 1.0)
					{
						scale = std::sqrt(scale);
						HLBFGS_DSCAL(N, scale, q);
					}
					CONJUGATE_GRADIENT_UPDATE(N, q, prev_q_update,
						prev_q_first_stage, INFO);
					cg_dginit = -HLBFGS_DDOT(N, q, q);
					if (M > 0 && INFO[2] > 0 && scale != 1.0)
						HLBFGS_DSCAL(N, scale, q);
				}
				else
				{
					if (M > 0 && INFO[2] > 0)
					{
						//use prev_g as temporary array
						for (int i = 0; i < N; i++)
						{
							prev_g[i] = std::sqrt(diag[i]);
							q[i] *= prev_g[i];
						}
					}
					CONJUGATE_GRADIENT_UPDATE(N, q, prev_q_update,
						prev_q_first_stage, INFO);
					cg_dginit = -HLBFGS_DDOT(N, q, q);
					if (M > 0 && INFO[2] > 0)
					{
#ifdef USE_OPENMP
#pragma omp parallel for private(i)
#endif
						for (i = 0; i < N; i++)
						{
							q[i] *= prev_g[i];
						}
					}

				}
				INFO[12] = 1;
			}
		}

		if (INFO[2] > 0 && M > 0)
		{
			HLBFGS_UPDATE_Second_Step(N, M, q, s, y, rho, alpha, bound,
				cur_pos, INFO[2]);

			cur_pos = (cur_pos + 1) % M;
		}

		//store g and x
		std::copy(&x[0], &x[N], prev_x_vec.begin());
		std::copy(&g[0], &g[N], prev_g_vec.begin());
		prev_f = f;
		//linesearch, find new x
		bool blinesearch = true;
		if (INFO[2] == 0)
		{
			gnorm = HLBFGS_DNRM2(N, g);
			//if(gnorm > 1)
			stp = 1.0 / gnorm;
			//else
			//	stp = 1;
		}
		else
		{
			stp = 1;
		}

		info = 0;

		do
		{
			MCSRCH(&N, x, &f, g, q, &stp, &ftol, &gtol, &xtol, &stpmin,
				&stpmax, &maxfev, &info, &nfev, wa, keep, rkeep, INFO[10]
			== 0 ? 0 : (&cg_dginit));
			blinesearch = (info == -1);
			if (blinesearch)
			{
				evalfunc(N, x, prev_x, &f, g);
				INFO[1]++;
			}

			if (INFO[9] == 1 && prev_f > f) //modify line search to avoid too many function calls
			{
				info = 1;
				break;
			}

		} while (blinesearch);

		gnorm = HLBFGS_DNRM2(N, g);
		INFO[2]++;
		newiteration(INFO[2], INFO[1], x, &f, g, &gnorm);
		double xnorm = HLBFGS_DNRM2(N, x);
		xnorm = 1 > xnorm ? 1 : xnorm;
		rkeep[2] = gnorm;
		rkeep[8] = xnorm;

		if (info != 1)
		{
			HLBFGS_MESSAGE(INFO[5] != 0, 1, PARAMETERS);
			break;
		}
		if (gnorm / xnorm <= PARAMETERS[5])
		{
			HLBFGS_MESSAGE(INFO[5] != 0, 2, PARAMETERS);
			break;
		}
		if (gnorm < PARAMETERS[6])
		{
			HLBFGS_MESSAGE(INFO[5] != 0, 3, PARAMETERS);
			break;
		}
		if (stp < stpmin || stp > stpmax)
		{
			HLBFGS_MESSAGE(INFO[5] != 0, 4, PARAMETERS);
			break;
		}
		if (INFO[2] > INFO[4])
		{
			HLBFGS_MESSAGE(INFO[5] != 0, 5, PARAMETERS);
			break;
		}

	} while (true);
}

#endif
