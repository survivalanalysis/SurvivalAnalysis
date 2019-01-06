#include <stdio.h>
#include <stdlib.h>
/**#include <math.h>**/
#include <Rcpp.h>
#include <cmath>
#include <math.h>
#if defined(_OPENMP)
#include <omp.h>
#endif
/**
#include "frain.h"
**/
using namespace Rcpp;
using namespace std;
//#define E 2.7182818284590452353602874713527


#define likely(x) __builtin_expect ((x), 1)
#define unlikely(x) __builtin_expect ((x), 0)

#define Y(times,i,t) times[i] >= t

//#define MY_ADD_LN(x,y) (isinf(x) && isinf(y)? y: (x > y? x + log(1+exp(y-x)): y + log(1 + exp(x-y))))
double MY_ADD_LN(double x,double y) {
	if (isinf(x) && isinf(y)) {
		return y;
	} else {
		if (x > y) {
			return x + log(1+exp(y-x));
		}
		else {
			return y + log(1 + exp(x-y));
		}
	}
}

double* R2C_vec(NumericVector a) {
	int n = a.size();
	double* a_c = (double*)malloc(sizeof(double)*n);
	for(int i = 0; i < n; ++i) {
        	a_c[i] = a[i];
  	}
	return a_c;
}

double* R2C_mat(NumericMatrix a) {
	int n = a.nrow();
	int m = a.ncol();
	double* a_c = (double*)malloc(sizeof(double)*n*m);
	for(int i = 0; i < n; ++i) {
		for(int j = 0; j < m; ++j) {
        		a_c[i*m + j] = a(i,j);
		}
  	}
	return a_c;
}


/**#define N(times,delta,i,t) (times[i] < t && delta[i] != 0)? 1:0 
**/

int N(double* times, double* delta, int i, double t) {
	if (times[i] < t and delta[i] > 0) {
		return 1;
	}
	return 0;
}

inline
double exp4(register double x) {
  x = 1.0 + x / 8192;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  x *= x;
  return x;
}

#define my_exp(x) exp4(x)


inline int find_m_pos(double f, double* __restrict__ arr, int size) {
	int n = size/2;	
	while (!(arr[n] <= f and f <= arr[n + 1])) {
		if (f > arr[n]) n = size + (size-n) / 2;
		else n = n/2;
	}

	return n;
}


int binarySearch(double* __restrict__ arr, int l, int r, double x)
{
   if (r >= l)
   {
        int mid = l + (r - l)/2;
 
        // If the element is present at the middle itself
        if (arr[mid] == x)  return mid;
 
        // If element is smaller than mid, then it can only be present
        // in left subarray
        if (arr[mid] > x) return binarySearch(arr, l, mid-1, x);
 
        // Else the element can only be present in right subarray
        return binarySearch(arr, mid+1, r, x);
   }
 
   // We reach here when element is not present in array
   return r;
}

/**
int sum_e_matrix(double* __restrict__ a, int n, double* __restrict__  b,int m, double* times, double* __restrict__  sortedTimes, double* __restrict__ res)
{
	sum_e_matrix2(a,n,b,m,times,sortedTimes,res);
	return 0;
}**/


int sum_e_matrix(double* __restrict__ a, int n, double* __restrict__  b,int m, double* times, double* __restrict__  sortedTimes, double* __restrict__ res)
{
        for (int i=0; i < n; i++) {
                res[i] = 0;
        }

        double diff[m];

        for (int i=1; i < m; i++) {
                diff[i] = b[i]-b[i-1];
        }
        diff[0] = b[0];

        for (int i=0; i < n; i++) {
                register double sum = 0;
                register double cur_a = a[i];
                register double cur_t = times[i];
                if (unlikely(cur_t < sortedTimes[0])) {
                        sum += (cur_t/sortedTimes[0])*diff[0];
                } else {
                        sum += diff[0];
                }
                res[i] = sum;
        }

        #pragma omp parallel for
        for (int i=0; i < n; i++) {
                register double sum = 0;
                register double cur_a = a[i];
                register double cur_t = times[i];
                if (cur_t < sortedTimes[0]) {
                        continue;
                }
                for (register int j=1;  j< m; j++) {
                        if (unlikely(cur_t <= sortedTimes[j])) {
                                sum += exp(b[j - 1]*cur_a)*diff[j]*((cur_t-sortedTimes[j-1])/(sortedTimes[j]-sortedTimes[j-1]));
                                break;
                        }
                        sum += exp(b[j - 1]*cur_a)*diff[j];
                }
                /**if (cur_t > sortedTimes[m-1]) {
                        sum += exp(b[m-1]*cur_a)*diff[m-1]*((cur_t-sortedTimes[m-1])/(sortedTimes[m-1]-sortedTimes[m-2]));
                }**/

                res[i] += sum;
        }

        return 0;
}

/**
int ln_sum_e_matrix(double* __restrict__ a, int n, double* __restrict__  b,int m, double* times, double* __restrict__  sortedTimes, double* __restrict__ res)
{
	ln_sum_e_matrix2(a,n,b,m,times,sortedTimes,res);
	return 0;
}**/

int ln_sum_e_matrix(double* __restrict__ a, int n, double* __restrict__  b,int m, double* times, double* __restrict__  sortedTimes, double* __restrict__ res)
{
        for (int i=0; i < n; i++) {
                res[i] = 0;
        }

        double diff[m];

        for (int i=1; i < m; i++) {
                diff[i] = b[i]-b[i-1];
        }
        diff[0] = b[0];

        for (int i=0; i < n; i++) {
                register double sum = 0;
                register double cur_a = a[i];
                register double cur_t = times[i];
                if (unlikely(cur_t < sortedTimes[0])) {
                        sum += (cur_t/sortedTimes[0])*diff[0];
                } else {
                        sum += diff[0];
                }
                res[i] = log(sum);
        }

        #pragma omp parallel for
        for (int i=0; i < n; i++) {
                register double sum = res[i];
                register double cur_a = a[i];
                register double cur_t = times[i];
                if (cur_t < sortedTimes[0]) {
                        continue;
                }
                for (register int j=1;  j< m; j++) {
                        if (unlikely(cur_t <= sortedTimes[j])) {
								register double tmp = b[j - 1]*cur_a+ log(diff[j]*((cur_t-sortedTimes[j-1])/(sortedTimes[j]-sortedTimes[j-1])));
                            if (sum > tmp) {
								sum = sum + log(1 + exp(tmp-sum));
							}
							else {
								sum = tmp + log(1 + exp(sum-tmp));
							}
							break;
                        }
					register double tmp = b[j - 1]*cur_a + log(diff[j]);
					if (sum > tmp) {
						sum = sum + log(1 + exp(tmp-sum));
					}
					else {
						sum = tmp + log(1 + exp(sum-tmp));
					}
                }

		/** Decided not to do chupchik after last event time**/
                /**if (cur_t > sortedTimes[m-1]) {
                        sum += exp(b[m-1]*cur_a)*diff[m-1]*((cur_t-sortedTimes[m-1])/(sortedTimes[m-1]-sortedTimes[m-2]));
                }**/

                res[i] = sum;
        }

        return 0;
}





int get_gz_vec(double* __restrict__ z,  double* __restrict__  g, double* __restrict__  res,int n_rows, int n_variables) {
	for (int i=0; i < n_rows; i++) {
		res[i] = 0;
		for (int j=0; j < n_variables; j++) {
			res[i] += z[i*n_variables + j] * g[j];
		}
	}
	return 0;	
}


int get_exp_vec(double* __restrict__ in,  double* __restrict__  out,int len) {
	for (int i=0; i < len; i++) {
		out[i] = exp(in[i]);
	}
	return 0;	
}
/**

int derive_A_combined_c(
			double theta,
			double* __restrict__ z,
			double* __restrict__  g,
			double* __restrict__  H_times,
			double* __restrict__ times,
			double* __restrict__  sortedTimes,
			double* __restrict__ res,
			int n_rows,
			int n_columns,
			int n_variables)
{
	cuda_derive_A_combined_c(theta,
			z,
			g,
			H_times,
			times,
			sortedTimes,
			res,
			n_rows,
			n_columns,
			n_variables);

	return 0;
}
**/


// Derive A matrix by theta and all gamma variales 
// Res expected to be a n*p+1 matrix where p is the number of explenatory variables.
// The first column of res will be the derivative by theta
// The next columns will be the derivative by gamma
int derive_A_combined_c(
			double theta,
			double* __restrict__ z,
			double* __restrict__  g,
			double* __restrict__  H_times,
			double* __restrict__ times,
			double* __restrict__  sortedTimes,
			double* __restrict__ res,
			int n_rows,
			int n_columns,
			int n_variables)
{
	for (int i=0; i < n_rows; i++) {
		for (int j=0; j <= n_variables + 1; j++){
			res[i*(n_variables + 2) + j] = 0;
		}
	}


	double diff[n_columns];
	double gz_vec[n_rows];
	double gz_vec_exp[n_rows];

        for (int i=1; i < n_columns; i++) {
                diff[i] = H_times[i]-H_times[i-1];
        }
        diff[0] = H_times[0];

	get_gz_vec(z,g,gz_vec,n_rows,n_variables);
	get_exp_vec(gz_vec,gz_vec_exp,n_rows);

        for (int i=0; i < n_rows; i++) {
                register double val = 0;
                register double cur_t = times[i];
                if (unlikely(cur_t < sortedTimes[0])) {
                        val += exp(gz_vec[i])*(cur_t/sortedTimes[0])*diff[0];
                } else {
                        val += exp(gz_vec[i])*diff[0];
                }
		res[i*(n_variables + 2)] = log(val);
		int tmp = 0;
		res[i*(n_variables + 2) + 1] = log(tmp);
		for (register int k=2; k <= n_variables + 1; k++) {
			res[i*(n_variables + 2) + k] = log(val*(abs(z[i*n_variables + k-2]))); 
		}

        }


        #pragma omp parallel for
        for (int i=0; i < n_rows; i++) {
		 register double val = 0;
                register double cur_a = theta*gz_vec_exp[i];
                register double cur_t = times[i];
                if (cur_t < sortedTimes[0]) {
                        continue;
                }
                for (register int j=1;  j< n_columns; j++) {
                        if (unlikely(cur_t <= sortedTimes[j])) {
							val = gz_vec[i] + (H_times[j - 1]*cur_a)  + log(diff[j]*((cur_t-sortedTimes[j-1])/(sortedTimes[j]-sortedTimes[j-1])));
                        } else {
                        	val = gz_vec[i] + (H_times[j - 1]*cur_a) + log(diff[j]);
			}
			res[i*(n_variables + 2)] = MY_ADD_LN(res[i*(n_variables + 2)], val);			
			res[i*(n_variables + 2) + 1] = MY_ADD_LN(res[i*(n_variables + 2) + 1], val + log(gz_vec_exp[i]*H_times[j-1]));
						
			for (register int k=2; k <= n_variables + 1; k++) {
				res[i*(n_variables + 2) + k] = MY_ADD_LN(res[i*(n_variables + 2) + k], val + log((abs(z[i*n_variables + k-2]) + theta*abs(z[i*n_variables + k-2])*gz_vec_exp[i]*H_times[j-1]))); 
			}
			
			if (unlikely(cur_t <= sortedTimes[j])) break;
                }
                /**if (cur_t > sortedTimes[n_columns-1]) {
                        val = exp(gz_vec[i])*exp(H_times[n_columns-1]*cur_a)*diff[n_columns-1]*((cur_t-sortedTimes[n_columns-1])/(sortedTimes[n_columns-1]-sortedTimes[n_columns-2]));
			res[i*(n_variables + 2)] += val;
			res[i*(n_variables + 2) + 1] += val*gz_vec_exp[i]*H_times[n_columns-1];
			for (register int k=2; k <= n_variables + 1; k++) {
				res[i*(n_variables + 2) + k] += val*(z[i*n_variables + k-2] + theta*z[i*n_variables + k-2]*gz_vec_exp[i]*H_times[n_columns-1]); 
			}
                }**/



        }
	return 0;
}

// [[Rcpp::export]]
NumericMatrix derive_A_combined(NumericVector theta_hat, NumericMatrix z, NumericVector g, NumericVector H_times, NumericVector times, NumericVector sortedTimes) {
  int n_rows = z.nrow();
  int n_columns = H_times.size();
  int n_variables = z.ncol();
  int tmp = g.size();
  //printf("*** %d\n", tmp);

  double* z_c = R2C_mat(z);
  double* g_c = R2C_vec(g);
 
  double* H_times_c = R2C_vec(H_times);
  double* times_c =  R2C_vec(times);
  double* sortedTimes_c = R2C_vec(sortedTimes);
  //double out_c[n_rows][n_variables + 1];
  double out_c[n_rows * (n_variables + 2)];
  double theta = theta_hat[0];

  NumericMatrix out(n_rows,n_variables + 2);
  //printf("1\n");
  derive_A_combined_c(theta,z_c,g_c,H_times_c,times_c,sortedTimes_c,(double*)out_c,n_rows,n_columns,n_variables);
  //printf("2\n");

  free(z_c);
  free(g_c);
  free(H_times_c);
  free(times_c);
  free(sortedTimes_c);

  for(int i = 0; i < n_rows; ++i) {
	for(int j = 0; j <= n_variables + 1; ++j) {
    		out(i,j) = out_c[i*(n_variables + 2) + j];
	}
  }

  return out;
}





inline void a_star_at_t(int t_pos,
			 double* __restrict__ gz,
                         double theta_hat,
                         double* __restrict__ Htimes,
                         double* __restrict__ res ,
                         int length) {
		int i = 0;
		for (;i < length; i++) {
			res[i] = exp(gz[i]);
			res[i] *= exp(theta_hat*res[i]*Htimes[t_pos]);
		}
	}

inline void diff_vec(double* __restrict__ vec , double* __restrict__ ret, int size) {
	for (int i = 0; i < size; i++) {
		ret[i] = vec[i + 1] - vec[i];
	}
}

inline void diff_vec_compete(double* __restrict__ vec , double* __restrict__ ret, int size) {
	for (int i = 0; i < size - 1; i++) {
		ret[i] = vec[i+1] - vec[i];
	}
}

inline void copy_num_vector(NumericVector vec , double* __restrict__ ret) {
	int size = vec.size();
	for (int i = 0; i < size; i++) {
		ret[i] = vec[i];
	}
}


struct competeA {
	int pos;
	double t;
	double val;
};

inline void update_compete_A(competeA* ret, int n_ev_compete, double gz_exp, double theta_hat, double* __restrict__ estimated_H_compete, double* diff_compete, double* compete_sorted_time, double t) {
	int pos = ret->pos;
	if (pos == n_ev_compete-1) {
		//ret->val += gz_exp * exp(theta_hat*gz_exp*estimated_H_compete[pos])* diff_compete[pos-1] * (t - ret->t) / (compete_sorted_time[pos] - compete_sorted_time[pos-1]);
		ret->t = t;
		return;
	}
	if (t < compete_sorted_time[pos + 1]){
		ret->val += gz_exp * exp(theta_hat*gz_exp*estimated_H_compete[pos])* diff_compete[pos] * (t - ret->t) / (compete_sorted_time[pos+1] - compete_sorted_time[pos]);
		ret->t = t;
		return;
	}
	ret->val += gz_exp * exp(theta_hat*gz_exp*estimated_H_compete[pos])* diff_compete[pos] * (compete_sorted_time[pos + 1] - ret->t) / (compete_sorted_time[pos+1] - compete_sorted_time[pos]);

	pos += 1;
	while (t > compete_sorted_time[pos + 1] & pos != n_ev_compete - 1) {
		ret->val += gz_exp * exp(theta_hat*gz_exp*estimated_H_compete[pos])* diff_compete[pos];
		pos += 1;
	}
	if (pos == n_ev_compete -1) {
		//ret->val += gz_exp * exp(theta_hat*gz_exp*estimated_H_compete[pos])* diff_compete[pos-1] * (t - compete_sorted_time[pos]) / (compete_sorted_time[pos] - compete_sorted_time[pos-1]);
	} else {	
		ret->val += gz_exp * exp(theta_hat*gz_exp*estimated_H_compete[pos])* diff_compete[pos] * (t - compete_sorted_time[pos]) / (compete_sorted_time[pos+1] - compete_sorted_time[pos]);
	}
	ret->pos = pos;
	ret->t = t;

}


void c_estimate_H_JK(double* __restrict__ gz,
			  double* __restrict__ gz_compete,
                          double theta_hat,
                          double* __restrict__ estimated_H,
			  double* __restrict__ estimated_H_compete,
			  double* __restrict__ sorted_ev_times,
			  double* __restrict__ sorted_ev_times_compete,
			  int n_ev,
			  int n_ev_compete,
			  double* __restrict__ times,
			  double* __restrict__ delta,
			  int n_subjects,
			  double* weights_subjects,
			  double* weights_events,
			  double* __restrict__ ret) {

	double A[n_subjects];
	double gz_vec_exp[n_subjects];
	double gz_vec_exp_compete[n_subjects];

	double diff[n_ev], diff_compete[n_ev_compete];

	struct competeA compA[n_subjects];
	for (int i = 0; i < n_subjects; i++) {
		compA[i].pos = 0;
		compA[i].val = 0;
		compA[i].t = 0;
		A[i] = 0;
	}
		
	double E_val = 0;
	register int event_itr = 0;
	register int subjet_itr = 0;


	diff_vec(estimated_H,diff,n_ev);
	diff_vec_compete(estimated_H_compete,diff_compete,n_ev_compete);
	get_exp_vec(gz,gz_vec_exp,n_subjects);
	get_exp_vec(gz_compete,gz_vec_exp_compete,n_subjects);

	for (event_itr = 0; event_itr < n_ev; event_itr++ ) {
		ret[event_itr] = 0;
	}

	//In order to work in parallel we calculate each row of A and add the contributon of this sample to ret (we sum)
	//But we must keep diffrent vectors for ret for diffrent threads so as to not cause collison vetween the threads when they try to update the ret of the same even time
	//So we create a new array(only if using parralization) the size of events*number of threads each thread will work on a shifted temp array
	//In the end a single thread will merge them.
	#if defined(_OPENMP)
	   double* S_private = NULL;	
	#endif	
	
	#pragma omp parallel private(event_itr,subjet_itr,E_val)
	{
		#if defined(_OPENMP)
		    const int nthreads = omp_get_num_threads();
		    const int ithread = omp_get_thread_num();

		    #pragma omp single 
		    {
			S_private = new double[n_ev*nthreads];
			for(int i=0; i<(n_ev*nthreads); i++) S_private[i] = 0;
		    }
		#endif	


		#pragma omp for 
		for (subjet_itr = 0; subjet_itr < n_subjects; subjet_itr++) {
			for (event_itr = 0; event_itr < n_ev; event_itr++ ) {
				register double astar = -1;
				if (Y(times,subjet_itr,sorted_ev_times[event_itr])) {
					if (astar == -1) {astar = gz_vec_exp[subjet_itr] * exp(theta_hat*gz_vec_exp[subjet_itr]*estimated_H[event_itr]);}
					double tmp = compA[subjet_itr].val;
					update_compete_A(&(compA[subjet_itr]),n_ev_compete, gz_vec_exp_compete[subjet_itr], theta_hat , estimated_H_compete ,diff_compete, sorted_ev_times_compete, sorted_ev_times[event_itr]);				
					A[subjet_itr] += (compA[subjet_itr].val - tmp);
					A[subjet_itr] += astar* diff[event_itr];

				} else {
					break;
				}
				if (Y(times,subjet_itr,sorted_ev_times[event_itr])) {
					//if (astar == -1) {astar = gz_vec_exp[subjet_itr] * exp(theta_hat*gz_vec_exp[subjet_itr]*estimated_H[event_itr]);}
					astar = gz_vec_exp[subjet_itr] * exp(theta_hat*gz_vec_exp[subjet_itr]*estimated_H[event_itr]);		
					E_val = (1.0/theta_hat + N(times, delta ,subjet_itr, sorted_ev_times[event_itr])) / (1.0/theta_hat + A[subjet_itr]);
					#if defined(_OPENMP)		
					S_private[ithread*n_ev +event_itr] += E_val*astar*weights_subjects[subjet_itr];
					#else
					ret[event_itr] += E_val*astar*weights_subjects[subjet_itr];
					#endif				
				}			
			}
		}
	
	    #if defined(_OPENMP)
	    #pragma omp single
	    {	
	    	for(int i=0; i<n_ev; i++) {
			for(int t=0; t<nthreads; t++) {
			    ret[i] += S_private[n_ev*t + i];
			}
	    	}
	    	delete S_private;
	    }
	    #endif

	}

	for (event_itr = 0; event_itr < n_ev; event_itr++ ) {
		//printf("%f\n",ret[event_itr]);
		ret[event_itr] = weights_events[event_itr] / ret[event_itr];
	}
}


// [[Rcpp::export]]
NumericVector estimate_H_JK(
			NumericVector gz,
			NumericVector gz_competing,
			NumericVector theta_hat,
			NumericVector sorted_ev_times,
			NumericVector sorted_ev_times_compete,
			NumericVector estimated_H,
			NumericVector competing_estimated_H,
			NumericVector times,
			NumericVector delta,
			NumericVector weights_subjects,
			NumericVector weights_events
) {

  int n_ev = estimated_H.size() - 1;
  int n_ev_compete = competing_estimated_H.size();
  int n_times = times.size();

  double gz_c[n_times];
  double gz_competing_c[n_times];

  double sorted_ev_times_c[n_ev];
  double sorted_ev_times_compete_c[n_ev_compete];

  double estimated_H_c[n_ev + 1];
  double competing_estimated_H_c[n_ev_compete];

  double times_c[n_times];
  double delta_c[n_times];

  double theta = theta_hat[0];

  double weights_subjects_c[n_times];
  double weights_events_c[n_ev];



  double ret[n_ev];

  NumericVector out(n_ev);



  copy_num_vector(gz,gz_c);
  copy_num_vector(gz_competing,gz_competing_c);
  copy_num_vector(sorted_ev_times,sorted_ev_times_c);
  copy_num_vector(sorted_ev_times_compete,sorted_ev_times_compete_c);
  copy_num_vector(estimated_H,estimated_H_c);
  copy_num_vector(competing_estimated_H,competing_estimated_H_c);
  copy_num_vector(times,times_c);
  copy_num_vector(delta,delta_c);
  copy_num_vector(weights_subjects,weights_subjects_c);
  copy_num_vector(weights_events,weights_events_c);

  c_estimate_H_JK(gz_c,
		  gz_competing_c,
		  theta,
		  estimated_H_c,
		  competing_estimated_H_c,
		  sorted_ev_times_c,
		  sorted_ev_times_compete_c,
		  n_ev,
		  n_ev_compete,
		  times_c,
		  delta_c,
		  n_times,
		  weights_subjects_c,
		  weights_events_c,
		  ret);


  for(int i = 0; i < n_ev; i++) {
    out[i] = ret[i];
  }

  return out;
}


#define MAX(a,b) ((a) > (b) ? a : b)

void c_estimate_H_23(double* __restrict__ gz,
                          double theta_hat,
                          double* __restrict__ estimated_H,
			  double* __restrict__ sorted_ev_times,
			  int n_ev,
			  double* __restrict__ times,
			  double* __restrict__ V_times,
			  double* __restrict__ delta,
			  double* __restrict__ N_1_tau,
			  double* __restrict__ A_1_tau,
			  int n_subjects,
			  double* weights_subjects,
			  double* weights_events,
			  double* __restrict__ ret) {

	double a_star[n_subjects];
	double A[n_subjects];
	double gz_vec_exp[n_subjects];
	double diff[n_ev];
		
	double E_val = 0;
	int event_itr = 0;
	int subjet_itr = 0;


	diff_vec(estimated_H,diff,n_ev);
	get_exp_vec(gz,gz_vec_exp,n_subjects);
	
	for (event_itr = 0; event_itr < n_ev; event_itr++ ) {
		ret[event_itr] = 0;
	}


	ret[0] = 0;
	for (subjet_itr = 0; subjet_itr < n_subjects; subjet_itr++) {
		A[subjet_itr] = A_1_tau[subjet_itr];
	}


	//In order to work in parallel we calculate each row of A and add the contributon of this sample to ret (we sum)
	//But we must keep diffrent vectors for ret for diffrent threads so as to not cause collison vetween the threads when they try to update the ret of the same even time
	//So we create a new array(only if using parralization) the size of events*number of threads each thread will work on a shifted temp array
	//In the end a single thread will merge them.
	#if defined(_OPENMP)
	   double* S_private = NULL;	
	#endif	
	
	#pragma omp parallel private(event_itr,subjet_itr,E_val)
	{
		#if defined(_OPENMP)
		    const int nthreads = omp_get_num_threads();
		    const int ithread = omp_get_thread_num();

		    #pragma omp single 
		    {
			S_private = new double[n_ev*nthreads];
			for(int i=0; i<(n_ev*nthreads); i++) S_private[i] = 0;
		    }
		#endif	

		#pragma omp for 
		for (subjet_itr = 0; subjet_itr < n_subjects; subjet_itr++) {
			for (event_itr = 0; event_itr < n_ev; event_itr++ ) {
				register double astar = -1;
				if (Y(times,subjet_itr,sorted_ev_times[event_itr])) {
					if (astar == -1) {astar = gz_vec_exp[subjet_itr] * exp(theta_hat*gz_vec_exp[subjet_itr]*estimated_H[event_itr]);}	
					A[subjet_itr] += astar * diff[event_itr];
				}
				else {
					break;
				}
				if (Y(times,subjet_itr,sorted_ev_times[event_itr]) && (sorted_ev_times[event_itr] >= V_times[subjet_itr])) {
					astar = gz_vec_exp[subjet_itr] * exp(theta_hat*gz_vec_exp[subjet_itr]*estimated_H[event_itr]);
					E_val = (1.0/theta_hat + N_1_tau[subjet_itr] + N(times, delta ,subjet_itr, sorted_ev_times[event_itr])) / (1.0/theta_hat + A[subjet_itr]);
					//E_val = (1.0/theta_hat + N(times, delta ,subjet_itr, sorted_ev_times[event_itr])) / (1.0/theta_hat + A[subjet_itr]);
					#if defined(_OPENMP)		
					S_private[ithread*n_ev +event_itr] += E_val*astar*weights_subjects[subjet_itr];
					#else
					ret[event_itr] += E_val*astar*weights_subjects[subjet_itr];
					#endif	
				}
			}
		}


	
	    #if defined(_OPENMP)
	    #pragma omp single
	    {	
	    	for(int i=0; i<n_ev; i++) {
			for(int t=0; t<nthreads; t++) {
			    ret[i] += S_private[n_ev*t + i];
			}
	    	}
	    	delete S_private;
	    }
	    #endif

	}

	for (event_itr = 0; event_itr < n_ev; event_itr++ ) {
		ret[event_itr] = weights_events[event_itr] / ret[event_itr];
	}
}



// [[Rcpp::export]]
NumericVector estimate_H_23(
			NumericVector gz,
			NumericVector theta_hat,
			NumericVector sorted_ev_times,
			NumericVector estimated_H,
			NumericVector times,
		        NumericVector V_times,
			NumericVector delta,
			NumericVector N_1_tau,
			NumericVector A_1_tau,			
			NumericVector weights_subjects,
			NumericVector weights_events
) {

  int n_ev = estimated_H.size() - 1;
  int n_times = times.size();

  double gz_c[n_times];

  double sorted_ev_times_c[n_ev];

  double estimated_H_c[n_ev + 1];
  double competing_estimated_H_c[n_ev];

  double times_c[n_times];
  double V_times_c[n_times];
  double delta_c[n_times];

  double theta = theta_hat[0];

  double ret[n_ev];



  double weights_subjects_c[n_times];
  double weights_events_c[n_ev];

  double N_1_tau_c[n_times];
  double A_1_tau_c[n_times];

  NumericVector out(n_ev);


  copy_num_vector(gz,gz_c);
  copy_num_vector(sorted_ev_times,sorted_ev_times_c);
  copy_num_vector(estimated_H,estimated_H_c);
  copy_num_vector(times,times_c);
  copy_num_vector(V_times,V_times_c);
  copy_num_vector(delta,delta_c);
  copy_num_vector(N_1_tau,N_1_tau_c);
  copy_num_vector(A_1_tau,A_1_tau_c);
  copy_num_vector(weights_subjects,weights_subjects_c);
  copy_num_vector(weights_events,weights_events_c);

  c_estimate_H_23(gz_c,
		  theta,
		  estimated_H_c,
		  sorted_ev_times_c,
		  n_ev,
		  times_c,
		  V_times_c,
		  delta_c,
		  N_1_tau_c,
		  A_1_tau_c,
		  n_times,
		  weights_subjects_c,
		  weights_events_c,
		  ret);

  for(int i = 0; i < n_ev; ++i) {
    out[i] = ret[i];
  }

  return out;
}



// [[Rcpp::export]]
NumericVector sum_e_sqr(NumericVector a, NumericVector b, NumericVector times, NumericVector sortedTimes) {
  int m = a.size();
  int n = b.size();

  double* a_c = R2C_vec(a);
  double* b_c = R2C_vec(b);
  double* t_c =  R2C_vec(times);
  double* st_c = R2C_vec(sortedTimes);
  double* out_c = (double* )malloc(sizeof(double)*m);

  NumericVector out(m);

  sum_e_matrix(a_c,m,b_c,n,t_c,st_c, out_c);



  free(a_c);
  free(b_c);
  free(t_c);
  free(st_c);

  for(int i = 0; i < m; ++i) {
    out[i] = out_c[i];
  }

  free(out_c);

  return out;
}


// [[Rcpp::export]]
NumericVector ln_sum_e_sqr(NumericVector a, NumericVector b, NumericVector times, NumericVector sortedTimes) {
  int m = a.size();
  int n = b.size();

  double* a_c = R2C_vec(a);
  double* b_c = R2C_vec(b);
  double* t_c =  R2C_vec(times);
  double* st_c = R2C_vec(sortedTimes);
  double* out_c = (double* )malloc(sizeof(double)*m);

  NumericVector out(m);

  ln_sum_e_matrix(a_c,m,b_c,n,t_c,st_c, out_c);



  free(a_c);
  free(b_c);
  free(t_c);
  free(st_c);

  for(int i = 0; i < m; ++i) {
    out[i] = out_c[i];
  }

  free(out_c);

  return out;
}



