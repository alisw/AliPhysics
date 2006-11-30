/*@(#) $Id$*/

#ifndef AliHLTFitUtilities
#define AliHLTFitUtilities

/*This we do because this file is read both with c and c++ compiler, 
  and extern "C" is needed only in case of c++. */
#ifdef __cplusplus
extern "C" 
#endif
void f2gauss5( double, double *, double *,double *,int );

#ifdef __cplusplus 
extern "C" 
#endif
int lev_marq_fit( double x[], double y[], double sig[], int NPT, double a[], int ia[], double dev[], int MA,
		  double *chisq_p, void (*funcs)(double, double [], double *, double [], int) );

#define max(x, y) (((x) > (y)) ? (x) : (y))
#define min(x, y) (((x) < (y)) ? (x) : (y))
#define mabs(x)   (((x) > 0) ? (x) : (-(x)))
#define nint(x)   ((int)((x) < 0 ? (x)-0.5 : (x)+0.5))
#define DBL(x)    ((double)(x))
#define veclen2(x,y)  (DBL(x)*DBL(x) + DBL(y)*DBL(y))
#define samesign(x,y)   ((((x)>=0 && (y)>=0) || ((x)<0&&(y)<0)) ? TRUE : FALSE )

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
#define SQR(x)       ((x)*(x))

#define NR_END 1
#define FREE_ARG char*
#define EPSILON             1.0E-12
#define FIT_PTS     2000
#define  FIT_MAXPAR   41
#define NUM_PARS 5

/*--- fitting 2-dimensional cluster --------------------------*/
struct DPOINT {
  double u;
  double v;
};
typedef struct DPOINT DPOINT;

extern  DPOINT *plane; 

typedef struct { 
					long   float_size;
					long   steps;
					double dmin;
					double dmax;
					double step_size;
				} exp_header_t;
typedef float FLOAT_SIZE;

#endif
