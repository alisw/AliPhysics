
#include "AliL3FitUtilities.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <setjmp.h>


jmp_buf env;

DPOINT *plane;

#ifndef PB
int build_exp_table()
/******************************************************************************/
{
	exp_header_t tab_h;
	int i;
	FLOAT_SIZE *exp_val;
	FILE  *fp;
	
	if( !(fp = fopen( "exp_table", "wb" )))
	  printf("build_exp_table: I/O error\n");
		
	tab_h.float_size = sizeof( FLOAT_SIZE );
	tab_h.steps     = 100000;
	tab_h.dmin      = -5;
	tab_h.dmax      = 1;
	tab_h.step_size = (tab_h.dmax - tab_h.dmin)/tab_h.steps;
	if( !(exp_val = (FLOAT_SIZE *)malloc( tab_h.steps * sizeof( FLOAT_SIZE ) ) ) )
	  printf("build_exp_table: malloc error\n");
	fwrite((const void *)&tab_h, (size_t)sizeof(exp_header_t), (size_t)1, fp );
	for( i=0; i<tab_h.steps; ++i ) {
		exp_val[i] = exp(tab_h.dmin + i*tab_h.step_size);
	}
	fwrite((const void *)exp_val, (size_t)sizeof(FLOAT_SIZE), (size_t)tab_h.steps, fp );
	
	free( (void *)exp_val );
	return 0;
}

#ifdef TEST_TABLE
test_table()
/**********************************************************/
{
	double i, delta, dmin, dmax, exp_tab();
	FILE   *fp;
	
	fp = fopen( "exp_test", "w" );
	
	delta = 6.0/10000.;
	dmin  = -5.1;
	dmax  =  1.1;
	for( i=dmin; i<dmax; i+=delta ){
		fprintf( fp, "%lf\t%lf\t%lf\n", i, exp(i), exp_tab(i) );
	}
	return 0;
}
#endif
double exp_tab(double in)
{
	return exp( in );
}
#else
double exp_tab(in)
/**********************************************************/
double in;
{
	static exp_header_t tab_h;
	static FLOAT_SIZE   *exp_val=NULL, slope;
	FILE                *fp=NULL;
	
#ifdef HP
	if( isnan(in) ) {
		if( (num_nan /100)*100 == num_nan )
			fprintf( stderr, "exp_tab: NaN %d\n", num_nan );
		++num_nan;
		return 1;
	}
#endif
#ifdef MAC
		return exp( in );
#else
	if( !exp_val ) {
		if( !(fp = fopen( "exp_table", "rb" ))) {
			build_exp_table();
		}
		if( !(fp = fopen( "exp_table", "rb" )))
			printf("exp_tab: I/O error\n");
		fread(&tab_h, (size_t)sizeof(exp_header_t), (size_t)1, fp );
		if( tab_h.float_size != sizeof( FLOAT_SIZE ) )
			build_exp_table();

		if( !(exp_val = (FLOAT_SIZE *)malloc( tab_h.steps * sizeof( FLOAT_SIZE ) ) ) )
			printf("exp_tab: malloc error\n");
		fread(exp_val, (size_t)sizeof(FLOAT_SIZE), (size_t)tab_h.steps, fp );
		slope = tab_h.steps / (tab_h.dmax - tab_h.dmin);
	}
	if( in < tab_h.dmin || in > tab_h.dmax )
		return exp( in );
	else
		return exp_val[(int)((in-tab_h.dmin) * slope)];
#endif
}
#endif

void f2gauss5(double x,double a[],double *y,double dyda[],int na)
{
  int i,index;
  double fac,fac1,fac2,ex,ex1,ex2,arg1,arg2,u,v;
  
	
	index = nint(x);
	if( index < 0 || index >=FIT_PTS )
	{
		fprintf( stderr, "ff2gauss: wrong index %ld\n", index );
		return;
	}
	u     = plane[index].u;
	v     = plane[index].v;
	
	*y=0.0;
	for (i=1;i<=na-1;i+=5)
	{
		arg1      = (u-a[i+1])/a[i+2];
		arg2      = (v-a[i+3])/a[i+4];
		ex1       = exp_tab(-arg1*arg1);
		ex2       = exp_tab(-arg2*arg2);
		ex        = ex1*ex2;
		fac       = a[i]*ex*2.0;
		fac1      = fac * arg1;
		fac2      = fac * arg2;
		*y       += a[i] * ex;
		dyda[i]   = ex;
		dyda[i+1] = fac1/a[i+2];
		dyda[i+2] = fac1*arg1/a[i+2];
		dyda[i+3] = fac2/a[i+4];
		dyda[i+4] = fac2*arg2/a[i+4];
	}
}

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
  printf("%s\n",error_text);
  exit(1);
}

void free_vector(double *v, long nl, long nh)
/* free a double vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

double *vector(long nl,long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((unsigned int) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

double **matrix(long nrl,long nrh,long ncl,long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((unsigned int)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((unsigned int)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void gaussj(double **a, int n, double **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv,swap;

	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1) {
						free_ivector(ipiv,1,n);
						free_ivector(indxr,1,n);
						free_ivector(indxc,1,n);
						nrerror("gaussj: Singular Matrix-1");
					}
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0){
			free_ivector(ipiv,1,n);
			free_ivector(indxr,1,n);
			free_ivector(indxc,1,n);
			nrerror("gaussj: Singular Matrix-2");
		}
		if (mabs(a[icol][icol]) < EPSILON ) 
			nrerror("gaussj: a[icol][icol] == 0");
		pivinv=1.0/a[icol][icol];

		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}

void covsrt(double **covar, int ma, int ia[], int mfit)
{
	int i,j,k;
	double swap;

	for (i=mfit+1;i<=ma;i++)
		for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
	k=mfit;
	for (j=ma;j>=1;j--) {
		if (ia[j]) {
			for (i=1;i<=ma;i++) SWAP(covar[i][k],covar[i][j])
			for (i=1;i<=ma;i++) SWAP(covar[k][i],covar[j][i])
			k--;
		}
	}
}

void mrqcof(double x[], double y[], double sig[], int ndata, double a[], int ia[],
	int ma, double **alpha, double beta[], double *chisq,
	void (*funcs)(double, double [], double *, double [], int))
{
	int i,j,k,l,m,mfit=0;
	double ymod,wt,sig2i,dy,*dyda;

	dyda=vector(1,ma);

	for (j=1;j<=ma;j++)
		if (ia[j]) mfit++;
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=j;k++) alpha[j][k]=0.0;
		beta[j]=0.0;
	}
	*chisq=0.0;
	for (i=1;i<=ndata;i++) {
		(*funcs)(x[i],a,&ymod,dyda,ma);
		sig2i=1.0/(sig[i]*sig[i]);
		dy=y[i]-ymod;
		for (j=0,l=1;l<=ma;l++) {
			if (ia[l]) {
				wt=dyda[l]*sig2i;
				for (j++,k=0,m=1;m<=l;m++)
					if (ia[m]) alpha[j][++k] += wt*dyda[m];
				beta[j] += dy*wt;
			}
		}
		*chisq += dy*dy*sig2i;
	}
	for (j=2;j<=mfit;j++)
		for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];
	free_vector(dyda,1,ma);
}

void mrqmin(double x[], double y[], double sig[], int ndata, double a[], int ia[],
	int ma, double **covar, double **alpha, double *chisq,
	void (*funcs)(double, double [], double *, double [], int), double *alamda)
{
	void covsrt(double **covar, int ma, int ia[], int mfit);
	void gaussj(double **a, int n, double **b, int m);
	void mrqcof(double x[], double y[], double sig[], int ndata, double a[],
		int ia[], int ma, double **alpha, double beta[], double *chisq,
		void (*funcs)(double, double [], double *, double [], int));
	int j,k,l,m;
	static int mfit;
	static double ochisq,*atry,*beta,*da,**oneda;

	if (*alamda < 0.0) {
		atry=vector(1,ma);
		beta=vector(1,ma);
		da=vector(1,ma);
		for (mfit=0,j=1;j<=ma;j++)
			if (ia[j]) mfit++;
		oneda=matrix(1,mfit,1,1);
		*alamda=0.001;
		mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,chisq,funcs);
		ochisq=(*chisq);
		for (j=1;j<=ma;j++) atry[j]=a[j];
	}
	for (j=0,l=1;l<=ma;l++) {
		if (ia[l]) {
			for (j++,k=0,m=1;m<=ma;m++) {
				if (ia[m]) {
					k++;
					covar[j][k]=alpha[j][k];
				}
			}
			covar[j][j]=alpha[j][j]*(1.0+(*alamda));
			oneda[j][1]=beta[j];
		}
	}
	gaussj(covar,mfit,oneda,1);
	for (j=1;j<=mfit;j++) da[j]=oneda[j][1];
	if (*alamda == 0.0) {
		covsrt(covar,ma,ia,mfit);
		free_matrix(oneda,1,mfit,1,1);
		free_vector(da,1,ma);
		free_vector(beta,1,ma);
		free_vector(atry,1,ma);
		return;
	}
	for (j=0,l=1;l<=ma;l++)
		if (ia[l]) atry[l]=a[l]+da[++j];
	mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,chisq,funcs);
	if (*chisq < ochisq) {
		*alamda *= 0.1;
		ochisq=(*chisq);
		for (j=0,l=1;l<=ma;l++) {
			if (ia[l]) {
				for (j++,k=0,m=1;m<=ma;m++) {
					if (ia[m]) {
						k++;
						alpha[j][k]=covar[j][k];
					}
				}
				beta[j]=da[j];
				a[l]=atry[l];
			}
		}
	} else {
		*alamda *= 10.0;
		*chisq=ochisq;
	}
}


int lev_marq_fit( double x[], double y[], double sig[], int NPT, double a[], int ia[], double dev[], int MA,
                  double *chisq_p, void (*funcs)(double, double [], double *, double [], int) )
{
	int     i,itst,k;
	double   alamda,chisq,ochisq,**covar,**alpha;
	
	covar=matrix(1,MA,1,MA);
	alpha=matrix(1,MA,1,MA);

	if( setjmp(env) == 0 ) {	
		alamda = -1;
		mrqmin(x,y,sig,NPT,a,ia,MA,covar,alpha,&chisq,funcs,&alamda);
		k=1;
		itst=0;
		for (;;) {
			k++;
			ochisq=chisq;
			mrqmin(x,y,sig,NPT,a,ia,MA,covar,alpha,&chisq,funcs,&alamda);
			if (chisq > ochisq)
				itst=0;
			else if (fabs(ochisq-chisq) < 0.1)
				itst++;
			if (itst < 4) continue;
			alamda=0.0;
			mrqmin(x,y,sig,NPT,a,ia,MA,covar,alpha,&chisq,funcs,&alamda);
			*chisq_p = chisq;
			for (i=1;i<=MA;i++) 
				dev[i] = sqrt(covar[i][i]);
			break;
		}
		free_matrix(alpha,1,MA,1,MA);
		free_matrix(covar,1,MA,1,MA);
		return 0;
	}
	else {
	  //if( control_g.print_fit_errors==2 )
	  fprintf( stderr, " runtime error\n" );

		free_matrix(alpha,1,MA,1,MA);
		free_matrix(covar,1,MA,1,MA);
		return -1;
	}
}

