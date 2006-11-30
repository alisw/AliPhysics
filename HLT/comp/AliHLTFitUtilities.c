/* $Id$ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <setjmp.h>
#include "AliHLTFitUtilities.h"

jmp_buf env;

DPOINT *plane;

void f2gauss5(double x,double a[],double *y,double dyda[],int na)
{
  /*Function describing a sum of 2D gaussians with 5 parameters.
    There number of gaussians is na/5. */
  
  int i,index;
  double fac,fac1,fac2,ex,ex1,ex2,arg1,arg2,u,v;
  
  /*printf("fitting na %d, with pad %f time %f amplitude %f padwidth %f timewidth %f\n",na,a[2],a[4],a[1],a[3],a[5]);*/
        index = nint(x);
	if( index < 0 || index >=FIT_PTS )
	{
		fprintf( stderr, "ff2gauss: wrong index %d\n", index );
		return;
	}
	u     = plane[index].u;
	v     = plane[index].v;
	/*printf("u %f v %f\n",u,v);*/
	*y=0.0;
	for (i=1;i<=na-1;i+=5)
	{
		arg1      = (u-a[i+1])/a[i+2];
		arg2      = (v-a[i+3])/a[i+4];
		ex1       = exp(-0.5*arg1*arg1);
		ex2       = exp(-0.5*arg2*arg2);
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
  /*printf("%s\n",error_text);
    exit(1);*/
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

int gaussj(double **a, int n, double **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol=0,irow=0,j,k,l,ll;
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
						return -1;
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
			return -1;
		}
		if (mabs(a[icol][icol]) < EPSILON ) 
		  {
		    nrerror("gaussj: a[icol][icol] == 0");
		    return -1;
		  }
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
	return 0;
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

int mrqmin(double x[], double y[], double sig[], int ndata, double a[], int ia[],
	int ma, double **covar, double **alpha, double *chisq,
	void (*funcs)(double, double [], double *, double [], int), double *alamda)
{
	void covsrt(double **covar, int ma, int ia[], int mfit);
	int gaussj(double **a, int n, double **b, int m);
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
	if(gaussj(covar,mfit,oneda,1)<0)
	  {
	    free_matrix(oneda,1,mfit,1,1);
	    free_vector(da,1,ma);
	    free_vector(beta,1,ma);
	    free_vector(atry,1,ma);
	    return -1;
	  }
	for (j=1;j<=mfit;j++) da[j]=oneda[j][1];
	if (*alamda == 0.0) {
	        covsrt(covar,ma,ia,mfit);
		free_matrix(oneda,1,mfit,1,1);
		free_vector(da,1,ma);
		free_vector(beta,1,ma);
		free_vector(atry,1,ma);
		return 0;
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
	return 0;
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
		if(mrqmin(x,y,sig,NPT,a,ia,MA,covar,alpha,&chisq,funcs,&alamda)<0)
		  {
		    free_matrix(alpha,1,MA,1,MA);
		    free_matrix(covar,1,MA,1,MA);
		    return -1;
		  }
		k=1;
		itst=0;
		for (;;) {
			k++;
			ochisq=chisq;
			if(mrqmin(x,y,sig,NPT,a,ia,MA,covar,alpha,&chisq,funcs,&alamda)<0)
			  {
			    free_matrix(alpha,1,MA,1,MA);
			    free_matrix(covar,1,MA,1,MA);
			    return -1;
			  }
			if (chisq > ochisq)
				itst=0;
			else if (fabs(ochisq-chisq) < 0.1)
				itst++;
			if (itst < 4) continue;
			alamda=0.0;
			if(mrqmin(x,y,sig,NPT,a,ia,MA,covar,alpha,&chisq,funcs,&alamda)<0)
			  {
			    free_matrix(alpha,1,MA,1,MA);
			    free_matrix(covar,1,MA,1,MA);
			    return -1;
			  }
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
	  /*if( control_g.print_fit_errors==2 )*/
	  fprintf( stderr, " runtime error\n" );

	        free_matrix(alpha,1,MA,1,MA);
		free_matrix(covar,1,MA,1,MA);
		return -1;
	}
}

