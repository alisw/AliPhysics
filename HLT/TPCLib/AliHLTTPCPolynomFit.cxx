// **************************************************************************
// This file is property of and copyright by the ALICE HLT Project          *
// ALICE Experiment at CERN, All rights reserved.                           *
//                                                                          *
// Primary Authors: Vito Nordloh <vito.nordloh@vitonordloh.de>              *
//                  Sergey Gorbunov <sergey.gorbunov@fias.uni-frankfurt.de> *
//                  for The ALICE HLT Project.                              *
//                                                                          *
// Permission to use, copy, modify and distribute this software and its     *
// documentation strictly for non-commercial purposes is hereby granted     *
// without fee, provided that the above copyright notice appears in all     *
// copies and that both the copyright notice and this permission notice     *
// appear in the supporting documentation. The authors make no claims       *
// about the suitability of this software for any purpose. It is            *
// provided "as is" without express or implied warranty.                    *
//                                                                          *
//***************************************************************************

#include "AliHLTTPCPolynomFit.h"
#include "TMath.h"

// #include "TMatrixDSym.h"
// #include "TMatrixD.h"
// #include "TVectorD.h"


void AliHLTTPCPolynomFit::Reset( int nCoefficients )
{
  if( nCoefficients>=0 ){
    fN = nCoefficients;
    delete[] fA;
    delete[] fB;
    fA = new long double [fN*(fN+1)/2];
    fB = new long double [fN];  
  }
  
  for( int i=0; i<fN; i++){
    for( int j=0; j<=i; j++) fA[i*(i+1)/2+j] = 0;
    fB[i] = 0;
  }
}

void AliHLTTPCPolynomFit::AddMeasurement( float f[], float m)
{
  for( int i=0; i<fN; i++){
    for( int j=0; j<=i; j++){      
      fA[i*(i+1)/2+j]+= ((long double)f[i])*((long double)f[j]);
    }
    fB[i]+=m*f[i];
  }
}


int AliHLTTPCPolynomFit::Fit( float Coefficients[] )
{
  int ret = invS(fA,fN);  
  for( int i=0; i<fN; i++){
    long double s = 0;
    for( int j=0; j<=i; j++)
      s+= fA[i*(i+1)/2+j]*fB[j];
    for( int j=i+1; j<fN; j++)
      s+= fA[j*(j+1)/2+i]*fB[j];
    
    Coefficients[i] = s;
  }
  return ret;  
}


int AliHLTTPCPolynomFit::invS( long double A[], int N )
{
  int ret = 0;
  
  const long double ZERO = 1.E-20;
  
  // input: simmetric > 0 NxN matrix A = {a11,a21,a22,a31..a33,..}  
  // output: inverse A, in case of problems fill zero and return 1
  
  // A->low triangular Anew : A = Anew x Anew^T
  // method:
  // for(j=1,N) for(i=j,N) Aij=(Aii-sum_{k=1}^{j-1}Aik*Ajk )/Ajj
  //   
  
  {
    long double *j1 = A, *jj = A;
    for( int j=1; j<=N; j1+=j++, jj+=j ){
      long double *ik = j1, x = 0;
      while( ik!=jj ){
	x -= (*ik) * (*ik);
	ik++;
      }
      x += *ik;
      if( x > ZERO ){
	x = sqrt(x);
	*ik = x;
	ik++;
	x = 1 / x;
	for( int step=1; step<=N-j; ik+=++step ){ // ik==Ai1
	  long double sum = 0;
	  for( long double *jk=j1; jk!=jj; sum += (*(jk++)) * (*(ik++)) ) ;
	  *ik = (*ik - sum) * x; // ik == Aij
	}
      }else{
	long double *ji=jj;
	for( int i=j; i<N; i++ ) *(ji+=i) = 0.;
	ret = -1;
      }   
    }
  }
    
  // A -> Ainv
  // method : 
  // for(i=1,N){ 
  //   Aii = 1/Aii; 
  //   for(j=1,i-1) Aij=-(sum_{k=j}^{i-1} Aik * Akj) / Aii ;
  // }
  
  {
    long double *ii=A,*ij=A;
    for( int i = 1; i<=N; ij=ii+1, ii+=++i ){
      if( *ii > ZERO ){
	long double x = -(*ii = 1./ *ii);
	{ 
	  long double *jj = A;
	  for( int j=1; j<i; jj+=++j, ij++ ){
	    long double *ik = ij, *kj = jj, sum = 0.;
	    for( int k=j; ik!=ii; kj+=k++, ik++ ){
	      sum += *ik * *kj;
	    }
	    *kj = sum * x;
	  }
	}
      }else{      
	for( long double *ik = ij; ik!=ii+1; ik++ ){
	  *ik = 0.;
	}
	ret = -1;
      }
    }
  }
  
  // A -> A^T x A
  // method: 
  // Aij = sum_{k=i}^N Aki * Akj
  
  {
    long double *ii=A, *ij=A;
    for( int i=1; i<=N; ii+=++i ){
      do{ 
	long double *ki = ii, *kj = ij, sum = 0.;
	for( int k=i; k<=N; ki+=k, kj+=k++ ) sum += (*ki) * (*kj);
	*ij = sum;
      }while( (ij++)!=ii );
    }    
  }
  return ret;    
}
 
