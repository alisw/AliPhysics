/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

#include <TF1.h>
#include <TMath.h>

#include "AliSignalProcesor.h"


ClassImp(AliSignalProcesor)


Double_t asymgauss(Double_t* x, Double_t* par)
{
  // par[0] = normalization
  // par[1] = mean
  // par[2] = sigma
  // norm0  = 1
  // par[3] = lambda0
  // par[4] = norm1
  // par[5] = lambda1
  //
  
  Double_t par1save = par[1];    
  Double_t par2save = par[2];
  Double_t par3save = par[3];
  Double_t par5save = par[5];
  Double_t dx   = x[0]-par1save;
  //
  //
  Double_t  sigma2  = par2save*par2save;
  Double_t  sqrt2   = TMath::Sqrt(2.);
  if (-par3save*(dx-0.5*par3save*sigma2)>100) return 0;   // avoid overflow
  if (-par5save*(dx-0.5*par5save*sigma2)>100) return 0;   // avoid overflow 
  if (TMath::Abs(par[4])>1) return 0;
  Double_t  exp1    = par3save*TMath::Exp(-par3save*(dx-0.5*par3save*sigma2))
    *(1-TMath::Erf((par3save*sigma2-dx)/(sqrt2*par2save)));

  Double_t  exp2    = par5save*TMath::Exp(-par5save*(dx-0.5*par5save*sigma2))
    *(1-TMath::Erf((par5save*sigma2-dx)/(sqrt2*par2save)));


  return par[0]*(exp1+par[4]*exp2);
}

Double_t asymgaussN(Double_t* x, Double_t* par)
{
  // par[0] = normalization
  // par[1] = mean
  // par[2] = sigma
  // norm0  = 1
  // par[3] = lambda0
  // par[4] = norm1
  // par[5] = lambda1
  //
  
  Double_t par1save = par[1];    
  Double_t par2save = par[2];
  Double_t par3save = par[3];
  Double_t par5save = par[5];
  Double_t dx   = x[0]-par1save;
  //
  //
  Double_t  sigma2  = par2save*par2save;
  Double_t  sqrt2   = TMath::Sqrt(2.);
  if (-par3save*(dx-0.5*par3save*sigma2)>100) return 0;   // avoid overflow
  if (-par5save*(dx-0.5*par5save*sigma2)>100) return 0;   // avoid overflow 
  if (TMath::Abs(par[4])>=1) return 0;
  
  Double_t  exp1    = par3save*TMath::Exp(-par3save*(dx-0.5*par3save*sigma2))
    *0.5*(1-TMath::Erf((par3save*sigma2-dx)/(sqrt2*par2save)));

  Double_t  exp2    = par5save*TMath::Exp(-par5save*(dx-0.5*par5save*sigma2))
    *0.5*(1-TMath::Erf((par5save*sigma2-dx)/(sqrt2*par2save)));


  return par[0]*(1.*exp1+par[4]*exp2)/(1.+par[4]);
}


TF1 * AliSignalProcesor::GetAsymGauss()
{
  TF1 * f1 = new TF1("asymg",asymgaussN,-10,40,6);
  return f1;
}



void AliSignalProcesor::SplineSmoother(Double_t *ampin, Double_t *ampout, Int_t n)
{  
  //
  //
  Float_t in[10000];
  Float_t out[10000];
  in[0] = ampin[0];
  in[1] = (ampin[0]+ampin[1])*0.5;
  in[2*(n-1)]    = ampin[n-1];
  in[2*(n-1)+1]  = ampin[n-1];
  //
  // add charge to the end
  for (Int_t i=0;i<10;i++){
    in[2*(n-1)+i]=ampin[n-1];
  }

  //
  for (Int_t i=1;i<n-1;i++){
    in[2*i]    = ampin[i];
    in[2*i+1]  = (9.*(ampin[i]+ampin[i+1])-(ampin[i-1]+ampin[i+2]))/16.;    
  }
  //
  out[0] = in[0];
  for (Int_t i=1;i<=2*n;i++){
    out[i]  = (9.*(in[i]+in[i+1])-(in[i-1]+in[i+2]))/16.;    
  }
  //
  //
  for (int i=0;i<n;i++){
    ampout[i]      = out[2*i+1]; 
  }
}




void AliSignalProcesor::TailCancelationALTRO(Double_t *ampin, Double_t *ampout, Float_t K, Float_t L, 
		      Int_t n)
{
  //
  // ALTRO
  Float_t temp;
  ampout[0]  = ampin[0];
  temp = ampin[0];
  for (int i=1;i<n;i++){
    ampout[i]   = ampin[i]   + (K-L)*temp;
    temp        = ampin[i]   +  K*temp;
  }
}

//
//
void AliSignalProcesor::TailCancelationTRD(Double_t *ampin, Double_t *ampout, Float_t r, Float_t c, 
		      Int_t n)
{
  //TRD
  //
  Double_t reminder=0;
  //
  for (Int_t i=0; i<n; i++){
    ampout[i] = ampin[i]-reminder;
    //
    reminder  = r*(reminder+c*ampout[i]);
  }
  
}

void AliSignalProcesor::TailMaker(Double_t *ampin, Double_t *ampout, Float_t lambda, 
		      Int_t n)
{
  
  Double_t l = TMath::Exp(-lambda);
  //
  Float_t temp=0;
  for (Int_t i=n-1; i>0; i--){
    ampout[i] = ampin[i]+temp;
    //
    temp  = l*(temp+ampin[i]);
  }
}

void AliSignalProcesor::TailCancelationALTRO1(Double_t *ampin, Double_t *ampout, Float_t norm, 
					    Float_t lambda, Int_t n)
{
  
  Double_t l = TMath::Exp(-lambda);
  Double_t k = l*(1.-norm*lambda);

  return TailCancelationALTRO(ampin,ampout,k,l,n);
}


void AliSignalProcesor::TailCancelationTRD1(Double_t *ampin, Double_t *ampout, Float_t norm, 
					  Float_t lambda, Int_t n)
{
  //
  //
  Double_t r = TMath::Exp(-lambda);
  Double_t c = norm*lambda;
  return TailCancelationTRD(ampin,ampout,r,c,n);
}




void AliSignalProcesor::TailCancelationMI(Double_t *ampin, Double_t *ampout, Float_t norm, 
					Float_t lambda, Int_t n)
{
  
  Double_t L = TMath::Exp(-lambda*0.5);
  Double_t K = L*(1.-norm*lambda*0.5);
  //
  //
  Float_t in[10000];
  Float_t out[10000];
  for (Int_t i=0;i<n*2+20;i++){
    in[i] = 0;
    out[i]= 0;
  }
  in[0] = ampin[0];
  in[1] = (ampin[0]+ampin[1])*0.5;
  in[2*(n-1)]    = ampin[n-1];
  in[2*(n-1)+1]  = ampin[n-1];
  //
  for (Int_t i=1;i<n-2;i++){
    in[2*i]    = ampin[i];
    in[2*i+1]  = (9.*(ampin[i]+ampin[i+1])-(ampin[i-1]+ampin[i+2]))/16;    
  }
  //
  Float_t temp;
  out[0]     = in[0];
  temp       = in[0];
  for (int i=1;i<=2*n;i++){
    out[i]      = in[i]   + (K-L)*temp;
    temp        = in[i]   +  K*temp;
  }
  //
  //
  for (int i=0;i<n;i++){
    ampout[i]      = out[2*i+1];
  }
}





void AliSignalProcesor::TailMakerSpline(Double_t *ampin, Double_t *ampout, Float_t lambda, 
		      Int_t n)
{
  
  Double_t l = TMath::Exp(-lambda*0.5);
  //
  //
  Float_t in[10000];
  Float_t out[10000];
  for (Int_t i=0;i<n*2+20;i++){
    in[i] = 0;
    out[i]= 0;
  }
  in[0] = ampin[0];
  in[1] = (ampin[0]+ampin[1])*0.5;
  in[2*(n-1)]    = ampin[n-1];
  in[2*(n-1)+1]  = ampin[n-1];
  //
  // add charge to the end
  for (Int_t i=0;i<10;i++){
    in[2*(n-1)+i]=ampin[n-1];
  }

  //
  for (Int_t i=1;i<n-2;i++){
    in[2*i]    = ampin[i];
    in[2*i+1]  = (9.*(ampin[i]+ampin[i+1])-(ampin[i-1]+ampin[i+2]))/16;    
  }
  //
  //
  Float_t temp;
  out[2*n]     = in[2*n];
  temp         = 0;
  for (int i=2*n;i>=0;i--){
    out[i]      = in[i]   + temp;
    temp        = l*(temp+in[i]);
  }
  //
  //
  for (int i=0;i<n;i++){
    ampout[i]      = out[2*i+1]; 
  }
}
