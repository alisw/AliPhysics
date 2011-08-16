/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

#include <TMath.h>
#include "AliEventPlaneResolution.h"

/* $Id: AliITSOnlineSDDInjectors.cxx 48265 2011-03-09 23:36:06Z prino $ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Implementation of the class to compute event plane resolution //
// for flow analyses                                             //
// Origin: F.Prino, Torino, prino@to.infn.it                     //
//                                                               //
///////////////////////////////////////////////////////////////////

ClassImp(AliEventPlaneResolution)

//______________________________________________________________________
AliEventPlaneResolution::AliEventPlaneResolution():TObject(),
  fK(1),
  fSubRes(1.)
{
  // Default contructor
}


//______________________________________________________________________
AliEventPlaneResolution::AliEventPlaneResolution(Int_t k):
  TObject(),
  fK(k),
  fSubRes(1.)
{
  // Standard constructor
}


//______________________________________________________________________
Double_t AliEventPlaneResolution::Pol(Double_t x, Int_t k){
  // compute chi from polynomial approximation
  Double_t c[5];
  if(k==1){ 
    c[0]=0.626657; c[1]=0.; c[2]=-0.09694; c[3]=0.02754; c[4]=-0.002283;
  }
  else if(k==2){
    c[0]=0.; c[1]=0.25; c[2]=-0.011414; c[3]=-0.034726; c[4]=0.006815;
  }
  else return -1;
  return c[0]*x+c[1]*x*x+c[2]*x*x*x+c[3]*x*x*x*x+c[4]*x*x*x*x*x;
}

//______________________________________________________________________
Double_t AliEventPlaneResolution:: ResolK1(Double_t x){
  return TMath::Sqrt(TMath::Pi()/8)*x*TMath::Exp(-x*x/4)*(TMath::BesselI0(x*x/4)+TMath::BesselI1(x*x/4));
}


//______________________________________________________________________
Double_t AliEventPlaneResolution::FindChi(Double_t res,  Int_t k){
  // compute chi variable (=v2*sqrt(N)) from external values

  Double_t x1=0;
  Double_t x2=15;
  Double_t y1,y2;
  if(k==1){
    y1=ResolK1(x1)-res;
    y2=ResolK1(x2)-res;
  }
  else if(k==2){
    y1=Pol(x1,2)-res;
    y2=Pol(x2,2)-res;
  }
  else return -1;

  if(y1*y2>0) return -1;
  if(y1==0) return y1;
  if(y2==0) return y2;
  Double_t xmed,ymed;
  Int_t jiter=0;
  while((x2-x1)>0.0001){
    xmed=0.5*(x1+x2);
    if(k==1){
      y1=ResolK1(x1)-res;
      y2=ResolK1(x2)-res;
      ymed=ResolK1(xmed)-res;
    }
    else if(k==2){
      y1=Pol(x1,2)-res;
      y2=Pol(x2,2)-res;
      ymed=Pol(xmed,2)-res;
    }
    else return -1;
    if((y1*ymed)<0) x2=xmed;
    if((y2*ymed)<0)x1=xmed;
    if(ymed==0) return xmed;
    jiter++;
  }
  return 0.5*(x1+x2);
}

//______________________________________________________________________
Double_t AliEventPlaneResolution::GetFullEvResol(Double_t resSub, Int_t k){
  // computes event plane resolution starting from sub event resolution
  Double_t chisub=FindChi(resSub,k);
  Double_t chifull=chisub*TMath::Sqrt(2);
  if(k==1) return ResolK1(chifull);
  else if(k==2) return Pol(chifull,2);
  else{
    printf("k should be <=2\n");
    return 1.;
  }
}

//______________________________________________________________________
Double_t AliEventPlaneResolution::GetFullEvResol(const TH1F* hSubEvCorr, Int_t k){
  // computes event plane resolution starting from sub event correlation histogram
  if(!hSubEvCorr) return 1.;
  Double_t resSub=GetSubEvResol(hSubEvCorr);
  return GetFullEvResol(resSub,k);
}
