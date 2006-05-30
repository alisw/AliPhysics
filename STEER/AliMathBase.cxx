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


///////////////////////////////////////////////////////////////////////////
// Class AliMathBase
// 
// Subset of  matheamtical functions  not included in the TMath
//

///////////////////////////////////////////////////////////////////////////
#include "TMath.h"
#include "AliMathBase.h"
#include "Riostream.h"
 
ClassImp(AliMathBase) // Class implementation to enable ROOT I/O
 
AliMathBase::AliMathBase() : TObject()
{
// Default constructor
}
///////////////////////////////////////////////////////////////////////////
AliMathBase::~AliMathBase()
{
// Destructor
}


//_____________________________________________________________________________
void AliMathBase::EvaluateUni(Int_t nvectors, Double_t *data, Double_t &mean
                           , Double_t &sigma, Int_t hh)
{
  //
  // Robust estimator in 1D case MI version - (faster than ROOT version)
  //
  // For the univariate case
  // estimates of location and scatter are returned in mean and sigma parameters
  // the algorithm works on the same principle as in multivariate case -
  // it finds a subset of size hh with smallest sigma, and then returns mean and
  // sigma of this subset
  //

  if (hh==0)
    hh=(nvectors+2)/2;
  Double_t faclts[]={2.6477,2.5092,2.3826,2.2662,2.1587,2.0589,1.9660,1.879,1.7973,1.7203,1.6473};
  Int_t *index=new Int_t[nvectors];
  TMath::Sort(nvectors, data, index, kFALSE);
  
  Int_t    nquant = TMath::Min(Int_t(Double_t(((hh*1./nvectors)-0.5)*40))+1, 11);
  Double_t factor = faclts[nquant-1];
  
  Double_t sumx  =0;
  Double_t sumx2 =0;
  Int_t    bestindex = -1;
  Double_t bestmean  = 0; 
  Double_t bestsigma = data[index[nvectors-1]]-data[index[0]];   // maximal possible sigma
  for (Int_t i=0; i<hh; i++){
    sumx  += data[index[i]];
    sumx2 += data[index[i]]*data[index[i]];
  }
  
  Double_t norm = 1./Double_t(hh);
  Double_t norm2 = 1./Double_t(hh-1);
  for (Int_t i=hh; i<nvectors; i++){
    Double_t cmean  = sumx*norm;
    Double_t csigma = (sumx2 - hh*cmean*cmean)*norm2;
    if (csigma<bestsigma){
      bestmean  = cmean;
      bestsigma = csigma;
      bestindex = i-hh;
    }
    
    sumx  += data[index[i]]-data[index[i-hh]];
    sumx2 += data[index[i]]*data[index[i]]-data[index[i-hh]]*data[index[i-hh]];
  }
  
  Double_t bstd=factor*TMath::Sqrt(TMath::Abs(bestsigma));
  mean  = bestmean;
  sigma = bstd;
  delete [] index;

}



void AliMathBase::EvaluateUniExternal(Int_t nvectors, Double_t *data, Double_t &mean, Double_t &sigma, Int_t hh,  Float_t externalfactor)
{
  // Modified version of ROOT robust EvaluateUni
  // robust estimator in 1D case MI version
  // added external factor to include precision of external measurement
  // 

  if (hh==0)
    hh=(nvectors+2)/2;
  Double_t faclts[]={2.6477,2.5092,2.3826,2.2662,2.1587,2.0589,1.9660,1.879,1.7973,1.7203,1.6473};
  Int_t *index=new Int_t[nvectors];
  TMath::Sort(nvectors, data, index, kFALSE);
  //
  Int_t    nquant = TMath::Min(Int_t(Double_t(((hh*1./nvectors)-0.5)*40))+1, 11);
  Double_t factor = faclts[0];
  if (nquant>0){
    // fix proper normalization - Anja
    factor = faclts[nquant-1];
  }

  //
  //
  Double_t sumx  =0;
  Double_t sumx2 =0;
  Int_t    bestindex = -1;
  Double_t bestmean  = 0; 
  Double_t bestsigma = -1;
  for (Int_t i=0; i<hh; i++){
    sumx  += data[index[i]];
    sumx2 += data[index[i]]*data[index[i]];
  }
  //   
  Double_t kfactor = 2.*externalfactor - externalfactor*externalfactor;
  Double_t norm = 1./Double_t(hh);
  for (Int_t i=hh; i<nvectors; i++){
    Double_t cmean  = sumx*norm;
    Double_t csigma = (sumx2*norm - cmean*cmean*kfactor);
    if (csigma<bestsigma ||  bestsigma<0){
      bestmean  = cmean;
      bestsigma = csigma;
      bestindex = i-hh;
    }
    //
    //
    sumx  += data[index[i]]-data[index[i-hh]];
    sumx2 += data[index[i]]*data[index[i]]-data[index[i-hh]]*data[index[i-hh]];
  }
  
  Double_t bstd=factor*TMath::Sqrt(TMath::Abs(bestsigma));
  mean  = bestmean;
  sigma = bstd;
  delete [] index;
}


//_____________________________________________________________________________
Int_t AliMathBase::Freq(Int_t n, const Int_t *inlist
                        , Int_t *outlist, Bool_t down)
{    
  //
  //  Sort eleements according occurancy 
  //  The size of output array has is 2*n 
  //

  Int_t * sindexS = new Int_t[n];     // temp array for sorting
  Int_t * sindexF = new Int_t[2*n];   
  for (Int_t i=0;i<n;i++) sindexF[i]=0;
  //
  TMath::Sort(n,inlist, sindexS, down);  
  Int_t last      = inlist[sindexS[0]];
  Int_t val       = last;
  sindexF[0]      = 1;
  sindexF[0+n]    = last;
  Int_t countPos  = 0;
  //
  //  find frequency
  for(Int_t i=1;i<n; i++){
    val = inlist[sindexS[i]];
    if (last == val)   sindexF[countPos]++;
    else{      
      countPos++;
      sindexF[countPos+n] = val;
      sindexF[countPos]++;
      last =val;
    }
  }
  if (last==val) countPos++;
  // sort according frequency
  TMath::Sort(countPos, sindexF, sindexS, kTRUE);
  for (Int_t i=0;i<countPos;i++){
    outlist[2*i  ] = sindexF[sindexS[i]+n];
    outlist[2*i+1] = sindexF[sindexS[i]];
  }
  delete [] sindexS;
  delete [] sindexF;
  
  return countPos;

}
