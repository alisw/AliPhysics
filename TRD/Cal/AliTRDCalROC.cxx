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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Calibration base class for a single ROC                                  //
//  Contains one UShort_t value per pad                                      //
//  However, values are set and get as float, there are stored internally as //
//  (UShort_t) value * 10000                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TStyle.h>
#include <TMath.h>
#include <TH2F.h>
#include <TH1F.h>

#include "AliMathBase.h"

#include "AliTRDCalROC.h"

ClassImp(AliTRDCalROC)

//_____________________________________________________________________________
AliTRDCalROC::AliTRDCalROC()
  :TObject()
  ,fPla(0)
  ,fCha(0)
  ,fNrows(0)
  ,fNcols(0)
  ,fNchannels(0)
  ,fData(0)
{
  //
  // Default constructor
  //

}

//_____________________________________________________________________________
AliTRDCalROC::AliTRDCalROC(Int_t p, Int_t c)
  :TObject()
  ,fPla(p)
  ,fCha(c)
  ,fNrows(0)
  ,fNcols(144)
  ,fNchannels(0)
  ,fData(0)
{
  //
  // Constructor that initializes a given pad plane type
  //

  //
  // The pad plane parameter
  //
  switch (p) {
  case 0:
    if (c == 2) {
      // L0C0 type
      fNrows        =  12;
    }
    else {
      // L0C1 type
      fNrows        =  16;
    }
    break;
  case 1:
    if (c == 2) {
      // L1C0 type
      fNrows        =  12;
    }
    else {
      // L1C1 type
      fNrows        =  16;
    }
    break;
  case 2:
    if (c == 2) {
      // L2C0 type
      fNrows        =  12;
    }
    else {
      // L2C1 type
      fNrows        =  16;
    }
    break;
  case 3:
    if (c == 2) {
      // L3C0 type
      fNrows        =  12;
    }
    else {
      // L3C1 type
      fNrows        =  16;
    }
    break;
  case 4:
    if (c == 2) {
      // L4C0 type
      fNrows        =  12;
    }
    else {
      // L4C1 type
      fNrows        =  16;
    }
    break;
  case 5:
    if (c == 2) {
      // L5C0 type
      fNrows        =  12;
    }
    else {
      // L5C1 type
      fNrows        =  16;
    }
    break;
  };

  fNchannels = fNrows * fNcols;
  if (fNchannels != 0) {
    fData = new UShort_t[fNchannels];
  }

  for (Int_t i=0; i<fNchannels; ++i) {
    fData[i] = 0;
  }

}

//_____________________________________________________________________________
AliTRDCalROC::AliTRDCalROC(const AliTRDCalROC &c)
  :TObject(c)
  ,fPla(c.fPla)
  ,fCha(c.fCha)
  ,fNrows(c.fNrows)
  ,fNcols(c.fNcols)
  ,fNchannels(c.fNchannels)
  ,fData(0)
{
  //
  // AliTRDCalROC copy constructor
  //

  Int_t iBin = 0;

  fData = new UShort_t[fNchannels];
  for (iBin = 0; iBin < fNchannels; iBin++) {
    fData[iBin] = ((AliTRDCalROC &) c).fData[iBin];
  }

}

//_____________________________________________________________________________
AliTRDCalROC::~AliTRDCalROC()
{
  //
  // AliTRDCalROC destructor
  //

  if (fData) {
    delete [] fData;
    fData = 0;
  }

}

//_____________________________________________________________________________
AliTRDCalROC &AliTRDCalROC::operator=(const AliTRDCalROC &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliTRDCalROC &) c).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDCalROC::Copy(TObject &c) const
{
  //
  // Copy function
  //

  ((AliTRDCalROC &) c).fPla          = fPla;
  ((AliTRDCalROC &) c).fCha          = fCha;

  ((AliTRDCalROC &) c).fNrows        = fNrows;
  ((AliTRDCalROC &) c).fNcols        = fNcols;

  Int_t iBin = 0;

  ((AliTRDCalROC &) c).fNchannels = fNchannels;

  if (((AliTRDCalROC &) c).fData) delete [] ((AliTRDCalROC &) c).fData;
  ((AliTRDCalROC &) c).fData = new UShort_t[fNchannels];
  for (iBin = 0; iBin < fNchannels; iBin++) {
    ((AliTRDCalROC &) c).fData[iBin] = fData[iBin];
  }
  
  TObject::Copy(c);

}

//___________________________________________________________________________________
Double_t AliTRDCalROC::GetMean(AliTRDCalROC* const outlierROC) const
{
  //
  // Calculate the mean
  //

   Double_t *ddata = new Double_t[fNchannels];
   Int_t nPoints = 0;
   for (Int_t i=0;i<fNchannels;i++) {
      if ((!outlierROC) || (!(outlierROC->GetValue(i)))) {
	//if(fData[i] > 0.000000000000001){
         ddata[nPoints]= (Double_t) fData[i]/10000;
         nPoints++;
	 //}
      }
   }
   Double_t mean = TMath::Mean(nPoints,ddata);
   delete [] ddata;
   return mean;
}

//_______________________________________________________________________________________
Double_t AliTRDCalROC::GetMedian(AliTRDCalROC* const outlierROC) const
{
  //
  // Calculate the median
  //

  Double_t *ddata = new Double_t[fNchannels];
   Int_t nPoints = 0;
   for (Int_t i=0;i<fNchannels;i++) {
       if ((!outlierROC) || (!(outlierROC->GetValue(i)))) {
	 if(fData[i] > 0.000000000000001){         
	   ddata[nPoints]= (Double_t) fData[i]/10000;
	   nPoints++;
	 }
       }
   }
   Double_t mean = TMath::Median(nPoints,ddata);
   delete [] ddata;
   return mean;
}

//____________________________________________________________________________________________
Double_t AliTRDCalROC::GetRMS(AliTRDCalROC* const outlierROC) const
{
  //
  // Calculate the RMS
  //

  Double_t *ddata = new Double_t[fNchannels];
  Int_t nPoints = 0;
  for (Int_t i=0;i<fNchannels;i++) {
    if ((!outlierROC) || (!(outlierROC->GetValue(i)))) {
      //if(fData[i] > 0.000000000000001){
         ddata[nPoints]= (Double_t)fData[i]/10000;
         nPoints++;
	 //}
    }
  }
  Double_t mean = TMath::RMS(nPoints,ddata);
  delete [] ddata;
  return mean;
}

//______________________________________________________________________________________________
Double_t AliTRDCalROC::GetLTM(Double_t *sigma, Double_t fraction, AliTRDCalROC* const outlierROC)
{
  //
  //  Calculate LTM mean and sigma
  //

  Double_t *ddata = new Double_t[fNchannels];
  Double_t mean=0, lsigma=0;
  UInt_t nPoints = 0;
  for (Int_t i=0;i<fNchannels;i++) {
     if (!outlierROC || !(outlierROC->GetValue(i))) {
       if(fData[i] > 0.000000000000001){
	 ddata[nPoints]= (Double_t) fData[i]/10000;
	 nPoints++;
       }
     }
  }
  Int_t hh = TMath::Min(TMath::Nint(fraction *nPoints), Int_t(nPoints));
  AliMathBase::EvaluateUni(nPoints,ddata, mean, lsigma, hh);
  if (sigma) *sigma=lsigma;
  delete [] ddata;
  return mean;
}

//___________________________________________________________________________________
Bool_t AliTRDCalROC::Add(Float_t c1)
{
  //
  // add constant
  //

  Bool_t result = kTRUE;
  for (Int_t  idata = 0; idata< fNchannels; idata++) {
    if(((GetValue(idata)+c1) <= 6.5535) && ((GetValue(idata)+c1) >= 0.0)) SetValue(idata,GetValue(idata)+c1);
    else {
      SetValue(idata,0.0);
      result = kFALSE;
    }
  }
  return result;
}

//_______________________________________________________________________________________
Bool_t AliTRDCalROC::Multiply(Float_t c1)
{
  //
  // multiply constant
  //

  Bool_t result = kTRUE;
  if(c1 < 0) return kFALSE;
  for (Int_t  idata = 0; idata< fNchannels; idata++) {
    if((GetValue(idata)*c1) <= 6.5535) SetValue(idata,GetValue(idata)*c1);
    else {
      SetValue(idata,0.0);
      result = kFALSE;
    }
  }
  return result;
}

//____________________________________________________________________________________________
Bool_t AliTRDCalROC::Add(const AliTRDCalROC * roc, Double_t c1)
{
  //
  // add values 
  //

  Bool_t result = kTRUE;
  for (Int_t  idata = 0; idata< fNchannels; idata++){
    if(((GetValue(idata)+roc->GetValue(idata)*c1) <= 6.5535) && 
       ((GetValue(idata)+roc->GetValue(idata)*c1) >= 0.0)) 
          SetValue(idata,GetValue(idata)+roc->GetValue(idata)*c1);
    else {
      SetValue(idata,0.0);
      result = kFALSE;
    }
  }
  return result;
}

//____________________________________________________________________________________________
Bool_t AliTRDCalROC::Multiply(const AliTRDCalROC*  roc) 
{
  //
  // multiply values - per by pad
  //

  Bool_t result = kTRUE;
  for (Int_t  idata = 0; idata< fNchannels; idata++){
    if((GetValue(idata)*roc->GetValue(idata)) <= 6.5535) 
      SetValue(idata,GetValue(idata)*roc->GetValue(idata));
    else {
      SetValue(idata,0.0);
      result = kFALSE;
    }
  }
  return result;
}

//______________________________________________________________________________________________
Bool_t AliTRDCalROC::Divide(const AliTRDCalROC*  roc) 
{
  //
  // divide values 
  //

  Bool_t result = kTRUE;
  Float_t kEpsilon=0.00000000000000001;
  for (Int_t  idata = 0; idata< fNchannels; idata++){
    if (TMath::Abs(roc->GetValue(idata))>kEpsilon){
      if((GetValue(idata)/roc->GetValue(idata)) <= 6.5535) 
        SetValue(idata,GetValue(idata)/roc->GetValue(idata));
      else {
	SetValue(idata,0.0);
	result = kFALSE;
      }
    }
    else result = kFALSE;
  }
  return result;
}
//______________________________________________________________________________________________
Bool_t AliTRDCalROC::Unfold() 
{
  //
  // Compute the mean value per pad col
  // Divide with this value each pad col
  // This is for the noise study
  // Return kFALSE if one or more of the pad col was not normalised
  //

  Bool_t result = kTRUE;
  Float_t kEpsilon=0.00000000000000001;
  
  // calcul the mean value per col
  for(Int_t icol = 0; icol < fNcols; icol++){

    Float_t mean = 0.0;
    Float_t nb   = 0.0;

    for(Int_t irow = 0; irow < fNrows; irow++){
      if((GetValue(icol,irow) > 0.06) && (GetValue(icol,irow) < 0.15)){
	mean += GetValue(icol,irow);
	nb += 1.0;
      }
    }
    
    if(nb > kEpsilon) {
      
      mean = mean/nb;

      if(mean > kEpsilon){
	for(Int_t irow = 0; irow < fNrows; irow++){
	  Float_t value = GetValue(icol,irow);
	  SetValue(icol,irow,(Float_t)(value/mean));
	}    
      }
      else result = kFALSE;
      
    }

    else result = kFALSE;

    
  }
 

  return result;
}

//__________________________________________________________________________________
TH2F * AliTRDCalROC::MakeHisto2D(Float_t min, Float_t max,Int_t type,  Float_t mu)
{
  //
  // make 2D histo
  // type -1 = user defined range
  //       0 = nsigma cut nsigma=min
  //       1 = delta cut around median delta=min

  Float_t kEpsilonr = 0.005;
  gStyle->SetPalette(1);
  
  if (type>=0){
    if (type==0){
      // nsigma range
      Float_t mean  = GetMean();
      Float_t sigma = GetRMS();
      Float_t nsigma = TMath::Abs(min);
      sigma *= nsigma;
      if(sigma < kEpsilonr) sigma = kEpsilonr;
      min = mean-sigma;
      max = mean+sigma;
    }
    if (type==1){
      // fixed range
      Float_t mean   = GetMedian();
      Float_t  delta = min;
      min = mean-delta;
      max = mean+delta;
    }
    if (type==2){
      Double_t sigma;
      Float_t mean  = GetLTM(&sigma,max);
      sigma*=min;
      if(sigma < kEpsilonr) sigma = kEpsilonr;
      min = mean-sigma;
      max = mean+sigma;

    }
  }
  
  char  name[1000];
  sprintf(name,"%s 2D Plane %d Chamber %d",GetTitle(),fPla, fCha);
  TH2F * his = new TH2F(name,name,fNrows,0, fNrows, fNcols, 0, fNcols);
  for (Int_t irow=0; irow<fNrows; irow++){
    for (Int_t icol=0; icol<fNcols; icol++){
      his->Fill(irow+0.5,icol+0.5,GetValue(icol,irow)*mu);
    }
  }
  his->SetStats(0);
  his->SetMaximum(max);
  his->SetMinimum(min);
  return his;
}

//_______________________________________________________________________________________
TH1F * AliTRDCalROC::MakeHisto1D(Float_t min, Float_t max,Int_t type,  Float_t mu)
{
  //
  // make 1D histo
  // type -1 = user defined range
  //       0 = nsigma cut nsigma=min
  //       1 = delta cut around median delta=min

  Float_t kEpsilonr = 0.005;

  if (type>=0){
    if (type==0){
      // nsigma range
      Float_t mean  = GetMean();
      Float_t sigma = GetRMS();
      Float_t nsigma = TMath::Abs(min);
      sigma *= nsigma;
      if(sigma < kEpsilonr) sigma = kEpsilonr;
      min = mean-sigma;
      max = mean+sigma;
    }
    if (type==1){
      // fixed range
      Float_t mean   = GetMedian();
      Float_t  delta = min;
      min = mean-delta;
      max = mean+delta;
    }
    if (type==2){
      //
      // LTM mean +- nsigma
      //
      Double_t sigma;
      Float_t mean  = GetLTM(&sigma,max);
      sigma*=min;
      if(sigma < kEpsilonr) sigma = kEpsilonr;
      min = mean-sigma;
      max = mean+sigma;
    }
  }
  char  name[1000];
  sprintf(name,"%s 1D Plane %d Chamber %d",GetTitle(),fPla, fCha);
  TH1F * his = new TH1F(name,name,100, min,max);
  for (Int_t irow=0; irow<fNrows; irow++){
    for (Int_t icol=0; icol<fNcols; icol++){
      his->Fill(GetValue(icol,irow)*mu);
    }
  }
  return his;
}
