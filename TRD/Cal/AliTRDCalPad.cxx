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
//  TRD calibration class for parameters which saved per pad                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TStyle.h>

//#include "AliMathBase.h"

#include "AliTRDCalPad.h"
#include "AliTRDCalROC.h"
#include "AliTRDCalDet.h"
#include "AliTRDpadPlane.h"
#include "AliTRDgeometry.h"

ClassImp(AliTRDCalPad)

//_____________________________________________________________________________
AliTRDCalPad::AliTRDCalPad():TNamed()
{
  //
  // AliTRDCalPad default constructor
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    fROC[idet] = 0;
  }

}

//_____________________________________________________________________________
AliTRDCalPad::AliTRDCalPad(const Text_t *name, const Text_t *title)
                :TNamed(name,title)
{
  //
  // AliTRDCalPad constructor
  //

  for (Int_t isec = 0; isec < kNsect; isec++) {
    for (Int_t ipla = 0; ipla < kNplan; ipla++) {
      for (Int_t icha = 0; icha < kNcham; icha++) {
        Int_t idet = GetDet(ipla,icha,isec);
        fROC[idet] = new AliTRDCalROC(ipla,icha);
      }
    }
  }

}

//_____________________________________________________________________________
AliTRDCalPad::AliTRDCalPad(const AliTRDCalPad &c)
  :TNamed(c)
{
  //
  // AliTRDCalPad copy constructor
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    fROC[idet] = new AliTRDCalROC(*((AliTRDCalPad &) c).fROC[idet]);
  }

}

//_____________________________________________________________________________
AliTRDCalPad::~AliTRDCalPad()
{
  //
  // AliTRDCalPad destructor
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    if (fROC[idet]) {
      delete fROC[idet];
      fROC[idet] = 0;
    }
  }

}

//_____________________________________________________________________________
AliTRDCalPad &AliTRDCalPad::operator=(const AliTRDCalPad &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliTRDCalPad &) c).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDCalPad::Copy(TObject &c) const
{
  //
  // Copy function
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    if (((AliTRDCalPad &) c).fROC[idet]) {
      delete ((AliTRDCalPad &) c).fROC[idet];
    }
    ((AliTRDCalPad &) c).fROC[idet] = new AliTRDCalROC();
    if (fROC[idet]) {
      fROC[idet]->Copy(*((AliTRDCalPad &) c).fROC[idet]);
    }
  }

  TObject::Copy(c);

}

//_____________________________________________________________________________
Bool_t AliTRDCalPad::ScaleROCs(const AliTRDCalDet* values)
{
  // 
  // Scales ROCs of this class with the values from the class <values>
  // Is used if an AliTRDCalPad object defines local variations of a parameter
  // defined per detector using a AliTRDCalDet class
  //
  
  if (!values)
    return kFALSE;
  Bool_t result = kTRUE;
  for (Int_t idet = 0; idet < kNdet; idet++) {
    if (fROC[idet]) { 
      if(!fROC[idet]->Multiply(values->GetValue(idet))) result = kFALSE;
    }
  }
  return result;
}

//_____________________________________________________________________________
void AliTRDCalPad::SetCalROC(Int_t det, AliTRDCalROC* calroc)
{
  // 
  // Set the AliTRDCalROC to this one
  //
  
  if (!calroc) return;
  if (fROC[det]) { 
    for(Int_t icol = 0; icol < calroc->GetNcols(); icol++){
      for(Int_t irow = 0; irow < calroc->GetNrows(); irow++){
	fROC[det]->SetValue(icol,irow,calroc->GetValue(icol,irow));
      }
    }
  }

}
//_____________________________________________________________________________
Double_t AliTRDCalPad::GetMeanRMS(Double_t &rms, const AliTRDCalDet *calDet, Int_t type)
{
    //
    // Calculate mean an RMS of all rocs
    // If calDet correct the CalROC from the detector coefficient
    // type == 0 for gain and vdrift
    // type == 1 for t0 
    //
  Double_t factor = 0.0;
  if(type == 0) factor = 1.0;
  Double_t sum = 0, sum2 = 0, n=0, val=0;
  for (Int_t idet = 0; idet < kNdet; idet++) {
    if(calDet) factor = calDet->GetValue(idet);
    AliTRDCalROC *calRoc = fROC[idet];
    if ( calRoc ){
      for (Int_t irow=0; irow<calRoc->GetNrows(); irow++){
	for (Int_t icol=0; icol<calRoc->GetNcols(); icol++){
	  if(type == 0) val = calRoc->GetValue(icol,irow)*factor;
	  else  val = calRoc->GetValue(icol,irow)+factor;
	  sum+=val;
	  sum2+=val*val;
	  n++;
	}
      }
    }
  }
  Double_t n1 = 1./n;
  Double_t mean = sum*n1;
  rms  = TMath::Sqrt(TMath::Abs(sum2*n1-mean*mean));
  return mean;
}

//_____________________________________________________________________________
Double_t AliTRDCalPad::GetMean(const AliTRDCalDet *calDet, Int_t type
                             , AliTRDCalPad* const outlierPad)
{
    //
    // return mean of the mean of all ROCs
    // If calDet correct the CalROC from the detector coefficient
    // type == 0 for gain and vdrift
    // type == 1 for t0 
    //
  Double_t factor = 0.0;
  if(type == 0) factor = 1.0;
  Double_t arr[kNdet];
  Int_t n=0;
  for (Int_t idet = 0; idet < kNdet; idet++) {
    if(calDet) factor = calDet->GetValue(idet);
    AliTRDCalROC *calRoc = fROC[idet];
    if ( calRoc ){
      AliTRDCalROC* outlierROC = 0;
      if (outlierPad) outlierROC = outlierPad->GetCalROC(idet);
      if(type == 0)  arr[n] = calRoc->GetMean(outlierROC)*factor;
      else arr[n] = calRoc->GetMean(outlierROC)+factor;
      n++;
    }
  }
  return TMath::Mean(n,arr);
}

//_____________________________________________________________________________
Double_t AliTRDCalPad::GetRMS(const AliTRDCalDet *calDet, Int_t type
                            , AliTRDCalPad* const outlierPad)
{
    //
    // return mean of the RMS of all ROCs
    // If calDet correct the CalROC from the detector coefficient
    // type == 0 for gain and vdrift
    // type == 1 for t0 
    //
  Double_t factor = 0.0;
  if(type == 0) factor = 1.0;
  Double_t arr[kNdet];
  Int_t n=0;
  for (Int_t idet = 0; idet < kNdet; idet++) {
    if(calDet) factor = calDet->GetValue(idet);
    AliTRDCalROC *calRoc = fROC[idet];
    if ( calRoc ){
      AliTRDCalROC* outlierROC = 0;
      if (outlierPad) outlierROC = outlierPad->GetCalROC(idet);
      if(type == 0)  arr[n] = calRoc->GetRMS(outlierROC)*factor;
      else  arr[n] = calRoc->GetRMS(outlierROC);
      n++;
    }
  }
    return TMath::Mean(n,arr);
}

//_____________________________________________________________________________
Double_t AliTRDCalPad::GetMedian(const AliTRDCalDet *calDet, Int_t type
                               , AliTRDCalPad* const outlierPad)
{
    //
    // return mean of the median of all ROCs
    // If AliTRDCalDet, the correct the CalROC from the detector coefficient
    // type == 0 for gain and vdrift
    // type == 1 for t0 
    //
  Double_t factor = 0.0;
  if(type == 0) factor = 1.0;
  Double_t arr[kNdet];
  Int_t n=0;
  for (Int_t idet = 0; idet < kNdet; idet++) {
    if(calDet) factor = calDet->GetValue(idet);
    AliTRDCalROC *calRoc = fROC[idet];
    if ( calRoc ){
      AliTRDCalROC* outlierROC = 0;
      if (outlierPad) outlierROC = outlierPad->GetCalROC(idet);
      if(type == 0)  arr[n] = calRoc->GetMedian(outlierROC)*factor;
      else  arr[n] = calRoc->GetMedian(outlierROC)+factor;
      n++;
    }
  }
  return TMath::Mean(n,arr);
}

//_____________________________________________________________________________
Double_t AliTRDCalPad::GetLTM(Double_t *sigma, Double_t fraction
                            , const AliTRDCalDet *calDet, Int_t type
                            , AliTRDCalPad* const outlierPad)
{
    //
    // return mean of the LTM and sigma of all ROCs
    // If calDet correct the CalROC from the detector coefficient
    // type == 0 for gain and vdrift
    // type == 1 for t0 
    //
  Double_t factor = 0.0;
  if(type == 0) factor = 1.0;
  Double_t arrm[kNdet];
  Double_t arrs[kNdet];
  Double_t *sTemp=0x0;
  Int_t n=0;
  
  for (Int_t idet = 0; idet < kNdet; idet++) {
    if(calDet) factor = calDet->GetValue(idet);
    AliTRDCalROC *calRoc = fROC[idet];
    if ( calRoc ){
      if ( sigma ) sTemp=arrs+n;
      AliTRDCalROC* outlierROC = 0;
      if (outlierPad) outlierROC = outlierPad->GetCalROC(idet);
      if(type == 0)  arrm[n] = calRoc->GetLTM(sTemp,fraction, outlierROC)*factor;
      else arrm[n] = calRoc->GetLTM(sTemp,fraction, outlierROC)+factor;
      n++;
    }
  }
  if ( sigma ) *sigma = TMath::Mean(n,arrs);
  return TMath::Mean(n,arrm);
}

//_____________________________________________________________________________
TH1F * AliTRDCalPad::MakeHisto1D(const AliTRDCalDet *calDet, Int_t typedet
                               , Float_t min, Float_t max,Int_t type)
{
  //
  // make 1D histo
  // type -1 = user defined range
  //       0 = nsigma cut nsigma=min
  // If calDet correct the CalROC from the detector coefficient
  // typedet == 0 for gain and vdrift
  // typedet == 1 for t0
  //
  Float_t kEpsilonr = 0.005;

  Double_t factor = 0.0;
  if(typedet == 0) factor = 1.0;
 
  if (type>=0){
    if (type==0){
      // nsigma range 
      Float_t mean  = GetMean(calDet,typedet);
      Float_t sigma = 0.0;
      if(GetRMS(calDet,typedet) > kEpsilonr) sigma = GetRMS(calDet,typedet);
      else {
	Double_t rms = 0.0;
	sigma = GetMeanRMS(rms,calDet,typedet);
      }
      Float_t nsigma = TMath::Abs(min);
      sigma *= nsigma;
      if(sigma < kEpsilonr) sigma = kEpsilonr;
      min = mean-sigma;
      max = mean+sigma;
    }
    if (type==1){
      // fixed range
      Float_t mean   = GetMedian(calDet,typedet);
      Float_t  delta = min;
      min = mean-delta;
      max = mean+delta;
    }
    if (type==2){
      //
      // LTM mean +- nsigma
      //
      Double_t sigma;
      Float_t mean  = GetLTM(&sigma,max,calDet,typedet);
      sigma*=min;
      if(sigma < kEpsilonr) sigma = kEpsilonr;
      min = mean-sigma;
      max = mean+sigma;
    }
  }
  char  name[1000];
  snprintf(name,1000,"%s Pad 1D",GetTitle());
  TH1F * his = new TH1F(name,name,100, min,max);
    for (Int_t idet = 0; idet < kNdet; idet++) {
      if(calDet) factor = calDet->GetValue(idet);
      if (fROC[idet]){
	for (Int_t irow=0; irow<fROC[idet]->GetNrows(); irow++){
	  for (Int_t icol=0; icol<fROC[idet]->GetNcols(); icol++){
	    if(typedet == 0) his->Fill(fROC[idet]->GetValue(irow,icol)*factor);
	    else his->Fill(fROC[idet]->GetValue(irow,icol)+factor);
	  }
	}
      }
    }
    return his;
}

//_____________________________________________________________________________
TH2F *AliTRDCalPad::MakeHisto2DSmPl(Int_t sm, Int_t pl, const AliTRDCalDet *calDet
                                  , Int_t typedet, Float_t min, Float_t max,Int_t type)
{
  //
  // Make 2D graph
  // sm    - supermodule number
  // pl    - plane number
  // If calDet correct the CalROC from the detector coefficient
  // typedet == 0 for gain and vdrift
  // typedet == 1 for t0
  //
  gStyle->SetPalette(1);
  Double_t factor = 0.0;
  if(typedet == 0) factor = 1.0;

  Float_t kEpsilon = 0.000000000001;
  Float_t kEpsilonr = 0.005;

  AliTRDgeometry *trdGeo = new AliTRDgeometry();

  if (type>=0){
    if (type==0){
      // nsigma range
      Float_t mean  = GetMean(calDet,typedet);
      Float_t sigma = 0.0;
      if(GetRMS(calDet,typedet) > kEpsilonr) sigma = GetRMS(calDet,typedet);
      else {
	Double_t rms = 0.0;
	sigma = GetMeanRMS(rms,calDet,typedet);
      }
      Float_t nsigma = TMath::Abs(min);
      sigma *= nsigma;
      if(sigma < kEpsilonr) sigma = kEpsilonr;
      min = mean-sigma;
      max = mean+sigma;
    }
    if (type==1){
      // fixed range
      Float_t mean   = GetMedian(calDet,typedet);
      Float_t  delta = min;
      min = mean-delta;
      max = mean+delta;
    }
    if (type==2){
      //
      // LTM mean +- nsigma
      //
      Double_t sigma;
      Float_t mean  = GetLTM(&sigma,max,calDet,typedet);
      sigma*=min;
      if(sigma < kEpsilonr) sigma = kEpsilonr;
      min = mean-sigma;
      max = mean+sigma;
    }
  }
  
  AliTRDpadPlane *padPlane0 = trdGeo->GetPadPlane(pl,0);
  Double_t row0    = padPlane0->GetRow0();
  Double_t col0    = padPlane0->GetCol0();

  char  name[1000];
  snprintf(name,1000,"%s Pad 2D sm %d pl %d",GetTitle(),sm,pl);
  TH2F * his = new TH2F( name, name, 76,-TMath::Abs(row0),TMath::Abs(row0),144,-TMath::Abs(col0),TMath::Abs(col0));

  // Where we begin
  Int_t offsetsmpl = 30*sm+pl;
  

  for (Int_t k = 0; k < kNcham; k++){
    Int_t det = offsetsmpl+k*6;
    if(calDet) factor = calDet->GetValue(det);
    if (fROC[det]){
      AliTRDCalROC * calRoc = fROC[det];
      for (Int_t irow=0; irow<calRoc->GetNrows(); irow++){
	for (Int_t icol=0; icol<calRoc->GetNcols(); icol++){
	  if (TMath::Abs(calRoc->GetValue(icol,irow))>kEpsilon){
	    Int_t binz     = 0;
	    Int_t kb       = kNcham-1-k;
	    Int_t krow     = calRoc->GetNrows()-1-irow;
	    Int_t kcol     = calRoc->GetNcols()-1-icol;
	    if(kb > 2) binz = 16*(kb-1)+12+krow+1;
	    else binz = 16*kb+krow+1; 
	    Int_t biny = kcol+1;
	    Float_t value = calRoc->GetValue(icol,irow);
	    if(typedet == 0) his->SetBinContent(binz,biny,value*factor);
	    else his->SetBinContent(binz,biny,value+factor);
	  }
	}
      }
    }
  }
  his->SetXTitle("z (cm)");
  his->SetYTitle("y (cm)");
  his->SetStats(0);
  his->SetMaximum(max);
  his->SetMinimum(min);
  delete trdGeo;
  return his;
}

//_____________________________________________________________________________
TH2F *AliTRDCalPad::MakeHisto2DCh(Int_t ch, const AliTRDCalDet *calDet, Int_t typedet
                                , Float_t min, Float_t max,Int_t type)
{
  //
  // Make 2D graph mean value in z direction
  // ch    - chamber
  // If calDet correct the CalROC from the detector coefficient
  // typedet == 0 for gain and vdrift
  // typedet == 1 for t0
  //
  gStyle->SetPalette(1);
  Double_t factor = 0.0;
  if(typedet == 0) factor = 1.0;

  Float_t kEpsilonr = 0.005;
  
  if (type>=0){
    if (type==0){
      // nsigma range
      Float_t mean  = GetMean(calDet,typedet);
      Float_t sigma = 0.0;
      if(GetRMS(calDet,typedet) > kEpsilonr) sigma = GetRMS(calDet,typedet);
      else {
	Double_t rms = 0.0;
	sigma = GetMeanRMS(rms,calDet,typedet);
      }
      Float_t nsigma = TMath::Abs(min);
      sigma *= nsigma;
      if(sigma < kEpsilonr) sigma = kEpsilonr;
      min = mean-sigma;
      max = mean+sigma;
    }
    if (type==1){
      // fixed range
      Float_t mean   = GetMedian(calDet,typedet);
      Float_t  delta = min;
      min = mean-delta;
      max = mean+delta;
    }
    if (type==2){
      //
      // LTM mean +- nsigma
      //
      Double_t sigma;
      Float_t mean  = GetLTM(&sigma,max,calDet,typedet);
      sigma*=min;
      if(sigma < kEpsilonr) sigma = kEpsilonr;
      min = mean-sigma;
      max = mean+sigma;
    }
  }

  AliTRDgeometry *trdGeo = new AliTRDgeometry();
      
  Float_t kEpsilon = 0.000000000001;

  Double_t poslocal[3]  = {0.0,0.0,0.0};
  Double_t posglobal[3] = {0.0,0.0,0.0};
  
  char  name[1000];
  snprintf(name,1000,"%s Pad 2D ch %d",GetTitle(),ch);
  TH2F * his = new TH2F( name, name, 400,-400.0,400.0,400,-400.0,400.0);

  // Where we begin
  Int_t offsetch = 6*ch;
  

  for (Int_t isec = 0; isec < kNsect; isec++){
    for(Int_t ipl = 0; ipl < kNplan; ipl++){
      Int_t det = offsetch+isec*30+ipl;
      if(calDet) factor = calDet->GetValue(det);
      if (fROC[det]){
	AliTRDCalROC * calRoc = fROC[det];
	AliTRDpadPlane *padPlane = trdGeo->GetPadPlane(ipl,ch);
	for (Int_t icol=0; icol<calRoc->GetNcols(); icol++){
	  poslocal[0] = trdGeo->GetTime0(ipl);
	  poslocal[2] = padPlane->GetRowPos(0);
	  poslocal[1] = padPlane->GetColPos(icol);
	  trdGeo->RotateBack(det,poslocal,posglobal);
	  Int_t binx = 1+TMath::Nint((posglobal[0]+400.0)*0.5);
	  Int_t biny = 1+TMath::Nint((posglobal[1]+400.0)*0.5);
	  Float_t value = 0.0;
	  Int_t nb = 0;
	  for (Int_t irow=0; irow<calRoc->GetNrows(); irow++){
	    if (TMath::Abs(calRoc->GetValue(icol,irow))>kEpsilon){
	      value += calRoc->GetValue(icol,irow);
	      nb++;	    
	    }
	  }
          if (nb > 0) {
	    value = value/nb;
	  }
	  if(typedet == 0) his->SetBinContent(binx,biny,value*factor);
	  else his->SetBinContent(binx,biny,value+factor);
	}
      }
    }    
  }
  his->SetXTitle("x (cm)");
  his->SetYTitle("y (cm)");
  his->SetStats(0);
  his->SetMaximum(max);
  his->SetMinimum(min);
  delete trdGeo;
  return his;
}

//_____________________________________________________________________________
Bool_t AliTRDCalPad::Add(Float_t c1)
{
    //
    // add constant for all channels of all ROCs
    //

  Bool_t result = kTRUE;
  for (Int_t idet = 0; idet < kNdet; idet++) {
    if (fROC[idet]){
      if(!fROC[idet]->Add(c1)) result = kFALSE;
    }
  }
  return result;
}

//_____________________________________________________________________________
Bool_t AliTRDCalPad::Multiply(Float_t c1)
{
    //
    // multiply constant for all channels of all ROCs
    //
  Bool_t result = kTRUE;
  for (Int_t idet = 0; idet < kNdet; idet++) {
    if (fROC[idet]){
      if(!fROC[idet]->Multiply(c1)) result = kFALSE;
    }
  }
  return result;
}

//_____________________________________________________________________________
Bool_t AliTRDCalPad::Add(const AliTRDCalPad * pad, Double_t c1
                       , const AliTRDCalDet * calDet1, const AliTRDCalDet *calDet2
                       , Int_t type)
{
    //
    // add calpad channel by channel multiplied by c1 - all ROCs
    // If calDet1 and calDet2, the correct the CalROC from the detector coefficient
    // then you have calDet1 and the calPad together
    // calDet2 and pad together
    // typedet == 0 for gain and vdrift
    // typedet == 1 for t0
    //
  Float_t kEpsilon = 0.000000000001;
 
  Double_t factor1 = 0.0;
  Double_t factor2 = 0.0;
  if(type == 0) {
    factor1 = 1.0;
    factor2 = 1.0;
  }
  Bool_t result = kTRUE;
  for (Int_t idet = 0; idet < kNdet; idet++) {
    if(calDet1) factor1 = calDet1->GetValue(idet);
    if(calDet2) factor2 = calDet2->GetValue(idet);
    if (fROC[idet]){
      if(type == 0){
	if(TMath::Abs(factor1) > kEpsilon){
	  if(!fROC[idet]->Add(pad->GetCalROC(idet),c1*factor2/factor1)) result = kFALSE;
	}
	else result = kFALSE;
      }
      else{
	AliTRDCalROC *croc = new AliTRDCalROC((const AliTRDCalROC) *pad->GetCalROC(idet));
	if(!croc->Add(factor2)) result = kFALSE;
	if(!fROC[idet]->Add(croc,c1)) result = kFALSE;
      }
    }	 
  }
  return result;
}

//_____________________________________________________________________________
Bool_t AliTRDCalPad::Multiply(const AliTRDCalPad * pad, const AliTRDCalDet * calDet1
                            , const AliTRDCalDet *calDet2, Int_t type)
{
    //
    // multiply calpad channel by channel - all ROCs
    // If calDet1 and calDet2, the correct the CalROC from the detector coefficient
    // then you have calDet1 and the calPad together
    // calDet2 and pad together
    // typedet == 0 for gain and vdrift
    // typedet == 1 for t0
    //
  Float_t kEpsilon = 0.000000000001;
  Bool_t result = kTRUE;
  Double_t factor1 = 0.0;
  Double_t factor2 = 0.0;
  if(type == 0){
    factor1 = 1.0;
    factor2 = 1.0;
  }
  for (Int_t idet = 0; idet < kNdet; idet++) {
    if(calDet1) factor1 = calDet1->GetValue(idet);
    if(calDet2) factor2 = calDet2->GetValue(idet);
    if (fROC[idet]){
      if(type == 0){
	if(TMath::Abs(factor1) > kEpsilon){
	  AliTRDCalROC *croc = new AliTRDCalROC((const AliTRDCalROC) *pad->GetCalROC(idet));
	  if(!croc->Multiply(factor2)) result = kFALSE;
	  if(!fROC[idet]->Multiply(croc)) result = kFALSE;
	}
	else result = kFALSE;
      }
      else{
	AliTRDCalROC *croc2 = new AliTRDCalROC((const AliTRDCalROC) *pad->GetCalROC(idet));
	if(!croc2->Add(factor2)) result = kFALSE;
	if(!fROC[idet]->Add(factor1)) result = kFALSE;
	if(!fROC[idet]->Multiply(croc2)) result = kFALSE;
	if(!fROC[idet]->Add(-factor1)) result = kFALSE;
      }
    }
  }
  return result;
}

//_____________________________________________________________________________
Bool_t AliTRDCalPad::Divide(const AliTRDCalPad * pad, const AliTRDCalDet * calDet1
                          , const AliTRDCalDet *calDet2, Int_t type)
{
    //
    // divide calpad channel by channel - all ROCs
    // If calDet1 and calDet2, the correct the CalROC from the detector coefficient
    // then you have calDet1 and the calPad together
    // calDet2 and pad together
    // typedet == 0 for gain and vdrift
    // typedet == 1 for t0
    //
  Float_t kEpsilon = 0.000000000001;
  Bool_t result = kTRUE;
  Double_t factor1 = 0.0;
  Double_t factor2 = 0.0;
  if(type == 0){
    factor1 = 1.0;
    factor2 = 1.0;
  }
  for (Int_t idet = 0; idet < kNdet; idet++) {
    if(calDet1) factor1 = calDet1->GetValue(idet);
    if(calDet2) factor2 = calDet2->GetValue(idet);
    if (fROC[idet]){
      if(type == 0){
	if(TMath::Abs(factor1) > kEpsilon){
	  AliTRDCalROC *croc = new AliTRDCalROC((const AliTRDCalROC) *pad->GetCalROC(idet));
	  if(!croc->Multiply(factor2)) result = kFALSE;
	  if(!fROC[idet]->Divide(croc)) result = kFALSE;
	}
	else result = kFALSE;
      }
      else{
	AliTRDCalROC *croc2 = new AliTRDCalROC((const AliTRDCalROC) *pad->GetCalROC(idet));
	if(!croc2->Add(factor2)) result = kFALSE;
	if(!fROC[idet]->Add(factor1)) result = kFALSE;
	if(!fROC[idet]->Divide(croc2)) result = kFALSE;
	if(!fROC[idet]->Add(-factor1)) result = kFALSE;
      }
    }
  }
  return result;
}
