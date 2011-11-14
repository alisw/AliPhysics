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
//  TRD calibration class for parameters which saved per detector            //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TStyle.h>

#include "AliTRDCalDet.h"
#include "AliTRDgeometry.h"
#include "AliMathBase.h"
#include "AliTRDpadPlane.h"

ClassImp(AliTRDCalDet)

//_____________________________________________________________________________
AliTRDCalDet::AliTRDCalDet():TNamed()
{
  //
  // AliTRDCalDet default constructor
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    fData[idet] = 0;
  }

}

//_____________________________________________________________________________
AliTRDCalDet::AliTRDCalDet(const Text_t *name, const Text_t *title)
                :TNamed(name,title)
{
  //
  // AliTRDCalDet constructor
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    fData[idet] = 0;
  }

}

//_____________________________________________________________________________
AliTRDCalDet::AliTRDCalDet(const AliTRDCalDet &c):TNamed(c)
{
  //
  // AliTRDCalDet copy constructor
  //

  ((AliTRDCalDet &) c).Copy(*this);

}

///_____________________________________________________________________________
AliTRDCalDet::~AliTRDCalDet()
{
  //
  // AliTRDCalDet destructor
  //

}

//_____________________________________________________________________________
AliTRDCalDet &AliTRDCalDet::operator=(const AliTRDCalDet &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliTRDCalDet &) c).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDCalDet::Copy(TObject &c) const
{
  //
  // Copy function
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    ((AliTRDCalDet &) c).fData[idet] = fData[idet];
  }

  TObject::Copy(c);

}

//___________________________________________________________________________________
Double_t AliTRDCalDet::GetMean(AliTRDCalDet * const outlierDet) const
{
  //
  // Calculate the mean
  //

  if (!outlierDet) return TMath::Mean(kNdet,fData);
  Double_t *ddata = new Double_t[kNdet];
  Int_t nPoints = 0;
  for (Int_t i=0;i<kNdet;i++) {
    if (!(outlierDet->GetValue(i))) {
      ddata[nPoints]= fData[nPoints];
      nPoints++;
    }
  }
  Double_t mean = TMath::Mean(nPoints,ddata);
  delete [] ddata;
  return mean;
}

//_______________________________________________________________________________________
Double_t AliTRDCalDet::GetMedian(AliTRDCalDet * const outlierDet) const
{
  //
  // Calculate the median
  //

  if (!outlierDet) return (Double_t) TMath::Median(kNdet,fData);
  Double_t *ddata = new Double_t[kNdet];
  Int_t nPoints = 0;
  for (Int_t i=0;i<kNdet;i++) {
    if (!(outlierDet->GetValue(i))) {
      ddata[nPoints]= fData[nPoints];
      nPoints++;
    }
  }
  Double_t mean = TMath::Median(nPoints,ddata);
  delete [] ddata;
  return mean;

}

//____________________________________________________________________________________________
Double_t AliTRDCalDet::GetRMS(AliTRDCalDet * const outlierDet) const
{
  //
  // Calculate the RMS
  //

  if (!outlierDet) return TMath::RMS(kNdet,fData);
  Double_t *ddata = new Double_t[kNdet];
  Int_t nPoints = 0;
  for (Int_t i=0;i<kNdet;i++) {
    if (!(outlierDet->GetValue(i))) {
      ddata[nPoints]= fData[nPoints];
      nPoints++;
    }
  }
  Double_t mean = TMath::RMS(nPoints,ddata);
  delete [] ddata;
  return mean;
}

//____________________________________________________________________________________________
Double_t AliTRDCalDet::GetRMSRobust(Double_t robust) const
{
  //
  // Calculate the robust RMS
  //

  // sorted
  Int_t *index = new Int_t[kNdet];
  TMath::Sort((Int_t)kNdet,fData,index);
 
  // reject
  Double_t reject = (Int_t) (kNdet*(1-robust)/2.0);
  if(reject <= 0.0) reject = 0.0;
  if(reject >= kNdet) reject = 0.0;
  //printf("Rejecter %f\n",reject);

  Double_t *ddata = new Double_t[kNdet];
  Int_t nPoints = 0;
  for (Int_t i=0;i<kNdet;i++) {
    Bool_t rej = kFALSE;
    for(Int_t k = 0; k < reject; k++){
      if(i==index[k]) rej = kTRUE;
      if(i==index[kNdet-(k+1)]) rej  = kTRUE;
    }
    if(!rej){
      ddata[nPoints]= fData[i];
      nPoints++;
    }
  }
  //printf("Number of points %d\n",nPoints);
  Double_t mean = TMath::RMS(nPoints,ddata);
  delete [] ddata;
  delete [] index;
  return mean;
}

//______________________________________________________________________________________________
Double_t AliTRDCalDet::GetLTM(Double_t *sigma
                            , Double_t fraction
                            , AliTRDCalDet * const outlierDet)
{
  //
  //  Calculate LTM mean and sigma
  //

  Double_t *ddata = new Double_t[kNdet];
  Double_t mean=0, lsigma=0;
  UInt_t nPoints = 0;
  for (Int_t i=0;i<kNdet;i++) {
     if (!outlierDet || !(outlierDet->GetValue(i))) {
        ddata[nPoints]= fData[nPoints];
        nPoints++;
     }
  }
  Int_t hh = TMath::Min(TMath::Nint(fraction *nPoints), Int_t(nPoints));
  AliMathBase::EvaluateUni(nPoints,ddata, mean, lsigma, hh);
  if (sigma) *sigma=lsigma;
  delete [] ddata;
  return mean;
}

//_________________________________________________________________________
TH1F * AliTRDCalDet::MakeHisto1Distribution(Float_t min, Float_t max,Int_t type)
{
  //
  // make 1D histo
  // type -1 = user defined range
  //       0 = nsigma cut nsigma=min
  //       1 = delta cut around median delta=min
  //

  if (type>=0){
    if (type==0){
      // nsigma range
      Float_t mean  = GetMean();
      Float_t sigma = GetRMS();
      Float_t nsigma = TMath::Abs(min);
      min = mean-nsigma*sigma;
      max = mean+nsigma*sigma;
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
      min = mean-sigma;
      max = mean+sigma;
    }
  }
  char  name[1000];
  snprintf(name,1000,"%s CalDet 1Distribution",GetTitle());
  TH1F * his = new TH1F(name,name,100, min,max);
  for (Int_t idet=0; idet<kNdet; idet++){
    his->Fill(GetValue(idet));
  }
  return his;
}

//________________________________________________________________________________
TH1F * AliTRDCalDet::MakeHisto1DAsFunctionOfDet(Float_t min, Float_t max,Int_t type)
{
  //
  // make 1D histo
  // type -1 = user defined range
  //       0 = nsigma cut nsigma=min
  //       1 = delta cut around median delta=min
  //

  if (type>=0){
    if (type==0){
      // nsigma range
      Float_t mean  = GetMean();
      Float_t sigma = GetRMS();
      Float_t nsigma = TMath::Abs(min);
      min = mean-nsigma*sigma;
      max = mean+nsigma*sigma;
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
      min = mean-sigma;
      max = mean+sigma;

    }
  }
 
  char  name[1000];
  snprintf(name,1000,"%s CalDet as function of det",GetTitle());
  TH1F * his = new TH1F(name,name,kNdet, 0, kNdet);
  for(Int_t det = 0; det< kNdet; det++){
    his->Fill(det+0.5,GetValue(det));
  }

  his->SetMaximum(max);
  his->SetMinimum(min);
  return his;
}

//_____________________________________________________________________________
TH2F *AliTRDCalDet::MakeHisto2DCh(Int_t ch, Float_t min, Float_t max, Int_t type)
{
  //
  // Make 2D graph
  // ch    - chamber
  // type -1 = user defined range
  //       0 = nsigma cut nsigma=min
  //       1 = delta cut around median delta=min
  //

  gStyle->SetPalette(1);
  if (type>=0){
    if (type==0){
      // nsigma range
      Float_t mean  = GetMean();
      Float_t sigma = GetRMS();
      Float_t nsigma = TMath::Abs(min);
      min = mean-nsigma*sigma;
      max = mean+nsigma*sigma;
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
      min = mean-sigma;
      max = mean+sigma;

    }
  }
    
  AliTRDgeometry *trdGeo = new AliTRDgeometry();

  Double_t poslocal[3]  = {0.0,0.0,0.0};
  Double_t posglobal[3] = {0.0,0.0,0.0};
  
  char  name[1000];
  snprintf(name,1000,"%s CalDet 2D ch %d",GetTitle(),ch);
  TH2F * his = new TH2F(name, name, 400,-400.0,400.0,400,-400.0,400.0);

  // Where we begin
  Int_t offsetch = 6*ch;
  

  for (Int_t isec = 0; isec < kNsect; isec++){
    for(Int_t ipl = 0; ipl < kNplan; ipl++){
      Int_t det   = offsetch+isec*30+ipl;
      AliTRDpadPlane *padPlane = trdGeo->GetPadPlane(ipl,ch);
      for (Int_t icol=0; icol<padPlane->GetNcols(); icol++){
	poslocal[0] = trdGeo->GetTime0(ipl);
	poslocal[2] = padPlane->GetRowPos(0);
	poslocal[1] = padPlane->GetColPos(icol);
        trdGeo->RotateBack(det,poslocal,posglobal);
	Int_t binx = 1+TMath::Nint((posglobal[0]+400.0)*0.5);
	Int_t biny = 1+TMath::Nint((posglobal[1]+400.0)*0.5);
	his->SetBinContent(binx,biny,fData[det]);
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
TH2F *AliTRDCalDet::MakeHisto2DSmPl(Int_t sm, Int_t pl, Float_t min, Float_t max, Int_t type)
{
  //
  // Make 2D graph
  // sm    - supermodule number
  // pl    - plane number
  // type -1 = user defined range
  //       0 = nsigma cut nsigma=min
  //       1 = delta cut around median delta=min
  //

  gStyle->SetPalette(1);
  if (type>=0){
    if (type==0){
      // nsigma range
      Float_t mean  = GetMean();
      Float_t sigma = GetRMS();
      Float_t nsigma = TMath::Abs(min);
      min = mean-nsigma*sigma;
      max = mean+nsigma*sigma;
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
      min = mean-sigma;
      max = mean+sigma;

    }
  }
     
  AliTRDgeometry *trdGeo = new AliTRDgeometry();
  AliTRDpadPlane *padPlane0 = trdGeo->GetPadPlane(pl,0);
  Double_t row0    = padPlane0->GetRow0();
  Double_t col0    = padPlane0->GetCol0();

  char  name[1000];
  snprintf(name,1000,"%s CalDet 2D sm %d and pl %d",GetTitle(),sm,pl);
  TH2F * his = new TH2F( name, name, 5,  -TMath::Abs(row0),  TMath::Abs(row0)
                                   , 4,-2*TMath::Abs(col0),2*TMath::Abs(col0));

  // Where we begin
  Int_t offsetsmpl = 30*sm+pl;
  

  for (Int_t k = 0; k < kNcham; k++){
    Int_t det = offsetsmpl+k*6;
    Int_t kb  = kNcham-1-k;
    his->SetBinContent(kb+1,2,fData[det]);
    his->SetBinContent(kb+1,3,fData[det]);
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
void AliTRDCalDet::Add(Float_t c1)
{
  //
  // Add constant for all detectors
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
     fData[idet] += c1;
  }
}

//_____________________________________________________________________________
void AliTRDCalDet::Multiply(Float_t c1)
{
    //
    // multiply constant for all detectors
    //
    for (Int_t idet = 0; idet < kNdet; idet++) {
      fData[idet] *= c1;
    }
}

//_____________________________________________________________________________
void AliTRDCalDet::Add(const AliTRDCalDet * calDet, Double_t c1)
{
    //
    // add caldet channel by channel multiplied by c1
    //
    for (Int_t idet = 0; idet < kNdet; idet++) {
      fData[idet] += calDet->GetValue(idet)*c1;
    }
}

//_____________________________________________________________________________
void AliTRDCalDet::Multiply(const AliTRDCalDet * calDet)
{
    //
    // multiply caldet channel by channel 
    //
    for (Int_t idet = 0; idet < kNdet; idet++) {
      fData[idet] *= calDet->GetValue(idet);
    }
}

//_____________________________________________________________________________
void AliTRDCalDet::Divide(const AliTRDCalDet * calDet)
{
    //
    // divide caldet channel by channel 
    //
 Float_t kEpsilon=0.00000000000000001;
    for (Int_t idet = 0; idet < kNdet; idet++) {
      if (TMath::Abs(calDet->GetValue(idet))>kEpsilon){
	fData[idet] /= calDet->GetValue(idet);
	}
    }
}
//_____________________________________________________________________________
Double_t AliTRDCalDet::CalcMean(Bool_t wghtPads, Int_t &calib)
{
  // Calculate the mean value after rejection of the chambers not calibrated
  // wghPads = kTRUE weighted with the number of pads in case of a AliTRDCalPad created (t0)
  // calib = number of used chambers for the mean calculation

  Int_t iSM;
  Double_t sum = 0.0;
  Int_t ndet = 0;
  Double_t meanALL = 0.0;
  Double_t meanWP = 0.0;
  Double_t pads=0.0;
  Double_t padsALL=(144*16*24+144*12*6)*18;
  Double_t *meanSM = new Double_t[18];
  Double_t *meanSMWP = new Double_t[18];
  
  for (Int_t i = 0; i < 18; i++) {
    meanSM[i]=0.0;
    meanSMWP[i]=0.0;
  }

  Int_t det = 0;
  while(det < 540) {
    Float_t val= fData[det];
    iSM = (Int_t)(det / (6*5));
    pads=(((Int_t) (det % (6 * 5)) / 6) == 2) ? 144*12 : 144*16;
    meanALL+=val/540.;
    meanSM[iSM]+=val/30.;
    meanWP+=val*(pads/padsALL);
    meanSMWP[iSM]+=val*(pads/(padsALL/18.));
    
    /*    
	  printf(" det %d  val %.3f meanALL %.5f meanWP %.5f meanSM[%d] %.5f meanSMWP[%d] %.5f \n",
	  det,
	  val,
	  meanALL,
	  meanWP,
	  iSM,
	  meanSM[iSM],
	  iSM,
	  meanSMWP[iSM]);
    */
   
    det++;
  }

  // debug
  /*
    printf(" ALL %.5f \n",meanALL);
    printf(" WP %.5f \n",meanWP);
    for(Int_t i=0; i<18; i++) printf(" SM %02d %.5f \n",i,meanSM[i]);
    for(Int_t i=0; i<18; i++) printf(" SM %02d %.5f \n",i,meanSMWP[i]);
  */

  det=0;
  while(det < 540) {
    Float_t val= fData[det];
    if (( (!wghtPads) &&
          (TMath::Abs(val - meanALL) > 0.0001) &&
          (TMath::Abs(val - meanSM[(Int_t)(det / (6*5))]) > 0.0001) ) ||
        ( (wghtPads) &&
          (TMath::Abs(val - meanWP) > 0.0001) &&
          (TMath::Abs(val - meanSMWP[(Int_t)(det / (6*5) )]) > 0.0001) )
        ) {
      sum+=val;
      ndet++;
    }
    det++;
  }

  delete []meanSM;
  delete []meanSMWP;

  calib=ndet;
  return (sum!=0.0 ? sum/ndet : -1);
}
//_____________________________________________________________________________
Double_t AliTRDCalDet::CalcMean(Bool_t wghtPads)
{
  // Calculate the mean value after rejection of the chambers not calibrated
  // wghPads = kTRUE weighted with the number of pads in case of a AliTRDCalPad created (t0)
  // calib = number of used chambers for the mean calculation

  Int_t iSM;
  Double_t sum = 0.0;
  Int_t ndet = 0;
  Double_t meanALL = 0.0;
  Double_t meanWP = 0.0;
  Double_t pads=0.0;
  Double_t padsALL=(144*16*24+144*12*6)*18;
  Double_t *meanSM = new Double_t[18];
  Double_t *meanSMWP = new Double_t[18];
  
  for (Int_t i = 0; i < 18; i++) {
    meanSM[i]=0.0;
    meanSMWP[i]=0.0;
  }

  Int_t det = 0;
  while(det < 540) {
    Float_t val= fData[det];
    iSM = (Int_t)(det / (6*5));
    pads=(((Int_t) (det % (6 * 5)) / 6) == 2) ? 144*12 : 144*16;
    meanALL+=val/540.;
    meanSM[iSM]+=val/30.;
    meanWP+=val*(pads/padsALL);
    meanSMWP[iSM]+=val*(pads/(padsALL/18.));
    
    /*    
	  printf(" det %d  val %.3f meanALL %.5f meanWP %.5f meanSM[%d] %.5f meanSMWP[%d] %.5f \n",
	  det,
	  val,
	  meanALL,
	  meanWP,
	  iSM,
	  meanSM[iSM],
	  iSM,
	  meanSMWP[iSM]);
    */
   
    det++;
  }

  // debug
  /*
    printf(" ALL %.5f \n",meanALL);
    printf(" WP %.5f \n",meanWP);
    for(Int_t i=0; i<18; i++) printf(" SM %02d %.5f \n",i,meanSM[i]);
    for(Int_t i=0; i<18; i++) printf(" SM %02d %.5f \n",i,meanSMWP[i]);
  */

  det=0;
  while(det < 540) {
    Float_t val= fData[det];
    if (( (!wghtPads) &&
          (TMath::Abs(val - meanALL) > 0.0001) &&
          (TMath::Abs(val - meanSM[(Int_t)(det / (6*5))]) > 0.0001) ) ||
        ( (wghtPads) &&
          (TMath::Abs(val - meanWP) > 0.0001) &&
          (TMath::Abs(val - meanSMWP[(Int_t)(det / (6*5) )]) > 0.0001) )
        ) {
      sum+=val;
      ndet++;
    }
    det++;
  }

  delete []meanSM;
  delete []meanSMWP;

  return (sum!=0.0 ? sum/ndet : -1);
}
//_____________________________________________________________________________
Double_t AliTRDCalDet::CalcRMS(Bool_t wghtPads, Int_t &calib)
{
  // Calculate the RMS value after rejection of the chambers not calibrated
  // wghPads = kTRUE weighted with the number of pads in case of a AliTRDCalPad created (t0)
  // calib = number of used chambers for the mean calculation
  
  Int_t iSM;
  Double_t sum = 0.0;
  Int_t ndet = 0;
  Double_t meanALL = 0.0;
  Double_t meanWP = 0.0;
  Double_t pads=0.0;
  Double_t padsALL=(144*16*24+144*12*6)*18;
  Double_t *meanSM = new Double_t[18];
  Double_t *meanSMWP = new Double_t[18];
  
  for (Int_t i = 0; i < 18; i++) {
    meanSM[i]=0.0;
    meanSMWP[i]=0.0;
  }
  
  Int_t det = 0;
  while(det < 540) {
    Float_t val= fData[det];
    iSM = (Int_t)(det / (6*5));
    pads=(((Int_t) (det % (6 * 5)) / 6) == 2) ? 144*12 : 144*16;
    meanALL+=val/540.;
    meanSM[iSM]+=val/30.;
    meanWP+=val*(pads/padsALL);
    meanSMWP[iSM]+=val*(pads/(padsALL/18.));
    det++;
  }
  
  Double_t mean=0.0;
  if(!wghtPads) mean= meanALL;
  if(wghtPads) mean= meanWP;
  
  det=0;
  while(det < 540) {
    Float_t val= fData[det];
    if (( (!wghtPads) &&
		 (TMath::Abs(val - meanALL) > 0.0001) &&
		 (TMath::Abs(val - meanSM[(Int_t)(det / (6*5))]) > 0.0001) ) ||
        ( (wghtPads) &&
		 (TMath::Abs(val - meanWP) > 0.0001) &&
		 (TMath::Abs(val - meanSMWP[(Int_t)(det / (6*5) )]) > 0.0001) )
        ) {
      sum+=(val-mean)*(val-mean);
      ndet++;
    }
    det++;
  }
  
  delete []meanSM;
  delete []meanSMWP;
  
  calib=ndet;
  return (sum!=0.0 ? TMath::Sqrt(sum/ndet) : -1);
}
//_____________________________________________________________________________
Double_t AliTRDCalDet::CalcRMS(Bool_t wghtPads)
{
  // Calculate the RMS value after rejection of the chambers not calibrated
  // wghPads = kTRUE weighted with the number of pads in case of a AliTRDCalPad created (t0)
  // calib = number of used chambers for the mean calculation
  
  Int_t iSM;
  Double_t sum = 0.0;
  Int_t ndet = 0;
  Double_t meanALL = 0.0;
  Double_t meanWP = 0.0;
  Double_t pads=0.0;
  Double_t padsALL=(144*16*24+144*12*6)*18;
  Double_t *meanSM = new Double_t[18];
  Double_t *meanSMWP = new Double_t[18];
  
  for (Int_t i = 0; i < 18; i++) {
    meanSM[i]=0.0;
    meanSMWP[i]=0.0;
  }
  
  Int_t det = 0;
  while(det < 540) {
    Float_t val= fData[det];
    iSM = (Int_t)(det / (6*5));
    pads=(((Int_t) (det % (6 * 5)) / 6) == 2) ? 144*12 : 144*16;
    meanALL+=val/540.;
    meanSM[iSM]+=val/30.;
    meanWP+=val*(pads/padsALL);
    meanSMWP[iSM]+=val*(pads/(padsALL/18.));
    det++;
  }
  
  Double_t mean=0.0;
  if(!wghtPads) mean= meanALL;
  if(wghtPads) mean= meanWP;
  
  det=0;
  while(det < 540) {
    Float_t val= fData[det];
    if (( (!wghtPads) &&
		 (TMath::Abs(val - meanALL) > 0.0001) &&
		 (TMath::Abs(val - meanSM[(Int_t)(det / (6*5))]) > 0.0001) ) ||
        ( (wghtPads) &&
		 (TMath::Abs(val - meanWP) > 0.0001) &&
		 (TMath::Abs(val - meanSMWP[(Int_t)(det / (6*5) )]) > 0.0001) )
        ) {
      sum+=(val-mean)*(val-mean);
      ndet++;
    }
    det++;
  }
  
  delete []meanSM;
  delete []meanSMWP;
  
  return (sum!=0.0 ? TMath::Sqrt(sum/ndet) : -1);
}
