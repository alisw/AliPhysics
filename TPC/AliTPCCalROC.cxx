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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Calibration base class for a single ROC                                  //
//  Contains one float value per pad                                         //
//     mapping of the pads taken form AliTPCROC                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTPCCalROC.h"
#include "TMath.h"
#include "TClass.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "AliMathBase.h"
ClassImp(AliTPCCalROC)


//_____________________________________________________________________________
AliTPCCalROC::AliTPCCalROC()
             :TObject(),
	      fSector(0),
	      fNChannels(0),
	      fNRows(0),
	      fIndexes(0),
	      fData(0)
{
  //
  // Default constructor
  //

}

//_____________________________________________________________________________
AliTPCCalROC::AliTPCCalROC(UInt_t  sector)
             :TObject(),
	      fSector(0),
	      fNChannels(0),
	      fNRows(0),
	      fIndexes(0),
	      fData(0)
{
  //
  // Constructor that initializes a given sector
  //
  fSector = sector;
  fNChannels    =  AliTPCROC::Instance()->GetNChannels(fSector);
  fNRows        =  AliTPCROC::Instance()->GetNRows(fSector);
  fIndexes      =  AliTPCROC::Instance()->GetRowIndexes(fSector);
  fData = new Float_t[fNChannels];
  for (UInt_t  idata = 0; idata< fNChannels; idata++) fData[idata] = 0.;
}

//_____________________________________________________________________________
AliTPCCalROC::AliTPCCalROC(const AliTPCCalROC &c)
             :TObject(c),
	      fSector(0),
	      fNChannels(0),
	      fNRows(0),
	      fIndexes(0),
	      fData(0)
{
  //
  // AliTPCCalROC copy constructor
  //
  fSector = c.fSector;
  fNChannels    =  AliTPCROC::Instance()->GetNChannels(fSector);
  fNRows        =  AliTPCROC::Instance()->GetNRows(fSector);
  fIndexes      =  AliTPCROC::Instance()->GetRowIndexes(fSector);
  //
  fData   = new Float_t[fNChannels];
  for (UInt_t  idata = 0; idata< fNChannels; idata++) fData[idata] = c.fData[idata];
}
//____________________________________________________________________________
AliTPCCalROC & AliTPCCalROC::operator =(const AliTPCCalROC & param)
{
  //
  // assignment operator - dummy
  //
  fData=param.fData;
  return (*this);
}


//_____________________________________________________________________________
AliTPCCalROC::~AliTPCCalROC()
{
  //
  // AliTPCCalROC destructor
  //
  if (fData) {
    delete [] fData;
    fData = 0;
  }
}



void AliTPCCalROC::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliTPCCalROC.
   if (R__b.IsReading()) {
      AliTPCCalROC::Class()->ReadBuffer(R__b, this);
      fIndexes =  AliTPCROC::Instance()->GetRowIndexes(fSector);
   } else {
      AliTPCCalROC::Class()->WriteBuffer(R__b,this);
   }
}


Double_t AliTPCCalROC::GetLTM(Double_t *sigma, Double_t fraction){
  //
  //  Calculate LTM mean and sigma
  //
  Double_t *ddata = new Double_t[fNChannels];
  Double_t mean=0, lsigma=0;
  Int_t hh = TMath::Min(TMath::Nint(fraction *fNChannels), Int_t(fNChannels));
  for (UInt_t i=0;i<fNChannels;i++) ddata[i]= fData[i];
  AliMathBase::EvaluateUni(UInt_t(fNChannels),ddata, mean, lsigma, hh);
  if (sigma) *sigma=lsigma;
  delete [] ddata;
  return mean;
}

TH1F * AliTPCCalROC::MakeHisto1D(Float_t min, Float_t max,Int_t type){
  //
  // make 1D histo
  // type -1 = user defined range
  //       0 = nsigma cut nsigma=min
  //       1 = delta cut around median delta=min
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
  sprintf(name,"%s ROC 1D%d",GetTitle(),fSector);
  TH1F * his = new TH1F(name,name,100, min,max);
  for (UInt_t irow=0; irow<fNRows; irow++){
    UInt_t npads = (Int_t)GetNPads(irow);
    for (UInt_t ipad=0; ipad<=npads; ipad++){
      his->Fill(GetValue(irow,ipad));
    }
  }
  return his;
}



TH2F * AliTPCCalROC::MakeHisto2D(Float_t min, Float_t max,Int_t type){
  //
  // make 2D histo
  // type -1 = user defined range
  //       0 = nsigma cut nsigma=min
  //       1 = delta cut around median delta=min
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
  UInt_t maxPad = 0;
  for (UInt_t irow=0; irow<fNRows; irow++){
    if (GetNPads(irow)>maxPad) maxPad = GetNPads(irow);
  }
  char  name[1000];
  sprintf(name,"%s ROC%d",GetTitle(),fSector);
  TH2F * his = new TH2F(name,name,fNRows+10,-5, fNRows+5, maxPad+10, -(Int_t(maxPad/2))-5, maxPad/2+5);
  for (UInt_t irow=0; irow<fNRows; irow++){
    UInt_t npads = (Int_t)GetNPads(irow);
    for (UInt_t ipad=0; ipad<=npads; ipad++){
      his->Fill(irow+0.5,Int_t(ipad)-Int_t(npads/2)+0.5,GetValue(irow,ipad));
    }
  }
  his->SetMaximum(max);
  his->SetMinimum(min);
  return his;
}

TH2F * AliTPCCalROC::MakeHistoOutliers(Float_t delta, Float_t fraction, Int_t type){
  //
  // Make Histogram with outliers
  // mode = 0 - sigma cut used
  // mode = 1 - absolute cut used
  // fraction - fraction of values used to define sigma
  // delta - in mode 0 - nsigma cut
  //            mode 1 - delta cut
  Double_t sigma;
  Float_t mean  = GetLTM(&sigma,fraction);  
  if (type==0) delta*=sigma; 
  UInt_t maxPad = 0;
  for (UInt_t irow=0; irow<fNRows; irow++){
    if (GetNPads(irow)>maxPad) maxPad = GetNPads(irow);
  }

  char  name[1000];
  sprintf(name,"%s ROC Outliers%d",GetTitle(),fSector);
  TH2F * his = new TH2F(name,name,fNRows+10,-5, fNRows+5, maxPad+10, -(Int_t(maxPad/2))-5, maxPad/2+5);
  for (UInt_t irow=0; irow<fNRows; irow++){
    UInt_t npads = (Int_t)GetNPads(irow);
    for (UInt_t ipad=0; ipad<=npads; ipad++){
      if (TMath::Abs(GetValue(irow,ipad)-mean)>delta)
	his->Fill(irow+0.5,Int_t(ipad)-Int_t(npads/2)+0.5,1);
    }
  }
  return his;
}



void AliTPCCalROC::Draw(Option_t* opt){
  //
  // create histogram with values and draw it
  //
  TH1 * his=0; 
  TString option=opt;
  option.ToUpper();
  if (option.Contains("1D")){
    his = MakeHisto1D();
  }
  else{
    his = MakeHisto2D();
  }
  his->Draw(option);
}





void AliTPCCalROC::Test(){
  //
  // example function to show functionality and tes AliTPCCalROC
  //
  AliTPCCalROC  roc0(0);  
  for (UInt_t irow = 0; irow <roc0.GetNrows(); irow++){
    for (UInt_t ipad = 0; ipad <roc0.GetNPads(irow); ipad++){
      Float_t value  = irow+ipad/1000.;
      roc0.SetValue(irow,ipad,value);
    }
  }
  //
  AliTPCCalROC roc1(roc0);
  for (UInt_t irow = 0; irow <roc1.GetNrows(); irow++){
    for (UInt_t ipad = 0; ipad <roc1.GetNPads(irow); ipad++){
      Float_t value  = irow+ipad/1000.;
      if (roc1.GetValue(irow,ipad)!=value){
	printf("Read/Write error\trow=%d\tpad=%d\n",irow,ipad);
      }
    }
  }  
  TFile f("calcTest.root","recreate");
  roc0.Write("Roc0");
  AliTPCCalROC * roc2 = (AliTPCCalROC*)f.Get("Roc0");
  f.Close();
  //
  for (UInt_t irow = 0; irow <roc0.GetNrows(); irow++){
    if (roc0.GetNPads(irow)!=roc2->GetNPads(irow))
      printf("NPads - Read/Write error\trow=%d\n",irow);
    for (UInt_t ipad = 0; ipad <roc1.GetNPads(irow); ipad++){
      Float_t value  = irow+ipad/1000.;
      if (roc2->GetValue(irow,ipad)!=value){
	printf("Read/Write error\trow=%d\tpad=%d\n",irow,ipad);
      }
    }
  }   
}

