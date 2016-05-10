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


/// \class AliTPCCalROC
///
///     Calibration base class for a single ROC
///     Contains one float value per pad
///     mapping of the pads taken form AliTPCROC
///

// ROOT includes 
#include "TMath.h"
#include "TClass.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TArrayI.h"

//
//
#include "AliTPCCalROC.h"
#include "AliMathBase.h"

#include "TRandom3.h"      // only needed by the AliTPCCalROCTest() method

/// \cond CLASSIMP
ClassImp(AliTPCCalROC)
/// \endcond


//_____________________________________________________________________________
AliTPCCalROC::AliTPCCalROC()
             :TNamed(),
	      fSector(0),
	      fNChannels(0),
	      fNRows(0),
	      fkIndexes(0),
	      fData(0)
{
  //
  // Default constructor
  //

}

//_____________________________________________________________________________
AliTPCCalROC::AliTPCCalROC(UInt_t  sector)
             :TNamed(),
	      fSector(0),
	      fNChannels(0),
	      fNRows(0),
	      fkIndexes(0),
	      fData(0)
{
  /// Constructor that initializes a given sector

  fSector = sector;
  fNChannels    =  AliTPCROC::Instance()->GetNChannels(fSector);
  fNRows        =  AliTPCROC::Instance()->GetNRows(fSector);
  fkIndexes      =  AliTPCROC::Instance()->GetRowIndexes(fSector);
  fData = new Float_t[fNChannels];
  Reset();

}

//_____________________________________________________________________________
AliTPCCalROC::AliTPCCalROC(const AliTPCCalROC &c)
             :TNamed(c),
	      fSector(0),
	      fNChannels(0),
	      fNRows(0),
	      fkIndexes(0),
	      fData(0)
{
  /// AliTPCCalROC copy constructor

  fSector = c.fSector;
  fNChannels    =  AliTPCROC::Instance()->GetNChannels(fSector);
  fNRows        =  AliTPCROC::Instance()->GetNRows(fSector);
  fkIndexes      =  AliTPCROC::Instance()->GetRowIndexes(fSector);
  //
  fData   = new Float_t[fNChannels];
  for (UInt_t  idata = 0; idata< fNChannels; idata++) fData[idata] = c.fData[idata];
}
//____________________________________________________________________________
AliTPCCalROC & AliTPCCalROC::operator =(const AliTPCCalROC & param)
{
  /// assignment operator - dummy

  if (this == &param) return (*this);
  fSector       = param.fSector;
  fNChannels    =  AliTPCROC::Instance()->GetNChannels(fSector);
  fNRows        =  AliTPCROC::Instance()->GetNRows(fSector);
  fkIndexes     =  AliTPCROC::Instance()->GetRowIndexes(fSector);
  //
  if (fData) delete [] fData;
  fData         = new Float_t[fNChannels];
  for (UInt_t  idata = 0; idata< fNChannels; idata++) fData[idata] = param.fData[idata];
  return (*this);
}


//_____________________________________________________________________________
AliTPCCalROC::~AliTPCCalROC()
{
  /// AliTPCCalROC destructor

  if (fData) {
    delete [] fData;
    fData = 0;
  }
}



void AliTPCCalROC::Streamer(TBuffer &R__b)
{
   /// Stream an object of class AliTPCCalROC.

   if (R__b.IsReading()) {
      AliTPCCalROC::Class()->ReadBuffer(R__b, this);
      fkIndexes =  AliTPCROC::Instance()->GetRowIndexes(fSector);
   } else {
      AliTPCCalROC::Class()->WriteBuffer(R__b,this);
   }
}

void AliTPCCalROC::Print(Option_t* /*option*/) const{
  //
  // Print summary content
  //
  printf("|ROC%d|\t%f|\t%f|\t%f|\n",fSector, GetMean(),GetMedian(),GetRMS());
}



//

Bool_t AliTPCCalROC::MedianFilter(Int_t deltaRow, Int_t deltaPad, AliTPCCalROC* outlierROC,  Bool_t doEdge){
  ///   Modify content of the object - raplace value by median in neighorhood

  Float_t *newBuffer=new Float_t[fNChannels] ;
  Double_t *cacheBuffer=new Double_t[fNChannels];
  //
  for (Int_t iRow=0; iRow< Int_t(fNRows); iRow++){
    Int_t nPads=GetNPads(iRow); // number of rows in current row
    for (Int_t iPad=0; iPad<nPads; iPad++){
      Double_t value=GetValue(iRow,iPad);
      Int_t counter=0;
      //
      for (Int_t dRow=-deltaRow; dRow<=deltaRow; dRow++){
	Int_t jRow=iRow+dRow;  //take the row - mirror if ouside of range
	Float_t sign0=1.;
	if (jRow<0) sign0=-1.;
	if (UInt_t(jRow)>=fNRows) sign0=-1.;	
	Int_t jPads= GetNPads(iRow+sign0*dRow);
	Int_t offset=(nPads-jPads)/2.;
	//
	for (Int_t dPad=-deltaPad; dPad<=deltaPad; dPad++){
	  Float_t sign=sign0;
	  Int_t jPad=iPad+dPad+offset;  //take the pad - mirror if ouside of range
	  Int_t kRow=jRow;
	  if (jPad<0) sign=-1;
	  if (jPad>=jPads) sign=-1;
	  if (sign<0){
	    kRow=iRow-dRow;
	    jPad=iPad-dPad+offset;
	    if (!doEdge) continue;
	  }	  
	  if (IsInRange(UInt_t(kRow),jPad)){
	    Bool_t isOutlier=(outlierROC==NULL)?kFALSE:outlierROC->GetValue(kRow,jPad)>0;
	    if (!isOutlier){
	      cacheBuffer[counter]=sign*(GetValue(kRow,jPad)-value);
	      counter++;
	    }
	  }
	}
      }
      newBuffer[fkIndexes[iRow]+iPad]=0.;
      if (counter>1) newBuffer[fkIndexes[iRow]+iPad] = TMath::Median(counter, cacheBuffer)+value;
    }
  }
  memcpy(fData, newBuffer,GetNchannels()*sizeof(Float_t));
  delete []newBuffer;
  delete []cacheBuffer;
  return kTRUE;
}



Bool_t AliTPCCalROC::LTMFilter(Int_t deltaRow, Int_t deltaPad, Float_t fraction, Int_t type, AliTPCCalROC* outlierROC,  Bool_t doEdge){
  ///   Modify content of the class
  ///   write LTM mean or median

  if (fraction<0 || fraction>1) return kFALSE;
  Float_t *newBuffer=new Float_t[fNChannels] ;
  Double_t *cacheBuffer=new Double_t[fNChannels];
  //
  for (Int_t iRow=0; iRow< Int_t(fNRows); iRow++){
    Int_t nPads=GetNPads(iRow); // number of rows in current row
    for (Int_t iPad=0; iPad<nPads; iPad++){
      Double_t value=GetValue(iRow,iPad);
      Int_t counter=0;
      //
      for (Int_t dRow=-deltaRow; dRow<=deltaRow; dRow++){
	Int_t jRow=iRow+dRow;  //take the row - mirror if ouside of range
	Float_t sign0=1.;
	if (jRow<0) sign0=-1.;
	if (UInt_t(jRow)>=fNRows) sign0=-1.;	
	Int_t jPads= GetNPads(iRow+sign0*dRow);
	Int_t offset=(nPads-jPads)/2.;
	//
	for (Int_t dPad=-deltaPad; dPad<=deltaPad; dPad++){
	  Float_t sign=sign0;
	  Int_t jPad=iPad+dPad+offset;  //take the pad - mirror if ouside of range
	  Int_t kRow=jRow;
	  if (jPad<0) sign=-1;
	  if (jPad>=jPads) sign=-1;
	  if (sign<0){
	    if (!doEdge) continue;
	    kRow=iRow-dRow;
	    jPad=iPad-dPad+offset;
	  } 
	  if (IsInRange(UInt_t(kRow),jPad)){
	    Bool_t isOutlier=(outlierROC==NULL)?kFALSE:outlierROC->GetValue(kRow,jPad)>0;
	    if (!isOutlier){
	      cacheBuffer[counter]=sign*(GetValue(kRow,jPad)-value);
	      counter++;
	    }
	  }
	}
      }
      Double_t mean=0,rms=0;
      if (TMath::Nint(fraction*Double_t(counter))>1 ) AliMathBase::EvaluateUni(counter,cacheBuffer,mean,rms,TMath::Nint(fraction*Double_t(counter)));
      mean+=value;
      newBuffer[fkIndexes[iRow]+iPad] = (type==0)? mean:rms;
    }
  }
  memcpy(fData, newBuffer,GetNchannels()*sizeof(Float_t));
  delete []newBuffer;
  delete []cacheBuffer;
  return kTRUE;
}

Bool_t AliTPCCalROC::Convolute(Double_t sigmaPad, Double_t sigmaRow, AliTPCCalROC*outlierROC, TF1 */*fpad*/, TF1 */*frow*/){
  /// convolute the calibration with function fpad,frow
  /// in range +-4 sigma

  Float_t *newBuffer=new Float_t[fNChannels] ;
  //
  for (Int_t iRow=0; iRow< Int_t(fNRows); iRow++){
    Int_t nPads=GetNPads(iRow); // number of rows in current row
    for (Int_t iPad=0; iPad<nPads; iPad++){
      Int_t jRow0=TMath::Max(TMath::Nint(iRow-sigmaRow*4.),0);
      Int_t jRow1=TMath::Min(TMath::Nint(iRow+sigmaRow*4.),Int_t(fNRows));
      Int_t jPad0=TMath::Max(TMath::Nint(iPad-sigmaPad*4.),0);
      Int_t jPad1=TMath::Min(TMath::Nint(iRow+sigmaPad*4.),Int_t(nPads));
      //
      Double_t sumW=0;
      Double_t sumCal=0;
      for (Int_t jRow=jRow0; jRow<=jRow1; jRow++){
	for (Int_t jPad=jPad0; jPad<=jPad1; jPad++){
	  if (!IsInRange(jRow,jPad)) continue;
	  Bool_t isOutlier=(outlierROC==NULL)?kFALSE:outlierROC->GetValue(jRow,jPad)>0;
	  if (!isOutlier){
	    Double_t weight= TMath::Gaus(jPad,iPad,sigmaPad)*TMath::Gaus(jRow,iRow,sigmaRow);
	    sumCal+=weight*GetValue(jRow,jPad);
	    sumW+=weight;
	  }
	}
      }
      if (sumW>0){
	sumCal/=sumW;
	newBuffer[fkIndexes[iRow]+iPad]=sumCal;
      }
    }
  }
  memcpy(fData, newBuffer,GetNchannels()*sizeof(Float_t));
  delete []newBuffer;
  return kTRUE;
}


//





// algebra fuctions:

void AliTPCCalROC::Add(Float_t c1){
  /// add c1 to each channel of the ROC  

  for (UInt_t  idata = 0; idata< fNChannels; idata++) fData[idata]+=c1;
}


void AliTPCCalROC::Multiply(Float_t c1){
  /// multiply each channel of the ROC with c1

  for (UInt_t  idata = 0; idata< fNChannels; idata++) fData[idata]*=c1;
}


void AliTPCCalROC::Add(const AliTPCCalROC * roc, Double_t c1){
  /// multiply AliTPCCalROC roc by c1 and add each channel to the coresponing channel in the ROC
  ///  - pad by pad 

  if (!roc) return;
  for (UInt_t  idata = 0; idata< fNChannels; idata++){
    fData[idata]+=roc->fData[idata]*c1;
  }
}


void AliTPCCalROC::Multiply(const AliTPCCalROC*  roc) {
  /// multiply each channel of the ROC with the coresponding channel of 'roc'
  ///     - pad by pad -

  if (!roc) return;
  for (UInt_t  idata = 0; idata< fNChannels; idata++){
    fData[idata]*=roc->fData[idata];
  }
}


void AliTPCCalROC::Divide(const AliTPCCalROC*  roc) {
  /// divide each channel of the ROC by the coresponding channel of 'roc'
  ///     - pad by pad -

  if (!roc) return;
  Float_t kEpsilon=0.00000000000000001;
  for (UInt_t  idata = 0; idata< fNChannels; idata++){
    if (TMath::Abs(roc->fData[idata])>kEpsilon)
      fData[idata]/=roc->fData[idata];
  }
}

void AliTPCCalROC::Reset()
{
  /// reset to ZERO
  
  memset(fData,0,sizeof(Float_t)*fNChannels); // set all values to 0
}

//____________________________________________________________
Double_t AliTPCCalROC::GetStats(EStatType statType, AliTPCCalROC *const outlierROC, EPadType padType) const
{
  ///  returns the mean or median or RMS value depending on the statType 
  ///  of the ROC or medium or long pad region
  ///  pads with value != 0 in outlierROC are not used for the calculation
  ///  return 0 if no data is accepted by the outlier cuts 

  if(fSector<36) padType=kAll;
  if (!outlierROC) {
    if(padType==kAll) {
      if     (statType==kMean)       return TMath::Mean      (fNChannels, fData);
      else if(statType==kMedian)     return TMath::Median    (fNChannels, fData);
      else if(statType==kRMS)        return TMath::RMS       (fNChannels, fData);
      else if(statType==kMinElement) return TMath::MinElement(fNChannels, fData);
      else if(statType==kMaxElement) return TMath::MaxElement(fNChannels, fData);
    }
    else if(padType==kOROCmedium) {
      if     (statType==kMean)       return TMath::Mean      (fkIndexes[64], fData);
      else if(statType==kMedian)     return TMath::Median    (fkIndexes[64], fData);
      else if(statType==kRMS)        return TMath::RMS       (fkIndexes[64], fData);
      else if(statType==kMinElement) return TMath::MinElement(fkIndexes[64], fData);
      else if(statType==kMaxElement) return TMath::MaxElement(fkIndexes[64], fData);
    }
    else if(padType==kOROClong) {
      const Int_t offset=fkIndexes[64];
      if     (statType==kMean)       return TMath::Mean      (fNChannels-offset, fData+offset);
      else if(statType==kMedian)     return TMath::Median    (fNChannels-offset, fData+offset);
      else if(statType==kRMS)        return TMath::RMS       (fNChannels-offset, fData+offset);
      else if(statType==kMinElement) return TMath::MinElement(fNChannels-offset, fData+offset);
      else if(statType==kMaxElement) return TMath::MaxElement(fNChannels-offset, fData+offset);
    }
  }//end of no outlierROC

  Float_t ddata[fNChannels];

  Int_t nPoints = 0;
  UInt_t indexMin = 0;
  UInt_t indexMax = fNChannels;
  if(padType==kOROCmedium) indexMax=fkIndexes[64];
  else if(padType==kOROClong) indexMin=fkIndexes[64];

  for (UInt_t i=indexMin;i<indexMax;i++) {
    if (outlierROC->GetValue(i)>1e-20) continue;
    ddata[nPoints]= fData[i];
    nPoints++;
  }
  Double_t value = 0;
  if(nPoints>0){
    if     (statType==kMean)       value = TMath::Mean      (nPoints, ddata);
    else if(statType==kMedian)     value = TMath::Median    (nPoints, ddata);
    else if(statType==kRMS)        value = TMath::RMS       (nPoints, ddata);
    else if(statType==kMinElement) value = TMath::MinElement(nPoints, ddata);
    else if(statType==kMaxElement) value = TMath::MaxElement(nPoints, ddata);
  }
  return value;
}

//_____________________________________________________________________________
Double_t AliTPCCalROC::GetMean(AliTPCCalROC *const outlierROC,EPadType padType) const
{
  ///  returns the mean value of the ROC
  ///  pads with value != 0 in outlierROC are not used for the calculation
  ///  return 0 if no data is accepted by the outlier cuts 

  return GetStats(kMean,outlierROC,padType);
}

//_____________________________________________________________________________
Double_t AliTPCCalROC::GetMedian(AliTPCCalROC *const outlierROC,EPadType padType) const
{
  ///  returns the median value of the ROC
  ///  pads with value != 0 in outlierROC are not used for the calculation
  ///  return 0 if no data is accepted by the outlier cuts 

  return GetStats(kMedian,outlierROC,padType);
}

//_____________________________________________________________________________
Double_t AliTPCCalROC::GetRMS(AliTPCCalROC *const outlierROC,EPadType padType) const
{
   ///  returns the RMS value of the ROC
   ///  pads with value != 0 in outlierROC are not used for the calculation
   ///  return 0 if no data is accepted by the outlier cuts 

  return GetStats(kRMS,outlierROC,padType);
}

//_____________________________________________________________________________
Double_t AliTPCCalROC::GetMinElement(AliTPCCalROC *const outlierROC,EPadType padType) const
{
  ///  returns the MinElement value of the ROC
  ///  pads with value != 0 in outlierROC are not used for the calculation
  ///  return 0 if no data is accepted by the outlier cuts

  return GetStats(kMinElement,outlierROC,padType);
}

//_____________________________________________________________________________
Double_t AliTPCCalROC::GetMaxElement(AliTPCCalROC *const outlierROC,EPadType padType) const
{
  ///  returns the MaxElement value of the ROC
  ///  pads with value != 0 in outlierROC are not used for the calculation
  ///  return 0 if no data is accepted by the outlier cuts

  return GetStats(kMaxElement,outlierROC,padType);
}

//_____________________________________________________________________________
Double_t AliTPCCalROC::GetLTM(Double_t *const sigma, Double_t fraction, AliTPCCalROC *const outlierROC, EPadType padType){
  ///  returns the LTM and sigma
  ///  pads with value != 0 in outlierROC are not used for the calculation
  ///  return 0 if no data is accepted by the outlier cuts 
  ///  LTM for different padType

  if(fSector<36) padType=kAll;
  Double_t ddata[fNChannels];
  UInt_t nPoints = 0;
  UInt_t indexMin = 0;
  UInt_t indexMax = fNChannels;
  if(padType==kOROCmedium) indexMax=fkIndexes[64];
  else if(padType==kOROClong) indexMin=fkIndexes[64];
  for (UInt_t i=indexMin;i<indexMax;i++) {
     if (outlierROC && (outlierROC->GetValue(i) >1e-20)) continue;
     ddata[nPoints]= fData[i];
     nPoints++;
  }

  Double_t ltm =0, lsigma=0;
  if(nPoints>0) {
    Int_t hh = TMath::Min(TMath::Nint(fraction *nPoints), Int_t(nPoints));
    AliMathBase::EvaluateUni(nPoints,ddata, ltm, lsigma, hh);
    if (sigma) *sigma=lsigma;
  }
  
  return ltm;
}

//_____________________________________________________________________________
TH1F * AliTPCCalROC::MakeHisto1D(Float_t min, Float_t max,Int_t type){
  /// make 1D histo
  /// type -1 = user defined range
  ///       0 = nsigma cut nsigma=min
  ///       1 = delta cut around median delta=min

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
  TString name=Form("%s ROC 1D%d",GetTitle(),fSector);
  TH1F * his = new TH1F(name.Data(),name.Data(),100, min,max);
  for (UInt_t irow=0; irow<fNRows; irow++){
    UInt_t npads = (Int_t)GetNPads(irow);
    for (UInt_t ipad=0; ipad<npads; ipad++){
      his->Fill(GetValue(irow,ipad));
    }
  }
  return his;
}



TH2F * AliTPCCalROC::MakeHisto2D(Float_t min, Float_t max,Int_t type){
  /// make 2D histo
  /// type -1 = user defined range
  ///       0 = nsigma cut nsigma=min
  ///       1 = delta cut around median delta=min

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

  TString name=Form("%s ROC%d",GetTitle(),fSector);
  TH2F * his = new TH2F(name.Data(),name.Data(),fNRows+10,-5, fNRows+5, maxPad+10, -(Int_t(maxPad/2))-5, maxPad/2+5);
  for (UInt_t irow=0; irow<fNRows; irow++){
    UInt_t npads = (Int_t)GetNPads(irow);
    for (UInt_t ipad=0; ipad<npads; ipad++){
      his->Fill(irow+0.5,Int_t(ipad)-Int_t(npads/2)+0.5,GetValue(irow,ipad));
    }
  }
  his->SetMaximum(max);
  his->SetMinimum(min);
  return his;
}

TH2F * AliTPCCalROC::MakeHistoOutliers(Float_t delta, Float_t fraction, Int_t type){
  /// Make Histogram with outliers
  /// mode = 0 - sigma cut used
  /// mode = 1 - absolute cut used
  /// fraction - fraction of values used to define sigma
  /// delta - in mode 0 - nsigma cut
  ///            mode 1 - delta cut

  Double_t sigma;
  Float_t mean  = GetLTM(&sigma,fraction);  
  if (type==0) delta*=sigma; 
  UInt_t maxPad = 0;
  for (UInt_t irow=0; irow<fNRows; irow++){
    if (GetNPads(irow)>maxPad) maxPad = GetNPads(irow);
  }

  TString name=Form("%s ROC Outliers%d",GetTitle(),fSector);
  TH2F * his = new TH2F(name.Data(),name.Data(),fNRows+10,-5, fNRows+5, maxPad+10, -(Int_t(maxPad/2))-5, maxPad/2+5);
  for (UInt_t irow=0; irow<fNRows; irow++){
    UInt_t npads = (Int_t)GetNPads(irow);
    for (UInt_t ipad=0; ipad<npads; ipad++){
      if (TMath::Abs(GetValue(irow,ipad)-mean)>delta)
	his->Fill(irow+0.5,Int_t(ipad)-Int_t(npads/2)+0.5,1);
    }
  }
  return his;
}



void AliTPCCalROC::Draw(Option_t* opt){
  /// create histogram with values and draw it

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





void AliTPCCalROC::Test() {
  /// example function to show functionality and test AliTPCCalROC

  Float_t kEpsilon=0.00001;
  
  // create CalROC with defined values
  AliTPCCalROC  roc0(0);  
  for (UInt_t irow = 0; irow <roc0.GetNrows(); irow++){
    for (UInt_t ipad = 0; ipad <roc0.GetNPads(irow); ipad++){
      Float_t value  = irow+ipad/1000.;
      roc0.SetValue(irow,ipad,value);
    }
  }
  
  // copy CalROC, readout values and compare to original
  AliTPCCalROC roc1(roc0);
  for (UInt_t irow = 0; irow <roc1.GetNrows(); irow++){
    for (UInt_t ipad = 0; ipad <roc1.GetNPads(irow); ipad++){
      Float_t value  = irow+ipad/1000.;
      if (roc1.GetValue(irow,ipad)!=value){
        printf("Read/Write error\trow=%d\tpad=%d\n",irow,ipad);
      }
    }
  }

  // write original CalROC to file
  TFile f("calcTest.root","recreate");
  roc0.Write("Roc0");
  AliTPCCalROC * roc2 = (AliTPCCalROC*)f.Get("Roc0");
  f.Close();
  
  // read CalROC from file and compare to original CalROC
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

  // 
  // Algebra Tests
  //
  
  // Add constant
  AliTPCCalROC roc3(roc0);
  roc3.Add(1.5);
  for (UInt_t irow = 0; irow <roc3.GetNrows(); irow++){
    for (UInt_t ipad = 0; ipad <roc3.GetNPads(irow); ipad++){
      Float_t value  = irow+ipad/1000. + 1.5;
      if (TMath::Abs(roc3.GetValue(irow,ipad)-value) > kEpsilon){
        printf("Add constant - error\trow=%d\tpad=%d\n",irow,ipad);
      }
    }
  }

  // Add another CalROC
  AliTPCCalROC roc4(roc0);
  roc4.Add(&roc0, -1.5);
  for (UInt_t irow = 0; irow <roc4.GetNrows(); irow++){
    for (UInt_t ipad = 0; ipad <roc4.GetNPads(irow); ipad++){
      Float_t value  = irow+ipad/1000. - 1.5 * (irow+ipad/1000.);
      if (TMath::Abs(roc4.GetValue(irow,ipad)-value) > kEpsilon){
        printf("Add CalROC - error\trow=%d\tpad=%d\n",irow,ipad);
      }
    }
  }
  
  // Multiply with constant
  AliTPCCalROC roc5(roc0);
  roc5.Multiply(-1.4);
  for (UInt_t irow = 0; irow <roc5.GetNrows(); irow++){
    for (UInt_t ipad = 0; ipad <roc5.GetNPads(irow); ipad++){
      Float_t value  = (irow+ipad/1000.) * (-1.4);
      if (TMath::Abs(roc5.GetValue(irow,ipad)-value) > kEpsilon){
        printf("Multiply with constant - error\trow=%d\tpad=%d\n",irow,ipad);
      }
    }
  }

  // Multiply another CalROC
  AliTPCCalROC roc6(roc0);
  roc6.Multiply(&roc0);
  for (UInt_t irow = 0; irow <roc6.GetNrows(); irow++){
    for (UInt_t ipad = 0; ipad <roc6.GetNPads(irow); ipad++){
      Float_t value  = (irow+ipad/1000.) * (irow+ipad/1000.);
      if (TMath::Abs(roc6.GetValue(irow,ipad)-value) > kEpsilon*100){
        printf("Multiply with CalROC - error\trow=%d\tpad=%d\n",irow,ipad);
      }
    }
  }


  // Divide by CalROC
  AliTPCCalROC roc7(roc0);
  roc7.Divide(&roc0);
  for (UInt_t irow = 0; irow <roc7.GetNrows(); irow++){
    for (UInt_t ipad = 0; ipad <roc7.GetNPads(irow); ipad++){
      Float_t value  = 1.;
      if (irow+ipad == 0) value = roc0.GetValue(irow,ipad);
      if (TMath::Abs(roc7.GetValue(irow,ipad)-value) > kEpsilon){
        printf("Multiply with CalROC - error\trow=%d\tpad=%d\n",irow,ipad);
      }
    }
  }

  //
  // Statistics Test
  //
  
  // create CalROC with defined values
  TRandom3 rnd(0);
  AliTPCCalROC sroc0(0);
  for (UInt_t ichannel = 0; ichannel < sroc0.GetNchannels(); ichannel++){
    Float_t value  = rnd.Gaus(10., 2.);
    sroc0.SetValue(ichannel,value);
  }

  printf("Mean   (should be close to 10): %f\n", sroc0.GetMean());
  printf("RMS    (should be close to  2): %f\n", sroc0.GetRMS());
  printf("Median (should be close to 10): %f\n", sroc0.GetMedian());
  printf("LTM    (should be close to 10): %f\n", sroc0.GetLTM());

  //AliTPCCalROC* sroc1 = sroc0.LocalFit(4, 8);
  
  //delete sroc1;

//        std::cout << TMath::Abs(roc5.GetValue(irow,ipad)-value) << std::endl;
}


AliTPCCalROC * AliTPCCalROC::LocalFit(Int_t rowRadius, Int_t padRadius, AliTPCCalROC* ROCoutliers, Bool_t robust, Double_t chi2Threshold, Double_t robustFraction) {
  /// MakeLocalFit - smoothing
  /// returns a AliTPCCalROC with smoothed data
  /// rowRadius and padRadius specifies a window around a given pad. 
  /// The data of this window are fitted with a parabolic function. 
  /// This function is evaluated at the pad's position.
  /// At the edges the window is shifted, so that the specified pad is not anymore in the center of the window. 
  /// rowRadius  -  radius - rows to be used for smoothing
  /// padradius  -  radius - pads to be used for smoothing
  /// ROCoutlier -  map of outliers - pads not to be used for local smoothing
  /// robust     -  robust method of fitting  - (much slower)
  /// chi2Threshold: Threshold for chi2 when EvalRobust is called
  /// robustFraction: Fraction of data that will be used in EvalRobust

  AliTPCCalROC * xROCfitted = new AliTPCCalROC(fSector);
  TLinearFitter fitterQ(6,"hyp5");
   // TLinearFitter fitterQ(6,"x0++x1++x2++x3++x4++x5");
  fitterQ.StoreData(kTRUE);
  for (UInt_t row=0; row < GetNrows(); row++) {
    //std::cout << "Entering row " << row << " of " << GetNrows() << " @ sector "<< fSector << " for local fitting... "<< std::endl;
    for (UInt_t pad=0; pad < GetNPads(row); pad++)
      xROCfitted->SetValue(row, pad, GetNeighbourhoodValue(&fitterQ, row, pad, rowRadius, padRadius, ROCoutliers, robust, chi2Threshold, robustFraction));
  }
  return xROCfitted;
}


Double_t AliTPCCalROC::GetNeighbourhoodValue(TLinearFitter* fitterQ, Int_t row, Int_t pad, Int_t rRadius, Int_t pRadius, AliTPCCalROC *const ROCoutliers, Bool_t robust, Double_t chi2Threshold, Double_t robustFraction) {
  /// AliTPCCalROC::GetNeighbourhoodValue - smoothing - PRIVATE
  /// in this function the fit for LocalFit is done

  fitterQ->ClearPoints();
  TVectorD fitParam(6);
  Int_t    npoints=0;
  Double_t xx[6];
  Float_t dlx, dly;
  Float_t lPad[3] = {0};
  Float_t localXY[3] = {0};
  
  AliTPCROC* tpcROCinstance = AliTPCROC::Instance();
  tpcROCinstance->GetPositionLocal(fSector, row, pad, lPad);  // calculate position lPad by pad and row number
  
  TArrayI *neighbourhoodRows = 0;
  TArrayI *neighbourhoodPads = 0;
  
  //std::cerr << "Trying to get neighbourhood for row " << row << ", pad " << pad << std::endl;
  GetNeighbourhood(neighbourhoodRows, neighbourhoodPads, row, pad, rRadius, pRadius);
  //std::cerr << "Got neighbourhood for row " << row << ", pad " << pad << std::endl;
  
  Int_t r, p;
  for (Int_t i=0; i < (2*rRadius+1)*(2*pRadius+1); i++) {
    r = neighbourhoodRows->At(i);
    p = neighbourhoodPads->At(i);
    if (r == -1 || p == -1) continue;   // window is bigger than ROC
    tpcROCinstance->GetPositionLocal(fSector, r, p, localXY);   // calculate position localXY by pad and row number
    dlx = lPad[0] - localXY[0];
    dly = lPad[1] - localXY[1];
    //xx[0] = 1;
    xx[1] = dlx;
    xx[2] = dly;
    xx[3] = dlx*dlx;
    xx[4] = dly*dly;
    xx[5] = dlx*dly;
    if (!ROCoutliers || ROCoutliers->GetValue(r,p) != 1) {
      fitterQ->AddPoint(&xx[1], GetValue(r, p), 1);
      npoints++;
    }
  }
  
  delete neighbourhoodRows;
  delete neighbourhoodPads;
  
  if (npoints < 0.5 * ((2*rRadius+1)*(2*pRadius+1)) ) {
    // std::cerr << "Too few data points for fitting @ row " << row << ", pad " << pad << " in sector " << fSector << std::endl;
    return 0.;  // for diagnostic
  }
  fitterQ->Eval();
  fitterQ->GetParameters(fitParam);
  Float_t chi2Q = 0;
  if (robust) chi2Q = fitterQ->GetChisquare()/(npoints-6.);
  //if (robust) chi2Q = fitterQ->GetChisquare()/(npoints-6.);
  if (robust && chi2Q > chi2Threshold) {
    //std::cout << "robust fitter called... " << std::endl;
    fitterQ->EvalRobust(robustFraction);
    fitterQ->GetParameters(fitParam);
  }
  Double_t value = fitParam[0];
  
  //if (value < 0) std::cerr << "negative fit-value " << value << " in sector "<< this->fSector << " @ row: " << row << " and pad: " << pad << ", with fitter Chi2 = " << chi2Q <<  std::endl;
  return value;
}




void AliTPCCalROC::GetNeighbourhood(TArrayI* &rowArray, TArrayI* &padArray, Int_t row, Int_t pad, Int_t rRadius, Int_t pRadius) {
  /// AliTPCCalROC::GetNeighbourhood - PRIVATE
  /// in this function the window for LocalFit is determined

  rowArray = new TArrayI((2*rRadius+1)*(2*pRadius+1));
  padArray = new TArrayI((2*rRadius+1)*(2*pRadius+1));
  
  Int_t rmin = row - rRadius;
  UInt_t rmax = row + rRadius;
  
  // if window goes out of ROC
  if (rmin < 0) {
    rmax = rmax - rmin;
    rmin = 0;
  }
  if (rmax >= GetNrows()) {
    rmin = rmin - (rmax - GetNrows()+1);
    rmax = GetNrows() - 1;
      if (rmin  < 0 ) rmin = 0; // if the window is bigger than the ROC
  }
  
  Int_t pmin, pmax;
  Int_t i = 0;
  
  for (UInt_t r = rmin; r <= rmax; r++) {
    pmin = pad - pRadius;
    pmax = pad + pRadius;
    if (pmin < 0) {
      pmax = pmax - pmin;
      pmin = 0;
    }
    if (pmax >= (Int_t)GetNPads(r)) {
      pmin = pmin - (pmax - GetNPads(r)+1);
      pmax = GetNPads(r) - 1;
      if (pmin  < 0 ) pmin = 0; // if the window is bigger than the ROC
    }
    for (Int_t p = pmin; p <= pmax; p++) {
      (*rowArray)[i] = r;
      (*padArray)[i] = p;
      i++;
    }
  }
  for (Int_t j = i; j < rowArray->GetSize(); j++){  // unused padArray-entries, in the case that the window is bigger than the ROC
    //std::cout << "trying to write -1" << std::endl;
    (*rowArray)[j] = -1;
    (*padArray)[j] = -1;
    //std::cout << "writing -1" << std::endl;
  }
}



void AliTPCCalROC::GlobalFit(const AliTPCCalROC* ROCoutliers, Bool_t robust, TVectorD &fitParam, TMatrixD &covMatrix, Float_t & chi2, Int_t fitType, Double_t chi2Threshold, Double_t robustFraction, Double_t err, EPadType padType){
  /// Makes a  GlobalFit for the given sector and return fit-parameters, covariance and chi2
  /// update of GlobalFit in 2015 to have different gain in OROC medium and long pads
  /// The origin of the fit function is the center of the ROC!
  /// fitType == 0: fit plane function
  /// fitType == 1: fit parabolic function
  /// ROCoutliers - pads with value !=0 are not used in fitting procedure
  /// chi2Threshold: Threshold for chi2 when EvalRobust is called
  /// robustFraction: Fraction of data that will be used in EvalRobust
  /// err: error of the data points

  if(fSector<36) padType=kAll;

  TLinearFitter* fitterG = 0;
  Double_t xx[6];
  
  if (fitType  == 1) {
    fitterG = new TLinearFitter (6,"x0++x1++x2++x3++x4++x5");
    fitParam.ResizeTo(6);
    covMatrix.ResizeTo(6,6);
  } else {
    fitterG = new TLinearFitter(3,"x0++x1++x2");
    fitParam.ResizeTo(3);
    covMatrix.ResizeTo(3,3);
  }
  fitterG->StoreData(kTRUE);   
  fitterG->ClearPoints();
  Int_t    npoints=0;
  
  Float_t dlx, dly;
  Float_t centerPad[3] = {0};
  Float_t localXY[3] = {0};
  
  AliTPCROC* tpcROCinstance = AliTPCROC::Instance();
  if(padType==kAll)
    tpcROCinstance->GetPositionLocal(fSector, GetNrows()/2, GetNPads(GetNrows()/2)/2, centerPad);  // calculate center of ROC 
  else if(padType==kOROCmedium)  
    tpcROCinstance->GetPositionLocal(fSector, 64/2, GetNPads(64/2)/2, centerPad);  // calculate center of ROC medium pads
  else if(padType==kOROClong)
    tpcROCinstance->GetPositionLocal(fSector, (GetNrows()-64)/2+64, GetNPads((GetNrows()-64)/2+64)/2, centerPad);  // calculate center of ROC long pads

  UInt_t irowMin = 0;
  UInt_t irowMax = GetNrows();
  if(padType==kOROCmedium) irowMax = 64;
  else if(padType==kOROClong) irowMin = 64;

  // loop over all channels and read data into fitterG
  for (UInt_t irow = irowMin; irow < irowMax; irow++) {
    for (UInt_t ipad = 0; ipad < GetNPads(irow); ipad++) {
      // fill fitterG
      if (ROCoutliers && ROCoutliers->GetValue(irow, ipad) != 0) continue;
      tpcROCinstance->GetPositionLocal(fSector, irow, ipad, localXY);   // calculate position localXY by pad and row number
      dlx = localXY[0] - centerPad[0];
      dly = localXY[1] - centerPad[1];
      xx[0] = 1;
      xx[1] = dlx;
      xx[2] = dly;
      xx[3] = dlx*dlx;
      xx[4] = dly*dly;
      xx[5] = dlx*dly;
      npoints++;
      fitterG->AddPoint(xx, GetValue(irow, ipad), err);
    }
  }
  if(npoints>10) { // make sure there is something to fit
    fitterG->Eval();
    fitterG->GetParameters(fitParam);
    fitterG->GetCovarianceMatrix(covMatrix);
    if (fitType == 1)
      chi2 = fitterG->GetChisquare()/(npoints-6.);
    else chi2 = fitterG->GetChisquare()/(npoints-3.);
    if (robust && chi2 > chi2Threshold) {
      //    std::cout << "robust fitter called... " << std::endl;
      fitterG->EvalRobust(robustFraction);
      fitterG->GetParameters(fitParam);
    }
  } else {
    // set parameteres to 0
    Int_t nParameters = 3;
    if (fitType  == 1)
      nParameters = 6;

    for(Int_t i = 0; i < nParameters; i++)
      fitParam[i] = 0;
  }
  
  delete fitterG;
}

AliTPCCalROC* AliTPCCalROC::CreateGlobalFitCalROC(TVectorD &fitParam, Int_t sector, EPadType padType, AliTPCCalROC *oldTPCCalROC ){
  /// Create ROC with global fit parameters
  /// The origin of the fit function is the center of the ROC!
  /// loop over all channels, write fit values into new ROC and return it
  /// In 2015 different gain in medium and long pads 
  /// New CalROC is created for medium pads.
  /// For long pads oldTPCCalROC already contains values for medium pads.

  if(sector<36) padType=kAll;

  Float_t dlx, dly;
  Float_t centerPad[3] = {0};
  Float_t localXY[3] = {0};
  AliTPCCalROC * xROCfitted = NULL;
  if(oldTPCCalROC!=0 && padType!=kAll) xROCfitted = oldTPCCalROC;
  else xROCfitted = new AliTPCCalROC(sector);
  AliTPCROC* tpcROCinstance = AliTPCROC::Instance();
  if(padType==kAll)
    tpcROCinstance->GetPositionLocal(sector, xROCfitted->GetNrows()/2, xROCfitted->GetNPads(xROCfitted->GetNrows()/2)/2, centerPad);  // calculate center of ROC 
  else if(padType==kOROCmedium)
    tpcROCinstance->GetPositionLocal(sector, 64/2, xROCfitted->GetNPads(64/2)/2, centerPad);  // calculate center of ROC medium pads 
  else if(padType==kOROClong)
    tpcROCinstance->GetPositionLocal(sector, (xROCfitted->GetNrows()-64)/2+64, xROCfitted->GetNPads((xROCfitted->GetNrows()-64)/2+64)/2, centerPad);  // calculate center of ROC long pads 

  UInt_t irowMin = 0;
  UInt_t irowMax = xROCfitted->GetNrows();
  if(padType==kOROCmedium) irowMax = 64;
  else if(padType==kOROClong) irowMin = 64;

  Int_t fitType = 1;
  if (fitParam.GetNoElements() == 6) fitType = 1;
  else fitType = 0;
  Double_t value = 0;
  if (fitType == 1) { // parabolic fit
    for (UInt_t irow = irowMin; irow < irowMax; irow++) {
      for (UInt_t ipad = 0; ipad < xROCfitted->GetNPads(irow); ipad++) {
	tpcROCinstance->GetPositionLocal(sector, irow, ipad, localXY);   // calculate position localXY by pad and row number
	dlx = localXY[0] - centerPad[0];
	dly = localXY[1] - centerPad[1];
	value = fitParam[0] + fitParam[1]*dlx + fitParam[2]*dly + fitParam[3]*dlx*dlx + fitParam[4]*dly*dly + fitParam[5]*dlx*dly;
	xROCfitted->SetValue(irow, ipad, value);
      }
    }   
  }
  else {  // linear fit
    for (UInt_t irow = irowMin; irow < irowMax; irow++) {
      for (UInt_t ipad = 0; ipad < xROCfitted->GetNPads(irow); ipad++) {
	tpcROCinstance->GetPositionLocal(sector, irow, ipad, localXY);   // calculate position localXY by pad and row number
	dlx = localXY[0] - centerPad[0];
	dly = localXY[1] - centerPad[1];
	value = fitParam[0] + fitParam[1]*dlx + fitParam[2]*dly;
	xROCfitted->SetValue(irow, ipad, value);
      }
    }   
  }
  return xROCfitted;
}

