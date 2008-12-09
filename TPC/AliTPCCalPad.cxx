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
//  TPC calibration class for parameters which are saved per pad             //
//  Each AliTPCCalPad consists of 72 AliTPCCalROC-objects                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliTPCCalPad.h"
#include "AliTPCCalROC.h"
#include <TObjArray.h>
#include <TAxis.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TH2F.h>
#include "TTreeStream.h"
#include "TFile.h"
#include "TKey.h"
#include <iostream>

ClassImp(AliTPCCalPad)

//_____________________________________________________________________________
AliTPCCalPad::AliTPCCalPad():TNamed()
{
  //
  // AliTPCCalPad default constructor
  //

  for (Int_t isec = 0; isec < kNsec; isec++) {
    fROC[isec] = 0;
  }

}

//_____________________________________________________________________________
AliTPCCalPad::AliTPCCalPad(const Text_t *name, const Text_t *title)
                :TNamed(name,title)
{
  //
  // AliTPCCalPad constructor
  //
  for (Int_t isec = 0; isec < kNsec; isec++) {
    fROC[isec] = new AliTPCCalROC(isec);
  }
}


//_____________________________________________________________________________
AliTPCCalPad::AliTPCCalPad(const AliTPCCalPad &c):TNamed(c)
{
  //
  // AliTPCCalPad copy constructor
  //

  for (Int_t isec = 0; isec < kNsec; isec++) {
         fROC[isec] = 0;
     if (c.fROC[isec])
       fROC[isec] = new AliTPCCalROC(*(c.fROC[isec]));
  }
}

//_____________________________________________________________________________
AliTPCCalPad::AliTPCCalPad(TObjArray * array):TNamed()
{
  //
  // AliTPCCalPad default constructor
  //

  for (Int_t isec = 0; isec < kNsec; isec++) {
    fROC[isec] = (AliTPCCalROC *)array->At(isec);
  }

}


///_____________________________________________________________________________
AliTPCCalPad::~AliTPCCalPad()
{
  //
  // AliTPCCalPad destructor
  //

  for (Int_t isec = 0; isec < kNsec; isec++) {
    if (fROC[isec]) {
      delete fROC[isec];
      fROC[isec] = 0;
    }
  }

}

//_____________________________________________________________________________
AliTPCCalPad &AliTPCCalPad::operator=(const AliTPCCalPad &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliTPCCalPad &) c).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTPCCalPad::Copy(TObject &c) const
{
  //
  // Copy function
  //

  for (Int_t isec = 0; isec < kNsec; isec++) {
    if (fROC[isec]) {
      fROC[isec]->Copy(*((AliTPCCalPad &) c).fROC[isec]);
    }
  }
  TObject::Copy(c);
}


void AliTPCCalPad::SetCalROC(AliTPCCalROC* roc, Int_t sector){
   //
   // Set AliTPCCalROC copies values from 'roc'
   // if sector == -1 the sector specified in 'roc' is used
   // else sector specified in 'roc' is ignored and specified sector is filled
   //
   if (sector == -1) sector = roc->GetSector();
   if (!fROC[sector]) fROC[sector] = new AliTPCCalROC(sector);
   for (UInt_t ichannel = 0; ichannel < roc->GetNchannels(); ichannel++) 
      fROC[sector]->SetValue(ichannel, roc->GetValue(ichannel));
}



//_____________________________________________________________________________
void AliTPCCalPad::Add(Float_t c1)
{
    //
    // add constant c1 to all channels of all ROCs
    //

    for (Int_t isec = 0; isec < kNsec; isec++) {
	if (fROC[isec]){
	    fROC[isec]->Add(c1);
	}
    }
}

//_____________________________________________________________________________
void AliTPCCalPad::Multiply(Float_t c1)
{
  //
  // multiply each channel of all ROCs with c1
  //    
    for (Int_t isec = 0; isec < kNsec; isec++) {
	if (fROC[isec]){
	    fROC[isec]->Multiply(c1);
	}
    }
}

//_____________________________________________________________________________
void AliTPCCalPad::Add(const AliTPCCalPad * pad, Double_t c1)
{
  //
  // multiply AliTPCCalPad 'pad' by c1 and add each channel to the coresponing channel in all ROCs
  //  - pad by pad -
  //
    for (Int_t isec = 0; isec < kNsec; isec++) {
	if (fROC[isec] && pad->GetCalROC(isec)){
	    fROC[isec]->Add(pad->GetCalROC(isec),c1);
	}
    }
}

//_____________________________________________________________________________
void AliTPCCalPad::Multiply(const AliTPCCalPad * pad)
{
  //
  // multiply each channel of all ROCs with the coresponding channel of 'pad'
  //     - pad by pad -
  //
   for (Int_t isec = 0; isec < kNsec; isec++) {
	if (fROC[isec]){
	    fROC[isec]->Multiply(pad->GetCalROC(isec));
	}
    }
}

//_____________________________________________________________________________
void AliTPCCalPad::Divide(const AliTPCCalPad * pad)
{
  //
  // divide each channel of all ROCs by the coresponding channel of 'pad'
  //     - pad by pad -
  //    
    for (Int_t isec = 0; isec < kNsec; isec++) {
	if (fROC[isec]){
	    fROC[isec]->Divide(pad->GetCalROC(isec));
	}
    }
}

//_____________________________________________________________________________
TGraph  *  AliTPCCalPad::MakeGraph(Int_t type, Float_t ratio){
  //
  //   type=1 - mean
  //        2 - median
  //        3 - LTM
  //
  Int_t npoints = 0;
  for (Int_t i=0;i<72;i++) if (fROC[i]) npoints++;
  TGraph * graph = new TGraph(npoints);
  npoints=0;   
  for (Int_t isec=0;isec<72;isec++){
    if (!fROC[isec]) continue;
    if (type==0)  graph->SetPoint(npoints,isec,fROC[isec]->GetMean());      
    if (type==1)  graph->SetPoint(npoints,isec,fROC[isec]->GetMedian());
    if (type==2)  graph->SetPoint(npoints,isec,fROC[isec]->GetLTM(0,ratio));    
    npoints++;
  }

  graph->GetXaxis()->SetTitle("Sector"); 
  if (type==0) {
    graph->GetYaxis()->SetTitle("Mean");   
    graph->SetMarkerStyle(22);    
  }
  if (type==1) {
    graph->GetYaxis()->SetTitle("Median");   
    graph->SetMarkerStyle(22);    
  }
  if (type==2) {
      graph->GetYaxis()->SetTitle(Form("Mean%f",ratio));      
      graph->SetMarkerStyle(24);
  }

  return graph;
}

//_____________________________________________________________________________
Double_t AliTPCCalPad::GetMeanRMS(Double_t &rms)
{
    //
    // Calculates mean and RMS of all ROCs
    //
    Double_t sum = 0, sum2 = 0, n=0, val=0;
    for (Int_t isec = 0; isec < kNsec; isec++) {
        AliTPCCalROC *calRoc = fROC[isec];
	if ( calRoc ){
	    for (UInt_t irow=0; irow<calRoc->GetNrows(); irow++){
		for (UInt_t ipad=0; ipad<calRoc->GetNPads(irow); ipad++){
		    val = calRoc->GetValue(irow,ipad);
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
Double_t AliTPCCalPad::GetMean(AliTPCCalPad* outlierPad)
{
    //
    // return mean of the mean of all ROCs
    //
    Double_t arr[kNsec];
    Int_t n=0;
    for (Int_t isec = 0; isec < kNsec; isec++) {
       AliTPCCalROC *calRoc = fROC[isec];
       if ( calRoc ){
          AliTPCCalROC* outlierROC = 0;
          if (outlierPad) outlierROC = outlierPad->GetCalROC(isec);
	       arr[n] = calRoc->GetMean(outlierROC);
          n++;
       }
    }
    return TMath::Mean(n,arr);
}

//_____________________________________________________________________________
Double_t AliTPCCalPad::GetRMS(AliTPCCalPad* outlierPad)
{
    //
    // return mean of the RMS of all ROCs
    //
    Double_t arr[kNsec];
    Int_t n=0;
    for (Int_t isec = 0; isec < kNsec; isec++) {
       AliTPCCalROC *calRoc = fROC[isec];
       if ( calRoc ){
          AliTPCCalROC* outlierROC = 0;
          if (outlierPad) outlierROC = outlierPad->GetCalROC(isec);
          arr[n] = calRoc->GetRMS(outlierROC);
          n++;
       }
    }
    return TMath::Mean(n,arr);
}

//_____________________________________________________________________________
Double_t AliTPCCalPad::GetMedian(AliTPCCalPad* outlierPad)
{
    //
    // return mean of the median of all ROCs
    //
    Double_t arr[kNsec];
    Int_t n=0;
    for (Int_t isec = 0; isec < kNsec; isec++) {
       AliTPCCalROC *calRoc = fROC[isec];
       if ( calRoc ){
          AliTPCCalROC* outlierROC = 0;
          if (outlierPad) outlierROC = outlierPad->GetCalROC(isec);
          arr[n] = calRoc->GetMedian(outlierROC);
          n++;
       }
    }
    return TMath::Mean(n,arr);
}

//_____________________________________________________________________________
Double_t AliTPCCalPad::GetLTM(Double_t *sigma, Double_t fraction, AliTPCCalPad* outlierPad)
{
    //
    // return mean of the LTM and sigma of all ROCs
    //
    Double_t arrm[kNsec];
    Double_t arrs[kNsec];
    Double_t *sTemp=0x0;
    Int_t n=0;

    for (Int_t isec = 0; isec < kNsec; isec++) {
        AliTPCCalROC *calRoc = fROC[isec];
	if ( calRoc ){
	    if ( sigma ) sTemp=arrs+n;
       AliTPCCalROC* outlierROC = 0;
       if (outlierPad) outlierROC = outlierPad->GetCalROC(isec);
	    arrm[n] = calRoc->GetLTM(sTemp,fraction, outlierROC);
            n++;
	}
    }
    if ( sigma ) *sigma = TMath::Mean(n,arrs);
    return TMath::Mean(n,arrm);
}

//_____________________________________________________________________________
TH1F * AliTPCCalPad::MakeHisto1D(Float_t min, Float_t max,Int_t type){
  //
  // make 1D histo
  // type -1 = user defined range
  //       0 = nsigma cut nsigma=min
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
  sprintf(name,"%s Pad 1D",GetTitle());
  TH1F * his = new TH1F(name,name,100, min,max);
    for (Int_t isec = 0; isec < kNsec; isec++) {
	if (fROC[isec]){
	    for (UInt_t irow=0; irow<fROC[isec]->GetNrows(); irow++){
		UInt_t npads = (Int_t)fROC[isec]->GetNPads(irow);
		for (UInt_t ipad=0; ipad<npads; ipad++){
		    his->Fill(fROC[isec]->GetValue(irow,ipad));
		}
	    }
	}
    }
  return his;
}

//_____________________________________________________________________________
TH2F *AliTPCCalPad::MakeHisto2D(Int_t side){
  //
  // Make 2D graph
  // side  -  specify the side A = 0 C = 1
  // type  -  used types of determination of boundaries in z
  //
  Float_t kEpsilon = 0.000000000001;
  TH2F * his = new TH2F(GetName(), GetName(), 250,-250,250,250,-250,250);
  AliTPCROC * roc  = AliTPCROC::Instance(); 
  for (Int_t isec=0; isec<72; isec++){
    if (side==0 && isec%36>=18) continue;
    if (side>0 && isec%36<18) continue;
    if (fROC[isec]){
      AliTPCCalROC * calRoc = fROC[isec];
      for (UInt_t irow=0; irow<calRoc->GetNrows(); irow++)
	for (UInt_t ipad=0; ipad<calRoc->GetNPads(irow); ipad++)
	  if (TMath::Abs(calRoc->GetValue(irow,ipad))>kEpsilon){
	    Float_t xyz[3];
	    roc->GetPositionGlobal(isec,irow,ipad,xyz);
	    Int_t binx = 1+TMath::Nint((xyz[0]+250.)*0.5);
	    Int_t biny = 1+TMath::Nint((xyz[1]+250.)*0.5);
	    Float_t value = calRoc->GetValue(irow,ipad);	    
	    his->SetBinContent(binx,biny,value);
	  }
    }
  }
  his->SetXTitle("x (cm)");
  his->SetYTitle("y (cm)");
  return his;
}


AliTPCCalPad* AliTPCCalPad::LocalFit(const char* padName, Int_t rowRadius, Int_t padRadius, AliTPCCalPad* PadOutliers, Bool_t robust, Double_t chi2Threshold, Double_t robustFraction, Bool_t printCurrentSector) const {
   //
   // Loops over all AliTPCCalROCs and performs a localFit in each ROC
   // AliTPCCalPad with fit-data is returned
   // rowRadius and padRadius specifies a window around a given pad in one sector. 
   // The data of this window are fitted with a parabolic function. 
   // This function is evaluated at the pad's position.
   // At the edges the window is shifted, so that the specified pad is not anymore in the center of the window. 
   // rowRadius  -  radius - rows to be used for smoothing
   // padradius  -  radius - pads to be used for smoothing
   // ROCoutlier -  map of outliers - pads not to be used for local smoothing
   // robust     -  robust method of fitting  - (much slower)
   // chi2Threshold: Threshold for chi2 when EvalRobust is called
   // robustFraction: Fraction of data that will be used in EvalRobust
   //
   //
   AliTPCCalPad* pad = new AliTPCCalPad(padName, padName);
   for (Int_t isec = 0; isec < 72; isec++){
      if (printCurrentSector) std::cout << "LocalFit in sector " << isec << "\r" << std::flush;
      if (PadOutliers)
         pad->SetCalROC(GetCalROC(isec)->LocalFit(rowRadius, padRadius, PadOutliers->GetCalROC(isec), robust, chi2Threshold, robustFraction));
      else 
         pad->SetCalROC(GetCalROC(isec)->LocalFit(rowRadius, padRadius, 0, robust, chi2Threshold, robustFraction));
   }
   return pad;
}


AliTPCCalPad* AliTPCCalPad::GlobalFit(const char* padName, AliTPCCalPad* PadOutliers, Bool_t robust, Int_t fitType, Double_t chi2Threshold, Double_t robustFraction, Double_t err, TObjArray *fitParArr, TObjArray *fitCovArr){
   //
   // Loops over all AliTPCCalROCs and performs a globalFit in each ROC
   // AliTPCCalPad with fit-data is returned
   // chi2Threshold: Threshold for chi2 when EvalRobust is called
   // robustFraction: Fraction of data that will be used in EvalRobust
   // chi2Threshold: Threshold for chi2 when EvalRobust is called
   // robustFraction: Fraction of data that will be used in EvalRobust
   // err: error of the data points
   // if fitParArr and/or fitCovArr is given, write fitParameters and/or covariance Matrices into the array
   //
   AliTPCCalPad* pad = new AliTPCCalPad(padName, padName);
   TVectorD fitParam(0);
   TMatrixD covMatrix(0,0);
   Float_t chi2 = 0;
   for (Int_t isec = 0; isec < 72; isec++){
      if (PadOutliers)
         GetCalROC(isec)->GlobalFit(PadOutliers->GetCalROC(isec), robust, fitParam, covMatrix, chi2, fitType, chi2Threshold, robustFraction, err);
      else 
         GetCalROC(isec)->GlobalFit(0, robust, fitParam, covMatrix, chi2, fitType, chi2Threshold, robustFraction, err);

      AliTPCCalROC *roc=AliTPCCalROC::CreateGlobalFitCalROC(fitParam, isec);
      pad->SetCalROC(roc);
      delete roc;
      if ( fitParArr ) fitParArr->AddAtAndExpand(new TVectorD(fitParam), isec);
      if ( fitCovArr ) fitCovArr->AddAtAndExpand(new TMatrixD(covMatrix), isec);
   }
   return pad;
}


void AliTPCCalPad::GlobalSidesFit(const AliTPCCalPad* PadOutliers, TVectorD &fitParamSideA, TVectorD &fitParamSideC,TMatrixD &covMatrixSideA, TMatrixD &covMatrixSideC, Float_t & chi2SideA, Float_t & chi2SideC, Int_t fitType, Bool_t robust, Double_t chi2Threshold, Double_t robustFraction){
  //
  // Makes a  GlobalFit over each side and return fit-parameters, covariance and chi2 for each side
  // fitType == 0: fit plane function
  // fitType == 1: fit parabolic function
  // PadOutliers - pads with value !=0 are not used in fitting procedure
  // chi2Threshold: Threshold for chi2 when EvalRobust is called
  // robustFraction: Fraction of data that will be used in EvalRobust
  //
  TLinearFitter* fitterGA = 0;
  TLinearFitter* fitterGC = 0;
  
  if (fitType  == 1) {
    fitterGA = new TLinearFitter (6,"x0++x1++x2++x3++x4++x5");
    fitterGC = new TLinearFitter (6,"x0++x1++x2++x3++x4++x5");
  }
  else {
    fitterGA = new TLinearFitter(3,"x0++x1++x2");
    fitterGC = new TLinearFitter(3,"x0++x1++x2");
  }
  fitterGA->StoreData(kTRUE);   
  fitterGC->StoreData(kTRUE);   
  fitterGA->ClearPoints();
  fitterGC->ClearPoints();
  Double_t xx[6];  
  Int_t    npointsA=0;
  Int_t    npointsC=0;
  
  Float_t localXY[3] = {0}; // pad's position, needed to get the pad's position
  Float_t lx, ly;  // pads position
  
  AliTPCROC* tpcROCinstance = AliTPCROC::Instance();  // to calculate the pad's position
  
  // loop over all sectors and pads and read data into fitterGA and fitterGC 
  if (fitType == 1) {  
  // parabolic fit
    fitParamSideA.ResizeTo(6);
    fitParamSideC.ResizeTo(6);
    covMatrixSideA.ResizeTo(6,6);
    covMatrixSideC.ResizeTo(6,6);
    for (UInt_t isec = 0; isec<72; isec++){
      for (UInt_t irow = 0; irow < GetCalROC(isec)->GetNrows(); irow++) {
         for (UInt_t ipad = 0; ipad < GetCalROC(isec)->GetNPads(irow); ipad++) {
            // fill fitterG
            tpcROCinstance->GetPositionLocal(isec, irow, ipad, localXY);   // calculate position localXY by sector, pad and row number
            lx = localXY[0];
            ly = localXY[1];
            xx[0] = 1;
            xx[1] = lx;
            xx[2] = ly;
            xx[3] = lx*lx;
            xx[4] = ly*ly;
            xx[5] = lx*ly;
            if (!PadOutliers || PadOutliers->GetCalROC(isec)->GetValue(irow, ipad) != 1) {
            // if given pad is no outlier, add it to TLinearFitter, decide to which of both
//                sector  0 - 17: IROC, A
//                sector 18 - 35: IROC, C
//                sector 36 - 53: OROC, A
//                sector 54 - 71: CROC, C
               if (isec <= 17 || (isec >= 36 && isec <= 53)) { // Side A
                  npointsA++;
                  fitterGA->AddPoint(xx, GetCalROC(isec)->GetValue(irow, ipad), 1);  
               }
               else { // side C
                  npointsC++;
                  fitterGC->AddPoint(xx, GetCalROC(isec)->GetValue(irow, ipad), 1);  
               }
            }
         }
      }
    }
  }
  else {   
  // linear fit
    fitParamSideA.ResizeTo(3);
    fitParamSideC.ResizeTo(3);
    covMatrixSideA.ResizeTo(3,3);
    covMatrixSideC.ResizeTo(3,3);
    
    for (UInt_t isec = 0; isec<72; isec++){
      for (UInt_t irow = 0; irow < GetCalROC(isec)->GetNrows(); irow++) {
         for (UInt_t ipad = 0; ipad < GetCalROC(isec)->GetNPads(irow); ipad++) {
            // fill fitterG
            tpcROCinstance->GetPositionLocal(isec, irow, ipad, localXY);   // calculate position localXY by sector, pad and row number
            lx = localXY[0];
            ly = localXY[1];
            xx[0] = 1;
            xx[1] = lx;
            xx[2] = ly;
            if (!PadOutliers || PadOutliers->GetCalROC(isec)->GetValue(irow, ipad) != 1) {
            // if given pad is no outlier, add it to TLinearFitter, decide to which of both
//                sector  0 - 17: IROC, A
//                sector 18 - 35: IROC, C
//                sector 36 - 53: OROC, A
//                sector 54 - 71: CROC, C
               if (isec <= 17 || (isec >= 36 && isec <= 53)) { 
               // Side A
                  npointsA++;
                  fitterGA->AddPoint(xx, GetCalROC(isec)->GetValue(irow, ipad), 1);  
               }
               else { 
               // side C
                  npointsC++;
                  fitterGC->AddPoint(xx, GetCalROC(isec)->GetValue(irow, ipad), 1);  
               }
            }
         }
      }
    }
  }    
  
  fitterGA->Eval();
  fitterGC->Eval();
  fitterGA->GetParameters(fitParamSideA);
  fitterGC->GetParameters(fitParamSideC);
  fitterGA->GetCovarianceMatrix(covMatrixSideA);
  fitterGC->GetCovarianceMatrix(covMatrixSideC);
  if (fitType == 1){
    chi2SideA = fitterGA->GetChisquare()/(npointsA-6.);
    chi2SideC = fitterGC->GetChisquare()/(npointsC-6.);
  }
  else {
   chi2SideA = fitterGA->GetChisquare()/(npointsA-3.);
   chi2SideC = fitterGC->GetChisquare()/(npointsC-3.);
  }
  if (robust && chi2SideA > chi2Threshold) {
    //    std::cout << "robust fitter called... " << std::endl;
    fitterGA->EvalRobust(robustFraction);
    fitterGA->GetParameters(fitParamSideA);
  }
  if (robust && chi2SideC > chi2Threshold) {
    //    std::cout << "robust fitter called... " << std::endl;
    fitterGC->EvalRobust(robustFraction);
    fitterGC->GetParameters(fitParamSideC);
  }
  delete fitterGA;
  delete fitterGC;
}


