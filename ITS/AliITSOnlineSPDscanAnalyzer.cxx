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

/* $Id$ */

////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                 //
// This class is used in the detector algorithm framework //
// to process the data stored in special container files  //
// (see AliITSOnlineSPDscan). For instance, minimum       //
// threshold values can be extracted.                     //
////////////////////////////////////////////////////////////

#include "AliITSOnlineSPDscanAnalyzer.h"
#include "AliITSOnlineSPDscan.h"
#include "AliITSOnlineSPDscanSingle.h"
#include "AliITSOnlineSPDscanMultiple.h"
#include "AliITSOnlineSPDscanMeanTh.h"
#include "AliITSOnlineCalibrationSPDhandler.h"
#include "AliITSRawStreamSPD.h"
#include <TStyle.h>
#include <TMath.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH2F.h>
#include <TError.h>
#include <iostream>
#include <fstream>

Double_t itsSpdErrorf(Double_t *x, Double_t *par){
  if (par[2]<0) par[2]=0;
  Double_t val = par[2]+(0.12*256*32-par[2])*(0.5+0.5*TMath::Erf((x[0]-par[0])/par[1]/sqrt(2.)));
  return val;
}
//Double_t itsSpdErrorfOrig(Double_t *x, Double_t *par){
//  return 0.5+0.5*TMath::Erf((x[0]-par[0])/par[1]/sqrt(2.));
//}


AliITSOnlineSPDscanAnalyzer::AliITSOnlineSPDscanAnalyzer(const Char_t *fileName, AliITSOnlineCalibrationSPDhandler *handler) :
  fType(99),fDacId(99),fFileName(fileName),fScanObj(NULL),fHandler(handler),fTriggers(NULL),
  fOverWrite(kFALSE),fNoiseThreshold(0.01),fNoiseMinimumEvents(100),
  fMinNrStepsBeforeIncrease(5),fMinIncreaseFromBaseLine(2),fStepDownDacSafe(5),fMaxBaseLineLevel(10)
{
  // constructor
  for (UInt_t chipNr=0; chipNr<11; chipNr++) {
    for (UInt_t hs=0; hs<6; hs++) {
      fMeanMultiplicity[hs][chipNr]=NULL;
      fHitEventEfficiency[hs][chipNr]=NULL;
    }
  }
  for (UInt_t mod=0; mod<240; mod++) {
    fbModuleScanned[mod]=kFALSE;
  }

  Init();
}

AliITSOnlineSPDscanAnalyzer::AliITSOnlineSPDscanAnalyzer(const AliITSOnlineSPDscanAnalyzer& handle) :
  fType(99),fDacId(99),fFileName("."),fScanObj(NULL),fHandler(NULL),fTriggers(NULL),
  fOverWrite(kFALSE),fNoiseThreshold(0.01),fNoiseMinimumEvents(100),
  fMinNrStepsBeforeIncrease(5),fMinIncreaseFromBaseLine(2),fStepDownDacSafe(5),fMaxBaseLineLevel(10)
{
  // copy constructor, only copies the filename and params (not the processed data)
  fFileName=handle.fFileName;
  fOverWrite=handle.fOverWrite;
  fNoiseThreshold=handle.fNoiseThreshold;
  fNoiseMinimumEvents=handle.fNoiseMinimumEvents;
  fMinNrStepsBeforeIncrease=handle.fMinNrStepsBeforeIncrease;
  fMinIncreaseFromBaseLine=handle.fMinIncreaseFromBaseLine;
  fStepDownDacSafe=handle.fStepDownDacSafe;
  fMaxBaseLineLevel=handle.fMaxBaseLineLevel;

  for (UInt_t chipNr=0; chipNr<11; chipNr++) {
    for (UInt_t hs=0; hs<6; hs++) {
      fMeanMultiplicity[hs][chipNr]=NULL;
      fHitEventEfficiency[hs][chipNr]=NULL;
    }
  }
  for (UInt_t mod=0; mod<240; mod++) {
    fbModuleScanned[mod]=kFALSE;
  }

  Init();
}

AliITSOnlineSPDscanAnalyzer::~AliITSOnlineSPDscanAnalyzer() {
  // destructor
  for (UInt_t hs=0; hs<6; hs++) {
    for (UInt_t chipNr=0; chipNr<11; chipNr++) {
      if (fMeanMultiplicity[hs][chipNr]!=NULL) {
	delete fMeanMultiplicity[hs][chipNr];
      }
      if (fHitEventEfficiency[hs][chipNr]!=NULL) {
	delete fHitEventEfficiency[hs][chipNr];
      }
    }
  }
  if (fTriggers!=NULL) delete fTriggers;
  if (fScanObj!=NULL) delete fScanObj;
}

AliITSOnlineSPDscanAnalyzer& AliITSOnlineSPDscanAnalyzer::operator=(const AliITSOnlineSPDscanAnalyzer& handle) {
  // assignment operator, only copies the filename and params (not the processed data)
  if (this!=&handle) {
    for (UInt_t hs=0; hs<6; hs++) {
      for (UInt_t chipNr=0; chipNr<11; chipNr++) {
	if (fMeanMultiplicity[hs][chipNr]!=NULL) {
	  delete fMeanMultiplicity[hs][chipNr];
	}
	if (fHitEventEfficiency[hs][chipNr]!=NULL) {
	  delete fHitEventEfficiency[hs][chipNr];
	}
      }
    }
    if (fTriggers!=NULL) delete fTriggers;
    if (fScanObj!=NULL) delete fScanObj;
   
    fFileName=handle.fFileName;
    fOverWrite=handle.fOverWrite;
    fNoiseThreshold=handle.fNoiseThreshold;
    fNoiseMinimumEvents=handle.fNoiseMinimumEvents;
    fMinNrStepsBeforeIncrease=handle.fMinNrStepsBeforeIncrease;
    fMinIncreaseFromBaseLine=handle.fMinIncreaseFromBaseLine;
    fStepDownDacSafe=handle.fStepDownDacSafe;
    fMaxBaseLineLevel=handle.fMaxBaseLineLevel;

    for (UInt_t chipNr=0; chipNr<11; chipNr++) {
      for (UInt_t hs=0; hs<6; hs++) {
	fMeanMultiplicity[hs][chipNr]=NULL;
	fHitEventEfficiency[hs][chipNr]=NULL;
      }
    }
    for (UInt_t mod=0; mod<240; mod++) {
      fbModuleScanned[mod]=kFALSE;
    }
    fTriggers=NULL;
    fHandler=NULL;
    
    fScanObj=NULL;
    fType=99;
    fDacId=99;

    Init();    
  }
  return *this;
}

void AliITSOnlineSPDscanAnalyzer::Init() {
  // first checks type of container and then initializes container obj
  FILE* fp0 = fopen(fFileName.Data(), "r");
  if (fp0 == NULL) {
    return;
  }
  else {
    fclose(fp0);
  }
  fScanObj = new AliITSOnlineSPDscan(fFileName.Data());
  fType = fScanObj->GetType();
  delete fScanObj;

  // init container
  switch(fType) {
  case kUNIMA:
  case kNOISE:
    fScanObj = new AliITSOnlineSPDscanSingle(fFileName.Data());
    break;
  case kMINTH:
  case kDAC:
  case kDELAY:
    fScanObj = new AliITSOnlineSPDscanMultiple(fFileName.Data());
    fDacId = ((AliITSOnlineSPDscanMultiple*)fScanObj)->GetDacId();
    break;
  case kMEANTH:
    fScanObj = new AliITSOnlineSPDscanMeanTh(fFileName.Data());
    fDacId = ((AliITSOnlineSPDscanMeanTh*)fScanObj)->GetDacId();
    break;
  default:
    Error("AliITSOnlineSPDscanAnalyzer::Init","Type %d not defined!",fType);
    fScanObj=NULL;
    return;
    break;
  }

}

void AliITSOnlineSPDscanAnalyzer::SetParam(const Char_t *pname, const Char_t *pval) {
  // set a parameter
  TString name = pname;
  TString val = pval;
  if (name.CompareTo("fOverWrite")==0) {
    if (val.CompareTo("YES")==0 || val.CompareTo("1")==0) {
      fOverWrite = kTRUE;
    }
    else fOverWrite = kFALSE;
  }
  else if (name.CompareTo("fNoiseThreshold")==0) {
    fNoiseThreshold = val.Atof();
  }
  else if (name.CompareTo("fNoiseMinimumEvents")==0) {
    fNoiseMinimumEvents = val.Atoi();
  }
  else if (name.CompareTo("fMinNrStepsBeforeIncrease")==0) {
    fMinNrStepsBeforeIncrease = val.Atoi();
  }
  else if (name.CompareTo("fMinIncreaseFromBaseLine")==0) {
    fMinIncreaseFromBaseLine = val.Atof();
  }
  else if (name.CompareTo("fStepDownDacSafe")==0) {
    fStepDownDacSafe = val.Atoi();
  }
  else if (name.CompareTo("fMaxBaseLineLevel")==0) {
    fMaxBaseLineLevel = val.Atof();
  }
  else {
    Error("AliITSOnlineSPDscanAnalyzer::SetParam","Parameter %s in configuration file unknown.",name.Data());
  }
}

void AliITSOnlineSPDscanAnalyzer::ReadParamsFromLocation(const Char_t *dirName) {
  // opens file (default name) in dir dirName and reads parameters from it
  TString paramsFileName = Form("%s/standal_params.txt",dirName);
  ifstream paramsFile;
  paramsFile.open(paramsFileName, ifstream::in);
  if (paramsFile.fail()) {
    printf("No config file (%s) present. Using default tuning parameters.\n",paramsFileName.Data());
  }
  else {
    while(1) {
      Char_t paramN[50];
      Char_t paramV[50];
      paramsFile >> paramN;
      if (paramsFile.eof()) break;
      paramsFile >> paramV;
      SetParam(paramN,paramV);
      if (paramsFile.eof()) break;
    }
    paramsFile.close();
  }
}

Bool_t AliITSOnlineSPDscanAnalyzer::IsChipPresent(UInt_t hs, UInt_t chipNr) {
  // is the chip present?
  if (fScanObj==NULL) {
    Warning("AliITSOnlineSPDscanAnalyzer::IsChipPresent","No data!");
    return kFALSE;
  }
  return fScanObj->GetChipPresent(hs,chipNr);
}

Bool_t AliITSOnlineSPDscanAnalyzer::ProcessDeadPixels(/*Char_t *oldcalibDir*/) {
  // process dead pixel data, for uniformity scan, 
  // NB: This will not be the general way of finding dead pixels.
  if (fScanObj==NULL) {
    Warning("AliITSOnlineSPDscanAnalyzer::ProcessDeadPixels","No data!");
    return kFALSE;
  }
  // should be type kUNIMA
  if (fType!=kUNIMA) {
    Warning("AliITSOnlineSPDscanAnalyzer::ProcessDeadPixels","Dead pixels only for scan type %d.",kUNIMA);
    return kFALSE;
  }
  // handler should be initialized
  if (fHandler==NULL) {
    Error("AliITSOnlineSPDscanAnalyzer::ProcessDeadPixels","Calibration handler is not initialized!");
    return kFALSE;
  }

  UInt_t routerNr = fScanObj->GetRouterNr();
  UInt_t rowStart = fScanObj->GetRowStart();
  UInt_t rowEnd   = fScanObj->GetRowEnd();
  for (UInt_t hs=0; hs<6; hs++) {
    for (UInt_t chipNr=0; chipNr<10; chipNr++) {
      if (fScanObj->GetChipPresent(hs,chipNr) && fScanObj->GetAverageMultiplicity(0,hs,chipNr)>0) { // check the status of the chippresent parameter in the mood header!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (fOverWrite) {fHandler->ResetDeadForChip(routerNr,hs,chipNr);}
	for (UInt_t col=0; col<32; col++) {
	  for (UInt_t row=rowStart; row<=rowEnd; row++) {
	    if (col!=1 && col!=9 && col!=17 && col!=25) { //exclude test columns!!!
	      if (fScanObj->GetHits(0,hs,chipNr,col,row)==0) {
		fHandler->SetDeadPixel(routerNr,hs,chipNr,col,row);
	      }
	    }
	  }
	}
      }
    }
  }
  return kTRUE;
}



Bool_t AliITSOnlineSPDscanAnalyzer::ProcessNoisyPixels(/*Char_t *oldcalibDir*/) {
  // process noisy pixel data
  if (fScanObj==NULL) {
    Warning("AliITSOnlineSPDscanAnalyzer::ProcessNoisyPixels","No data!");
    return kFALSE;
  }
  // should be type kNOISE
  if (fType != kNOISE) {
    Warning("AliITSOnlineSPDscanAnalyzer::ProcessNoisyPixels","Noisy pixels only for scan type %d.",kNOISE);
    return kFALSE;
  }
  // handler should be initialized
  if (fHandler==NULL) {
    Error("AliITSOnlineSPDscanAnalyzer::ProcessNoisyPixels","Calibration handler is not initialized!");
    return kFALSE;
  }
  // check if enough statistics
  if (fScanObj->GetTriggers(0)<fNoiseMinimumEvents) {
    Warning("AliITSOnlineSPDscanAnalyzer::ProcessNoisyPixels","Process noisy: Too few events.");
    return kFALSE;
  }

  UInt_t routerNr = fScanObj->GetRouterNr();
  for (UInt_t hs=0; hs<6; hs++) {
    for (UInt_t chipNr=0; chipNr<10; chipNr++) {
      if (fScanObj->GetChipPresent(hs,chipNr)) { // check the status of the chippresent parameter in the mood header!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (fOverWrite) {fHandler->ResetNoisyForChip(routerNr,hs,chipNr);}
	for (UInt_t col=0; col<32; col++) {
	  for (UInt_t row=0; row<256; row++) {
	    if (fScanObj->GetHitsEfficiency(0,hs,chipNr,col,row)>fNoiseThreshold) {
	      fHandler->SetNoisyPixel(routerNr,hs,chipNr,col,row);
	    }
	  }
	}
      }
    }
  }
  return kTRUE;
}

Int_t AliITSOnlineSPDscanAnalyzer::GetDelay(UInt_t hs, UInt_t chipNr) {
  // get delay
  if (hs>=6 || chipNr>10) return -1;
  if (fScanObj==NULL) {
    Warning("AliITSOnlineSPDscanAnalyzer::GetDelay","No data!");
    return -1;
  }
  // should be type kDELAY or kDAC with id 42 (delay_ctrl)
  if (fType!=kDELAY && (fType!=kDAC || fDacId!=42)) {
    Warning("AliITSOnlineSPDscanAnalyzer::GetDelay","Delay only for scan type %d or %d and dac_id 42.",kDELAY,kDAC);
    return -1;
  }
  if (fMeanMultiplicity[hs][chipNr]==NULL) {
    if (!ProcessMeanMultiplicity()) {
      return -1;
    }
  }

  UInt_t maxStep=0;
  Float_t maxVal=0;
  for (UInt_t step=0; step<fScanObj->GetNSteps(); step++) {
    Double_t thisDac;
    Double_t thisMult;
    fMeanMultiplicity[hs][chipNr]->GetPoint(step,thisDac,thisMult);
    if (thisMult > maxVal) {
      maxVal = thisMult;
      maxStep = step;
    }
  }

  if (maxVal>0) {
    return ((AliITSOnlineSPDscanMultiple*)fScanObj)->GetDacValue(maxStep);
  }
  else {
    return -1;
  }

}


Int_t AliITSOnlineSPDscanAnalyzer::GetNrNoisyUnima(UInt_t hs, UInt_t chipNr) {
  // in case of a uniformity scan, returns the nr of noisy pixels, (here > 200 hits)
  if (hs>=6 || chipNr>10) return -1;
  if (fScanObj==NULL) {
    Error("AliITSOnlineSPDscanAnalyzer::GetNrNoisyUnima","No data!");
    return kFALSE;
  }
  // should be type kUNIMA
  if (fType != kUNIMA) {
    Error("AliITSOnlineSPDscanAnalyzer::GetNrNoisyUnima","Noisy pixels Unima only for scan type %d.",kUNIMA);
    return kFALSE;
  }
  if (fScanObj->GetTriggers(0)!=25600) {
    Error("AliITSOnlineSPDscanAnalyzer::GetNrNoisyUnima","Process noisy unima: Incorrect number of events (!=25600.");
    return kFALSE;
  }

  Int_t nrNoisy=0;
  if (fScanObj->GetChipPresent(hs,chipNr)) { // check the status of the chippresent parameter in the mood header!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for (UInt_t col=0; col<32; col++) {
      for (UInt_t row=0; row<256; row++) {
	if (fScanObj->GetHits(0,hs,chipNr,col,row)>200) {
	  nrNoisy++;
	}
      }
    }
  }
  else {
    return -1;
  }
  return nrNoisy;
}

Int_t AliITSOnlineSPDscanAnalyzer::FindLastMinThDac(UInt_t hs, UInt_t chipNr) {
  // returns dac value where fMinIncreaseFromBaseLine reached
  if (hs>=6 || chipNr>10) return -1;
  if (fMeanMultiplicity[hs][chipNr]==NULL) {
    if (!ProcessMeanMultiplicity()) {
      return -1;
    }
  }
  Double_t firstVal, dummy1;
  fMeanMultiplicity[hs][chipNr]->GetPoint(0,dummy1,firstVal);
  UInt_t step=0;
  while (step<fScanObj->GetNSteps()-1) {
    Double_t graphVal, dummy2;
    fMeanMultiplicity[hs][chipNr]->GetPoint(step+1,dummy2,graphVal);
    if (graphVal>firstVal+fMinIncreaseFromBaseLine) break;
    step++;
  }
  if (step==fScanObj->GetNSteps()-1) return -1;
  return ((AliITSOnlineSPDscanMultiple*)fScanObj)->GetDacValue(step);
}

Int_t AliITSOnlineSPDscanAnalyzer::FindClosestLowerStep(Float_t dacValueInput) {
  // returns step closest (lower) to a dacvalue 
  UInt_t step=0;
  while (step<fScanObj->GetNSteps()-1) {
    Int_t dacVal = ((AliITSOnlineSPDscanMultiple*)fScanObj)->GetDacValue(step+1);
    if (dacVal>=dacValueInput) break;
    step++;
  }
  return step;
}

Float_t AliITSOnlineSPDscanAnalyzer::GetCompareLine(UInt_t step, UInt_t hs, UInt_t chipNr, Float_t basePar2) {
  // returns value to compare mean mult with (when finding min th)
  if (hs>=6 || chipNr>10) return -1;
  if (step<fMinNrStepsBeforeIncrease) return -1;
  Float_t baseLine = basePar2;
  if (baseLine<0) baseLine=0;
  Float_t baseAdd;
  Double_t baseM=0;
  Double_t baseS=0;
  Double_t d,m;
  for (UInt_t st=1;st<2*step/3;st++) { // skip first point...
    fMeanMultiplicity[hs][chipNr]->GetPoint(st,d,m);
    baseM+=m-baseLine;
    baseS+=(m-baseLine)*(m-baseLine);
  }
  baseAdd=2*sqrt( baseS/(2*step/3-1) - (baseM/(2*step/3-1))*(baseM/(2*step/3-1)) );
  baseAdd+=0.03; // magic number
  if (baseAdd>fMinIncreaseFromBaseLine) baseAdd=fMinIncreaseFromBaseLine;
  return baseLine + baseAdd;
}

Int_t AliITSOnlineSPDscanAnalyzer::GetMinTh(UInt_t hs, UInt_t chipNr) {
  // calculates and returns the minimum threshold
  if (hs>=6 || chipNr>10) return -1;
  if (fScanObj==NULL) {
    Error("AliITSOnlineSPDscanAnalyzer::GetMinTh","No data!");
    return -1;
  }
  // should be type  kMINTH  or  kDAC with id 39 (pre_vth)
  if (fType!=kMINTH && (fType!=kDAC || fDacId!=39)) {
    Error("AliITSOnlineSPDscanAnalyzer::GetMinTh","MinTh only for scan type %d OR %d with dac_id 39.",kMINTH,kDAC);
    return -1;
  }
  if (fMeanMultiplicity[hs][chipNr]==NULL) {
    if (!ProcessMeanMultiplicity()) {
      return -1;
    }
  }

  Int_t lastDac = FindLastMinThDac(hs,chipNr);
  if (lastDac==-1) {
    Warning("AliITSOnlineSPDscanAnalyzer::GetMinTh","HS%d, Chip%d: Increase of Mean Multiplicity by %1.2f never reached.",hs,chipNr,fMinIncreaseFromBaseLine);
    return -1;
  }

  Int_t minDac = ((AliITSOnlineSPDscanMultiple*)fScanObj)->GetDacValue(0);
  TString funcName = Form("Fit minth func HS%d CHIP%d",hs,chipNr);
  TF1 *minThFunc = new TF1(funcName.Data(),itsSpdErrorf,100,500,3);
  minThFunc->SetParameter(0,lastDac+10);
  minThFunc->SetParameter(1,2);
  minThFunc->SetParameter(2,0);
  minThFunc->SetParName(0,"Mean");
  minThFunc->SetParName(1,"Sigma");
  minThFunc->SetParName(2,"BaseLine");
  minThFunc->SetLineWidth(1);
  if (fMeanMultiplicity[hs][chipNr]==NULL) {
    if (!ProcessMeanMultiplicity()) {
      return -1;
    }
  }
  fMeanMultiplicity[hs][chipNr]->Fit(funcName,"Q0","",minDac,lastDac);

  //  Double_t mean = fMinThFunc[hs][chipNr]->GetParameter(0);
  //  Double_t sigma = fMinThFunc[hs][chipNr]->GetParameter(1);
  Double_t baseLine = minThFunc->GetParameter(2);
  delete minThFunc;

  if (baseLine>fMaxBaseLineLevel) {
    Warning("AliITSOnlineSPDscanAnalyzer::GetMinTh","HS%d, Chip%d: BaseLine too large (%1.2f>%1.2f).",hs,chipNr,baseLine,fMaxBaseLineLevel);
    return -1;
  }
  UInt_t step=FindClosestLowerStep(lastDac);
  Float_t compareLine = GetCompareLine(step,hs,chipNr,baseLine);
  if (compareLine==-1) {
    Warning("AliITSOnlineSPDscanAnalyzer::GetMinTh","HS%d, Chip%d: Not enough steps (%d<%d) before increase to get a compare line.",hs,chipNr,step,fMinNrStepsBeforeIncrease);
    return -1;
  }

  Double_t mult, dummy;
  mult=1000;
  while (mult > compareLine && step>0) {
    fMeanMultiplicity[hs][chipNr]->GetPoint(step,dummy,mult);
    step--;
  }
  Int_t minth = ((AliITSOnlineSPDscanMultiple*)fScanObj)->GetDacValue(step+1)-fStepDownDacSafe;

  if (step>0) {
    return minth;
  }
  else {
    Warning("AliITSOnlineSPDscanAnalyzer::GetMinTh","HS%d, Chip%d: Did not find a point below the compare line (%f).",hs,chipNr,compareLine);
    return -1;
  }
}



Bool_t AliITSOnlineSPDscanAnalyzer::ProcessMeanMultiplicity() {
  // process mean multiplicity data
  if (fScanObj==NULL) {
    Error("AliITSOnlineSPDscanAnalyzer::ProcessMeanMultiplicity","No data!");
    return kFALSE;
  }
  for (UInt_t step=0; step<fScanObj->GetNSteps(); step++) {
    for (UInt_t hs=0; hs<6; hs++) {
      for (UInt_t chipNr=0; chipNr<11; chipNr++) {
	//	if (fScanObj->GetChipPresent(hs,chipNr)) { // check the status of the chippresent parameter in the mood header!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (step==0) {
	  if (fMeanMultiplicity[hs][chipNr]!=NULL) {
	    delete fMeanMultiplicity[hs][chipNr];
	  }
	  fMeanMultiplicity[hs][chipNr] = new TGraph();
	}
	Float_t multiplMean=fScanObj->GetAverageMultiplicity(step,hs,chipNr);
	if (fType==kMINTH || fType==kMEANTH || fType==kDAC || fType==kDELAY) {
	  fMeanMultiplicity[hs][chipNr]->SetPoint(step,((AliITSOnlineSPDscanMultiple*)fScanObj)->GetDacValue(step),multiplMean);
	}
	else {
	  fMeanMultiplicity[hs][chipNr]->SetPoint(step,0,multiplMean);
	}
      }
      //      }
    }
  }
  return kTRUE;
}

TGraph* AliITSOnlineSPDscanAnalyzer::GetMeanMultiplicityG(UInt_t hs, UInt_t chipNr) {
  // returns mean multiplicity graph
  if (hs>=6 || chipNr>10) return NULL;
  if (fMeanMultiplicity[hs][chipNr]==NULL) {
    if (!ProcessMeanMultiplicity()) {
      return NULL;
    }
  }
  return fMeanMultiplicity[hs][chipNr];
}

Bool_t AliITSOnlineSPDscanAnalyzer::ProcessHitEventEfficiency() {
  // process hit event efficiency
  if (fScanObj==NULL) {
    Error("AliITSOnlineSPDscanAnalyzer::ProcessHitEventEfficiency","No data!");
    return kFALSE;
  }
  for (UInt_t step=0; step<fScanObj->GetNSteps(); step++) {
    for (UInt_t hs=0; hs<6; hs++) {
      for (UInt_t chipNr=0; chipNr<11; chipNr++) {
	//	if (fScanObj->GetChipPresent(hs,chipNr)) { // check the status of the chippresent parameter in the mood header!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if (step==0) {
	  if (fHitEventEfficiency[hs][chipNr]!=NULL) {
	    delete fHitEventEfficiency[hs][chipNr];
	  }
	  fHitEventEfficiency[hs][chipNr] = new TGraph();
	}
	Float_t efficiency=fScanObj->GetHitEventsEfficiency(step,hs,chipNr);
	if (fType==kMINTH || fType==kMEANTH || fType==kDAC || fType==kDELAY) {
	  fHitEventEfficiency[hs][chipNr]->SetPoint(step,((AliITSOnlineSPDscanMultiple*)fScanObj)->GetDacValue(step),efficiency);
	}
	else {
	  fHitEventEfficiency[hs][chipNr]->SetPoint(step,0,efficiency);
	}
      }
      //      }
    }
  }
  return kTRUE;
}

TGraph* AliITSOnlineSPDscanAnalyzer::GetHitEventEfficiencyG(UInt_t hs, UInt_t chipNr) {
  // returns hit event efficiency graph
  if (hs>=6 || chipNr>10) return NULL;
  if (fHitEventEfficiency[hs][chipNr]==NULL) {
    if (!ProcessHitEventEfficiency()) {
      return NULL;
    }
  }
  return fHitEventEfficiency[hs][chipNr];
}


Bool_t AliITSOnlineSPDscanAnalyzer::ProcessNrTriggers() {
  // process nr of triggers data
  if (fScanObj==NULL) {
    Error("AliITSOnlineSPDscanAnalyzer::ProcessNrTriggers","No data!");
    return kFALSE;
  }
  for (UInt_t step=0; step<fScanObj->GetNSteps(); step++) {
    if (step==0) {
      if (fTriggers!=NULL) {
	delete fTriggers;
      }
      fTriggers = new TGraph();
    }
    if (fType==kMINTH || fType==kMEANTH || fType==kDAC || fType==kDELAY) {
      fTriggers->SetPoint(step,((AliITSOnlineSPDscanMultiple*)fScanObj)->GetDacValue(step),fScanObj->GetTriggers(step));
    }
    else {
      fTriggers->SetPoint(step,0,fScanObj->GetTriggers(step));
    }
  }
  return kTRUE;
}

TGraph* AliITSOnlineSPDscanAnalyzer::GetNrTriggersG() {
  // returns nr of triggers graph
  if (fTriggers==NULL) {
    if (!ProcessNrTriggers()) {
      return NULL;
    }
  }
  return fTriggers;
}

Bool_t AliITSOnlineSPDscanAnalyzer::GetHalfStavePresent(UInt_t hs) {
  // returns half stave present info
  if (hs<6 && fScanObj!=NULL) {
    Int_t chipstatus=0;
    for (Int_t chip=0; chip<10; chip++) {
      chipstatus+=fScanObj->GetChipPresent(hs,chip);
    }
    if (chipstatus>0) return kTRUE;
  }
  return kFALSE;
}

UInt_t AliITSOnlineSPDscanAnalyzer::GetRouterNr() {
  // returns the router nr of scan obj
  if (fScanObj!=NULL) return fScanObj->GetRouterNr(); 
  else return 99;
}

TH2F* AliITSOnlineSPDscanAnalyzer::GetHitMapTot(UInt_t step) {
  // creates and returns a pointer to a hitmap histo (half sector style a la spdmood)
  if (fScanObj==NULL) {
    Error("AliITSOnlineSPDscanAnalyzer::GetHitMapTot","No data!");
    return NULL;
  }
  TString histoname;
  if (fType==kMINTH || fType==kMEANTH || fType==kDAC) {
    histoname = Form("Router %d , DAC %d",GetRouterNr(),((AliITSOnlineSPDscanMultiple*)fScanObj)->GetDacValue(step));
  }
  else {
    histoname = Form("Router %d ",GetRouterNr());
  }
  TH2F* fHitMapTot = new TH2F(histoname.Data(),histoname.Data(),32*10,-0.5,32*10-0.5,256*6,-0.5,256*6-0.5);
  fHitMapTot->SetNdivisions(-10,"X");
  fHitMapTot->SetNdivisions(-006,"Y");
  fHitMapTot->SetTickLength(0,"X");
  fHitMapTot->SetTickLength(0,"Y");
  fHitMapTot->GetXaxis()->SetLabelColor(gStyle->GetCanvasColor());
  fHitMapTot->GetYaxis()->SetLabelColor(gStyle->GetCanvasColor());
  for (UInt_t hs=0; hs<6; hs++) {
    for (UInt_t chipNr=0; chipNr<10; chipNr++) {
      for (UInt_t col=0; col<32; col++) {
	for (UInt_t row=0; row<256; row++) {
	  fHitMapTot->Fill(chipNr*32+col,(5-hs)*256+row,fScanObj->GetHits(step,hs,chipNr,col,row));
	}
      }
    }
  }
  return fHitMapTot;
}

TH2F* AliITSOnlineSPDscanAnalyzer::GetHitMapChip(UInt_t step, UInt_t hs, UInt_t chip) {
  // creates and returns a pointer to a hitmap histo (chip style a la spdmood)
  if (fScanObj==NULL) {
    Error("AliITSOnlineSPDscanAnalyzer::GetHitMapChip","No data!");
    return NULL;
  }

  TString histoName;
  TString histoTitle;
  histoName = Form("fChipHisto_%d_%d_%d", GetRouterNr(), hs, chip);
  histoTitle = Form("Eq ID %d, Half Stave %d, Chip %d", GetRouterNr(), hs, chip);

  TH2F *returnHisto = new TH2F(histoName.Data(), histoTitle.Data(), 32, -0.5, 31.5, 256, -0.5, 255.5);
  returnHisto->SetMinimum(0);
  for (UInt_t col=0; col<32; col++) {
    for (UInt_t row=0; row<256; row++) {
      returnHisto->Fill(col,row,fScanObj->GetHits(step,hs,chip,col,row));
    }
  }

  return returnHisto;
}
