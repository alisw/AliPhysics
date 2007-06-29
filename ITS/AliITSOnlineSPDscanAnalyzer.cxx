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

/* $Id$*/

////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                 //
// This class is used in the detector algorithm framework //
// to process the data stored in special container files  //
// (see AliITSOnlineSPDscan). For instance, minimum       //
// threshold values can be calculated.                    //
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

Double_t itsSpdErrorf(Double_t *x, Double_t *par){
  if (par[2]<0) par[2]=0;
  Double_t val = par[2]+(0.12*256*32-par[2])*(0.5+0.5*TMath::Erf((x[0]-par[0])/par[1]/sqrt(2.)));
  return val;
}
//Double_t itsSpdErrorfOrig(Double_t *x, Double_t *par){
//  return 0.5+0.5*TMath::Erf((x[0]-par[0])/par[1]/sqrt(2.));
//}


AliITSOnlineSPDscanAnalyzer::AliITSOnlineSPDscanAnalyzer(Char_t *fileName) :
  fType(99),fDacId(99),fScanObj(NULL),fTriggers(NULL),
  fOverWrite(kFALSE),fNoiseThreshold(0.01),fNoiseMinimumEvents(100),
  fMinNrStepsBeforeIncrease(5),fMinIncreaseFromBaseLine(2),fStepDownDacSafe(2),fMaxBaseLineLevel(10)
{
  // constructor
  sprintf(fFileName,"%s",fileName);
  for (UInt_t chipNr=0; chipNr<11; chipNr++) {
    for (UInt_t hs=0; hs<6; hs++) {
      fMeanMultiplicity[hs][chipNr]=NULL;
      fHitEventEfficiency[hs][chipNr]=NULL;
    }
  }
  for (Int_t module=0; module<240; module++) {
    fHandler[module]=NULL;
  }

  Init();
}

AliITSOnlineSPDscanAnalyzer::AliITSOnlineSPDscanAnalyzer(const AliITSOnlineSPDscanAnalyzer& handle) :
  fType(99),fDacId(99),fScanObj(NULL),fTriggers(NULL),
  fOverWrite(kFALSE),fNoiseThreshold(0.01),fNoiseMinimumEvents(100),
  fMinNrStepsBeforeIncrease(5),fMinIncreaseFromBaseLine(2),fStepDownDacSafe(2),fMaxBaseLineLevel(10)
{
  // copy constructor, only copies the filename (not the processed data)
  sprintf(fFileName,"%s",handle.fFileName);

  fScanObj=NULL;
  fType=99;
  fDacId=99;
  for (UInt_t chipNr=0; chipNr<11; chipNr++) {
    for (UInt_t hs=0; hs<6; hs++) {
      fMeanMultiplicity[hs][chipNr]=NULL;
      fHitEventEfficiency[hs][chipNr]=NULL;
    }
  }
  fTriggers=NULL;
  for (Int_t module=0; module<240; module++) {
    fHandler[module]=NULL;
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
  for (Int_t module=0; module<240; module++) {
    if (fHandler[module]!=NULL) {
      delete fHandler[module];
    }
  }
}

AliITSOnlineSPDscanAnalyzer& AliITSOnlineSPDscanAnalyzer::operator=(const AliITSOnlineSPDscanAnalyzer& handle) {
  // assignment operator, only copies the filename (not the processed data)
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
    for (Int_t module=0; module<240; module++) {
      if (fHandler[module]!=NULL) {
	delete fHandler[module];
      }
    }
   
    sprintf(fFileName,"%s",handle.fFileName);

    fScanObj=NULL;
    fType=99;
    fDacId=99;
    for (UInt_t chipNr=0; chipNr<11; chipNr++) {
      for (UInt_t hs=0; hs<6; hs++) {
	fMeanMultiplicity[hs][chipNr]=NULL;
	fHitEventEfficiency[hs][chipNr]=NULL;
      }
    }
    fTriggers=NULL;
    for (Int_t module=0; module<240; module++) {
      fHandler[module]=NULL;
    }
    
    Init();    
  }
  return *this;
}

void AliITSOnlineSPDscanAnalyzer::Init() {
  // first checks type of container and then initializes container obj
  FILE* fp0 = fopen(fFileName, "r");
  if (fp0 == NULL) {
    return;
  }
  else {
    fclose(fp0);
  }
  fScanObj = new AliITSOnlineSPDscan(fFileName);
  fType = fScanObj->GetType();
  delete fScanObj;

  // init container
  switch(fType) {
  case kUNIMA:
  case kNOISE:
    fScanObj = new AliITSOnlineSPDscanSingle(fFileName);
    break;
  case kMINTH:
  case kDAC:
  case kDELAY:
    fScanObj = new AliITSOnlineSPDscanMultiple(fFileName);
    fDacId = ((AliITSOnlineSPDscanMultiple*)fScanObj)->GetDacId();
    break;
  case kMEANTH:
    fScanObj = new AliITSOnlineSPDscanMeanTh(fFileName);
    fDacId = ((AliITSOnlineSPDscanMeanTh*)fScanObj)->GetDacId();
    break;
  default:
    printf("Type %d not defined!\n",fType);
    fScanObj=NULL;
    return;
    break;
  }

  // set some default values (these should later be read from text file)
  fOverWrite=kFALSE;
  fNoiseThreshold=0.01;
  fNoiseMinimumEvents=100;
  fMinNrStepsBeforeIncrease=6;
  fMinIncreaseFromBaseLine=2;
  fStepDownDacSafe=2;
  fMaxBaseLineLevel=10;

}


Bool_t AliITSOnlineSPDscanAnalyzer::ProcessDeadPixels(Char_t *oldcalibDir) {
  // process dead pixel data
  if (fScanObj==NULL) {
    printf("No data!\n");
    return kFALSE;
  }
  // should be type kUNIMA
  if (fType!=kUNIMA) {
    printf("Dead pixels only for scan type %d\n",kUNIMA);
    return kFALSE;
  }

  Int_t nrDead[240];
  for (Int_t i=0; i<240; i++) {
    nrDead[i]=0;
  }
  UInt_t routerNr = fScanObj->GetRouterNr();
  UInt_t rowStart = fScanObj->GetRowStart();
  UInt_t rowEnd   = fScanObj->GetRowEnd();
  for (UInt_t hs=0; hs<6; hs++) {
    for (UInt_t chipNr=0; chipNr<10; chipNr++) {
      if (fScanObj->GetChipPresent(hs,chipNr) && fScanObj->GetAverageMultiplicity(0,hs,chipNr)>0) { // check the status of the chippresent parameter in the mood header!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	UInt_t module = AliITSRawStreamSPD::GetModuleNumber(routerNr,hs*2+chipNr/5);
	if (!fHandler[module]) {
	  fHandler[module] = new AliITSOnlineCalibrationSPDhandler(module);
	}
	fHandler[module]->SetFileLocation(oldcalibDir);
	fHandler[module]->ReadFromFile();
	if (fOverWrite) {fHandler[module]->ResetDead();}
	for (UInt_t col=0; col<32; col++) {
	  for (UInt_t row=rowStart; row<=rowEnd; row++) {
	    if (col!=1 && col!=9 && col!=17 && col!=25) { //exclude test columns!!!
	      if (fScanObj->GetHits(0,hs,chipNr,col,row)==0) {
		Int_t newCol = 32*(chipNr%5) + col;
		if (fHandler[module]->SetDeadPixel(newCol,row)) {
		  nrDead[module]++;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  return kTRUE;
}



Bool_t AliITSOnlineSPDscanAnalyzer::ProcessNoisyPixels(Char_t *oldcalibDir) {
  // process noisy pixel data
  if (fScanObj==NULL) {
    printf("No data!\n");
    return kFALSE;
  }
  // should be type kNOISE
  if (fType != kNOISE) {
    printf("Noisy pixels only for scan type %d\n",kNOISE);
    return kFALSE;
  }

  Int_t nrNoisy[240];
  for (Int_t i=0; i<240; i++) {
    nrNoisy[i]=0;
  }
  if (fScanObj->GetTriggers(0)<fNoiseMinimumEvents) {
    printf("Error: Process noisy: Too few events.\n"); 
    return kFALSE;
  }
  UInt_t routerNr = fScanObj->GetRouterNr();
  for (UInt_t hs=0; hs<6; hs++) {
    for (UInt_t chipNr=0; chipNr<10; chipNr++) {
      if (fScanObj->GetChipPresent(hs,chipNr)) { // check the status of the chippresent parameter in the mood header!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	UInt_t module = AliITSRawStreamSPD::GetModuleNumber(routerNr,hs*2+chipNr/5);
	for (UInt_t col=0; col<32; col++) {
	  for (UInt_t row=0; row<256; row++) {
	    if (fScanObj->GetHitsEfficiency(0,hs,chipNr,col,row)>fNoiseThreshold) {
	      if (!fHandler[module]) {
		fHandler[module] = new AliITSOnlineCalibrationSPDhandler(module);
	      }
	      fHandler[module]->SetFileLocation(oldcalibDir);
	      fHandler[module]->ReadFromFile();
	      if (fOverWrite) {fHandler[module]->ResetNoisy();}
	      Int_t newCol = 32*(chipNr%5) + col;
	      if (fHandler[module]->SetNoisyPixel(newCol,row)) {
		nrNoisy[module]++;
	      }
	    }
	  }
	}
      }
    }
  }
  return kTRUE;
}

Bool_t AliITSOnlineSPDscanAnalyzer::SaveDeadNoisyPixels(UInt_t module, Char_t *calibDir) {
  // save dead and noisy pixels to file in dir calibDir
  if (fHandler[module]!=NULL) {
    fHandler[module]->SetFileLocation(calibDir);
    fHandler[module]->WriteToFile();
    return kTRUE;
  }
  return kFALSE;
}




Int_t AliITSOnlineSPDscanAnalyzer::GetDelay(UInt_t hs, UInt_t chipNr) {
  // get delay
  if (hs>=6 || chipNr>10) return -1;
  if (fScanObj==NULL) {
    printf("No data!\n");
    return -1;
  }
  // should be type kDELAY or kDAC with id 42 (delay_ctrl)
  if (fType!=kDELAY && (fType!=kDAC || fDacId!=42)) {
    printf("Delay only for scan type %d or %d and dac_id 42\n",kDELAY,kDAC);
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


// ********** old version of "GetDelay", which fits a gaussian to the meanmult graph
//Int_t AliITSOnlineSPDscanAnalyzer::GetDelay(UInt_t hs, UInt_t chipNr) {
//  // get delay
//  if (hs>=6 || chipNr>10) return -1;
//  if (fScanObj==NULL) {
//    printf("No data!\n");
//    return -1;
//  }
//  // should be type kDELAY or kDAC with id 42 (delay_ctrl)
//  if (fType!=kDELAY && (fType!=kDAC || fDacId!=42)) {
//    printf("Delay only for scan type %d or %d and dac_id 42\n",kDELAY,kDAC);
//    return -1;
//  }
//  if (fMeanMultiplicity[hs][chipNr]==NULL) {
//    if (!ProcessMeanMultiplicity()) {
//      return -1;
//    }
//  }
//
//  Char_t funcName[30];
//  sprintf(funcName,"Fit delay func HS%d CHIP%d",hs,chipNr);
//  Int_t minDac = ((AliITSOnlineSPDscanMultiple*)fScanObj)->GetDacValue(0);
//  Int_t maxDac = ((AliITSOnlineSPDscanMultiple*)fScanObj)->GetDacValue(fScanObj->GetNSteps()-1);
//  TF1* delayFunc = new TF1(funcName,"gaus",minDac,maxDac);
//  fMeanMultiplicity[hs][chipNr]->Fit(funcName,"Q0");
//  Double_t mean = delayFunc->GetParameter(1);
//  //  Double_t sigma = fDelayFunc[hs][chipNr]->GetParameter(2);
//  delete delayFunc;
//  if (mean>minDac && mean<maxDac) {
//    return (Int_t) (mean+0.5);
//  }
//  else {
//    return -1;
//  }
//}

Int_t AliITSOnlineSPDscanAnalyzer::GetNrNoisyUnima(UInt_t hs, UInt_t chipNr) {
  // in case of a uniformity scan, returns the nr of noisy pixels, (here > 200 hits)
  if (hs>=6 || chipNr>10) return -1;
  if (fScanObj==NULL) {
    printf("No data!\n");
    return kFALSE;
  }
  // should be type kUNIMA
  if (fType != kUNIMA) {
    printf("Noisy pixels Unima only for scan type %d\n",kUNIMA);
    return kFALSE;
  }
  if (fScanObj->GetTriggers(0)!=25600) {
    printf("Error: Process noisy unima: Incorrect number of events (!=25600.\n"); 
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
    printf("No data!\n");
    return -1;
  }
  // should be type  kMINTH  or  kDAC with id 39 (pre_vth)
  if (fType!=kMINTH && (fType!=kDAC || fDacId!=39)) {
    printf("MinTh only for scan type %d OR %d with dac_id 39\n",kMINTH,kDAC);
    return -1;
  }
  if (fMeanMultiplicity[hs][chipNr]==NULL) {
    if (!ProcessMeanMultiplicity()) {
      return -1;
    }
  }

  Int_t lastDac = FindLastMinThDac(hs,chipNr);
  if (lastDac==-1) {
    printf("GetMinTh: HS%d, Chip%d: Increase of Mean Multiplicity by %1.2f never reached.\n",hs,chipNr,fMinIncreaseFromBaseLine);
    return -1;
  }

  Int_t minDac = ((AliITSOnlineSPDscanMultiple*)fScanObj)->GetDacValue(0);
  Char_t funcName[30];
  sprintf(funcName,"Fit minth func HS%d CHIP%d",hs,chipNr);
  TF1 *minThFunc = new TF1(funcName,itsSpdErrorf,100,500,3);
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
    printf("GetMinTh: HS%d, Chip%d: BaseLine too large (%1.2f>%1.2f)\n",hs,chipNr,baseLine,fMaxBaseLineLevel);
    return -1;
  }
  UInt_t step=FindClosestLowerStep(lastDac);
  Float_t compareLine = GetCompareLine(step,hs,chipNr,baseLine);
  if (compareLine==-1) {
    printf("GetMinTh: HS%d, Chip%d: Not enough steps (%d<%d) before increase to get a compare line.\n",hs,chipNr,step,fMinNrStepsBeforeIncrease);
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
    printf("GetMinTh: HS%d, Chip%d: Did not find a point below the compare line (%f).\n",hs,chipNr,compareLine);
    return -1;
  }
}



Bool_t AliITSOnlineSPDscanAnalyzer::ProcessMeanMultiplicity() {
  // process mean multiplicity data
  if (fScanObj==NULL) {
    printf("No data!\n");
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
    printf("No data!\n");
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
    printf("No data!\n");
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

AliITSOnlineCalibrationSPDhandler* AliITSOnlineSPDscanAnalyzer::GetOnlineCalibrationHandler(UInt_t module) {
  // returns a pointer to the AliITSOnlineCalibrationSPDhandler
  if (module<240) return fHandler[module]; 
  else return NULL;
}

UInt_t AliITSOnlineSPDscanAnalyzer::GetRouterNr() {
  // returns the router nr of scan obj
  if (fScanObj!=NULL) return fScanObj->GetRouterNr(); 
  else return 99;
}

TH2F* AliITSOnlineSPDscanAnalyzer::GetHitMapTot(UInt_t step) {
  // creates and returns a pointer to a hitmap histo (half sector style a la spdmood)
  if (fScanObj==NULL) {
    printf("No data!\n");
    return NULL;
  }
  Char_t histoname[50];
  if (fType==kMINTH || fType==kMEANTH || fType==kDAC) {
    sprintf(histoname,"Router %d , DAC %d",GetRouterNr(),((AliITSOnlineSPDscanMultiple*)fScanObj)->GetDacValue(step));
  }
  else {
    sprintf(histoname,"Router %d ",GetRouterNr());
  }
  TH2F* fHitMapTot = new TH2F(histoname,histoname,32*10,-0.5,32*10-0.5,256*6,-0.5,256*6-0.5);
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
