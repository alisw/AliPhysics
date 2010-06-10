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

////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                 //
// This class is used in the detector algorithm framework //
// to process the data stored in special container files  //
// (see AliITSOnlineSPDphys).                             //
////////////////////////////////////////////////////////////

#include "AliITSOnlineSPDphysAnalyzer.h"
#include "AliITSOnlineSPDphys.h"
#include "AliITSOnlineCalibrationSPDhandler.h"
#include "AliITSRawStreamSPD.h"
#include "AliITSIntMap.h"
#include <TStyle.h>
#include <TMath.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH2F.h>
#include <TError.h>
#include <iostream>
#include <fstream>

AliITSOnlineSPDphysAnalyzer::AliITSOnlineSPDphysAnalyzer(const Char_t *fileName, AliITSOnlineCalibrationSPDhandler* handler, Bool_t readFromGridFile) :
  fFileName(fileName),fPhysObj(NULL),fHandler(handler),
  fNrEnoughStatChips(0),fNrDeadChips(0),fNrInefficientChips(0),
  fNrEqHits(0),fbDeadProcessed(kFALSE),
  fThreshNoisy(1e-9),fThreshDead(1e-9),
  fMinEventsForNoisy(10000),fMinEventsForDead(10000),
  fDefinitelyNoisyRatio(0.3),
  fMinNrEqHitsForDeadChips(60000),fRatioToMeanForInefficientChip(0.1)
{
  // constructor  
  Init(readFromGridFile);
}

AliITSOnlineSPDphysAnalyzer::AliITSOnlineSPDphysAnalyzer(AliITSOnlineSPDphys* physObj, AliITSOnlineCalibrationSPDhandler* handler) :
  fFileName("test.root"),fPhysObj(NULL),fHandler(handler),
  fNrEnoughStatChips(0),fNrDeadChips(0),fNrInefficientChips(0),
  fNrEqHits(0),fbDeadProcessed(kFALSE),
  fThreshNoisy(1e-9),fThreshDead(1e-9),
  fMinEventsForNoisy(10000),fMinEventsForDead(10000),
  fDefinitelyNoisyRatio(0.3),
  fMinNrEqHitsForDeadChips(60000),fRatioToMeanForInefficientChip(0.1)
{
  // alt constructor
  fPhysObj = physObj;
}

AliITSOnlineSPDphysAnalyzer::AliITSOnlineSPDphysAnalyzer(const AliITSOnlineSPDphysAnalyzer& handle) :
  fFileName("test.root"),fPhysObj(NULL),fHandler(NULL),
  fNrEnoughStatChips(0),fNrDeadChips(0),fNrInefficientChips(0),
  fNrEqHits(0),fbDeadProcessed(kFALSE),
  fThreshNoisy(1e-9),fThreshDead(1e-9),
  fMinEventsForNoisy(10000),fMinEventsForDead(10000),
  fDefinitelyNoisyRatio(0.3),
  fMinNrEqHitsForDeadChips(60000),fRatioToMeanForInefficientChip(0.1)
{
  // copy constructor, only copies the filename and params (not the processed data)
  fFileName=handle.fFileName;
  fThreshNoisy = handle.fThreshNoisy;
  fThreshDead = handle.fThreshDead;
  fMinEventsForNoisy = handle.fMinEventsForNoisy;
  fMinEventsForDead = handle.fMinEventsForDead;
  fDefinitelyNoisyRatio = handle.fDefinitelyNoisyRatio;
  fMinNrEqHitsForDeadChips = handle.fMinNrEqHitsForDeadChips;
  fRatioToMeanForInefficientChip = handle.fRatioToMeanForInefficientChip;

  Init();
}

AliITSOnlineSPDphysAnalyzer::~AliITSOnlineSPDphysAnalyzer() {
  // destructor
  if (fPhysObj!=NULL) delete fPhysObj;
}

AliITSOnlineSPDphysAnalyzer& AliITSOnlineSPDphysAnalyzer::operator=(const AliITSOnlineSPDphysAnalyzer& handle) {
  // assignment operator, only copies the filename and params (not the processed data)
  if (this!=&handle) {
    if (fPhysObj!=NULL) delete fPhysObj;
    
    fFileName=handle.fFileName;
    fThreshNoisy = handle.fThreshNoisy;
    fThreshDead = handle.fThreshDead;
    fMinEventsForNoisy = handle.fMinEventsForNoisy;
    fMinEventsForDead = handle.fMinEventsForDead;
    fDefinitelyNoisyRatio = handle.fDefinitelyNoisyRatio;
    fMinNrEqHitsForDeadChips = handle.fMinNrEqHitsForDeadChips;
    fRatioToMeanForInefficientChip = handle.fRatioToMeanForInefficientChip;

    fPhysObj = NULL;
    fHandler = NULL;
    fNrEnoughStatChips = 0;
    fNrDeadChips = 0;
    fNrInefficientChips = 0;
    fNrEqHits = 0;
    fbDeadProcessed = kFALSE;

    Init();    
  }
  return *this;
}

void AliITSOnlineSPDphysAnalyzer::Init(Bool_t readFromGridFile) {
  // initialize container obj
  if (!readFromGridFile) {
    FILE* fp0 = fopen(fFileName.Data(), "r");
    if (fp0 == NULL) {
      return;
    }
    else {
      fclose(fp0);
    }
  }
  fPhysObj = new AliITSOnlineSPDphys(fFileName.Data(), readFromGridFile);
}

void AliITSOnlineSPDphysAnalyzer::SetParam(const Char_t *pname, const Char_t *pval) {
  // set a parameter
  TString name = pname;
  TString val = pval;
  //  printf("Setting Param %s  to %s\n",name.Data(),val.Data());
  if (name.CompareTo("MistakeProbabilityNoisy")==0) {
    Double_t mistakeProbabilityNoisy = val.Atof();
    fThreshNoisy = mistakeProbabilityNoisy/(20*6*10*32*256);
  }
  else if (name.CompareTo("MistakeProbabilityDead")==0) {
    Double_t mistakeProbabilityDead = val.Atof();
    fThreshDead = mistakeProbabilityDead/(20*6*10*32*256);
  }
  else if (name.CompareTo("fMinEventsForNoisy")==0) {
    fMinEventsForNoisy = val.Atoi();
  }
  else if (name.CompareTo("fMinEventsForDead")==0) {
    fMinEventsForDead = val.Atoi();
  }
  else if (name.CompareTo("fDefinitelyNoisyRatio")==0) {
    fDefinitelyNoisyRatio = val.Atof();
  }
  else if (name.CompareTo("fMinNrEqHitsForDeadChips")==0) {
    fMinNrEqHitsForDeadChips = val.Atof();
  }
  else if (name.CompareTo("fRatioToMeanForInefficientChip")==0) {
    fRatioToMeanForInefficientChip = val.Atof();
  }
  else {
    Error("AliITSOnlineSPDphysAnalyzer::SetParam","Parameter %s in configuration file unknown.",name.Data());
  }
}

void AliITSOnlineSPDphysAnalyzer::ReadParamsFromLocation(const Char_t *dirName) {
  // opens file (default name) in dir dirName and reads parameters from it
  TString paramsFileName = Form("%s/physics_params.txt",dirName);
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


UInt_t AliITSOnlineSPDphysAnalyzer::ProcessNoisyPixels() {
  // process noisy pixel data , returns number of chips with enough statistics
  if (fPhysObj==NULL) {
    Warning("AliITSOnlineSPDphysAnalyzer::ProcessNoisyPixels","No data!");
    return 0;
  }
  // do we have enough events to even try the algorithm?
  if (GetNrEvents() < fMinEventsForNoisy) {
    Warning("AliITSOnlineSPDphysAnalyzer::ProcessNoisyPixels","Nr events (%d) < fMinEventsForNoisy (%d)!",GetNrEvents(),fMinEventsForNoisy);
    return 0;
  }
  // handler should be initialized
  if (fHandler==NULL) {
    Error("AliITSOnlineSPDphysAnalyzer::ProcessNoisyPixels","Calibration handler is not initialized!");
    return 0;
  }
  
  UInt_t nrEnoughStatChips = 0;

  for (UInt_t hs=0; hs<6; hs++) {
    for (UInt_t chip=0; chip<10; chip++) {

      UInt_t nrPixels = 0;
      UInt_t nrChipHits = 0;
      UInt_t nrMostHits = 0;
      for (UInt_t col=0; col<32; col++) {
	for (UInt_t row=0; row<256; row++) {
	  UInt_t nrHits = fPhysObj->GetHits(hs,chip,col,row);
	  nrChipHits += nrHits;
	  //	  if (nrHits>0) nrPixels++; // don't include pixels that might be dead
	  nrPixels++;
	  if (nrHits>fDefinitelyNoisyRatio*GetNrEvents()) {
	    fHandler->SetNoisyPixel(GetEqNr(),hs,chip,col,row);
	    nrPixels--;
	    nrChipHits-=nrHits;
	  }
	  else {
	    if (nrMostHits<nrHits) nrMostHits=nrHits;
	  }
	}
      }

      if (nrChipHits>0) { // otherwise there are for sure no noisy
	// Binomial with n events and probability p for pixel hit
	UInt_t n = GetNrEvents();
	if (nrPixels>0 && n>0) {

	  Double_t p = (Double_t)nrChipHits/nrPixels/n;

	  // Bin(n,k=0):
	  Double_t bin = pow((Double_t)(1-p),(Double_t)n);
	  // Bin(n,k)
	  UInt_t k=1;
	  while ((bin>fThreshNoisy || k<n*p) && k<=n) {
	    k++;
	    bin = bin*(n-k+1)/k*p/(1-p);
	  }
	  
	  // can we find noisy pixels...?
	  if (k<=n) {
	    //	    printf("eq %d , hs %d , chip %d : Noisy level = %d\n",GetEqNr(),hs,chip,k);
	    nrEnoughStatChips++;
	    // add noisy pixels to handler
	    UInt_t noiseLimit=k;
	    if (nrMostHits>=noiseLimit) {
	      for (UInt_t col=0; col<32; col++) {
		for (UInt_t row=0; row<256; row++) {
		  UInt_t nrHits = fPhysObj->GetHits(hs,chip,col,row);
		  if (nrHits >= noiseLimit) {
		    fHandler->SetNoisyPixel(GetEqNr(),hs,chip,col,row);
		  }
		}
	      }
	    }
	  }
	}

      }

    } // for chip
  } // for hs

  return nrEnoughStatChips;
}

UInt_t AliITSOnlineSPDphysAnalyzer::ProcessNoisyPixels(UInt_t eq, UInt_t nrEvts) {
  // process noisy pixel data , returns number of chips with enough statistics
  if (fPhysObj==NULL) {
    Warning("AliITSOnlineSPDphysAnalyzer::ProcessNoisyPixels","No data!");
    return 0;
  }
  // do we have enough events to even try the algorithm?
  if (nrEvts < fMinEventsForNoisy) {
    Warning("AliITSOnlineSPDphysAnalyzer::ProcessNoisyPixels","Nr events (%d) < fMinEventsForNoisy (%d)!",nrEvts,fMinEventsForNoisy);
    return 0;
  }
  // handler should be initialized
  if (fHandler==NULL) {
    Error("AliITSOnlineSPDphysAnalyzer::ProcessNoisyPixels","Calibration handler is not initialized!");
    return 0;
  }
  
  UInt_t nrEnoughStatChips = 0;

  for (UInt_t hs=0; hs<6; hs++) {
    for (UInt_t chip=0; chip<10; chip++) {

      UInt_t nrPixels = 0;
      UInt_t nrChipHits = 0;
      UInt_t nrMostHits = 0;
      for (UInt_t col=0; col<32; col++) {
	for (UInt_t row=0; row<256; row++) {
	  UInt_t nrHits = fPhysObj->GetHits(hs,chip,col,row);
	  nrChipHits += nrHits;
	  //	  if (nrHits>0) nrPixels++; // don't include pixels that might be dead
	  nrPixels++;
	  if (nrHits>fDefinitelyNoisyRatio*nrEvts) {
	    fHandler->SetNoisyPixel(eq,hs,chip,col,row);
	    nrPixels--;
	    nrChipHits-=nrHits;
	  }
	  else {
	    if (nrMostHits<nrHits) nrMostHits=nrHits;
	  }
	}
      }

      if (nrChipHits>0) { // otherwise there are for sure no noisy
	// Binomial with n events and probability p for pixel hit
	UInt_t n = nrEvts;
	if (nrPixels>0 && n>0) {

	  Double_t p = (Double_t)nrChipHits/nrPixels/n;

	  // Bin(n,k=0):
	  Double_t bin = pow((Double_t)(1-p),(Double_t)n);
	  // Bin(n,k)
	  UInt_t k=1;
	  while ((bin>fThreshNoisy || k<n*p) && k<=n) {
	    k++;
	    bin = bin*(n-k+1)/k*p/(1-p);
	  }
	  
	  // can we find noisy pixels...?
	  if (k<=n) {
	    //	    printf("eq %d , hs %d , chip %d : Noisy level = %d\n",GetEqNr(),hs,chip,k);
	    nrEnoughStatChips++;
	    // add noisy pixels to handler
	    UInt_t noiseLimit=k;
	    if (nrMostHits>=noiseLimit) {
	      for (UInt_t col=0; col<32; col++) {
		for (UInt_t row=0; row<256; row++) {
		  UInt_t nrHits = fPhysObj->GetHits(hs,chip,col,row);
		  if (nrHits >= noiseLimit) {
		    fHandler->SetNoisyPixel(eq,hs,chip,col,row);
		  }
		}
	      }
	    }
	  }
	}

      }

    } // for chip
  } // for hs

  return nrEnoughStatChips;
}


UInt_t AliITSOnlineSPDphysAnalyzer::ProcessDeadPixels() {
  // process dead pixel data , returns number of chips with enough statistics
  if (fPhysObj==NULL) {
    Warning("AliITSOnlineSPDphysAnalyzer::ProcessDeadPixels","No data!");
    return 0;
  }

  // do we have enough events to even try the algorithm?
  if (GetNrEvents() < fMinEventsForDead) {
    Warning("AliITSOnlineSPDphysAnalyzer::ProcessDeadPixels","Nr events (%d) < fMinEventsForDead (%d)!",GetNrEvents(),fMinEventsForDead);
    return 0;
  }
  // handler should be initialized
  if (fHandler==NULL) {
    Error("AliITSOnlineSPDphysAnalyzer::ProcessDeadPixels","Calibration handler is not initialized!");
    return 0;
  }

  AliITSIntMap* possiblyDead  = new AliITSIntMap();
  AliITSIntMap* possiblyIneff = new AliITSIntMap();

  fNrEnoughStatChips = 0;
  fNrDeadChips = 0;
  fNrInefficientChips = 0;
  UInt_t nrPossiblyDeadChips = 0;
  fNrEqHits = 0;


  for (UInt_t hs=0; hs<6; hs++) {
    if (!fHandler->IsActiveHS(GetEqNr(),hs)) {
      fNrDeadChips+=10;
    }
    else {
      for (UInt_t chip=0; chip<10; chip++) {
	if (!fHandler->IsActiveChip(GetEqNr(),hs,chip)) {
	  fNrDeadChips++;
	}
	else {
	  // perform search for individual dead pixels...
	  Bool_t good=kFALSE;

	  UInt_t nrPossiblyDeadPixels = 0;
	  UInt_t nrPixels = 0;
	  UInt_t nrChipHits = 0;
	  for (UInt_t col=0; col<32; col++) {
	    for (UInt_t row=0; row<256; row++) {
	      UInt_t nrHits = fPhysObj->GetHits(hs,chip,col,row);
	      nrChipHits += nrHits;
	      if (!fHandler->IsPixelNoisy(GetEqNr(),hs,chip,col,row)) {
		// don't include noisy pixels
		nrPixels++;
		if (nrHits==0) {
		  nrPossiblyDeadPixels++;
		}
		else {
		  fHandler->UnSetDeadPixel(GetEqNr(),hs,chip,col,row); // unset (no action unless dead before)
		}
	      }
	      else {
		nrChipHits -= nrHits; // this is needed when running offline (online nrHits should be 0 already)
	      }
	    }
	  }
	  fNrEqHits+=nrChipHits;

	  if (nrChipHits>0) {
	    // make sure the chip is not flagged as dead
	    fHandler->SetDeadChip(GetEqNr(),hs,chip,kFALSE);
	  }

	  if (nrPossiblyDeadPixels==0) {
	    // no need to see if we have enough statistics...
	    fNrEnoughStatChips++;
	    good=kTRUE;
	    //	printf("%3d",good);
	    //	if (chip==9) printf("\n");
	    continue;
	  }

	  if (nrChipHits==0) {
	    nrPossiblyDeadChips++;
	    possiblyDead->Insert(hs,chip);
	    good=kFALSE;
	    //	printf("%3d",good);
	    //	if (chip==9) printf("\n");
	    continue;
	  }

	  // Binomial with n events and probability p for pixel hit
	  UInt_t n = GetNrEvents();
	  if (nrPixels>0 && n>0) {

	    Double_t p = (Double_t)nrChipHits/nrPixels/n;

	    // probability of falsely assigning a dead pixel
	    Double_t falselyDeadProb = pow((Double_t)(1-p),(Double_t)n);
	    //	    printf("falselyprob=%e\n",falselyDeadProb);

	    // can we find dead pixels...?
	    if (falselyDeadProb<fThreshDead) {
	      fNrEnoughStatChips++;
	      good=kTRUE;
	      // add dead pixels to handler
	      for (UInt_t col=0; col<32; col++) {
		for (UInt_t row=0; row<256; row++) {
		  UInt_t nrHits = fPhysObj->GetHits(hs,chip,col,row);
		  if (nrHits==0) {
		    if (!fHandler->IsPixelNoisy(GetEqNr(),hs,chip,col,row)) {
		      // don't include noisy pixels
		      fHandler->SetDeadPixel(GetEqNr(),hs,chip,col,row);
		    }
		  }
		}
	      }
	    }
	    if (!good) {
	      // this might be an inefficient chip
	      possiblyIneff->Insert(hs*10+chip,nrChipHits);
	    }

	  }
	  else {
	    if (n>0) {
	      // this is a completely noisy chip... put in category enough stat
	      fNrEnoughStatChips++;
	      good=kTRUE;
	    }
	  }

	  //      printf("%3d",good);
	  //      if (chip==9) printf("\n");

	}
      } // for chip
    }
  } // for hs
 

  Int_t key,val;

  // dead chips?
  if (fNrEqHits>fMinNrEqHitsForDeadChips) {
    while (possiblyDead->Pop(key,val)) {
      fHandler->SetDeadChip(GetEqNr(),key,val,kFALSE);
    }
    fNrDeadChips+=nrPossiblyDeadChips;
  }
  delete possiblyDead;

  // inefficient chips?
  while (possiblyIneff->Pop(key,val)) {
    if (val<fNrEqHits/60*fRatioToMeanForInefficientChip) {
      fNrInefficientChips++;
    }
  }
  delete possiblyIneff;


  fbDeadProcessed = kTRUE;

  return fNrEnoughStatChips;
}


UInt_t AliITSOnlineSPDphysAnalyzer::GetNrEnoughStatChips() {
  // returns nr of enough stat chips
  if (!fbDeadProcessed) ProcessDeadPixels();
  return fNrEnoughStatChips;
}
UInt_t AliITSOnlineSPDphysAnalyzer::GetNrDeadChips() {
  // returns nr of dead chips
  if (!fbDeadProcessed) ProcessDeadPixels();
  return fNrDeadChips;
}
UInt_t AliITSOnlineSPDphysAnalyzer::GetNrInefficientChips() {
  // returns nr of inefficient chips
  if (!fbDeadProcessed) ProcessDeadPixels();
  return fNrInefficientChips;
}
UInt_t AliITSOnlineSPDphysAnalyzer::GetNrNeedsMoreStatChips() {
  // returns nr of needs more stat chips
  if (!fbDeadProcessed) ProcessDeadPixels();
  return 60-fNrEnoughStatChips-fNrDeadChips-fNrInefficientChips;
}

UInt_t AliITSOnlineSPDphysAnalyzer::GetEqNr() const {
  // returns the eq nr of phys obj
  if (fPhysObj!=NULL) return fPhysObj->GetEqNr(); 
  else return 999;
}

UInt_t AliITSOnlineSPDphysAnalyzer::GetNrEvents() const {
  // returns the nr of events of phys obj
  if (fPhysObj!=NULL) return fPhysObj->GetNrEvents(); 
  else return 0;
}

void AliITSOnlineSPDphysAnalyzer::Exponent(Double_t &val, Int_t &valExp) const {
  // put double in format with val and exp so that 1<val<10 - The actual value is val*10e(valExp)
  while (val>10) {
    val/=10;
    valExp++;
  }
  while (val<1) {
    val*=10;
    valExp--;
  }
}
//____________________________________________________________________________________________________
TH2F* AliITSOnlineSPDphysAnalyzer::GetPhysicalHitMapTot() {
  // creates and returns a pointer to a hitmap histo (equipment opened up)
  // physical representation of the half stave hitmap.

  if (fPhysObj==NULL) {
    Error("AliITSOnlineSPDphysAnalyzer::GetPhysicalHitMapTot","No data!");
    return NULL;
  }
  TString histoname = Form("Eq %d",GetEqNr());
  TH2F* hPhysicalHitMapTot = new TH2F(histoname.Data(),histoname.Data(),32*10,-0.5,32*10-0.5,256*6,-0.5,256*6-0.5);
  hPhysicalHitMapTot->SetNdivisions(-10,"X");
  hPhysicalHitMapTot->SetNdivisions(-006,"Y");
  hPhysicalHitMapTot->SetTickLength(0,"X");
  hPhysicalHitMapTot->SetTickLength(0,"Y");
  hPhysicalHitMapTot->GetXaxis()->SetLabelColor(gStyle->GetCanvasColor());
  hPhysicalHitMapTot->GetYaxis()->SetLabelColor(gStyle->GetCanvasColor());
  Int_t correctChip=-1;
  for (UInt_t hs=0; hs<6; hs++) {
    for (UInt_t chipNr=0; chipNr<10; chipNr++) {
      if(GetEqNr()<10) correctChip=9-chipNr;
      else correctChip=chipNr;
      for (UInt_t col=0; col<32; col++) {
        for (UInt_t row=0; row<256; row++) {
          if(hs>1) hPhysicalHitMapTot->Fill(correctChip*32+(31-col),(5-hs)*256+row,fPhysObj->GetHits(hs,chipNr,col,row));
          else hPhysicalHitMapTot->Fill(correctChip*32+col,(5-hs)*256+row,fPhysObj->GetHits(hs,chipNr,col,row));
        }
      }
    }
  }
  return hPhysicalHitMapTot;
}
//_____________________________________________________________________________________________________
TH2F* AliITSOnlineSPDphysAnalyzer::GetHitMapTot() {
  // creates and returns a pointer to a hitmap histo (half sector style a la spdmood)
  // This histogram  shown the read out numbering pattern, it is not the physical one.
  if (fPhysObj==NULL) {
    Error("AliITSOnlineSPDphysAnalyzer::GetHitMapTot","No data!");
    return NULL;
  }
  TString histoname = Form("Eq %d",GetEqNr());
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
	  fHitMapTot->Fill(chipNr*32+col,(5-hs)*256+row,fPhysObj->GetHits(hs,chipNr,col,row));
	}
      }
    }
  }
  return fHitMapTot;
}
//________________________________________________________________________________________________________
TH2F* AliITSOnlineSPDphysAnalyzer::GetHitMapChip(UInt_t hs, UInt_t chip) {
  // creates and returns a pointer to a hitmap histo (chip style a la spdmood)
  if (fPhysObj==NULL) {
    Error("AliITSOnlineSPDphysAnalyzer::GetHitMapChip","No data!");
    return NULL;
  }

  TString histoName;
  TString histoTitle;
  histoName = Form("fChipHisto_%d_%d_%d", GetEqNr(), hs, chip);
  histoTitle = Form("Eq ID %d, Half Stave %d, Chip %d", GetEqNr(), hs, chip);

  TH2F *returnHisto = new TH2F(histoName.Data(), histoTitle.Data(), 32, -0.5, 31.5, 256, -0.5, 255.5);
  returnHisto->SetMinimum(0);
  for (UInt_t col=0; col<32; col++) {
    for (UInt_t row=0; row<256; row++) {
      if(hs<2) returnHisto->Fill(31-col,row,fPhysObj->GetHits(hs,chip,col,row));
      else returnHisto->Fill(col,row,fPhysObj->GetHits(hs,chip,col,row));
    }
  }

  return returnHisto;
}

