#ifndef ALIANALYSISTASKZDC_H
#define ALIANALYSISTASKZDC_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*****************************************************
//   Class AliAnalysisTaskZDC
//   author: Chiara Oppedisano
//*****************************************************

#include "AliAnalysisTaskSE.h"

class TROOT;
class TSystem;
class TList;
class TFile;
class TH2F;
class TTree;

class AliAnalysisTaskZDC : public AliAnalysisTaskSE {

 public:

  AliAnalysisTaskZDC();
  AliAnalysisTaskZDC(const char *name);
  AliAnalysisTaskZDC& operator= (const AliAnalysisTaskZDC& ana);
  AliAnalysisTaskZDC(const AliAnalysisTaskZDC& c);
  virtual ~AliAnalysisTaskZDC();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  virtual void  SetDebugLevel(Int_t level) {fDebug = level;}
  void SetMCInput() {fIsMCInput = kTRUE;}
 
 private:

  Int_t    fDebug;		//  Debug flag
  Bool_t   fIsMCInput;  	// true when input is MC
  //
  TList   *fOutput;		//! list send on output slot 0
  TH1F    *fhTDCZNSum;  	//! TDC ZNC sum
  TH1F    *fhTDCZNDiff; 	//! TDC DIFF sum
  TH1F    *fhZNCSpectrum;	//! ZNC spectra
  TH1F    *fhZNASpectrum;	//! ZNA spectra
  TH1F    *fhZPCSpectrum;	//! ZPC spectra
  TH1F    *fhZPASpectrum;	//! ZPA spectra
  TH1F    *fhZEM1Spectrum;	//! ZEM1 spectra
  TH1F    *fhZEM2Spectrum;	//! ZEM2 spectra
  TH1F    *fhZNCpmc;		//! ZNC PMCs
  TH1F    *fhZNApmc;		//! ZNA PMCs
  TH1F    *fhZPCpmc;		//! ZPC PMCs
  TH1F    *fhZPApmc;		//! ZPA PMCs
  TH2F    *fhZNCCentroid;       //! ZNC centroid
  TH2F    *fhZNACentroid;       //! ZNA centroid
  TH1F    *fhZNCemd;	     	//! ZNC sum low gain chain
  TH1F    *fhZNAemd;	     	//! ZNA sum gain chain
  TH1F    *fhPMCZNCemd;      	//! ZNC PMC low gain chain
  TH1F    *fhPMCZNAemd;      	//! ZNA PMC low gain chain
  TH2F    *fDebunch;            //! TDC sum vs. diff
  
 
  ClassDef(AliAnalysisTaskZDC,1); 

};

#endif

