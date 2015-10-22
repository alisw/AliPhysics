#ifndef ALIANALYSISTASKZDCPP_H
#define ALIANALYSISTASKZDCPP_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*****************************************************
//   Class AliAnalysisTaskZDCpp
//   author: Chiara Oppedisano
//*****************************************************

#include "AliAnalysisTaskSE.h"

class TROOT;
class TSystem;
class TList;
class TFile;
class TH2F;
class TTree;

class AliAnalysisTaskZDCpp : public AliAnalysisTaskSE {

 public:

  AliAnalysisTaskZDCpp();
  AliAnalysisTaskZDCpp(const char *name);
  AliAnalysisTaskZDCpp& operator= (const AliAnalysisTaskZDCpp& ana);
  AliAnalysisTaskZDCpp(const AliAnalysisTaskZDCpp& c);
  virtual ~AliAnalysisTaskZDCpp();

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
  TH2F    *fDebunch;            //! TDC sum vs. diff
  TH1F    *fhTDCZNC;		//! ZNC TDC 
  TH1F    *fhTDCZNA;		//! ZNC TDC 
 
  ClassDef(AliAnalysisTaskZDCpp,1); 

};

#endif

