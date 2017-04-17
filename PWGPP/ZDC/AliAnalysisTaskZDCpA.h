#ifndef ALIANALYSISTASKZDCPA_H
#define ALIANALYSISTASKZDCPA_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*****************************************************
//   Class AliAnalysisTaskZDCpA
//   author: Chiara Oppedisano
//*****************************************************

#include "AliAnalysisTaskSE.h"

class TROOT;
class TSystem;
class TList;
class TFile;
class TH2F;
class TTree;

class AliAnalysisTaskZDCpA : public AliAnalysisTaskSE {

 public:

  AliAnalysisTaskZDCpA();
  AliAnalysisTaskZDCpA(const char *name);
  AliAnalysisTaskZDCpA& operator= (const AliAnalysisTaskZDCpA& ana);
  AliAnalysisTaskZDCpA(const AliAnalysisTaskZDCpA& c);
  virtual ~AliAnalysisTaskZDCpA();

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
  TH1F    *fhTDCZNC;  		//! TDC ZNC
  TH1F    *fhTDCZNA; 		//! TDC ZNA
  TH1F    *fhTDCZNSum;  	//! TDC ZNC sum
  TH1F    *fhTDCZNDiff; 	//! TDC DIFF sum
  TH1F    *fhZNCSumQ;	 	//! ZNC sum 4Q
  TH1F    *fhZNASumQ;		//! ZNA sum 4Q
  TH1F    *fhZPCSumQ;		//! ZPC sum 4Q
  TH1F    *fhZPASumQ;		//! ZPA sum 4Q
  TH1F    *fhZEM1Spectrum;	//! ZEM1 spectra
  TH1F    *fhZEM2Spectrum;	//! ZEM2 spectra
  TH1F    *fhZNCpmc;		//! ZNC PMCs
  TH1F    *fhZNApmc;		//! ZNA PMCs
  TH1F    *fhZPCpmc;		//! ZPC PMCs
  TH1F    *fhZPApmc;		//! ZPA PMCs
  TH1F    *fhZNCpmcUncalib;	//! uncalibrated ZNC PMCs
  TH1F    *fhZNApmcUncalib;	//! uncalibrated ZNA PMCs
  TH1F    *fhZPCpmcUncalib;	//! uncalibrated ZPC PMCs
  TH1F    *fhZPApmcUncalib;	//! uncalibrated ZPA PMCs
  TH2F    *fhZNCCentroid;       //! ZNC centroid
  TH2F    *fhZNACentroid;       //! ZNA centroid
  TH1F    *fhPMCZNCemdUncalib;  //! ZNC PMC low gain chain
  TH1F    *fhPMCZNAemdUncalib;  //! ZNA PMC low gain chain
  TH1F    *fhPMCZNCemd;      	//! ZNC PMC low gain chain
  TH1F    *fhPMCZNAemd;      	//! ZNA PMC low gain chain
  TH2F    *fDebunch;            //! TDC sum vs. diff
  TH1F    *fhTDCZNAcorr;      	//! ZNA corrected TDC
  TH1F    *fhTDCZNCcorr;      	//! ZNC corrected TDC


  ClassDef(AliAnalysisTaskZDCpA,4);

};

#endif
