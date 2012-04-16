#ifndef ALIANALYSISTASKZDCPBPB_H
#define ALIANALYSISTASKZDCPBPB_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*****************************************************
//   Class AliAnalysisTaskZDCPbPb
//   author: Chiara Oppedisano
//*****************************************************

#include "AliAnalysisTaskSE.h"

class TROOT;
class TSystem;
class TList;
class TFile;
class TH1F;
class TH2F;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskZDCPbPb : public AliAnalysisTaskSE {

 public:
  
  enum kAnalysisInput{kESD=1, kAOD=2}; 
  
  AliAnalysisTaskZDCPbPb();
  AliAnalysisTaskZDCPbPb(const char *name);
  AliAnalysisTaskZDCPbPb& operator= (const AliAnalysisTaskZDCPbPb& ana);
  AliAnalysisTaskZDCPbPb(const AliAnalysisTaskZDCPbPb& c);
  virtual ~AliAnalysisTaskZDCPbPb();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  virtual void  SetDebugLevel(Int_t level) {fDebug = level;}
  void SetInput(int input) {fAnalysisInput = input;}
  void SetMCInput() {fIsMCInput = kTRUE;}
  void SetCentralityRange(Float_t centrlow=0., Float_t centrup=100.) {fCentrLowLim=centrlow;
  fCentrUpLim=centrup;}
  void SetCentralityEstimator(TString centrest = "V0M") {fCentrEstimator=centrest;}              
 
 private:
  
  Int_t    fAnalysisInput;      // analysis input
  Bool_t   fIsMCInput;          // true when input is MC
  Float_t  fCentrLowLim;	// centrality lower limit
  Float_t  fCentrUpLim;		// centrality upper limit
  TString  fCentrEstimator;     // string for the centrality estimator   
  //
  TList   *fOutput;	   	//! list send on output slot 0
  //
  TH1F *fhZNCPM[5];		//! ZNC PM high res.
  TH1F *fhZNAPM[5];		//! ZNA PM high res.
  TH1F *fhZPCPM[5];		//! ZPC PM high res.
  TH1F *fhZPAPM[5];		//! ZPA PM high res.
  TH1F *fhZEM[2];		//! ZEM PM high res.
  TH1F *fhZNCPMlg[5];		//! ZNC PM low res.
  TH1F *fhZNAPMlg[5];		//! ZNA PM low res.
  TH1F *fhZPCPMlg[5];		//! ZPC PM low res.
  TH1F *fhZPAPMlg[5];		//! ZPA PM low res.
  TH1F *fhTDCraw[6];		//! raw TDC histos
  TH1F *fhTDC[6];		//! corrected TDC histos
  //
  TH2F *fhZNCvsZNA;		//! ZNC vs ZNA;
  TH2F *fhZPCvsZPA;		//! ZPC vs ZPA;
  TH2F *fhZDCCvsZDCCA;		//! ZDCC vs ZDCCA
  TH2F *fhZNvsZEM;		//! ZN vs ZEM;
  TH2F *fhZDCvsZEM;		//! ZDC vs ZEM;
  TH2F *fhZNvsVZERO;		//! ZN vs VZERO;
  TH2F *fhZDCvsVZERO;		//! ZDC vs VZERO;
  TH2F *fhZDCvsTracklets;	//! ZDC vs N_tracklets;
  TH2F *fhVZEROvsZEM;		//! VZERO vs ZEM;
  TH2F *fhDebunch;		//! Debunch;
  TH2F *fhZNCcentroid;		//! ZNC centroid
  TH2F *fhZNAcentroid;		//! ZNA centroid
  TH2F *fhPMCvsPMQ[4];		//! PMC vs sum PMQi
  //
  TH1F *fhAsymm;		//! ZN asymmetry
  TH2F *fhZNAvsAsymm;		//! ZNA vs asymmetry
  TH2F *fhZNCvsAsymm;		//! ZNC vs asymmetry
 
  ClassDef(AliAnalysisTaskZDCPbPb,1); 

};

#endif

