#ifndef ALIANALYSISTASKZDCPACALIB_H
#define ALIANALYSISTASKZDCPACALIB_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*****************************************************
//   Class AliAnalysisTaskZDCpAcalib
//   author: Chiara Oppedisano
//*****************************************************

#include "AliAnalysisTaskSE.h"

class TROOT;
class TSystem;
class TList;
class TFile;
class TTree;

class AliAnalysisTaskZDCpAcalib : public AliAnalysisTaskSE {

 public:

  AliAnalysisTaskZDCpAcalib();
  AliAnalysisTaskZDCpAcalib(const char *name);
  virtual ~AliAnalysisTaskZDCpAcalib();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  virtual void  SetDebugLevel(Int_t level) {fDebug = level;}
  void SetInput(const char* input) {fAnalysisInput = input;}
  void SetMCInput() {fIsMCInput = kTRUE;}
   void SetUseSpecialOutput(Bool_t v=kTRUE) {fUseSpecialOutput = v;}

 private:

  Int_t    fDebug;	       //  Debug flag
  TString  fAnalysisInput;     // "ESD", "AOD"
  Bool_t   fIsMCInput;         // true when input is MC
  Bool_t   fUseSpecialOutput;  // do we use special output instead of merging?
  //
  TList   *fOutput;	       //! list send on output slot 0
  //
  TTree   *fCentralityTree;    //! output tree
  //
  char     fTrigClass[100];    //  fired trigger classes
  //
  Bool_t   fIsEventSelected;   //  is physics selection on
  //
  Int_t    fIsV0ATriggered;    //  VOA decision
  Int_t    fIsV0CTriggered;    //  V0C decision
  //
  Float_t  fZNCtower[5];       //  ZNC 5 tower signals
  Float_t  fZNAtower[5];       //  ZNA 5 tower signals
  Float_t  fZNCtowerLG[5];     //  ZNC 5 tower signals
  Float_t  fZNAtowerLG[5];     //  ZPC 5 tower signals
  Float_t  fZPCtower[5];       //  ZNC 5 tower signals
  Float_t  fZPAtower[5];       //  ZNA 5 tower signals
  Float_t  fZPCtowerLG[5];     //  ZNC 5 tower signals
  Float_t  fZPAtowerLG[5];     //  ZPC 5 tower signals
  //
  Float_t    fZNCtdc[4];       // TDC raw values !ONLY FOR ESDs ANALYSIS!
  Float_t    fZNAtdc[4];       // TDC raw values !ONLY FOR ESDs ANALYSIS!
  Float_t    fZPCtdc[4];       // TDC raw values !ONLY FOR ESDs ANALYSIS!
  Float_t    fZPAtdc[4];       // TDC raw values !ONLY FOR ESDs ANALYSIS!
  
  AliAnalysisTaskZDCpAcalib& operator= (const AliAnalysisTaskZDCpAcalib& ana);
  AliAnalysisTaskZDCpAcalib(const AliAnalysisTaskZDCpAcalib& c);
  //
  ClassDef(AliAnalysisTaskZDCpAcalib,2); 

};

#endif

