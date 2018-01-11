#ifndef ALIANALYSISTASKZDCTREE_H
#define ALIANALYSISTASKZDCTREE_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//------------------------------------------------------
//   Based on class AliAnalysisTaskZDCTreeMaker by Chiara Oppedisano
//------------------------------------------------------

#include "AliAnalysisTaskSE.h"

class TROOT;
class TSystem;
class TList;
class TFile;
class TTree;

class AliAnalysisTaskZDCTree : public AliAnalysisTaskSE {

 public:

  AliAnalysisTaskZDCTree();
  AliAnalysisTaskZDCTree(const char *name);
  virtual ~AliAnalysisTaskZDCTree();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  virtual void  SetDebugLevel(Int_t level) {fDebug = level;}
  void SetInput(const char* input) {fAnalysisInput = input;}
  void SetMCInput() {fIsMCInput = kTRUE;}
   

 private:

  Int_t    fDebug;	   	//  Debug flag
  TString  fAnalysisInput; 	// "ESD", "AOD"
  Bool_t   fIsMCInput;          // true when input is MC
  //
  TList   *fOutput;	   	//! list send on output slot 0
  //
  TTree   *fZDCTree;     	// output tree
  //
  Bool_t fIsZEM1;		//
  Bool_t fIsZEM2;		//
  Bool_t fIsZNC;		//
  Bool_t fIsZPC;		//
  Bool_t fIsZNA;		//
  Bool_t fIsZPA;  		//
 
  Bool_t fIs1ZED;		//
  Bool_t fIsCTRUE; 		//

  Float_t  fZNCEnergy;		//  ZNC Energy
  Float_t  fZPCEnergy;		//  ZPC Energy
  Float_t  fZNAEnergy;		//  ZNA Energy
  Float_t  fZPAEnergy;		//  ZPA Energy
  Float_t  fZEM1Energy;		//  ZEM1 Energy
  Float_t  fZEM2Energy;		//  ZEM2 Energy
  
  Float_t  fZNCtower[5];	//  ZNC 5 tower signals
  Float_t  fZPCtower[5];	//  ZPC 5 tower signals
  Float_t  fZNAtower[5];	//  ZNA 5 tower signals
  Float_t  fZPAtower[5];	//  ZPA 5 tower signals

  Float_t  fZNCTDC[4];		//  ZNC TDC
  Float_t  fZPCTDC[4];		//  ZPC TDC
  Float_t  fZNATDC[4];		//  ZNA TDC
  Float_t  fZPATDC[4];		//  ZPA TDC

  AliAnalysisTaskZDCTree& operator= (const AliAnalysisTaskZDCTree& ana);
  AliAnalysisTaskZDCTree(const AliAnalysisTaskZDCTree& c);
  //
  ClassDef(AliAnalysisTaskZDCTree,1); 

};

#endif


