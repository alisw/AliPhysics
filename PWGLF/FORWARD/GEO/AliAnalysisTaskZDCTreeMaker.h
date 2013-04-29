#ifndef ALIANALYSISTASKZDCTREEMAKER_H
#define ALIANALYSISTASKZDCTREEMAKER_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*****************************************************
//   Class AliAnalysisTaskZDCTreeMaker
//   author: Chiara Oppedisano
//*****************************************************

#include "AliAnalysisTaskSE.h"

class TROOT;
class TSystem;
class TList;
class TFile;
class TTree;

class AliAnalysisTaskZDCTreeMaker : public AliAnalysisTaskSE {

 public:

  AliAnalysisTaskZDCTreeMaker();
  AliAnalysisTaskZDCTreeMaker(const char *name);
  virtual ~AliAnalysisTaskZDCTreeMaker();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

  virtual void  SetDebugLevel(Int_t level) {fDebug = level;}
  void SetInput(const char* input) {fAnalysisInput = input;}
  void SetMCInput() {fIsMCInput = kTRUE;}
   void SetUseSpecialOutput(Bool_t v=kTRUE) {fUseSpecialOutput = v;}

 private:

  Int_t    fDebug;	   	//  Debug flag
  TString  fAnalysisInput; 	// "ESD", "AOD"
  Bool_t   fIsMCInput;          // true when input is MC
  Bool_t   fUseSpecialOutput;   // do we use special output instead of merging?
  //
  TList   *fOutput;	   	//! list send on output slot 0
  //
  TTree   *fCentralityTree;     //! output tree
  //
  char     fTrigClass[100];	//  fired trigger classes
  //
  Bool_t   fIsEventSelected;    //  is physics selection on
  Bool_t   fIsPileupFromSPD;	//  is pilue up from SPD
  //
  Double_t fxVertex;		//  X vertex from ITS
  Double_t fyVertex;		//  Y vertex from ITS
  Double_t fzVertex;		//  Z vertex from ITS
  Bool_t   fVertexer3d;		//  Is vertex from 3d vertexer?
  //
  Int_t    fNTracklets;		//  no. tracklets
  Int_t    fNClusters[2];	//  no. clusters on SPD layers
  //
  Int_t    fIsV0ATriggered;	//  VOA decision
  Int_t    fIsV0CTriggered;	//  V0C decision
  Float_t  fMultV0A;		//  mult. V0A
  Float_t  fMultV0C;		//  mult. V0C
  //
  UInt_t   fESDFlag;		//  ZDC ESD flags
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
  Float_t  fZNCtowerLG[5];	//  ZNC 5 tower signals
  Float_t  fZPCtowerLG[5];	//  ZPC 5 tower signals
  Float_t  fZNAtowerLG[5];	//  ZNA 5 tower signals
  Float_t  fZPAtowerLG[5];	//  ZPA 5 tower signals
  //
  Int_t    fTDCZNC[4];	     // TDC raw values !ONLY FOR ESDs ANALYSIS!
  Int_t    fTDCZPC[4];	     // TDC raw values !ONLY FOR ESDs ANALYSIS!
  Int_t    fTDCZNA[4];	     // TDC raw values !ONLY FOR ESDs ANALYSIS!
  Int_t    fTDCZPA[4];	     // TDC raw values !ONLY FOR ESDs ANALYSIS!
//  Int_t    fTDCZEM1;	     // TDC raw values !ONLY FOR ESDs ANALYSIS!
//  Int_t    fTDCZEM2;	     // TDC raw values !ONLY FOR ESDs ANALYSIS!
  
  Float_t fCentralityV0M;       // Centrality from V0A+V0C
  Float_t fCentralityV0A;       // Centrality from V0A
  Float_t fCentralityV0C;       // Centrality from V0C
  Float_t fCentralityCL1;       // Centrality from Clusters in layer 1
  Float_t fCentralityZNA;       // Centrality from ZNA
  Float_t fCentralityZPA;       // Centrality from ZPA
  Float_t fCentralityZNC;       // Centrality from ZNC
  Float_t fCentralityZPC;       // Centrality from ZPC
  //
  AliAnalysisTaskZDCTreeMaker& operator= (const AliAnalysisTaskZDCTreeMaker& ana);
  AliAnalysisTaskZDCTreeMaker(const AliAnalysisTaskZDCTreeMaker& c);
  //
  ClassDef(AliAnalysisTaskZDCTreeMaker,3); 

};

#endif

