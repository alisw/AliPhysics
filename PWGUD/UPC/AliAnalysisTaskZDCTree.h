#ifndef AliAnalysisTaskZDCTree_H
#define AliAnalysisTaskZDCTree_H

/* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*****************************************************
//	Based on the class AliAnalysisTaskZDCTreeMaker by Chiara Oppedisano
//	Created by Uliana Dmitrieva uliana.dmitrieva@cern.ch on 11/01/2018
//*****************************************************

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

 private:

  Int_t			fDebug;			//  Debug flag
  
  AliESDEvent	*fESD;			//! input event
  TList			*fOutput;		//! list send on output slot 0
  TTree			*fZDCTree;		//! output tree
  //
  Int_t fRunNum;
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
  ClassDef(AliAnalysisTaskZDCTree,4); 

};

#endif


