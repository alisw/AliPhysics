/////////////////////////////////////////////////////
// AliAnalysisTaskFlowEventfoRP:
// analysis task to fill the flow event 
// and calculate the ReactionPlane for the AODheader
//////////////////////////////////////////////////////

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id: $ */

#ifndef AliAnalysisTaskFlowEventforRP_H
#define AliAnalysisTaskFlowEventforRP_H

//class AliESDEvent;
//class AliAODEvent;
class AliCFManager;
class AliFlowEventSimpleMaker;

#include "TString.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskFlowEventforRP : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskFlowEventforRP();
  AliAnalysisTaskFlowEventforRP(const char *name);
  virtual ~AliAnalysisTaskFlowEventforRP();
  
  //virtual void   ConnectInputData(Option_t *);
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetAnalysisType(TString type) { this->fAnalysisType = type; }
  TString GetAnalysisType() const    { return this->fAnalysisType; }

  void    SetMinMult(Int_t multmin)    {this->fMinMult = multmin; }
  Int_t   GetMinMult() const           {return this->fMinMult; }
  void    SetMaxMult(Int_t multmax)    {this->fMaxMult = multmax; }
  Int_t   GetMaxMult() const           {return this->fMaxMult; }
  
  void          SetCFManager1(AliCFManager* cfmgr) {this->fCFManager1 = cfmgr; } 
  AliCFManager* GetCFManager1()           {return this->fCFManager1; }
  void          SetCFManager2(AliCFManager* cfmgr) {this->fCFManager2 = cfmgr; } 
  AliCFManager* GetCFManager2()           {return this->fCFManager2; }
  
  
 private:

  AliAnalysisTaskFlowEventforRP(const AliAnalysisTaskFlowEventforRP& aAnalysisTask);
  AliAnalysisTaskFlowEventforRP& operator=(const AliAnalysisTaskFlowEventforRP& aAnalysisTask); 

  //TFile*        fOutputFile;              // temporary output file for testing
  //AliESDEvent*  fESD;                     // ESD object
  //AliAODEvent*  fAOD;                     // AOD object
  AliFlowEventSimpleMaker* fEventMaker;   // FlowEventSimple maker object
  TString       fAnalysisType;            // can be MC, ESD or AOD
  AliCFManager* fCFManager1;              // correction framework manager
  AliCFManager* fCFManager2;              // correction framework manager
  Int_t         fMinMult;                 // Minimum multiplicity from tracks selected using CORRFW
  Int_t         fMaxMult;                 // Maximum multiplicity from tracks selected using CORRFW 
    
  Double_t  fMCReactionPlaneAngle;   // the angle of the reaction plane from the MC truth
  
  ClassDef(AliAnalysisTaskFlowEventforRP, 1); // example of analysis
};

#endif

