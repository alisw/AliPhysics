/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id: $ */

#ifndef AliAnalysisTaskFlowEvent_H
#define AliAnalysisTaskFlowEvent_H

// AliAnalysisTaskFlowEvent:
// analysis task to fill the flow event and make it available to the methods


class AliESDEvent;
class AliAODEvent;
class AliCFManager;
class AliFlowEventSimpleMaker;
class TList;
class TRandom3;

#include "TString.h"
#include "AliAnalysisTask.h"

class AliAnalysisTaskFlowEvent : public AliAnalysisTask {
 public:
  AliAnalysisTaskFlowEvent();
  AliAnalysisTaskFlowEvent(const char *name,Bool_t QAon);
  AliAnalysisTaskFlowEvent(const char *name,Bool_t QAon,UInt_t);
  virtual ~AliAnalysisTaskFlowEvent();
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetAnalysisType(TString type) { this->fAnalysisType = type; }
  TString GetAnalysisType() const    { return this->fAnalysisType; }

  void          SetCFManager1(AliCFManager* cfmgr) {this->fCFManager1 = cfmgr; } 
  AliCFManager* GetCFManager1()           {return this->fCFManager1; }
  void          SetCFManager2(AliCFManager* cfmgr) {this->fCFManager2 = cfmgr; } 
  AliCFManager* GetCFManager2()           {return this->fCFManager2; }
  void          SetQAList1(TList* list)   {this->fQAInt = list; }
  TList*        GetQAList1()              {return this->fQAInt; }
  void          SetQAList2(TList* list)   {this->fQADiff = list; }
  TList*        GetQAList2()              {return this->fQADiff; }
  void          SetQAOn(Bool_t kt)        {this->fQA = kt; }
  Bool_t        GetQAOn()                 {return this->fQA; }
  //setters for adding by hand flow values (afterburner)
  void SetMCReactionPlaneAngle(Double_t fPhiRP)  { this->fMCReactionPlaneAngle = fPhiRP; }
  void SetNoOfLoops(Int_t noofl) {this->fNoOfLoops = noofl;}
  Int_t GetNoOfLoops() const {return this->fNoOfLoops;} 
  void SetEllipticFlowValue(Double_t elfv) {this->fEllipticFlowValue = elfv;}
  Double_t GetEllipticFlowValue() const {return this->fEllipticFlowValue;} 
  void SetSigmaEllipticFlowValue(Double_t sigelfv) {this->fSigmaEllipticFlowValue = sigelfv;}
  Double_t GetSigmaEllipticFlowValue() const {return this->fSigmaEllipticFlowValue;} 
  void SetMultiplicityOfEvent(Int_t multevnt) {this->fMultiplicityOfEvent = multevnt;}
  Int_t GetMultiplicityOfEvent() const {return this->fMultiplicityOfEvent;} 
  void SetSigmaMultiplicityOfEvent(Int_t sigmultevnt) {this->fSigmaMultiplicityOfEvent = sigmultevnt;}
  Int_t GetSigmaMultiplicityOfEvent() const {return this->fSigmaMultiplicityOfEvent;} 
  //end setters afterburner

 private:

  AliAnalysisTaskFlowEvent(const AliAnalysisTaskFlowEvent& aAnalysisTask);
  AliAnalysisTaskFlowEvent& operator=(const AliAnalysisTaskFlowEvent& aAnalysisTask); 

  //  TFile*        fOutputFile;              // temporary output file for testing
  AliESDEvent*  fESD;                     // ESD object
  AliAODEvent*  fAOD;                     // AOD object
  AliFlowEventSimpleMaker* fEventMaker;   // FlowEventSimple maker object
  TString       fAnalysisType;            // can be MC, ESD or AOD
  AliCFManager* fCFManager1;              // correction framework manager
  AliCFManager* fCFManager2;              // correction framework manager
  TList*        fQAInt;                   // QA histogram list
  TList*        fQADiff;                  // QA histogram list
  Bool_t fQA;                             // flag to set the filling of the QA hostograms
  // values afterburner
  Double_t  fMCReactionPlaneAngle;   // the angle of the reaction plane from the MC truth
  Int_t     fCount;   // counter for the number of events processed
  Int_t     fNoOfLoops; // number of times to use the same particle (nonflow) 
  Double_t  fEllipticFlowValue; // Add Flow. Must be in range [0,1].
  Double_t  fSigmaEllipticFlowValue; // Sigma Flow (Gaussian). Must be in range [0,1].
  Int_t     fMultiplicityOfEvent; // Set maximal multiplicity.
  Int_t     fSigmaMultiplicityOfEvent; // Sigma multiplicity (Gaussian).
    
  TRandom3* fMyTRandom3; // our TRandom3 generator
  // end afterburner
  
  ClassDef(AliAnalysisTaskFlowEvent, 1); // example of analysis
};

#endif

