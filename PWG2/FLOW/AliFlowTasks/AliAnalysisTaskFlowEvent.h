/////////////////////////////////////////////////////
// AliAnalysisTaskFlowEvent:
// analysis task to fill the flow event 
// and make it available to the flow analysis methods.
//////////////////////////////////////////////////////

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id: $ */

#ifndef AliAnalysisTaskFlowEvent_H
#define AliAnalysisTaskFlowEvent_H

class AliCFManager;
class AliFlowEventSimpleMaker;
class TList;
class TRandom3;
class AliAnalysisTaskSE;
class TString;

class AliAnalysisTaskFlowEvent : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskFlowEvent();
  AliAnalysisTaskFlowEvent(const char *name,TString RPtype,Bool_t QAon,UInt_t=666);
  virtual ~AliAnalysisTaskFlowEvent();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetAnalysisType(TString type) { this->fAnalysisType = type; }
  TString GetAnalysisType() const    { return this->fAnalysisType; }

  void SetRPType(TString rptype) { this->fRPType = rptype; }
  TString GetRPType() const    { return this->fRPType; }

  void    SetMinMult(Int_t multmin)    {this->fMinMult = multmin; }
  Int_t   GetMinMult() const           {return this->fMinMult; }
  void    SetMaxMult(Int_t multmax)    {this->fMaxMult = multmax; }
  Int_t   GetMaxMult() const           {return this->fMaxMult; }

  void SetSubeventEtaRange(Double_t minA, Double_t maxA, Double_t minB, Double_t maxB)
    {this->fMinA = minA; this->fMaxA = maxA; this->fMinB = minB; this->fMaxB = maxB; }
  Double_t GetMinA() const {return this->fMinA;}
  Double_t GetMaxA() const {return this->fMaxA;}
  Double_t GetMinB() const {return this->fMinB;}
  Double_t GetMaxB() const {return this->fMaxB;}
  
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

  // setters for common constants
  void SetNbinsMult( Int_t i ) { fNbinsMult = i; }
  void SetNbinsPt( Int_t i )   { fNbinsPt = i; }
  void SetNbinsPhi( Int_t i )  { fNbinsPhi = i; }
  void SetNbinsEta( Int_t i )  { fNbinsEta = i; }
  void SetNbinsQ( Int_t i )    { fNbinsQ = i; }
   
  void SetMultMin( Double_t i ) { fMultMin = i; }
  void SetMultMax( Double_t i ) { fMultMax = i; }
  void SetPtMin( Double_t i )   { fPtMin = i; }
  void SetPtMax( Double_t i )   { fPtMax = i; }
  void SetPhiMin( Double_t i )  { fPhiMin = i; }
  void SetPhiMax( Double_t i )  { fPhiMax = i; }
  void SetEtaMin( Double_t i )  { fEtaMin = i; }
  void SetEtaMax( Double_t i )  { fEtaMax = i; }
  void SetQMin( Double_t i )    { fQMin = i; }
  void SetQMax( Double_t i )    { fQMax = i; }
  // end setters common constants

  // setters for adding by hand flow values (afterburner)
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
  // end setters afterburner

 private:

  AliAnalysisTaskFlowEvent(const AliAnalysisTaskFlowEvent& aAnalysisTask);
  AliAnalysisTaskFlowEvent& operator=(const AliAnalysisTaskFlowEvent& aAnalysisTask); 

  //  TFile*        fOutputFile;    // temporary output file for testing
  //  AliESDEvent*  fESD;           // ESD object
  //  AliAODEvent*  fAOD;           // AOD object
  TString       fAnalysisType;      // can be MC, ESD or AOD
  TString       fRPType;            // can be Global or Tracklet or FMD
  AliCFManager* fCFManager1;        // correction framework manager
  AliCFManager* fCFManager2;        // correction framework manager
  TList*        fQAInt;             // QA histogram list
  TList*        fQADiff;            // QA histogram list
  Int_t         fMinMult;           // Minimum multiplicity from tracks selected using CORRFW
  Int_t         fMaxMult;           // Maximum multiplicity from tracks selected using CORRFW 
  Double_t      fMinA;              // Minimum of eta range for subevent A
  Double_t      fMaxA;              // Maximum of eta range for subevent A
  Double_t      fMinB;              // Minimum of eta range for subevent B
  Double_t      fMaxB;              // Maximum of eta range for subevent B

  Bool_t fQA;                       // flag to set the filling of the QA hostograms

  // setters for common constants
  //histogram sizes
  Int_t  fNbinsMult; // histogram size
  Int_t  fNbinsPt;   // histogram size
  Int_t  fNbinsPhi;  // histogram size
  Int_t  fNbinsEta;  // histogram size
  Int_t  fNbinsQ;    // histogram size
 
  // Histograms limits
  Double_t  fMultMin;  // histogram limit 
  Double_t  fMultMax;  // histogram limit
  Double_t  fPtMin;    // histogram limit
  Double_t  fPtMax;    // histogram limit
  Double_t  fPhiMin;   // histogram limit
  Double_t  fPhiMax;   // histogram limit
  Double_t  fEtaMin;   // histogram limit
  Double_t  fEtaMax;   // histogram limit
  Double_t  fQMin;     // histogram limit
  Double_t  fQMax;     // histogram limit
  // end common constants

  // values afterburner
  Double_t  fMCReactionPlaneAngle;     // the angle of the reaction plane from the MC truth
  Int_t     fCount;                    // counter for the number of events processed
  Int_t     fNoOfLoops;                // number of times to use the same particle (nonflow) 
  Double_t  fEllipticFlowValue;        // Add Flow. Must be in range [0,1].
  Double_t  fSigmaEllipticFlowValue;   // Sigma Flow (Gaussian). Must be in range [0,1].
  Int_t     fMultiplicityOfEvent;      // Set maximal multiplicity.
  Int_t     fSigmaMultiplicityOfEvent; // Sigma multiplicity (Gaussian).
    
  TRandom3* fMyTRandom3;     // TRandom3 generator
  Bool_t    fbAfterburnerOn;
  // end afterburner
  
  ClassDef(AliAnalysisTaskFlowEvent, 1); // example of analysis
};

#endif

