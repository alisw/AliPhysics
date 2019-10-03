/////////////////////////////////////////////////////
// AliAnalysisTaskFlowEvent:
// analysis task to fill the flow event 
// and make it available to the flow analysis methods.
//////////////////////////////////////////////////////

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id: $ */

#ifndef ALIANALYSISTASKFLOWEVENT_H
#define ALIANALYSISTASKFLOWEVENT_H

#include "AliFlowTrackSimple.h"

class AliCFManager;
class AliFlowEventCuts;
class AliFlowTrackCuts;
class AliFlowEventSimpleMaker;
class AliFlowEvent;
class TList;
class TF1;
class TRandom3;
class AliAnalysisTaskSE;
class TString;
class AliESDpid;

class AliAnalysisTaskFlowEvent : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskFlowEvent();
  AliAnalysisTaskFlowEvent(const char *name, TString RPtype = "", Bool_t QAon = kFALSE, UInt_t seed=666, Bool_t bCandidates=kFALSE);
  virtual ~AliAnalysisTaskFlowEvent();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual void   NotifyRun();

  void    SetAnalysisType(TString type) { this->fAnalysisType = type; }
  TString GetAnalysisType() const       { return this->fAnalysisType; }

  void    SetRPType(TString rptype) { this->fRPType = rptype; }
  TString GetRPType() const         { return this->fRPType; }

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
  
  void DefineDeadZone( Double_t etaMin, Double_t etaMax, Double_t phiMin, Double_t phiMax )
  {this->fExcludedEtaMin = etaMin; this->fExcludedEtaMax = etaMax; 
    this->fExcludedPhiMin = phiMin; this->fExcludedPhiMax = phiMax; }

  void          SetCutsEvent(AliFlowEventCuts* cutsEvent) {fCutsEvent=cutsEvent;}
  AliFlowEventCuts* GetCutsEvent() const {return fCutsEvent;}
  void          SetCutsRP(AliFlowTrackCuts* cutsRP) {fCutContainer->AddAt(cutsRP,0); fCutsRP=cutsRP; cutsRP->SetPOItype(0); }
  AliFlowTrackCuts* GetCutsRP() const {return fCutsRP;} //to be reimplemented
  void          SetCutsPOI(AliFlowTrackCuts* cutsPOI) {fCutContainer->AddAt(cutsPOI,1); fCutsPOI=cutsPOI; cutsPOI->SetPOItype(1); }
  AliFlowTrackCuts* GetCutsPOI() const {return fCutsPOI;} //to be reimplemented

  void          SetCFManager1(AliCFManager* cfmgr) {this->fCFManager1 = cfmgr; } 
  AliCFManager* GetCFManager1() const {return this->fCFManager1; }
  void          SetCFManager2(AliCFManager* cfmgr) {this->fCFManager2 = cfmgr; } 
  AliCFManager* GetCFManager2() const       {return this->fCFManager2; }
  TList*        GetQAList()      const      {return fQAList; }
  void          SetQAOn(Bool_t kt)        {fQAon = kt; }
  Bool_t        GetQAOn()   const         {return fQAon; }

  void          SetShuffleTracks(Bool_t b)  {fShuffleTracks=b;}

  void   SetPassMCeventToCutsObject(Bool_t passMC){this->fPassMCeventToCutsObject = passMC;}


  // setters for common constants
  void SetNbinsMult( Int_t i ) { fNbinsMult = i; }
  void SetNbinsPt( Int_t i )   { fNbinsPt = i; }
  void SetNbinsPhi( Int_t i )  { fNbinsPhi = i; }
  void SetNbinsEta( Int_t i )  { fNbinsEta = i; }
  void SetNbinsQ( Int_t i )    { fNbinsQ = i; }
  void SetNbinsMass( Int_t i ) { fNbinsMass = i; }
   
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
  void SetMassMin( Double_t i ) { fMassMin = i; }
  void SetMassMax( Double_t i ) { fMassMax = i; }
  void SetHistWeightvsPhiMin( Double_t i ) {fHistWeightvsPhiMin=i;}
  void SetHistWeightvsPhiMax( Double_t i ) {fHistWeightvsPhiMax=i;}
  // end setters common constants

  // setters for adding by hand flow values (afterburner)
  void SetAfterburnerOn(Bool_t b=kTRUE) {fAfterburnerOn=b;}
  void SetNonFlowNumberOfTrackClones(Int_t n) {fNonFlowNumberOfTrackClones=n;}
  void SetPtDifferentialV2( TF1 *gPtV2) {
    fDifferentialV2 = gPtV2;}
  void SetFlow( Double_t v1, Double_t v2, Double_t v3=0.0, Double_t v4=0.0, Double_t v5=0.0)
               {fV1=v1;fV2=v2;fV3=v3;fV4=v4;fV5=v5;}
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
  AliFlowEventCuts* fCutsEvent;     //event cuts
  AliFlowTrackCuts* fCutsRP;        //cuts for RPs
  AliFlowTrackCuts* fCutsPOI;       //cuts for POIs
  TList*            fCutContainer;  //contains the cut objects
  TList*        fQAList;             // QA histogram list
  Int_t         fMinMult;           // Minimum multiplicity from tracks selected using CORRFW
  Int_t         fMaxMult;           // Maximum multiplicity from tracks selected using CORRFW 
  Double_t      fMinA;              // Minimum of eta range for subevent A
  Double_t      fMaxA;              // Maximum of eta range for subevent A
  Double_t      fMinB;              // Minimum of eta range for subevent B
  Double_t      fMaxB;              // Maximum of eta range for subevent B

  Bool_t fQAon;                     // flag to set the filling of the QA hostograms
  Bool_t fLoadCandidates;           // true if reciving candidates collection
  Bool_t fPassMCeventToCutsObject;  // defaut: true 

  // setters for common constants
  //histogram sizes
  Int_t  fNbinsMult; // histogram size
  Int_t  fNbinsPt;   // histogram size
  Int_t  fNbinsPhi;  // histogram size
  Int_t  fNbinsEta;  // histogram size
  Int_t  fNbinsQ;    // histogram size
  Int_t  fNbinsMass; // histogram size
 
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
  Double_t  fMassMin;  // histogram limit
  Double_t  fMassMax;  // histogram limit
  Double_t fHistWeightvsPhiMin; //histogram limit
  Double_t fHistWeightvsPhiMax; //histogram limit
  // end common constants

  // Excluding a range
  Double_t  fExcludedEtaMin;  // excluded region limit 
  Double_t  fExcludedEtaMax;  // excluded region limit 
  Double_t  fExcludedPhiMin;  // excluded region limit 
  Double_t  fExcludedPhiMax;  // excluded region limit 
  // End of excluding a range

  // values afterburner
  Bool_t    fAfterburnerOn;              // do we afterburn?
  Int_t     fNonFlowNumberOfTrackClones; // number of times to clone the particles (nonflow) 
  Double_t  fV1;        // Add Flow. Must be in range [0,0.5].
  Double_t  fV2;        // Add Flow. Must be in range [0,0.5].
  Double_t  fV3;        // Add Flow. Must be in range [0,0.5].
  Double_t  fV4;        // Add Flow. Must be in range [0,0.5].
  Double_t  fV5;        // Add Flow. Must be in range [0,0.5].
  TF1 *fDifferentialV2; // pt-differential v2

  AliFlowEvent* fFlowEvent; //flowevent
  Bool_t fShuffleTracks;    //serve the tracks shuffled
    
  TRandom3* fMyTRandom3;     // TRandom3 generator
  // end afterburner
  
  ClassDef(AliAnalysisTaskFlowEvent, 1); // example of analysis
};

#endif

