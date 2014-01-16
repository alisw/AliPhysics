/////////////////////////////////////////////////////
// AliAnalysisTaskFlowCascade:
// Analysis task to select Xi and Omega candidates for flow analysis.
// Author: 
//////////////////////////////////////////////////////

/* Copyright(c) 1998-1999, ALICExperiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id: $ */

#ifndef AliAnalysisTaskFlowCascade_H
#define AliAnalysisTaskFlowCascade_H

#include "AliAnalysisTaskSE.h"

class AliESDtrackCuts;
class AliFlowEventCuts;
class TList;
class TH1D;
class TH2D;
class TH3D;
class AliFlowCandidateTrack;
class TVector3;
class AliPIDResponse;
class AliFlowEvent;
class AliAnalysisTaskFlowCascade : public AliAnalysisTaskSE {
  private:
  AliAnalysisTaskFlowCascade(const AliAnalysisTaskFlowCascade& analysisTask);
  AliAnalysisTaskFlowCascade& 
    operator=(const AliAnalysisTaskFlowCascade& analysisTask);

  //Progate to the primary vertex
  void Propagate(Double_t vv[3], Double_t x[3], Double_t p[3], Double_t bz, 
		 Short_t sign);

  void AddQAEvents();
  void AddQACandidates();
  void ReadFromESDv0(AliESDEvent *fESD);
  void ReadFromAODv0(AliAODEvent *fAOD);

  void AddCandidates();

  void MakeTrack( double mass, double pt, 
		  double phi, double eta, 
		  int iid, int jid, int kid );
  
  //  double fMinCent, fMaxCent;   //
  Int_t fSpecie; //0 for Xi and 1 for Omega
  Double_t fCascadeCuts[8]; // cuts for cascade selection
  Int_t fMassBins;   // to configure FLOWCOMMON                               
  Double_t fMinMass; // to configure FLOWCOMMON                               
  Double_t fMaxMass; // to configure FLOWCOMMON
  
  AliFlowEventCuts *fCutsEvent; // event cuts 
  AliFlowTrackCuts *fCutsRPTPC;    // cuts for RPs
  AliFlowTrackCuts *fCutsRPVZE;    // cuts for RPs
  AliFlowTrackCuts *fCutsPOI; // null cuts for POI
  AliFlowTrackCuts *fCutsDau; // cuts for daughters
  AliPIDResponse *fPIDResponse;
  AliFlowEvent  *fFlowEventTPC;  //flow event TPC                             
  AliFlowEvent  *fFlowEventVZE;  //flow event VZE 
  TObjArray *fCandidates; // Array of selected candidates
  TList *fQAList;               //! list for QA histos

 public:
  AliAnalysisTaskFlowCascade();
  AliAnalysisTaskFlowCascade(const char *name, 
			     AliFlowEventCuts *cutsEvent, 
			     AliFlowTrackCuts *cutsRPTPC,
			     AliFlowTrackCuts *cutsRPVZE,
			     /* AliESDtrackCuts */AliFlowTrackCuts *cutsDau);
  //void SetDebug() {fDebug = true;}
  void SetSpecie(int specie){fSpecie = specie;}
  void SetCommonConstants(Int_t massBins, Double_t minMass, Double_t maxMass);
  void SetCuts2010(int setOfCuts);
  virtual ~AliAnalysisTaskFlowCascade();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);
  virtual void NotifyRun();
  
  ClassDef(AliAnalysisTaskFlowCascade, 2);
};

#endif
