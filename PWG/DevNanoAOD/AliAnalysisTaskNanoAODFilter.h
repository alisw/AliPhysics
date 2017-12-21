/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
/* AliAnalysisTaskNanoAODFilter.h
 *
 * Template task producing a P_t spectrum and pseudorapidity distribution.
 * Includes explanations of physics and primary track selections
 *
 * Based on tutorial example from offline pages
 * Edited by Arvinder Palaha
 */
#ifndef ALIANALYSISTASKNANOAODESEFILTER_H
#define ALIANALYSISTASKNANOAODESEFILTER_H

class TH1F;
class TList;
class AliESDtrackCuts;
class AliAnalysisCuts;
class AliNanoAODReplicator;
class AliNanoAODCustomSetter;

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class AliAnalysisTaskNanoAODFilter : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskNanoAODFilter();
  AliAnalysisTaskNanoAODFilter(const char *name, Bool_t saveCutsFlag);
  virtual ~AliAnalysisTaskNanoAODFilter();
  
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void Terminate(Option_t *);
  virtual void FinishTaskOutput() ;

  Int_t GetMCMode() { return fMCMode; }
  void  SetMCMode (Int_t var) { fMCMode = var;}
  void AddFilteredAOD(const char* aodfilename, const char* title);

  AliAnalysisCuts *           GetEvtCuts() { return fEvtCuts; }
  AliAnalysisCuts *           GetTrkCuts() { return fTrkCuts; }
  AliNanoAODCustomSetter *    GetSetter()  { return fSetter; }
  TString                     GetVarList() { return fVarList; }
  TString                     GetVarListHead() { return fVarListHead; }
  Bool_t                      GetSaveCutsFlag() { return fSaveCutsFlag; }

  void  SetEvtCuts     (AliAnalysisCuts * var           ) { fEvtCuts = var;}
  void  SetTrkCuts     (AliAnalysisCuts * var           ) { fTrkCuts = var;}
  void  SetSetter      (AliNanoAODCustomSetter * var    ) { fSetter = var;}
  void  SetVarList     (TString var                     ) { fVarList = var;}
  void  SetVarListHead (TString var                     ) { fVarListHead = var;}
  void  SetVarFiredTriggerClasses (TString var          ) { fVarListHeader_fTC = var;}
  void  ReplicatorSaveVzero(Bool_t var ) {fSaveVzero=var;}
  void  ReplicatorSaveAODZDC(Bool_t var ) {fSaveAODZDC=var;}

  void SetInputArrayName(TString name) {fInputArrayName=name;}
  void SetOutputArrayName(TString name) {fOutputArrayName=name;}

private:
  Int_t fMCMode; // true if processing monte carlo. if > 1 not all MC particles are filtered
  AliNanoAODReplicator* fTrkrep       ; // ! replicator
  TString                  fVarList      ; // List of variables to be added to the special track
  TString fVarListHead; // List of variables to be added to the special header
  TString fVarListHeader_fTC; // List of fired trigger classes used in the NanoAOD generation if the fired Trigger Classes are safed
  AliAnalysisCuts * fEvtCuts; // Event cuts
  AliAnalysisCuts * fTrkCuts; // Track cuts

  AliNanoAODCustomSetter * fSetter; // setter for custom variables
  
  Bool_t fSaveCutsFlag; // If true, the event and track cuts are saved to disk. Can only be set in the constructor.
  Bool_t fSaveAODZDC;  // if kTRUE AliAODZDC will be saved in AliAODEvent
  Bool_t fSaveVzero; // if kTRUE AliAODVZERO will be saved in AliAODEvent

  TString fInputArrayName; // name of TObjectArray of Tracks
  TString fOutputArrayName; // name of TObjectArray of AliNanoAODTracks

  AliAnalysisTaskNanoAODFilter(const AliAnalysisTaskNanoAODFilter&); // not implemented
  AliAnalysisTaskNanoAODFilter& operator=(const AliAnalysisTaskNanoAODFilter&); // not implemented

  ClassDef(AliAnalysisTaskNanoAODFilter, 4); // example of analysis
};

#endif

