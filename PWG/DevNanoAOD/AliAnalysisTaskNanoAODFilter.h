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

class AliNanoAODCustomSetter;
class AliNanoFilterNormalisation;

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisCuts.h"
#include "AliNanoAODReplicator.h"
#include "AliEventCuts.h"
#include "AliNanoAODTrack.h"
#include "AliPID.h"
#include <list>

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
  void  AddFilteredAOD(const char* aodfilename, const char* title);

  void  AddEvtCuts     (AliAnalysisCuts * var           ) { fEvtCuts.push_back(var);}
  void  SetTrkCuts     (AliAnalysisCuts * var           ) { fReplicator->SetTrackCuts(var); if (fSaveCutsFlag) fQAOutput->Add(var);}
  void  AddSetter      (AliNanoAODCustomSetter * var    ) { fReplicator->AddCustomSetter(var); }

  void  SetVarListTrack(TString var                     ) { fReplicator->SetVarListTrack(var);}
  void  AddPIDField(AliNanoAODTrack::ENanoPIDResponse response, AliPID::EParticleType particle);
  void  SetVarListHeader(TString var                    ) { fReplicator->SetVarListHeader(var);}
  void  SetVarFiredTriggerClasses (TString var          ) { fReplicator->SetVarListHeaderTC(var);}
  void  SaveVzero(Bool_t var)                             { fReplicator->SetSaveVzero(var); }
  void  SaveZDC(Bool_t var)                               { fReplicator->SetSaveZDC(var); }
  void  SaveV0s(Bool_t var, AliAnalysisCuts* v0Cuts = 0)  { fReplicator->SetSaveV0s(var); fReplicator->SetV0Cuts(v0Cuts); if (fSaveCutsFlag && v0Cuts) fQAOutput->Add(v0Cuts); }
  void  SaveCascades(Bool_t var, AliAnalysisCuts* cuts = 0) { fReplicator->SetSaveCascades(var); fReplicator->SetCascadeCuts(cuts); if (fSaveCutsFlag && cuts) fQAOutput->Add(cuts); }
  void  SaveConversionPhotons(Bool_t var, AliAnalysisCuts* cuts = 0) { fReplicator->SetSaveConversionPhotons(var); fReplicator->SetConversionPhotonCuts(cuts); if (fSaveCutsFlag && cuts) fQAOutput->Add(cuts); }
  void  SaveConversionPhotonsFromDelta(Bool_t var, TString name, AliAnalysisCuts* cuts = 0) { fReplicator->SetSaveConversionPhotons(var); fReplicator->SetPhotonDeltaBranchName(name); fReplicator->SetConversionPhotonCuts(cuts); if (fSaveCutsFlag && cuts) fQAOutput->Add(cuts); }
  void  FilterMCStack(AliAnalysisCuts* cuts = nullptr) { fReplicator->SetMCParticleCuts(cuts); if (fSaveCutsFlag && cuts) fQAOutput->Add(cuts); }
  
  AliNanoAODReplicator* GetReplicator() { return fReplicator; }

  void SetInputArrayName(TString name) {fInputArrayName=name;}
  void SetOutputArrayName(TString name) {fOutputArrayName=name;}

  bool         fUseAliEventCuts;
  AliEventCuts fEventCuts;


protected:
  Int_t fMCMode; // true if processing monte carlo. if > 1 not all MC particles are filtered
  AliNanoAODReplicator* fReplicator; // replicator

  std::list<AliAnalysisCuts*> fEvtCuts; // Event cuts

  TList* fQAOutput;     //  Output list for cut objects and QA
  Bool_t fSaveCutsFlag; // If true, the event and track cuts are saved to disk. Can only be set in the constructor.

  TString fInputArrayName; // name of TObjectArray of Tracks
  TString fOutputArrayName; // name of TObjectArray of AliNanoAODTracks

  AliNanoFilterNormalisation* fNormalisation;          //!<! Normalisation object

  AliAnalysisTaskNanoAODFilter(const AliAnalysisTaskNanoAODFilter&); // not implemented
  AliAnalysisTaskNanoAODFilter& operator=(const AliAnalysisTaskNanoAODFilter&); // not implemented

  ClassDef(AliAnalysisTaskNanoAODFilter, 9); // Nano AOD Filter Task
};

#endif

