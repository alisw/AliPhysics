/*
 * **********************************************************
 * Virtual class for processing trees of AliReducedEventInfo
 * Authors: Ionut-Cristian Arsene, iarsene@cern.ch, i.c.arsene@fys.uio.no
 *                Jacobus Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
 * Creation date: 2015/10/01
 *********************************************************
 */


#ifndef ALIREDUCEDANALYSISTASKSE_H
#define ALIREDUCEDANALYSISTASKSE_H

#include <TNamed.h>
#include <TList.h>

#include "AliReducedVarManager.h"
#include "AliReducedInfoCut.h"

class AliReducedBaseEvent;
class AliHistogramManager;
class AliMixingHandler; 

//________________________________________________________________
class AliReducedAnalysisTaskSE : public TNamed {
  
public:
  AliReducedAnalysisTaskSE();
  AliReducedAnalysisTaskSE(const Char_t* name, const Char_t* title);
  virtual ~AliReducedAnalysisTaskSE();
  
  // Virtual functions, to be implemented in the inheriting analysis classes
  // initialization (typically called in AliAnalysisTask::UserCreateOutputObjects())
  virtual void Init();
  // process a given event (typically called in AliAnalysisTask::UserExec())
  virtual void Process();
  // finish, to be executed after all events were processed
  virtual void Finish();
  // add output objects;
  
  // setters
  void SetHistogramManager(AliHistogramManager* histos) {fHistosManager = histos;}
  void SetMixingHandler(AliMixingHandler* mixing) {fMixingHandler = mixing;}
  void SetEvent(AliReducedBaseEvent* event) {fEvent = event;}
  void AddEventCut(AliReducedInfoCut* cut) {fEventCuts.Add(cut);}
  void AddTrackCut(AliReducedInfoCut* cut) {fTrackCuts.Add(cut);}
  void AddPairCut(AliReducedInfoCut* cut) {fPairCuts.Add(cut);}

  //void SetInputObject(Int_t slot, TObject* obj) {fInputObjects[slot]=obj;}

  
  // getters
  AliHistogramManager* GetHistogramManager() const {return fHistosManager;}
  AliMixingHandler* GetMixingHandler() const {return fMixingHandler;}
  AliReducedBaseEvent* GetEvent() const {return fEvent;}
  
protected:
  AliReducedAnalysisTaskSE(const AliReducedAnalysisTaskSE& handler);             
  AliReducedAnalysisTaskSE& operator=(const AliReducedAnalysisTaskSE& handler);      
   
  AliHistogramManager* fHistosManager;   // Histogram manager
  AliMixingHandler* fMixingHandler;      // Mixed event handler
  
  AliReducedBaseEvent* fEvent;           //! current event to be processed
  Float_t fValues[AliReducedVarManager::kNVars];   // array of values to hold information for histograms
  
  TList fEventCuts;     // array of event cuts
  TList fTrackCuts;     // array of track cuts
  TList fPairCuts;      // array of pair cuts

  //TObject* fInputObjects[100];
  //TObject* fOutputObjects[100];
  
  virtual Bool_t IsEventSelected(AliReducedBaseEvent* event);
  virtual Bool_t IsTrackSelected(AliReducedBaseTrack* track);
  virtual Bool_t IsPairSelected(AliReducedBaseTrack* pair);
  
  ClassDef(AliReducedAnalysisTaskSE,3);
};

#endif
