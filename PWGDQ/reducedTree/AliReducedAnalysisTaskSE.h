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

#include <TObject.h>
#include <TString.h> 

#include "AliReducedVarManager.h"
#include "AliHistogramManager.h"
#include "AliReducedBaseEvent.h"

//________________________________________________________________
class AliReducedAnalysisTaskSE : public TObject {
  
public:
   enum ETreeWritingOptions {
      kBaseEventsWithBaseTracks=0,    // write basic event info and basic track info
      kBaseEventsWithFullTracks,           // write basic event info and full track info
      kFullEventsWithBaseTracks,           // write full event info and base track info
      kFullEventsWithFullTracks              // write full event info and full track info
   };
   
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
  
  void InitFilteredTree();
  
  // setters
  void SetEvent(AliReducedBaseEvent* event) {fEvent = event;}
  
  void SetFilteredTreeWritingOption(Int_t option)         {fFilteredTreeWritingOption = option;}
  void SetFilteredTreeActiveBranch(TString b)   {fActiveBranches+=b+";";}
  void SetFilteredTreeInactiveBranch(TString b) {fInactiveBranches+=b+";";}
  void SetProcessMC(Bool_t option=kTRUE) {fProcessMCInfo=option;}
  
  // getters
  virtual AliHistogramManager* GetHistogramManager() const = 0;
  AliReducedBaseEvent* GetEvent() const {return fEvent;}
  TTree* GetFilteredTree() {return fFilteredTree;}
  Int_t GetFilteredTreeWritingOption() const {return fFilteredTreeWritingOption;}
  Bool_t ProcessMC() const {return fProcessMCInfo;}
  
protected:
  AliReducedAnalysisTaskSE(const AliReducedAnalysisTaskSE& task);             
  AliReducedAnalysisTaskSE& operator=(const AliReducedAnalysisTaskSE& task);      
  
  TString fName;             // name
  TString fTitle;                // title
    
  AliReducedBaseEvent* fEvent;           //! current event to be processed
  Float_t fValues[AliReducedVarManager::kNVars];   // array of values to hold information for histograms
  
  TTree *fFilteredTree;                          //! tree to hold filtered reduced events
  TString fActiveBranches;                   // list of active output tree branches 
  TString fInactiveBranches;                // list of inactive output tree branches
  AliReducedBaseEvent *fFilteredEvent;     // filtered reduced event
  Int_t     fFilteredTreeWritingOption;     // one of the options described by ETreeWritingOptions
  Bool_t  fProcessMCInfo;
  
  ULong_t fEventCounter;   // event counter
  
  ClassDef(AliReducedAnalysisTaskSE, 4)
};

#endif
