#ifndef ALIDNDPT_H
#define ALIDNDPT_H

//------------------------------------------------------------------------------
// Abstract class for dNdPt analysis. All dNdPt components should derive from it.   
// 
// Author: J.Otwinowski 03/11/2008 
//------------------------------------------------------------------------------

class AliESDEvent; 
class AliMCEvent; 
class AliESDtrackCuts; 
class AlidNdPtEventCuts;
class AlidNdPtAcceptanceCuts;

#include "TNamed.h"
#include "TFolder.h"
#include "AliTriggerAnalysis.h"
#include "AlidNdPtHelper.h"

class AlidNdPt : public TNamed {
public:
  AlidNdPt(); 
  AlidNdPt(Char_t* name, Char_t* title);
  ~AlidNdPt();

  // Init data members
  virtual void Init() = 0;

  // Process events
  virtual void Process(AliESDEvent *const esdEvent=0, AliMCEvent *const mcEvent=0) = 0;

  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* list=0) = 0;

  // Analyse output histograms 
  virtual void Analyse() = 0;

  // Export analysed output objects to folder
  virtual TFolder *ExportToFolder(TObjArray * array=0) = 0;

  //
  AlidNdPt(const AlidNdPt&); // not implemented
  AlidNdPt& operator=(const AlidNdPt&); // not implemented

  //
  void SetEventCuts(AlidNdPtEventCuts* const cuts)              { fdNdPtEventCuts = cuts; }
  void SetAcceptanceCuts(AlidNdPtAcceptanceCuts* const cuts)    { fdNdPtAcceptanceCuts = cuts; }
  void SetTrackCuts(AliESDtrackCuts* const cuts)                { fEsdTrackCuts = cuts; }
  void SetUseMCInfo(const Bool_t info)                          { fUseMCInfo = info; }
  void SetAnalysisMode(const AlidNdPtHelper::AnalysisMode mode) { fAnalysisMode = mode; }
  void SetTrigger(const AliTriggerAnalysis::Trigger trigger)    { fTrigger = trigger; }
  void SetTriggerClass(const Char_t *triggerClass)              { fTriggerClass = triggerClass; }

  AlidNdPtEventCuts* GetEventCuts() const                       { return fdNdPtEventCuts; }
  AlidNdPtAcceptanceCuts* GetAcceptanceCuts() const             { return fdNdPtAcceptanceCuts; }
  AliESDtrackCuts* GetTrackCuts() const                         { return fEsdTrackCuts; }
  Bool_t IsUseMCInfo() const                                    { return fUseMCInfo; }
  AlidNdPtHelper::AnalysisMode GetAnalysisMode() const          { return fAnalysisMode; }
  AliTriggerAnalysis::Trigger GetTrigger() const                { return fTrigger; }
  const Char_t* GetTriggerClass() const                         { return fTriggerClass; }

private:

  AlidNdPtEventCuts      *fdNdPtEventCuts;      // event cuts
  AlidNdPtAcceptanceCuts *fdNdPtAcceptanceCuts; // acceptance cuts
  AliESDtrackCuts *fEsdTrackCuts;               // esd track cuts

  Bool_t fUseMCInfo;                            // use MC information
  AlidNdPtHelper::AnalysisMode fAnalysisMode;   // analysis mode TPC only, TPC + ITS
  AliTriggerAnalysis::Trigger fTrigger;         // trigger definition MB1, MB2 ...
  const Char_t * fTriggerClass;                 // trigger class

  ClassDef(AlidNdPt,1);
};

#endif
