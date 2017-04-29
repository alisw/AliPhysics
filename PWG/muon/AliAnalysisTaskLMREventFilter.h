#ifndef AliAnalysisTaskLMREventFilter_H
#define AliAnalysisTaskLMREventFilter_H

#include "AliAnalysisTaskSE.h"
#include "AliMuonTrackCuts.h"
#include "AliLMREvent.h"

class TClonesArray;
class AliAODEvent;
class AliVVertex;
class AliVParticle;

//====================================================================================================================================================

class  AliAnalysisTaskLMREventFilter : public AliAnalysisTaskSE {

public:
 
  AliAnalysisTaskLMREventFilter();
  AliAnalysisTaskLMREventFilter(const char *name, AliMuonTrackCuts *cuts);

  virtual ~AliAnalysisTaskLMREventFilter();
  
  
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);
  Bool_t IsSelectedTrigger(AliAODEvent *fAOD, Bool_t fillHisto,UInt_t &evtTrigSelect,UShort_t &physicsSelectionMask);
  virtual void NotifyRun();
private:
  
  AliAnalysisTaskLMREventFilter(const AliAnalysisTaskLMREventFilter&); // not implemented
  AliAnalysisTaskLMREventFilter& operator=(const AliAnalysisTaskLMREventFilter&);// not implemented
  

  AliMuonTrackCuts *fMuonTrackCuts; 
  TTree *fEventTree;   // Tree containing events
  TList *fOutputList;  // ! Output list

  AliLMREvent *fAliLMREvent;

  TH1D *fhTriggers;
  TH1D *fhL0TriggerInputMLL;
  TH1D *fhL0TriggerInputMUL;
  TH1D *fhL0TriggerInputMSL;
  TH2D *fhNMu;

  Int_t fL0TriggerInputMLL;
  Int_t fL0TriggerInputMUL;
  Int_t fL0TriggerInputMSL;
  Int_t fNTrigClass;
  TString fTriggerClasses[23];
  UInt_t fTriggerMask[23];

  ClassDef(AliAnalysisTaskLMREventFilter, 1) // example of analysis

};

//====================================================================================================================================================

#endif
