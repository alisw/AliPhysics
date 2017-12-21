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
  Bool_t IsSelectedTrigger(AliAODEvent *fAOD, Bool_t fillHisto,UShort_t &physicsSelectionMask,UShort_t &L0TriggerInput);
  virtual void NotifyRun();
private:
  
  AliAnalysisTaskLMREventFilter(const AliAnalysisTaskLMREventFilter&); // not implemented
  AliAnalysisTaskLMREventFilter& operator=(const AliAnalysisTaskLMREventFilter&);// not implemented
  

  AliMuonTrackCuts *fMuonTrackCuts; 
  TTree *fEventTree;   // Tree containing events
  TList *fOutputList;  // ! Output list

  AliLMREvent *fAliLMREvent;

  TH1D *fhTriggers;
  TH1D *fhBeamType;
  TH1D *fhL0TriggerInputMLL;
  TH1D *fhL0TriggerInputMUL;
  TH1D *fhL0TriggerInputMSL;
  TH1D *fhL0TriggerInputTVX;
  TH2D *fhNMu;

  Int_t fL0TriggerInputMLL;
  Int_t fL0TriggerInputMUL;
  Int_t fL0TriggerInputMSL;
  Int_t fL0TriggerInputTVX;
  Int_t fNTrigClass;
  TString fTriggerClasses[13];
  Int_t fminContributorsPileUp;
  ClassDef(AliAnalysisTaskLMREventFilter, 1) // example of analysis

};

//====================================================================================================================================================

#endif
