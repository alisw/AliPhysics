#ifndef ALIANALYSISTASKUPCFILTERSEMIFORWARD_H
#define ALIANALYSISTASKUPCFILTERSEMIFORWARD_H

// task for upc semiforward filter
// creates upc event from esd or aod
//
// jaroslav.adam@cern.ch

#include "AliAnalysisTaskSE.h"

class AliUPCEvent;
class AliAODMCHeader;
class AliMuonTrackCuts;
class AliESDtrackCuts;
class AliPIDResponse;
class AliTriggerAnalysis;

class AliAnalysisTaskUpcFilterSemiforward : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskUpcFilterSemiforward(const char *name="AliAnalysisTaskUpcFilterSemiforward");
  virtual ~AliAnalysisTaskUpcFilterSemiforward();

  void SetIsESD(Bool_t isESD) {fIsESD = isESD;}
  void SetIsMC(Bool_t isMC) {fIsMC = isMC;}
  virtual void UserCreateOutputObjects();
  virtual void NotifyRun();
  virtual void UserExec(Option_t *option);
  Bool_t RunAOD();
  void RunAODMC(TClonesArray *arrayMC, AliAODMCHeader *headerMC);
  Bool_t RunESD();
  void RunESDMC();
  virtual void Terminate(Option_t *);

 private:
  AliAnalysisTaskUpcFilterSemiforward(const AliAnalysisTaskUpcFilterSemiforward &o); // not implemented
  AliAnalysisTaskUpcFilterSemiforward &operator=(const AliAnalysisTaskUpcFilterSemiforward &o); // not implemented

  Bool_t fIsESD; // analysis type, ESD / AOD
  Bool_t fIsMC; // mc or data selection

  AliMuonTrackCuts *fMuonCuts; // class for muon track cuts, used for pDCA
  AliTriggerAnalysis *fTriggerAna; // class for trigger analysis, used for fired SPD FO
  AliESDtrackCuts **fCutsList; // array of pointers to filtering task for ESD tracks
  AliPIDResponse *fPIDResponse;     // PID response object

  TList *fHistList; // list of output histograms
  TH1I *fCounter; // analysis counter
  TH2I *fTriggerCounter; // counter of triggers per run
  AliUPCEvent *fUPCEvent; // output UPC event
  TTree *fUPCTree; // output tree

  ClassDef(AliAnalysisTaskUpcFilterSemiforward, 1); 
};

#endif









