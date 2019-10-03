#ifndef ALIANALYSISTASKUPCFILTER_H
#define ALIANALYSISTASKUPCFILTER_H

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

class AliAnalysisTaskUpcFilter : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskUpcFilter(const char *name="AliAnalysisTaskUpcFilter");
  virtual ~AliAnalysisTaskUpcFilter();

  void SetIsESD(Bool_t isESD) {fIsESD = isESD;}
  void SetIsMC(Bool_t isMC) {fIsMC = isMC;}
  void SetFillSPD(Bool_t fill=kTRUE) {fFillSPD = fill;}
  void SetAllTrg(Bool_t set);
  void SetTrgClass(Int_t idx, Bool_t set);
  void SetMuonTrackCutsPassName(const char *passname) {fMuonCutsPassName->Clear(); fMuonCutsPassName->SetString(passname);}
  virtual void UserCreateOutputObjects();
  virtual void NotifyRun();
  virtual void UserExec(Option_t *option);
  Bool_t RunAOD();
  void RunAODMC(TClonesArray *arrayMC, AliAODMCHeader *headerMC);
  Bool_t RunESD();
  void RunESDMC();
  virtual void Terminate(Option_t *);

 private:
  AliAnalysisTaskUpcFilter(const AliAnalysisTaskUpcFilter &o); // not implemented
  AliAnalysisTaskUpcFilter &operator=(const AliAnalysisTaskUpcFilter &o); // not implemented

  Bool_t fIsESD; // analysis type, ESD / AOD
  Bool_t fIsMC; // mc or data selection
  Bool_t fFillSPD; // fill SPD fired FO chips

  static const Int_t fgkNtrg = 51; // number of trigger classes

  Bool_t fTrgMask[fgkNtrg]; // allows to mask or un-mask the selected trigger classes

  AliMuonTrackCuts *fMuonCuts; // class for muon track cuts, used for pDCA
  AliTriggerAnalysis *fTriggerAna; // class for trigger analysis, used for fired SPD FO
  AliESDtrackCuts **fCutsList; // array of pointers to filtering task for ESD tracks
  AliPIDResponse *fPIDResponse;     // PID response object
  TObjString *fMuonCutsPassName; //-> name of pass used by AliMuonTrackCuts

  TList *fHistList; // list of output histograms
  TH1I *fCounter; // analysis counter
  TH2I *fTriggerCounter; // counter of triggers per run
  TH1I *fMuonCounter;
  AliUPCEvent *fUPCEvent; // output UPC event
  TTree *fUPCTree; // output tree

  enum EvtCount{ kAna=1, kTrg, kSpecific, kPass1, kPass2, kPassX, kWritten, kAOD, kMunTrack, kCenTrack, kESD, kPidErr };
  enum MuonCount{kMunAll=1, kMunRabs, kMunEta, kMunPDCA};

  ClassDef(AliAnalysisTaskUpcFilter, 1); 
};

#endif









