#ifndef AliAnalysisTaskSingleMuESD_cxx
#define AliAnalysisTaskSingleMuESD_cxx

/* $Id$ */ 

// analysis task for single muon analysis from the ESD events
// Authors: Bogdan Vulpescu, Nicole Bastid

class TNtuple;
class TH1F;
class AliESDEvent;
class TString;

#include "AliAnalysisTask.h"

class AliAnalysisTaskSingleMuESD : public AliAnalysisTask {
 public:
  AliAnalysisTaskSingleMuESD(const char *name = "AliAnalysisTaskSingleMuESD");
  virtual ~AliAnalysisTaskSingleMuESD() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  // force values for the trigger mask, for version where 
  // GetFiredTriggerClasses() does not work; "MUON" or "p-p"
  void SetTriggerType(const Char_t *trig);

 private:
  AliESDEvent *fESD;      // ESD object
  TNtuple     *fNtuple;   // ntuple with track variables
  TH1F        *fTrigger;  // ESD mask from CTP (trigger class)
  Int_t       fMaskTrig1MuL;  // trigger mask: single mu low  pt
  Int_t       fMaskTrig1MuH;  // trigger mask: single mu high pt
  Int_t       fMaskTrigUSL;   // trigger mask: un-like sign mu low  pt
  Int_t       fMaskTrigLSL;   // trigger mask: like sign mu low  pt
  Int_t       fMaskTrigUSH;   // trigger mask: un-like sign mu high pt
  Int_t       fMaskTrigLSH;   // trigger mask: like sign mu high pt
  Bool_t      fTriggerType;   // force the masks values

  AliAnalysisTaskSingleMuESD(const AliAnalysisTaskSingleMuESD&); // not implemented
  AliAnalysisTaskSingleMuESD& operator=(const AliAnalysisTaskSingleMuESD&); // not implemented
 
  void GetEffFitted(Double_t eta, Double_t pt, Double_t rEff[2]);
 
  ClassDef(AliAnalysisTaskSingleMuESD, 1); // single muon analysis from ESD
};

#endif
