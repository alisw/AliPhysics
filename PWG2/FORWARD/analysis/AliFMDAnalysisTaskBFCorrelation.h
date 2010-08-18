#ifndef ALIFMDANALYSISTASKBFCORRELATION_H
#define ALIFMDANALYSISTASKBFCORRELATION_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
#include "AliAnalysisTask.h"

#include "TObjArray.h"
#include "TObjString.h"
#include "TArrayI.h"
#include "TH1I.h"
#include "TH2.h"
#include "AliMCEvent.h"
#include "AliFMDFloatMap.h"
#include "TCanvas.h"

/**
 * @ingroup FMD_ana
 */
class AliFMDAnalysisTaskBFCorrelation : public AliAnalysisTask
{
 public:
    AliFMDAnalysisTaskBFCorrelation();
    AliFMDAnalysisTaskBFCorrelation(const char* name, Bool_t SE = kTRUE);
    virtual ~AliFMDAnalysisTaskBFCorrelation() {;}
  AliFMDAnalysisTaskBFCorrelation(const AliFMDAnalysisTaskBFCorrelation& o) : AliAnalysisTask(),
									      fDebug(o.fDebug),
									      fOutputList(0),
									      fInputList(0),
									      fInternalList(0),
									      fVertexString(o.fVertexString),
									      fStandalone(o.fStandalone),
									      fEvent(0),
									      fnBinsX(200),
									      fXmin(-6),
									      fXmax(6),
									      fnBinsY(20),
									      fYmin(2),
									      fYmax(2 * TMath::Pi()),
									      c(0),
									      debug0(0),
									      debug1(0)
 {}
  
  AliFMDAnalysisTaskBFCorrelation& operator=(const AliFMDAnalysisTaskBFCorrelation&) { return *this; }
  // Implementation of interface methods
  virtual void ConnectInputData(Option_t *option = "");
  virtual void CreateOutputObjects();
  virtual void Init() {}
  virtual void LocalInit() {Init();}
  virtual void Exec(Option_t *option);
  virtual void Terminate(Option_t *option);
  virtual void SetDebugLevel(Int_t level) {fDebug = level;}
  void SetInputList(TList* inputList) {fInputList = inputList;}
  void SetInputVertex(TObjString* vtxString) {fVertexString = vtxString;}
  void SetOutputList(TList* outputList) {fOutputList = outputList;}
  void CountESDHits();
  void CalculateValues(TString type);
  void SetBounds();
  void MergeEvent(TH2D *hMultTrVtxFull);
  
  void ProcessPrimary();
  
  TList* GetOutputList() {return fOutputList;}
   
 private:
  
  Int_t         fDebug;        //  Debug flag
  TList*        fOutputList;
  TList*        fInputList;
  TList*        fInternalList;
  TObjString*   fVertexString;
  Bool_t        fStandalone;
 
  Int_t   fEvent;
  Int_t   fnBinsX;
  Float_t fXmin;
  Float_t fXmax;
  Int_t   fnBinsY;
  Float_t fYmin;
  Float_t fYmax;

  TCanvas *c;

  Double_t debug0;
  Double_t debug1;

  ClassDef(AliFMDAnalysisTaskBFCorrelation, 0); // Analysis task for FMD analysis
};
 
#endif
// Local Variables:
//   mode: C++ 
// End:
