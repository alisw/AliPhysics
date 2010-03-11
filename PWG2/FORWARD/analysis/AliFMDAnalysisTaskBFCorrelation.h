#ifndef ALIFMDANALYSISTASKBFCORRELATION_H
#define ALIFMDANALYSISTASKBFCORRELATION_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
#include "AliAnalysisTask.h"

#include "TObjArray.h"
#include "TObjString.h"
#include "TArrayI.h"
#include "TH1I.h"
#include "AliMCEvent.h"
#include "AliFMDFloatMap.h"

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
									      fVertexString(o.fVertexString),
									      fStandalone(o.fStandalone)
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
  
  void ProcessPrimary();
  
  TList* GetOutputList() {return fOutputList;}
   
 private:
  Int_t         fDebug;        //  Debug flag
  TList*        fOutputList;
  TList*        fInputList;
  TObjString*   fVertexString;
  Bool_t        fStandalone;
  ClassDef(AliFMDAnalysisTaskBFCorrelation, 0); // Analysis task for FMD analysis
};
 
#endif
// Local Variables:
//   mode: C++ 
// End:
