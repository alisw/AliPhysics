#ifndef ALIFMDANALYSISTASKDNDETA_H
#define ALIFMDANALYSISTASKDNDETA_H
 
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
class AliFMDAnalysisTaskDndeta : public AliAnalysisTask
{
 public:
    AliFMDAnalysisTaskDndeta();
    AliFMDAnalysisTaskDndeta(const char* name, Bool_t SE = kTRUE);
    virtual ~AliFMDAnalysisTaskDndeta() {;}
 AliFMDAnalysisTaskDndeta(const AliFMDAnalysisTaskDndeta& o) : AliAnalysisTask(),
							       fDebug(o.fDebug),
							       fOutputList(0),
							       fInputList(0),
							       fVertexString(o.fVertexString),
							       fNevents(o.fNevents),
							       fNNSDevents(o.fNNSDevents),
							       fNMCevents(o.fNMCevents),
							       fNMCNSDevents(o.fNMCNSDevents),
							       fStandalone(o.fStandalone),
							       fLastTrackByStrip(o.fLastTrackByStrip),
							       fVtxEff(o.fVtxEff),
							       fVtxEffNSD(o.fVtxEffNSD)
  {}
  AliFMDAnalysisTaskDndeta& operator=(const AliFMDAnalysisTaskDndeta&) { return *this; }
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
  void SetVtxEfficiency(Float_t vtxeff) {fVtxEff = vtxeff;}
  void SetVtxEfficiencyNSD(Float_t vtxeff) {fVtxEffNSD = vtxeff;}
  
 private:
  Int_t          fDebug;        //  Debug flag
  TList*         fOutputList;
  TList*         fInputList;
  TObjString*    fVertexString;
  TH1I           fNevents;
  TH1I           fNNSDevents;
  TH1I           fNMCevents;
  TH1I           fNMCNSDevents;
  Bool_t         fStandalone;
  AliFMDFloatMap fLastTrackByStrip;
  Float_t        fVtxEff;             //Efficiency of vertex sel.
  Float_t        fVtxEffNSD;          //Efficiency of vertex sel., NSD
   
  ClassDef(AliFMDAnalysisTaskDndeta, 0); // Analysis task for FMD analysis
};
 
#endif
// Local Variables:
//   mode: C++ 
// End:
