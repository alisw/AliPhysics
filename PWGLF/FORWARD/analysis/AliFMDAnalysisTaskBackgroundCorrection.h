#ifndef ALIFMDANALYSISTASKBACKGROUNDCORRECTION_H
#define ALIFMDANALYSISTASKBACKGROUNDCORRECTION_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
#include "AliAnalysisTask.h"

#include "TObjArray.h"
#include "TObjString.h"
#include "TArrayI.h"
#include "TH1I.h"

/**
 * @ingroup FMD_ana
 * @brief Applu the background correction
 * particles. 
 * 
 */
class AliFMDAnalysisTaskBackgroundCorrection : public AliAnalysisTask
{
 public:
    AliFMDAnalysisTaskBackgroundCorrection();
    AliFMDAnalysisTaskBackgroundCorrection(const char* name, Bool_t SE = kTRUE);
    virtual ~AliFMDAnalysisTaskBackgroundCorrection() {;}
 AliFMDAnalysisTaskBackgroundCorrection(const AliFMDAnalysisTaskBackgroundCorrection& o) : AliAnalysisTask(),
      fDebug(o.fDebug),
      fOutputList(0),
      fInputList(0),
      fHitList(0),
      fVertexString(o.fVertexString),
      fNevents(o.fNevents),
      fStandalone(o.fStandalone), 
      fOutputVertexString(o.fOutputVertexString) {}
    AliFMDAnalysisTaskBackgroundCorrection& operator=(const AliFMDAnalysisTaskBackgroundCorrection&) { return *this; }
    // Implementation of interface methods
    virtual void ConnectInputData(Option_t *option = "");
    virtual void CreateOutputObjects();
    virtual void Init() {}
    virtual void LocalInit() {Init();}
    virtual void Exec(Option_t *option);
    virtual void Terminate(Option_t *option);
    virtual void SetDebugLevel(Int_t level) {fDebug = level;}
  void SetInputList(TList* inputList) {fInputList = inputList;}
  void SetOutputVertex(TObjString* vtxString) {fOutputVertexString = vtxString;}
  //void SetInputVtx(TObjString* vtxString) {fVertexString = vtxString;}
  void SetOutputList(TList* outputList) {fOutputList = outputList;}
  void SetHitList(TList* hitList) {fHitList = hitList;}
  void         CreatePerEventHistogram(Int_t vtxbin);
 private:
    Int_t         fDebug;        //  Debug flag
    TList*        fOutputList;
    TList*        fInputList;
    TList*        fHitList;
    TObjString*   fVertexString;
    TH1I          fNevents;
    Bool_t        fStandalone;
    TObjString*   fOutputVertexString;
    ClassDef(AliFMDAnalysisTaskBackgroundCorrection, 0); // Analysis task for FMD analysis
};
 
#endif
// Local Variables:
//   mode: C++
// End:
