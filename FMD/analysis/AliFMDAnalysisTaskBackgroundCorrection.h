#ifndef ALIFMDANALYSISTASKBACKGROUNDCORRECTION_H
#define ALIFMDANALYSISTASKBACKGROUNDCORRECTION_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
#include "AliAnalysisTask.h"

#include "TObjArray.h"
#include "TObjString.h"
#include "TArrayI.h"
#include "TH1I.h"

class AliFMDAnalysisTaskBackgroundCorrection : public AliAnalysisTask
{
 public:
    AliFMDAnalysisTaskBackgroundCorrection();
    AliFMDAnalysisTaskBackgroundCorrection(const char* name);
    virtual ~AliFMDAnalysisTaskBackgroundCorrection() {;}
 AliFMDAnalysisTaskBackgroundCorrection(const AliFMDAnalysisTaskBackgroundCorrection& o) : AliAnalysisTask(),
      fDebug(o.fDebug),
      fOutputList(),
      fArray(o.fArray),
      fInputArray(o.fInputArray),
      fVertexString(o.fVertexString),
      fNevents(o.fNevents)  {}
    AliFMDAnalysisTaskBackgroundCorrection& operator=(const AliFMDAnalysisTaskBackgroundCorrection&) { return *this; }
    // Implementation of interface methods
    virtual void ConnectInputData(Option_t *option = "");
    virtual void CreateOutputObjects();
    virtual void Init() {}
    virtual void LocalInit() {Init();}
    virtual void Exec(Option_t *option);
    virtual void Terminate(Option_t *option);
    virtual void SetDebugLevel(Int_t level) {fDebug = level;}
    
 private:
    Int_t         fDebug;        //  Debug flag
    TList         fOutputList;
    TObjArray     fArray;
    TObjArray*    fInputArray;
    TObjString*   fVertexString;
    TH1I          fNevents;
    ClassDef(AliFMDAnalysisTaskBackgroundCorrection, 0); // Analysis task for FMD analysis
};
 
#endif
