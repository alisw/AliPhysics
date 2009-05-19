#ifndef ALIFMDANALYSISTASKDENSITY_H
#define ALIFMDANALYSISTASKDENSITY_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
#include "AliAnalysisTask.h"

#include "TObjArray.h"
#include "AliESDFMD.h"
#include "AliESDVertex.h"
#include "TObjString.h"
#include "TTree.h"

class AliESDEvent;
class TChain;
class AliAODEvent;
class TF1;


class AliFMDAnalysisTaskDensity : public AliAnalysisTask
{
 public:
    AliFMDAnalysisTaskDensity();
    AliFMDAnalysisTaskDensity(const char* name, Bool_t SE = kTRUE);
    virtual ~AliFMDAnalysisTaskDensity() {;}
 AliFMDAnalysisTaskDensity(const AliFMDAnalysisTaskDensity& o) : AliAnalysisTask(),
      fDebug(o.fDebug),
      fOutputList(),
      fArray(o.fArray),
      fESD(o.fESD),
      fVertexString(o.fVertexString),
      fVertex(o.fVertex),
      fStandalone(o.fStandalone),
      fStatus(o.fStatus) {}
    AliFMDAnalysisTaskDensity& operator=(const AliFMDAnalysisTaskDensity&) { return *this; }
    // Implementation of interface methods
    virtual void ConnectInputData(Option_t *option);
    virtual void CreateOutputObjects();
    virtual void Init() {}
    virtual void LocalInit() {Init();}
    virtual void Exec(Option_t *option);
    virtual void Terminate(Option_t */*option*/) {}
    virtual void SetDebugLevel(Int_t level) {fDebug = level;}
    
    Float_t       GetAcceptanceCorrection(Char_t ring, UShort_t strip);
    void SetOutputList(TList* outlist) {fOutputList = outlist;}
    void SetInputESDFMD(AliESDFMD* esdfmd) {fESD = esdfmd;}
    void SetInputVertex(AliESDVertex* vertex) {fVertex = vertex;}
    Bool_t GetEventStatus() { return fStatus; }
    Float_t GetPhiFromSector(UShort_t det, Char_t ring, UShort_t sec);
 private:
    
    Int_t         fDebug;        //  Debug flag
    TList*        fOutputList;
    TObjArray     fArray;
    AliESDFMD*    fESD;
    TObjString    fVertexString;
    AliESDVertex* fVertex;
    Bool_t        fStandalone;
    Bool_t        fStatus;
    TF1*          fFuncPos;
    TF1*          fFuncNeg;
    
    ClassDef(AliFMDAnalysisTaskDensity, 0); // Analysis task for FMD analysis
};
 
#endif
