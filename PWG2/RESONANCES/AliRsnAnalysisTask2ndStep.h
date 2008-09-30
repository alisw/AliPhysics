//
// Class AliRsnAnalysisTask2ndStep
//
// AnalysisTask which collects an input of RSN events
// and produces histograms.
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//
#ifndef ALIRSNANALYSISTASK2NDSTEP_H
#define ALIRSNANALYSISTASK2NDSTEP_H

#include "AliRsnDaughter.h"
#include "AliAnalysisTaskSE.h"

class TObjArray;

class AliAODEvent;
class AliAODInputHandler;
class AliAnalysisManager;

class AliRsnEvent;
class AliRsnEventBuffer;
class AliRsnPairMgr;

class AliRsnAnalysisTask2ndStep : public AliAnalysisTaskSE
{
  public:
    AliRsnAnalysisTask2ndStep(const char *name = "AliRsnAnalysisTask2ndStep");
    virtual ~AliRsnAnalysisTask2ndStep() {/* Does nothing*/}

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *option);

    void AddPairMgr(AliRsnPairMgr *pairmgr);
    void CreateHandlers(AliAnalysisManager *am);
    void SetAnalysisMgr(AliAnalysisManager* theValue) { fAnalysisMgr = theValue; }
    AliAnalysisManager* GetAnalysisMgr() const { return fAnalysisMgr; }
    void ProcessEventAnalysis(AliRsnEvent *curEvent);
    void PostEventProcess();

  private:

    AliRsnAnalysisTask2ndStep(const AliRsnAnalysisTask2ndStep& copy) :
        AliAnalysisTaskSE(copy),fOutList(0x0),fPairMgrs(0x0),fEventBuffer(0x0),
        fRsnHandlerAOD(0x0),fAnalysisMgr(0x0) {}
    AliRsnAnalysisTask2ndStep& operator= (const AliRsnAnalysisTask2ndStep&) {return *this;}

    TList                       *fOutList;         // list of output histograms
    TObjArray                   *fPairMgrs;        // array if pair managers used
    AliRsnEventBuffer           *fEventBuffer;     // event buffer

    AliAODInputHandler          *fRsnHandlerAOD;   // AOD event handler
    AliAnalysisManager          *fAnalysisMgr;     // pointer to current AnalysisMgr

    ClassDef(AliRsnAnalysisTask2ndStep, 1)
};

#endif
