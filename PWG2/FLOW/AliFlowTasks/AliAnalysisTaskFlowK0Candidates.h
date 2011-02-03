/////////////////////////////////////////////////////
// AliAnalysisTaskFlowK0Candidates:
// Analysis task to select K0 candidates for flow analysis.
// Author: Carlos Perez (cperez@cern.ch)
//////////////////////////////////////////////////////

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id: $ */

#ifndef AliAnalysisTaskFlowK0Candidates_H
#define AliAnalysisTaskFlowK0Candidates_H

#include "AliAnalysisTaskSE.h"

class AliESDtrackCuts;
class AliFlowEventCuts;
class TList;
class TH1D;

class AliAnalysisTaskFlowK0Candidates : public AliAnalysisTaskSE {
  private:
    AliAnalysisTaskFlowK0Candidates(const AliAnalysisTaskFlowK0Candidates& analysisTask);
    AliAnalysisTaskFlowK0Candidates& operator=(const AliAnalysisTaskFlowK0Candidates& analysisTask);
    void AddQAEvents();
    void AddQACandidates();
    void ReadFromESD(AliESDEvent *fESD);
    void ReadFromAOD(AliAODEvent *fAOD);

    AliFlowEventCuts *fCutsEvent; // cuts for event
    AliESDtrackCuts  *fCuts;      // cuts for both pis
    TList *fQAList;               // list for QA histos (slot2)
    Double_t fMassMin, fMassMax;  // Mass cutting range
    Double_t fDLcut;
    TH1D *tEvent, *tMulti;
    TH1D *tMass[4], *tDCA[4], *tDL[4], *tCTP[4], *td0d0[4], *tPhi[4], *tEta[4], *tPt[4];
    TH1D *tAPhi[4], *tAEta[4], *tAPt[4], *tBPhi[4], *tBEta[4], *tBPt[4];

  public:
    AliAnalysisTaskFlowK0Candidates();
    AliAnalysisTaskFlowK0Candidates(const char *name, AliFlowEventCuts *cutsEvent, AliESDtrackCuts *cuts, Double_t MassMin, Double_t MassMax);
    virtual ~AliAnalysisTaskFlowK0Candidates();
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *);
    virtual void Terminate(Option_t *);
    virtual void NotifyRun();
    void SetDL( Double_t pMin ) { fDLcut = pMin; }


  ClassDef(AliAnalysisTaskFlowK0Candidates, 1);
};

#endif
