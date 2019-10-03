/////////////////////////////////////////////////////
// AliAnalysisTaskFlowK0Candidates:
// Analysis task to select K0 candidates for flow analysis.
// Author: Carlos Perez (cperez@cern.ch)
//////////////////////////////////////////////////////

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id: $ */

#ifndef ALIANALYSISTASKFLOWK0CANDIDATES_H
#define ALIANALYSISTASKFLOWK0CANDIDATES_H

#include "AliAnalysisTaskSE.h"

class AliESDtrackCuts;
class AliFlowEventCuts;
class TList;
class TH1D;

class AliAnalysisTaskFlowK0Candidates : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskFlowK0Candidates();
    AliAnalysisTaskFlowK0Candidates(const char *name, AliFlowEventCuts *cutsEvent, AliESDtrackCuts *cuts, Double_t MassMin, Double_t MassMax);
    virtual ~AliAnalysisTaskFlowK0Candidates();
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *);
    virtual void Terminate(Option_t *);
    virtual void NotifyRun();
    void SetDL( Double_t pMin ) { fDLcut = pMin; }

  private:
    AliAnalysisTaskFlowK0Candidates(const AliAnalysisTaskFlowK0Candidates& analysisTask);
    AliAnalysisTaskFlowK0Candidates& operator=(const AliAnalysisTaskFlowK0Candidates& analysisTask);
    void AddQAEvents();
    void AddQACandidates();
    void ReadFromESD(const AliESDEvent *fESD);
    void ReadFromAOD(const AliAODEvent *fAOD);

    AliFlowEventCuts *fCutsEvent; // cuts for event
    AliESDtrackCuts  *fCuts;      // cuts for both pis
    TList *fQAList;               // list for QA histos (slot2)
    Double_t fMassMin, fMassMax;  // Mass cutting range
    Double_t fDLcut;//
    TH1D *fEvent, *fMulti;//
    TH1D *fMass[4], *fDCA[4], *fDL[4], *fCTP[4], *fd0d0[4], *fPhi[4], *fEta[4], *fPt[4];//
    TH1D *fAPhi[4], *fAEta[4], *fAPt[4], *fBPhi[4], *fBEta[4], *fBPt[4];//

  ClassDef(AliAnalysisTaskFlowK0Candidates, 1);
};

#endif
