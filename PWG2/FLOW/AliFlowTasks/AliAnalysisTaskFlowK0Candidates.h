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

    AliFlowEventCuts *fCutsEvent; // cuts for event
    AliESDtrackCuts  *fCuts;      // cuts for both pis
    TList *fList;
    TH1D *fMass,  *fDCA,  *fDL,  *fCTP;
    TH1D *fMassF, *fDCAF, *fDLF, *fCTPF;
    Double_t fMassMin, fMassMax;     // Mass cutting range

  public:
    AliAnalysisTaskFlowK0Candidates();
    AliAnalysisTaskFlowK0Candidates(const char *name, AliFlowEventCuts *cutsEvent, AliESDtrackCuts *cuts, Double_t MassMin, Double_t MassMax);
    virtual ~AliAnalysisTaskFlowK0Candidates();
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *);
    virtual void Terminate(Option_t *);
    virtual void NotifyRun();

  ClassDef(AliAnalysisTaskFlowK0Candidates, 1);  // example of analysis
};

#endif
