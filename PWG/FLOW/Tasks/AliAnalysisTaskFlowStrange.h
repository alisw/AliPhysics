/////////////////////////////////////////////////////
// AliAnalysisTaskFlowStrange:
// Analysis task to select K0/Lambda candidates for flow analysis.
// Authors: Cristian Ivan (civan@cern.ch)
//          Carlos Perez (cperez@cern.ch)
//          Pawel Debski (pdebski@cern.ch)
//////////////////////////////////////////////////////

/* Copyright(c) 1998-1999, ALICExperiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id: $ */

#ifndef AliAnalysisTaskFlowStrange_H
#define AliAnalysisTaskFlowStrange_H

#include "AliAnalysisTaskSE.h"

class TList;
class TObjArray;
class AliESDtrackCuts;
class AliFlowEventCuts;
class AliPIDResponse;
class AliESDEvent;
class AliAODEvent;
class AliAODv0;
class AliESDv0;

class AliAnalysisTaskFlowStrange : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskFlowStrange();
    AliAnalysisTaskFlowStrange(const Char_t *name, AliFlowEventCuts *cutsEvent,
                               AliFlowTrackCuts *cutsRFPTPC, AliFlowTrackCuts *cutsRFPVZE,
			       AliESDtrackCuts *cutsDau);
    virtual ~AliAnalysisTaskFlowStrange();
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *);
    virtual void Terminate(Option_t *);
    void SetDebug() {fDebug = kTRUE;}
    void SetCuts2010(int setOfCuts);
    void SetK0L0(int specie) {fSpecie=specie;}
    void SetCommonConstants(Int_t massBins, Double_t minMass, Double_t maxMass);

  private:
    AliAnalysisTaskFlowStrange(const AliAnalysisTaskFlowStrange& analysisTask);
    AliAnalysisTaskFlowStrange& operator=(const AliAnalysisTaskFlowStrange& analysisTask);
    void AddQAEvents();
    void AddQACandidates();

    void ReadFromESDv0(AliESDEvent *tESD);
    void ReadFromAODv0(AliAODEvent *tAOD);
    Int_t PassesESDCuts(AliESDv0 *myV0, AliESDEvent *tESD);
    Int_t PassesAODCuts(AliAODv0 *myV0, AliAODEvent *tAOD);
    void AddCandidates();
    void MakeTrack( Double_t mass, Double_t pt, Double_t phi,
		    Double_t eta, Int_t iid, Int_t jid );
    AliPIDResponse *fPIDResponse;     //! PID response object
    Bool_t fDebug; // true if we want to inspect the main steps of the task
    Int_t fSpecie; // K0=>0 L0=>1
    Double_t fV0Cuts[9]; // v0 cuts: dl dca ctp d0 d0d0 qt minEta maxEta PID
    Int_t fMassBins;   // to configure FLOWCOMMON
    Double_t fMinMass; // to configure FLOWCOMMON
    Double_t fMaxMass; // to configure FLOWCOMMON
    AliFlowEventCuts *fCutsEvent; // event cuts
    AliFlowTrackCuts *fCutsRFPTPC; // RFP cuts
    AliFlowTrackCuts *fCutsRFPVZE; // RFP cuts
    AliFlowTrackCuts* fCutsPOI; // POI cuts
    AliESDtrackCuts  *fCutsDau; // Track quality cuts for the daughters
    AliFlowEvent  *fFlowEventTPC;  //flow event TPC
    AliFlowEvent  *fFlowEventVZE;  //flow event VZE
    TObjArray *fCandidates; // Array of selected candidates

    TList *fQAList; // stores the final list of output histograms

  ClassDef(AliAnalysisTaskFlowStrange, 1);
};

#endif
