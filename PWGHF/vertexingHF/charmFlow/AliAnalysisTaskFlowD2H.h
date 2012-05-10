/* Copyright(c) 1998-1999, ALICExperiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIANALYSISTASKFLOWD2H_H
#define ALIANALYSISTASKFLOWD2H_H

#include "AliAnalysisTaskSE.h"

//==============================================================================
// FlowD2H main task:
// >> Make flowEvent with RPcuts and POIcuts given in constructor and passes it
//    to the daughter tasks.
// >> The POIcuts are polymorphic based on the AliRDHFCuts class allowing the 
//    use of all charmed candidates reconstructed in the central barrel.
// Author: Carlos Perez (cperez@cern.ch)
//==============================================================================

class AliAODEvent;
class AliRDHFCuts;
class AliRDHFCutsD0toKpi;
class AliFlowEventCuts;
class AliFlowEvent;
class AliFlowCandidateTrack;
class AliFlowTrackCuts;
class TList;
class TH1D;

class AliAnalysisTaskFlowD2H : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskFlowD2H();
    AliAnalysisTaskFlowD2H( const Char_t *name, 
			    AliFlowTrackCuts *cutsTPC,
			    AliFlowTrackCuts *cutsVZE,
			    AliRDHFCuts *cutsPOIs, 
			    Int_t specie );
    void SetDebug() {fDebugV2 = true;}
    void SetSwapAsumption() {fSwap = true;}
    virtual ~AliAnalysisTaskFlowD2H();
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *);
    void SetCommonConstants(Int_t massBins, Double_t minMass, Double_t maxMass, Int_t ptWidth);

  private:
    AliAnalysisTaskFlowD2H(const AliAnalysisTaskFlowD2H& analysisTask);
    AliAnalysisTaskFlowD2H& operator=(const AliAnalysisTaskFlowD2H& analysisTask);
    void AddHistograms();
    void FillD0toKpi( const AliAODEvent *aod );
    void FillD0toKpipipi( const AliAODEvent *aod );
    void FillDStartoKpipi( const AliAODEvent *aod );
    void FillDplustoKpipi( const AliAODEvent *aod );
    void FillDstoKKpi( const AliAODEvent *aod );
    void FillJpsitoee( const AliAODEvent *aod );
    void FillLctoV0( const AliAODEvent *aod );
    void FillLctopKpi( const AliAODEvent *aod );
    void MakeTrack( Double_t mass, Double_t pt,
		    Double_t phi, Double_t eta,
		    Int_t nDaughters, const Int_t *iID );

    AliFlowEvent *fTPCEvent;
    AliFlowEvent *fVZEEvent;
    AliFlowTrackCuts *fCutsTPC;
    AliFlowTrackCuts *fCutsVZE;
    AliFlowTrackCuts *fNoPOIs;

    AliRDHFCuts *fCutsPOI; // cuts for POIs
    Int_t  fSource; // AliRDHFCuts::ESele
    Bool_t fDebugV2; // fully talkative task
    Bool_t fSwap;   // swap assumption (for neutral)

    Int_t fMassBins; // configures mass bins for the analysis
    Double_t fMinMass; // configures mass range for the analysis
    Double_t fMaxMass; // configures mass range for the analysis
    Int_t fPtBinWidth; // configures pt bin width for the analysis

    TList *fHList; // List for histos
    TH1D *fEvent;  // Event counter
    TH1D *fCC;     // CC histogram
    TH1D *fRFPMTPC; // Multiplicity RFPTPC
    TH1D *fRFPPhiTPC; // Phi RFPTPC

    TObjArray *fCandidates; // Array of selected candidates

  ClassDef(AliAnalysisTaskFlowD2H, 4);
};

#endif
