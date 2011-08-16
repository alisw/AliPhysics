/* Copyright(c) 1998-1999, ALICExperiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id: $ */

#ifndef AliAnalysisTaskFlowD2H_H
#define AliAnalysisTaskFlowD2H_H

#include "AliAnalysisTaskSE.h"

//==============================================================================
// FlowD2H main task:
// >> Make flowEvent with RPcuts and POIcuts given in constructor and passes it
//    to the daughter tasks.
// >> The POIcuts are polymorphic based on the AliRDHFCuts class allowing the 
//    use of all charmed candidates reconstructed in the central barrel.
// Author: Carlos Perez (cperez@cern.ch)
//==============================================================================

class TList;
class TH2D;
class TProfile;

class AliAODEvent;

class AliFlowEvent;
class AliFlowCandidateTrack;
class AliFlowTrackCuts;

class AliRDHFCuts;
class AliRDHFCutsD0toKpi;

class AliAnalysisTaskFlowD2H : public AliAnalysisTaskSE {
  public:
    AliAnalysisTaskFlowD2H();
    AliAnalysisTaskFlowD2H( const Char_t *name, AliFlowTrackCuts *cutsRPs, 
                            AliRDHFCuts *cutsPOIs, Int_t specie );
    void SetDebug() {fDebug = true;}
    void SetPOIEtaRange( Double_t minEta, Double_t maxEta )
                          { fPOIEta[0] = minEta; fPOIEta[1] = maxEta; }
    void SetFlowEtaRangeAB( Double_t minA, Double_t maxA, Double_t minB, Double_t maxB )
                          { fFlowEta[0] = minA; fFlowEta[1] = maxA;
                            fFlowEta[2] = minB; fFlowEta[3] = maxB; }
    void SetFlowPtRange( Int_t minPt, Int_t maxPt )
                       { fFlowPts[0] = minPt; fFlowPts[1] = maxPt; }
    void SetFlowBandRange( Int_t band, Double_t minMass, Double_t maxMass )
                         { fFlowBands[0][band] = minMass;
                           fFlowBands[1][band] = maxMass; } //TODO
    virtual ~AliAnalysisTaskFlowD2H();
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *);
    virtual void Terminate(Option_t *);
    virtual void NotifyRun();

  private:
    AliAnalysisTaskFlowD2H(const AliAnalysisTaskFlowD2H& analysisTask);
    AliAnalysisTaskFlowD2H& operator=(const AliAnalysisTaskFlowD2H& analysisTask);
    void AddHistograms();
    void FillD0toKpi(      AliAODEvent *aod, AliFlowEvent *mb[5]);
    void FillD0toKpipipi(  AliAODEvent *aod, AliFlowEvent *mb[5]);
    void FillDStartoKpipi( AliAODEvent *aod, AliFlowEvent *mb[5]);
    void FillDplustoKpipi( AliAODEvent *aod, AliFlowEvent *mb[5]);
    void FillDstoKKpi(     AliAODEvent *aod, AliFlowEvent *mb[5]);
    void FillJpsitoee(     AliAODEvent *aod, AliFlowEvent *mb[5]);
    void FillLctoV0(       AliAODEvent *aod, AliFlowEvent *mb[5]);
    void FillLctopKpi(     AliAODEvent *aod, AliFlowEvent *mb[5]);
    AliFlowCandidateTrack* MakeTrack( Double_t mass, Double_t pt,
                                      Double_t phi, Double_t eta,
                                      Int_t nDaughters, Int_t *iID );
    AliFlowTrackCuts *fCutsRP;  // cuts for RPs
    AliRDHFCuts      *fCutsPOI; // cuts for POIs
    Int_t  fSource;             // AliRDHFCuts::ESele
    Bool_t fDebug;              // fully talkative task
    TList *fHList;    // List for histos
    TProfile *fAnaCuts; // store analysis related cuts
    TH2D  *fEvent[2]; // Events histogram
    TH2D  *fMass[2];  // Mass spectra
    Double_t fPOIEta[2];        // Eta cut for POI
    Double_t fFlowEta[4];       // SP subEvents
    Int_t fFlowPts[2];          // Pt range
    Double_t fFlowBands[2][5];  // Mass bands TODO

  ClassDef(AliAnalysisTaskFlowD2H, 1);
};

#endif
