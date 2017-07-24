/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
//
// Modified version of AliAnalysisTaskCheckCascade.h
// Used bits of code from AliAnalysisTaskCheckPerformanceStrange
//
// --- David Dobrigkeit Chinellato
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#ifndef AliAnalysisTaskWeakDecayVertexer_H
#define AliAnalysisTaskWeakDecayVertexer_H

class TList;
class TH1F;

class AliESDpid;
class AliESDEvent;
class AliPhysicsSelection;

#include "AliEventCuts.h"

class AliAnalysisTaskWeakDecayVertexer : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskWeakDecayVertexer();
    AliAnalysisTaskWeakDecayVertexer(const char *name, TString lExtraOptions = "");
    virtual ~AliAnalysisTaskWeakDecayVertexer();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);

    void SetPreselectDedx (Bool_t lPreselectDedx= kTRUE ) {
        fkPreselectDedx   = lPreselectDedx;
    }
    void SetPreselectDedxLambda (Bool_t lPreselectDedx= kTRUE ) {
        fkPreselectDedxLambda   = lPreselectDedx;
    }
    void SetUseOnTheFlyV0Cascading( Bool_t lUseOnTheFlyV0Cascading = kTRUE ){
        //Highly experimental, use with care!
        fkUseOnTheFlyV0Cascading = lUseOnTheFlyV0Cascading;
    }
    void SetDoImprovedCascadeVertexFinding( Bool_t lOpt = kTRUE ){
        //Highly experimental, use with care!
        fkDoImprovedCascadeVertexFinding = lOpt;
    }
    void SetDoImprovedCascadePosition( Bool_t lOpt = kTRUE ){
        //Highly experimental, use with care!
        fkDoImprovedCascadePosition = lOpt;
    }
    void SetIfImprovedPerformInitialLinearPropag( Bool_t lOpt = kTRUE ){
        //Highly experimental, use with care!
        fkIfImprovedPerformInitialLinearPropag = lOpt;
    }
    void SetIfImprovedExtraPrecisionFactor( Double_t lOpt ){
        //Highly experimental, use with care!
        fkIfImprovedExtraPrecisionFactor = lOpt;
    }
    
//---------------------------------------------------------------------------------------
    //Task Configuration: trigger selection
    void SetSelectedTriggerClass(AliVEvent::EOfflineTriggerTypes trigType) { fTrigType = trigType;}
//---------------------------------------------------------------------------------------
    //Task Configuration: Meant to enable quick re-execution of vertexer if needed
    void SetRunV0Vertexer ( Bool_t lRunVertexer = kTRUE) {
        fkRunV0Vertexer = lRunVertexer;
    }
    void SetRunCascadeVertexer ( Bool_t lRunVertexer = kTRUE) {
        fkRunCascadeVertexer = lRunVertexer;
    }
    void SetUseUncheckedChargeCascadeVertexer ( Bool_t lOpt = kTRUE) {
        //WARNING: Experimental vertexer which disregards bachelor charge when creating candidates!
        //         The user has to take care... 
        fkUseUncheckedChargeCascadeVertexer = lOpt;
    }
    void SetDoV0Refit ( Bool_t lDoV0Refit = kTRUE) {
        fkDoV0Refit = lDoV0Refit;
    }
    void SetExtraCleanup ( Bool_t lExtraCleanup = kTRUE) {
        fkExtraCleanup = lExtraCleanup;
    }
//---------------------------------------------------------------------------------------
    void SetUseExtraEvSels ( Bool_t lUseExtraEvSels = kTRUE) {
        fkDoExtraEvSels = lUseExtraEvSels;
    }
//---------------------------------------------------------------------------------------
//Setters for the V0 Vertexer Parameters
    void SetV0VertexerMaxChisquare   ( Double_t lParameter ) {
        fV0VertexerSels[0] = lParameter;
    }
    void SetV0VertexerDCAFirstToPV   ( Double_t lParameter ) {
        fV0VertexerSels[1] = lParameter;
    }
    void SetV0VertexerDCASecondtoPV  ( Double_t lParameter ) {
        fV0VertexerSels[2] = lParameter;
    }
    void SetV0VertexerDCAV0Daughters ( Double_t lParameter ) {
        fV0VertexerSels[3] = lParameter;
    }
    void SetV0VertexerCosinePA       ( Double_t lParameter ) {
        fV0VertexerSels[4] = lParameter;
    }
    void SetV0VertexerMinRadius      ( Double_t lParameter ) {
        fV0VertexerSels[5] = lParameter;
    }
    void SetV0VertexerMaxRadius      ( Double_t lParameter ) {
        fV0VertexerSels[6] = lParameter;
    }
//---------------------------------------------------------------------------------------
//Setters for the Cascade Vertexer Parameters
    void SetCascVertexerMaxChisquare         ( Double_t lParameter ) {
        fCascadeVertexerSels[0] = lParameter;
    }
    void SetCascVertexerMinV0ImpactParameter ( Double_t lParameter ) {
        fCascadeVertexerSels[1] = lParameter;
    }
    void SetCascVertexerV0MassWindow         ( Double_t lParameter ) {
        fCascadeVertexerSels[2] = lParameter;
    }
    void SetCascVertexerDCABachToPV          ( Double_t lParameter ) {
        fCascadeVertexerSels[3] = lParameter;
    }
    void SetCascVertexerDCACascadeDaughters  ( Double_t lParameter ) {
        fCascadeVertexerSels[4] = lParameter;
    }
    void SetCascVertexerCascadeCosinePA      ( Double_t lParameter ) {
        fCascadeVertexerSels[5] = lParameter;
    }
    void SetCascVertexerCascadeMinRadius     ( Double_t lParameter ) {
        fCascadeVertexerSels[6] = lParameter;
    }
    void SetCascVertexerCascadeMaxRadius     ( Double_t lParameter ) {
        fCascadeVertexerSels[7] = lParameter;
    }
//---------------------------------------------------------------------------------------
    void SetMinPtCascade     ( Float_t lMinPt ) {
        fMinPtCascade = lMinPt;
    }
    void SetMaxPtCascade     ( Float_t lMaxPt ) {
        fMaxPtCascade = lMaxPt;
    }
    void SetCentralityInterval     ( Float_t lMinCent, Float_t lMaxCent ) {
        fMinCentrality = lMinCent;
        fMaxCentrality = lMaxCent;
    }
    void SetMassWindowAroundCascade     ( Double_t lMassWin ) {
        fMassWindowAroundCascade = lMassWin;
    }
//---------------------------------------------------------------------------------------
    //Functions for analysis Bookkeepinp
    // 1- Configure standard vertexing
    void SetupStandardVertexing();
    void SetupLooseVertexing();
//---------------------------------------------------------------------------------------
    //Re-vertex V0s
    Long_t Tracks2V0vertices(AliESDEvent *event);
    //Re-vertex Cascades
    Long_t V0sTracks2CascadeVertices(AliESDEvent *event);
    //Re-vertex Cascades without checking bachelor charge - V0 Mass hypo correspondence
    Long_t V0sTracks2CascadeVerticesUncheckedCharges(AliESDEvent *event);
    //Helper functions
    Double_t Det(Double_t a00, Double_t a01, Double_t a10, Double_t a11) const;
    Double_t Det(Double_t a00,Double_t a01,Double_t a02,
                 Double_t a10,Double_t a11,Double_t a12,
                 Double_t a20,Double_t a21,Double_t a22) const;
    Double_t PropagateToDCA(AliESDv0 *vtx,AliExternalTrackParam *trk, AliESDEvent *event, Double_t b);
    void Evaluate(const Double_t *h, Double_t t,
                  Double_t r[3],  //radius vector
                  Double_t g[3],  //first defivatives
                  Double_t gg[3]); //second derivatives
    void CheckChargeV0(AliESDv0 *v0);
    //---------------------------------------------------------------------------------------

private:
    // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
    // your data member object is created on the worker nodes and streaming is not needed.
    // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
    TList  *fListHist;      //! List of Cascade histograms

    AliPIDResponse *fPIDResponse;     // PID response object

    //Implementation of event selection utility
    AliEventCuts fEventCuts; /// Event cuts class

    //Objects Controlling Task Behaviour
    Bool_t fkPreselectDedx;
    Bool_t fkPreselectDedxLambda;
    Bool_t fkUseOnTheFlyV0Cascading;
    Bool_t fkDoImprovedCascadeVertexFinding;
    Bool_t fkDoImprovedCascadePosition;
    Bool_t fkIfImprovedPerformInitialLinearPropag;
    Double_t fkIfImprovedExtraPrecisionFactor;
    Bool_t fkDoExtraEvSels; //if true, rely on AliEventCuts

    //Objects Controlling Task Behaviour: has to be streamed!
    Bool_t    fkRunV0Vertexer;           // if true, re-run V0 vertexer
    Bool_t    fkRunCascadeVertexer;      // if true, re-run cascade vertexer
    Bool_t    fkUseUncheckedChargeCascadeVertexer; //if true, use cascade vertexer that does not check bachelor charge
    Bool_t    fkDoV0Refit;              // if true, will invoke AliESDv0::Refit in the vertexing procedure
    Bool_t    fkExtraCleanup;           //if true, perform pre-rejection of useless candidates before going through configs

    AliVEvent::EOfflineTriggerTypes fTrigType; // trigger type

    Double_t  fV0VertexerSels[7];        // Array to store the 7 values for the different selections V0 related
    Double_t  fCascadeVertexerSels[8];   // Array to store the 8 values for the different selections Casc. related

    //Min/Max pT for cascades
    Float_t fMinPtCascade; //minimum pt above which we keep candidates in TTree output
    Float_t fMaxPtCascade; //maximum pt below which we keep candidates in TTree output
    
    //Min/Max Centrality
    Float_t fMinCentrality; //centrality interval to actually regenerate candidates
    Float_t fMaxCentrality; //centrality interval to actually regenerate candidates
    
    //Mass Window around masses of interest
    Double_t fMassWindowAroundCascade; 

    //===========================================================================================
//   Histograms
//===========================================================================================

    TH1D *fHistEventCounter; //!
    TH1D *fHistCentrality; //!
    TH1D *fHistNumberOfCandidates; //!
    
    TH1D *fHistV0ToBachelorPropagationStatus; //! 

    AliAnalysisTaskWeakDecayVertexer(const AliAnalysisTaskWeakDecayVertexer&);            // not implemented
    AliAnalysisTaskWeakDecayVertexer& operator=(const AliAnalysisTaskWeakDecayVertexer&); // not implemented

    ClassDef(AliAnalysisTaskWeakDecayVertexer, 1);
    //1: first implementation
};

#endif
