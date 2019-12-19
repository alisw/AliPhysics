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
// For questions, comments, etc, please write to:
//      david.dobrigkeit.chinellato@cern.ch
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

#ifndef AliAnalysisTaskWeakDecayVertexer_H
#define AliAnalysisTaskWeakDecayVertexer_H

class TList;
class TH1F;

class AliV0HypSel;
class AliESDpid;
class AliESDEvent;
class AliPhysicsSelection;

#include "AliEventCuts.h"
//For mapping functionality
#include <map>

using namespace std;

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
    void SetResetInitialPositions( Bool_t lOpt = kTRUE ){
        //Highly experimental, use with care!
        fkResetInitialPositions = lOpt;
    }
    void SetDoImprovedDCAV0DauPropagation( Bool_t lOpt = kTRUE ){
        //Highly experimental, use with care!
        fkDoImprovedDCAV0DauPropagation = lOpt;
    }
    void SetDoMaterialCorrections( Bool_t lOpt = kTRUE ){
        //Highly experimental, use with care!
        //fkDoMaterialCorrection = lOpt;
        std::cout<<"THIS OPTION DOES NOTHING"<<std::endl; 
        //NEVER USE THIS
    }
    void SetXYCase1Preoptimization( Bool_t lOpt = kTRUE ){
        //Highly experimental, use with care!
        fkXYCase1 = lOpt;
    }
    void SetXYCase2Preoptimization( Bool_t lOpt = kTRUE ){
        //Highly experimental, use with care!
        fkXYCase2 = lOpt;
    }
    
    void SetDoImprovedDCACascDauPropagation( Bool_t lOpt = kTRUE ){
        //Highly experimental, use with care!
        fkDoImprovedDCACascDauPropagation = lOpt;
    }
    
    void SetDoXYPlanePreOptCascade( Bool_t lOpt = kTRUE ){
        //Highly experimental, use with care!
        fkDoXYPlanePreOptCascade = lOpt;
    }
    
    void SetDoPureGeometricMinimization( Bool_t lOpt = kTRUE ){
        //Highly experimental, use with care!
        fkDoPureGeometricMinimization = lOpt;
    }
    void SetDoV0Refit ( Bool_t lDoV0Refit = kTRUE) {
        fkDoV0Refit = lDoV0Refit;
    }
    void SetDoCascadeRefit ( Bool_t lDoCascadeRefit = kTRUE) {
        fkDoCascadeRefit = lDoCascadeRefit;
        //WARNING: Requires V0 refit for covariance matrix
        if( lDoCascadeRefit && !fkDoV0Refit ) fkDoV0Refit = kTRUE;
    }
    void SetMaxIterations (Long_t lMaxIter = 100){
        fMaxIterationsWhenMinimizing = lMaxIter;
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

    void SetExtraCleanup ( Bool_t lExtraCleanup = kTRUE) {
        fkExtraCleanup = lExtraCleanup;
    }
//---------------------------------------------------------------------------------------
    void SetRevertexAllEvents     ( Bool_t lOpt ) {
        fkRevertexAllEvents = lOpt;
    }
    void SetUseExtraEvSels ( Bool_t lUseExtraEvSels = kTRUE) {
        fkDoExtraEvSels = lUseExtraEvSels;
    }
    void SetUseStrictPileupCuts () {
        //This will enable Ionut's pileup rejection in AliEventCuts
        fEventCuts.fUseStrongVarCorrelationCut = true;
        fEventCuts.fUseVariablesCorrelationCuts = true;
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
    void SetMinPtV0     ( Float_t lMinPt ) {
        fMinPtV0 = lMinPt;
    }
    void SetMaxPtV0     ( Float_t lMaxPt ) {
        fMaxPtV0 = lMaxPt;
    }
    void SetForceResetV0s     ( Bool_t lOpt ) {
        fkForceResetV0s = lOpt;
    }
    void SetForceResetCascades     ( Bool_t lOpt ) {
        fkForceResetCascades = lOpt;
    }
    void SetCentralityInterval     ( Float_t lMinCent, Float_t lMaxCent ) {
        fMinCentrality = lMinCent;
        fMaxCentrality = lMaxCent;
    }
    
    //Modifications for V0 mass window selection (from Ruben) 
    void SetV0HypSel(TObjArray* selArr);
    const TObjArray* GetV0HypSelArray() const {return fV0HypSelArray;}
    void AddV0HypSel(const AliV0HypSel& sel);
    void AddStandardV0HypSel();
    
    void SetMassWindowAroundCascade     ( Double_t lMassWin ) {
        fMassWindowAroundCascade = lMassWin;
    }
    void SetMinXforXY     ( Double_t lMinX ) {
        fMinXforXYtest = lMinX;
    }
    void SetPreselectX     ( Bool_t lOpt ) {
        fkPreselectX = lOpt;
    }
    void SetSkipLargeXYDCA( Bool_t lOpt = kTRUE) {
        fkSkipLargeXYDCA=lOpt;
    }
    void SetUseMonteCarloAssociation( Bool_t lOpt = kTRUE) {
        fkMonteCarlo=lOpt;
    }
//---------------------------------------------------------------------------------------
    void SetUseImprovedFinding(){
        fkRunV0Vertexer = kTRUE;
        fkRunCascadeVertexer = kTRUE;
        fkDoImprovedDCAV0DauPropagation = kTRUE;
        fkDoImprovedDCACascDauPropagation = kTRUE;
        fkDoPureGeometricMinimization = kTRUE;
        fkDoV0Refit = kTRUE;
        fkDoCascadeRefit = kTRUE;
        fkXYCase1 = kTRUE;
        fkXYCase2 = kTRUE;
    }
//---------------------------------------------------------------------------------------
    void SetUseDefaultFinding(){
        fkRunV0Vertexer = kTRUE;
        fkRunCascadeVertexer = kTRUE;
        fkDoImprovedDCAV0DauPropagation = kFALSE;
        fkDoImprovedDCACascDauPropagation = kFALSE;
        fkDoPureGeometricMinimization = kFALSE;
        fkDoV0Refit = kFALSE;
        fkDoCascadeRefit = kFALSE;
        fkXYCase1 = kFALSE;
        fkXYCase2 = kFALSE;
    }
//---------------------------------------------------------------------------------------
    //Functions for analysis Bookkeepinp
    // 1- Configure standard vertexing
    void SetupStandardVertexing();
    void SetupLooseVertexing();
//---------------------------------------------------------------------------------------
    //Re-vertex V0s
    Long_t Tracks2V0vertices(AliESDEvent *event);

    //======================================================================
    //Re-vertex V0s based solely on perfect MC V0s
    //Warning: this cannot be called in real data without troubles,
    //         but may be particularly useful for performance studies
    //         as even very loose cuts will not incur in prohibitive
    //         performance costs
    Long_t Tracks2V0verticesMC(AliESDEvent *event);
    //======================================================================
    
    //Re-vertex Cascades
    Long_t V0sTracks2CascadeVertices(AliESDEvent *event);
    Long_t V0sTracks2CascadeVerticesMC(AliESDEvent *event);
    //Re-vertex Cascades without checking bachelor charge - V0 Mass hypo correspondence
    Long_t V0sTracks2CascadeVerticesUncheckedCharges(AliESDEvent *event);
    //Helper functions
    Double_t Det(Double_t a00, Double_t a01, Double_t a10, Double_t a11) const;
    Double_t Det(Double_t a00,Double_t a01,Double_t a02,
                 Double_t a10,Double_t a11,Double_t a12,
                 Double_t a20,Double_t a21,Double_t a22) const;
    Double_t PropagateToDCA(AliESDv0 *vtx,AliExternalTrackParam *trk, AliESDEvent *event, Double_t b, Double_t lBachMassForTracking=0.139);
    void Evaluate(const Double_t *h, Double_t t,
                  Double_t r[3],  //radius vector
                  Double_t g[3],  //first defivatives
                  Double_t gg[3]); //second derivatives
    void CheckChargeV0(AliESDv0 *v0);
    //---------------------------------------------------------------------------------------
    //Improved DCA V0 Dau
    Double_t GetDCAV0Dau ( AliExternalTrackParam *pt, AliExternalTrackParam *nt, Double_t &xp, Double_t &xn, Double_t b, Double_t lNegMassForTracking=0.139, Double_t lPosMassForTracking=0.139);
    void GetHelixCenter(const AliExternalTrackParam *track,Double_t center[2], Double_t b);
    //---------------------------------------------------------------------------------------
    
    //---------------------------------------------------------------------------------------
    // changes to enable AliExternalTrackParam inheritance from on-the-fly finder
    //selective reset: go over list of V0s and delete offline (0) or on-the-fly (1) V0s
    void SelectiveResetV0s(AliESDEvent *event, Int_t lType = 0);
    //Master switch
    void SetUseOptimalTrackParams (Bool_t lOpt){
        fkUseOptimalTrackParams = lOpt;
    }
    void SetUseOptimalTrackParamsBachelor (Bool_t lOpt){
        fkUseOptimalTrackParamsBachelor = lOpt;
    }
    //---------------------------------------------------------------------------------------
    

private:
    // Note : In ROOT, "//!" means "do not stream the data from Master node to Worker node" ...
    // your data member object is created on the worker nodes and streaming is not needed.
    // http://root.cern.ch/download/doc/11InputOutput.pdf, page 14
    TList  *fListHist;      //! List of Cascade histograms

    AliPIDResponse *fPIDResponse;     // PID response object

    //Implementation of event selection utility
    AliEventCuts fEventCuts; /// Event cuts class
    
    AliVEvent::EOfflineTriggerTypes fTrigType; // trigger type
    Bool_t fkDoExtraEvSels; //if true, rely on AliEventCuts
    //Min/Max Centrality
    Bool_t fkForceResetV0s;
    Bool_t fkForceResetCascades;
    Float_t fMinCentrality; //centrality interval to actually regenerate candidates
    Float_t fMaxCentrality; //centrality interval to actually regenerate candidates
    
    //Objects Controlling Task Behaviour
    Bool_t fkRevertexAllEvents; //Don't be smart. Re-vertex every single event 
    Bool_t fkPreselectDedx;
    Bool_t fkPreselectDedxLambda;
    Bool_t fkExtraCleanup;           //if true, perform pre-rejection of useless candidates before going through configs
    
    //Objects Controlling Task Behaviour: has to be streamed!
    Bool_t fkRunV0Vertexer;           // if true, re-run V0 vertexer
    Bool_t fkDoV0Refit;
    Bool_t fkXYCase1; //Circles-far-away case pre-optimization switch
    Bool_t fkXYCase2; //Circles-touch case pre-optimization switch (cowboy/sailor duality resolution)
    Bool_t fkResetInitialPositions; 
    Bool_t fkDoImprovedDCAV0DauPropagation;
    Bool_t fkDoMaterialCorrection; //Replace AliExternalTrackParam::PropagateTo with AliTrackerBase::PropagateTrackTo
    Int_t fRunNumber; //keep track of run number, needed to load geometry + invoke AliTrackerBase 
    
    Bool_t fkRunCascadeVertexer;      // if true, re-run cascade vertexer
    Bool_t fkUseUncheckedChargeCascadeVertexer; //if true, use cascade vertexer that does not check bachelor charge
    Bool_t fkUseOnTheFlyV0Cascading;
    Bool_t fkDoImprovedDCACascDauPropagation;
    Bool_t fkDoXYPlanePreOptCascade; 
    Bool_t fkDoPureGeometricMinimization;
    Bool_t fkDoCascadeRefit; //WARNING: needs DoV0Refit!
    Long_t fMaxIterationsWhenMinimizing;
    Bool_t fkPreselectX;
    Bool_t fkSkipLargeXYDCA;
    
    //Master MC switch
    Bool_t fkMonteCarlo; //do MC association in vertexing
    
    //Bool_t to conrtol the use of on-the-fly AliExternalTrackParams
    Bool_t fkUseOptimalTrackParams; //if true, use better track estimates from OTF V0s
    Bool_t fkUseOptimalTrackParamsBachelor; //if true, use better track estimates from OTF V0s
    
    //Min/Max pT for cascades
    Float_t fMinPtV0; //minimum pt above which we keep candidates in TTree output
    Float_t fMaxPtV0; //maximum pt below which we keep candidates in TTree output
    Float_t fMinPtCascade; //minimum pt above which we keep candidates in TTree output
    Float_t fMaxPtCascade; //maximum pt below which we keep candidates in TTree output

    //Mass window for V0s
    TObjArray* fV0HypSelArray; // array of V0 hypothesis to select
    
    //Mass Window around masses of interest (cascades)
    Double_t fMassWindowAroundCascade;
    
    Double_t fMinXforXYtest; //min X allowed for XY-plane preopt test
    
    Double_t  fV0VertexerSels[7];        // Array to store the 7 values for the different selections V0 related
    Double_t  fCascadeVertexerSels[8];   // Array to store the 8 values for the different selections Casc. related
    
    
    
    //(pair) -> (OTF index) map
    std::map<std::pair<int, int>, int> fOTFMap; //std::map to store index pair <-> OTF index equiv
    
//===========================================================================================
//   Histograms
//===========================================================================================

    TH1D *fHistEventCounter; //!
    TH1D *fHistCentrality; //!
    TH1D *fHistNumberOfCandidates; //!
    TH1D *fHistV0ToBachelorPropagationStatus; //!
    TH1D *fHistV0OptimalTrackParamUse; //!
    TH1D *fHistV0OptimalTrackParamUseBachelor; //!
    
    //V0 statistics
    TH1D *fHistV0Statistics; //! 

    AliAnalysisTaskWeakDecayVertexer(const AliAnalysisTaskWeakDecayVertexer&);            // not implemented
    AliAnalysisTaskWeakDecayVertexer& operator=(const AliAnalysisTaskWeakDecayVertexer&); // not implemented

    ClassDef(AliAnalysisTaskWeakDecayVertexer, 1);
    //1: first implementation
};

#endif
