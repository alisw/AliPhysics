//Constructed from AliAnalysisTaskStrangenessVsMultiplicityMCRun2.cxx
//Ejiro Umaka Mar 15 2018


class TTree;
class TParticle;
class TVector3;

//class AliMCEventHandler;
//class AliMCEvent;
//class AliStack;

class AliESDVertex;
class AliAODVertex;
class AliESDv0;
class AliAODv0;

#include <Riostream.h>
#include "TList.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TRandom3.h"

#include "AliESDEvent.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliCentrality.h"
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"

#include "AliCFContainer.h"
#include "AliMultiplicity.h"
#include "AliAODMCParticle.h"

#include "AliESDUtils.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTaskNetLambdaMCTrad.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskNetLambdaMCTrad)

AliAnalysisTaskNetLambdaMCTrad::AliAnalysisTaskNetLambdaMCTrad()
    : AliAnalysisTaskSE(), fListHist(0),fTreeEvent(0), fTreeV0(0), fPIDResponse(0), fESDtrackCuts(0), fUtils(0), fRand(0),

//---> Flags controlling Event Tree output
fkSaveEventTree    ( kTRUE ), //no downscaling in this tree so far

//---> Flags controlling V0 TTree output
fkSaveV0Tree       ( kTRUE ),
fkDownScaleV0      ( kFALSE  ),
fDownScaleFactorV0 ( 0  ),
fPtArray(NULL),
fNptBins(19),
fkPreselectDedx ( kFALSE ),
fkPreselectPID  ( kTRUE  ),
fkUseOnTheFlyV0Cascading( kFALSE ),
fkDebugWrongPIDForTracking ( kFALSE ),
fkDebugBump( kFALSE ),
fIsMC(kTRUE),


//---> Flag controlling trigger selection
fTrigType(AliVEvent::kINT7),

//---> Variables for fTreeEvent
fCentrality(0),
fMVPileupFlag(kFALSE),

//---> Variables for fTreeV0
fTreeVariableChi2V0(0),
fTreeVariableDcaV0Daughters(0),
fTreeVariableDcaV0ToPrimVertex(0),
fTreeVariableDcaPosToPrimVertex(0),
fTreeVariableDcaNegToPrimVertex(0),
fTreeVariableV0CosineOfPointingAngle(0),
fTreeVariableV0Radius(0),
fTreeVariablePt(0),
fTreeVariablePtMC(0),
fTreeVariableRapK0Short(0),
fTreeVariableRapLambda(0),
fTreeVariableRapMC(0),
fTreeVariableInvMassK0s(0),
fTreeVariableInvMassLambda(0),
fTreeVariableInvMassAntiLambda(0),
fTreeVariableAlphaV0(0),
fTreeVariablePtArmV0(0),
fTreeVariableNegEta(0),
fTreeVariablePosEta(0),
fTreeVariableMaxChi2PerCluster(0),
fTreeVariableMinTrackLength(0),

fTreeVariableNSigmasPosProton(0),
fTreeVariableNSigmasPosPion(0),
fTreeVariableNSigmasNegProton(0),
fTreeVariableNSigmasNegPion(0),

fTreeVariablePosPIDForTracking(-1),
fTreeVariableNegPIDForTracking(-1),
fTreeVariablePosdEdx(-1),
fTreeVariableNegdEdx(-1),
fTreeVariablePosInnerP(-1),
fTreeVariableNegInnerP(-1),
fTreeVariableNegTrackStatus(0),
fTreeVariablePosTrackStatus(0),
fTreeVariableNegDCAz(-1),
fTreeVariablePosDCAz(-1),

fTreeVariableDistOverTotMom(0),
fTreeVariableLeastNbrCrossedRows(0),
fTreeVariableLeastRatioCrossedRowsOverFindable(0),

fTreeVariableCentrality(0),
fTreeVariableMVPileupFlag(kFALSE),
//MC Variables
fTreeVariablePtMother(0),
fTreeVariableRapMother(0),
fTreeVariablePID(0),
fTreeVariablePIDPositive(0),
fTreeVariablePIDNegative(0),
fTreeVariablePIDMother(0),
fTreeVariablePrimaryStatus(0),
fTreeVariablePrimaryStatusMother(0),

//thnsparse
fPtBinNplusNminusCh(NULL),
fPtBinNplusNminusChTruth(NULL),


//Histos
fHistEventCounter(0),
fHistCentrality(0),
//V0s
fHistGeneratedPtVsYVsCentralityK0Short(0),
fHistGeneratedPtVsYVsCentralityLambda(0),
fHistGeneratedPtVsYVsCentralityAntiLambda(0),
//Added
f3dHistInvMassVsPtVsCentLambda(0),
f3dHistInvMassVsPtVsCentAntiLambda(0),
fhistLambdaRecPT(0),
fhistAntiLambdaRecPT(0),
fhistLambdaRecpri(0),
fhistAntiLambdaRecpri(0),
fHistGeneratedCentralityVsLambda(0),
fHistGeneratedCentralityVsAntiLambda(0)


//------------------------------------------------
// Tree Variables
{

}

AliAnalysisTaskNetLambdaMCTrad::AliAnalysisTaskNetLambdaMCTrad(Bool_t lSaveEventTree, Bool_t lSaveV0Tree, const char *name, TString lExtraOptions)
    : AliAnalysisTaskSE(name),  fListHist(0),  fTreeEvent(0), fTreeV0(0),  fPIDResponse(0), fESDtrackCuts(0), fUtils(0), fRand(0),

//---> Flags controlling Event Tree output
fkSaveEventTree    ( kFALSE ), //no downscaling in this tree so far

//---> Flags controlling V0 TTree output
fkSaveV0Tree       ( kTRUE ),
fkDownScaleV0      ( kFALSE  ),
fDownScaleFactorV0 ( 0  ),
fPtArray(NULL),
fNptBins(19),
fkPreselectDedx ( kFALSE ),
fkPreselectPID  ( kTRUE  ),
fkUseOnTheFlyV0Cascading( kFALSE ), 
fkDebugWrongPIDForTracking ( kFALSE ), //also for cascades...
fkDebugBump( kFALSE ),
fIsMC(kTRUE),


//---> Flag controlling trigger selection
fTrigType(AliVEvent::kINT7),

//---> Variables for fTreeEvent
fCentrality(0),
fMVPileupFlag(kFALSE),

//---> Variables for fTreeV0
fTreeVariableChi2V0(0),
fTreeVariableDcaV0Daughters(0),
fTreeVariableDcaV0ToPrimVertex(0),
fTreeVariableDcaPosToPrimVertex(0),
fTreeVariableDcaNegToPrimVertex(0),
fTreeVariableV0CosineOfPointingAngle(0),
fTreeVariableV0Radius(0),
fTreeVariablePt(0),
fTreeVariablePtMC(0),
fTreeVariableRapK0Short(0),
fTreeVariableRapLambda(0),
fTreeVariableRapMC(0),
fTreeVariableInvMassK0s(0),
fTreeVariableInvMassLambda(0),
fTreeVariableInvMassAntiLambda(0),
fTreeVariableAlphaV0(0),
fTreeVariablePtArmV0(0),
fTreeVariableNegEta(0),
fTreeVariablePosEta(0),
fTreeVariableMaxChi2PerCluster(0),
fTreeVariableMinTrackLength(0),

fTreeVariableNSigmasPosProton(0),
fTreeVariableNSigmasPosPion(0),
fTreeVariableNSigmasNegProton(0),
fTreeVariableNSigmasNegPion(0),

fTreeVariablePosPIDForTracking(-1),
fTreeVariableNegPIDForTracking(-1),
fTreeVariablePosdEdx(-1),
fTreeVariableNegdEdx(-1),
fTreeVariablePosInnerP(-1),
fTreeVariableNegInnerP(-1),
fTreeVariableNegTrackStatus(0),
fTreeVariablePosTrackStatus(0),
fTreeVariableNegDCAz(-1),
fTreeVariablePosDCAz(-1),

fTreeVariableDistOverTotMom(0),
fTreeVariableLeastNbrCrossedRows(0),
fTreeVariableLeastRatioCrossedRowsOverFindable(0),

fTreeVariableCentrality(0),
fTreeVariableMVPileupFlag(kFALSE),

//MC Variables
fTreeVariablePtMother(0),
fTreeVariableRapMother(0),
fTreeVariablePID(0),
fTreeVariablePIDPositive(0),
fTreeVariablePIDNegative(0),
fTreeVariablePIDMother(0),
fTreeVariablePrimaryStatus(0),
fTreeVariablePrimaryStatusMother(0),

//thnsparse
fPtBinNplusNminusCh(NULL),
fPtBinNplusNminusChTruth(NULL),

//Histos
fHistEventCounter(0),
fHistCentrality(0),
//V0s
fHistGeneratedPtVsYVsCentralityK0Short(0),
fHistGeneratedPtVsYVsCentralityLambda(0),
fHistGeneratedPtVsYVsCentralityAntiLambda(0),

//Added
f3dHistInvMassVsPtVsCentLambda(0),
f3dHistInvMassVsPtVsCentAntiLambda(0),
fhistLambdaRecPT(0),
fhistAntiLambdaRecPT(0),
fhistLambdaRecpri(0),
fhistAntiLambdaRecpri(0),
fHistGeneratedCentralityVsLambda(0),
fHistGeneratedCentralityVsAntiLambda(0)

{

    
    fkSaveEventTree    = lSaveEventTree;
    fkSaveV0Tree       = lSaveV0Tree;
    DefineOutput(1, TList::Class()); // Basic Histograms
//     DefineOutput(2, TTree::Class()); // Event Tree output
//     DefineOutput(3, TTree::Class()); // V0 Tree output

    
    
    //Special Debug Options (more to be added as needed)
    // A - Study Wrong PID for tracking bug
    if ( lExtraOptions.Contains("A") ) fkDebugWrongPIDForTracking = kTRUE;
    if ( lExtraOptions.Contains("B") ) fkDebugBump                = kTRUE;
}

AliAnalysisTaskNetLambdaMCTrad::~AliAnalysisTaskNetLambdaMCTrad()
{
    //------------------------------------------------
    // DESTRUCTOR
    //------------------------------------------------

    //Destroy output objects if present
    if (fListHist) {
        delete fListHist;
        fListHist = 0x0;
    }
    
   
    if (fTreeEvent) {
        delete fTreeEvent;
        fTreeEvent = 0x0;
    }
    if (fTreeV0) {
        delete fTreeV0;
        fTreeV0 = 0x0;
    }
    
    if (fUtils) {
        delete fUtils;
        fUtils = 0x0;
    }
    if (fRand) {
        delete fRand;
        fRand = 0x0;
    }
    delete [] fPtArray;

}

//________________________________________________________________________
void AliAnalysisTaskNetLambdaMCTrad::UserCreateOutputObjects()
{
    //------------------------------------------------
    // fTreeEvent: EbyE information
    //------------------------------------------------
    if(fkSaveEventTree){
        fTreeEvent = new TTree("fTreeEvent","Event");
        //Branch Definitions
        fTreeEvent->Branch("fCentrality",&fCentrality,"fCentrality/F");
        fTreeEvent->Branch("fMVPileupFlag",&fMVPileupFlag,"fMVPileupFlag/O");
    }
    
    //------------------------------------------------
    // fTreeV0: V0 Candidate Information
    //------------------------------------------------
    if(fkSaveV0Tree){
        //Create Basic V0 Output Tree
        fTreeV0 = new TTree( "fTreeV0", "V0 Candidates");
        //-----------BASIC-INFO---------------------------
        fTreeV0->Branch("fTreeVariableChi2V0",&fTreeVariableChi2V0,"fTreeVariableChi2V0/F");
        fTreeV0->Branch("fTreeVariableDcaV0Daughters",&fTreeVariableDcaV0Daughters,"fTreeVariableDcaV0Daughters/F");
        fTreeV0->Branch("fTreeVariableDcaV0ToPrimVertex",&fTreeVariableDcaV0ToPrimVertex,"fTreeVariableDcaV0ToPrimVertex/F");
        fTreeV0->Branch("fTreeVariableDcaPosToPrimVertex",&fTreeVariableDcaPosToPrimVertex,"fTreeVariableDcaPosToPrimVertex/F");
        fTreeV0->Branch("fTreeVariableDcaNegToPrimVertex",&fTreeVariableDcaNegToPrimVertex,"fTreeVariableDcaNegToPrimVertex/F");
        fTreeV0->Branch("fTreeVariableV0Radius",&fTreeVariableV0Radius,"fTreeVariableV0Radius/F");
        fTreeV0->Branch("fTreeVariablePt",&fTreeVariablePt,"fTreeVariablePt/F");
        fTreeV0->Branch("fTreeVariablePtMC",&fTreeVariablePtMC,"fTreeVariablePtMC/F");
        fTreeV0->Branch("fTreeVariableRapK0Short",&fTreeVariableRapK0Short,"fTreeVariableRapK0Short/F");
        fTreeV0->Branch("fTreeVariableRapLambda",&fTreeVariableRapLambda,"fTreeVariableRapLambda/F");
        fTreeV0->Branch("fTreeVariableRapMC",&fTreeVariableRapMC,"fTreeVariableRapMC/F");
        fTreeV0->Branch("fTreeVariableInvMassK0s",&fTreeVariableInvMassK0s,"fTreeVariableInvMassK0s/F");
        fTreeV0->Branch("fTreeVariableInvMassLambda",&fTreeVariableInvMassLambda,"fTreeVariableInvMassLambda/F");
        fTreeV0->Branch("fTreeVariableInvMassAntiLambda",&fTreeVariableInvMassAntiLambda,"fTreeVariableInvMassAntiLambda/F");
        fTreeV0->Branch("fTreeVariableV0CosineOfPointingAngle",&fTreeVariableV0CosineOfPointingAngle,"fTreeVariableV0CosineOfPointingAngle/F");
        fTreeV0->Branch("fTreeVariableAlphaV0",&fTreeVariableAlphaV0,"fTreeVariableAlphaV0/F");
        fTreeV0->Branch("fTreeVariablePtArmV0",&fTreeVariablePtArmV0,"fTreeVariablePtArmV0/F");
        fTreeV0->Branch("fTreeVariableLeastNbrCrossedRows",&fTreeVariableLeastNbrCrossedRows,"fTreeVariableLeastNbrCrossedRows/I");
        fTreeV0->Branch("fTreeVariableLeastRatioCrossedRowsOverFindable",&fTreeVariableLeastRatioCrossedRowsOverFindable,"fTreeVariableLeastRatioCrossedRowsOverFindable/F");
        fTreeV0->Branch("fTreeVariableMaxChi2PerCluster",&fTreeVariableMaxChi2PerCluster,"fTreeVariableMaxChi2PerCluster/F");
        fTreeV0->Branch("fTreeVariableMinTrackLength",&fTreeVariableMinTrackLength,"fTreeVariableMinTrackLength/F");
        fTreeV0->Branch("fTreeVariableDistOverTotMom",&fTreeVariableDistOverTotMom,"fTreeVariableDistOverTotMom/F");
        fTreeV0->Branch("fTreeVariableNSigmasPosProton",&fTreeVariableNSigmasPosProton,"fTreeVariableNSigmasPosProton/F");
        fTreeV0->Branch("fTreeVariableNSigmasPosPion",&fTreeVariableNSigmasPosPion,"fTreeVariableNSigmasPosPion/F");
        fTreeV0->Branch("fTreeVariableNSigmasNegProton",&fTreeVariableNSigmasNegProton,"fTreeVariableNSigmasNegProton/F");
        fTreeV0->Branch("fTreeVariableNSigmasNegPion",&fTreeVariableNSigmasNegPion,"fTreeVariableNSigmasNegPion/F");
        fTreeV0->Branch("fTreeVariableNegEta",&fTreeVariableNegEta,"fTreeVariableNegEta/F");
        fTreeV0->Branch("fTreeVariablePosEta",&fTreeVariablePosEta,"fTreeVariablePosEta/F");
        //-----------MULTIPLICITY-INFO--------------------
        fTreeV0->Branch("fTreeVariableCentrality",&fTreeVariableCentrality,"fTreeVariableCentrality/F");
        fTreeV0->Branch("fTreeVariableMVPileupFlag",&fTreeVariableMVPileupFlag,"fTreeVariableMVPileupFlag/O");
        //------------------------------------------------
        if ( fkDebugWrongPIDForTracking ){
            fTreeV0->Branch("fTreeVariablePosPIDForTracking",&fTreeVariablePosPIDForTracking,"fTreeVariablePosPIDForTracking/I");
            fTreeV0->Branch("fTreeVariableNegPIDForTracking",&fTreeVariableNegPIDForTracking,"fTreeVariableNegPIDForTracking/I");
            fTreeV0->Branch("fTreeVariablePosdEdx",&fTreeVariablePosdEdx,"fTreeVariablePosdEdx/F");
            fTreeV0->Branch("fTreeVariableNegdEdx",&fTreeVariableNegdEdx,"fTreeVariableNegdEdx/F");
            fTreeV0->Branch("fTreeVariablePosInnerP",&fTreeVariablePosInnerP,"fTreeVariablePosInnerP/F");
            fTreeV0->Branch("fTreeVariableNegInnerP",&fTreeVariableNegInnerP,"fTreeVariableNegInnerP/F");
            fTreeV0->Branch("fTreeVariableNegTrackStatus",&fTreeVariableNegTrackStatus,"fTreeVariableNegTrackStatus/l");
            fTreeV0->Branch("fTreeVariablePosTrackStatus",&fTreeVariablePosTrackStatus,"fTreeVariablePosTrackStatus/l");
            fTreeV0->Branch("fTreeVariableNegDCAz",&fTreeVariableNegDCAz,"fTreeVariableNegDCAz/F");
            fTreeV0->Branch("fTreeVariablePosDCAz",&fTreeVariablePosDCAz,"fTreeVariablePosDCAz/F");
        }
        //-----------MC Exclusive info--------------------
        fTreeV0->Branch("fTreeVariablePtMother",&fTreeVariablePtMother,"fTreeVariablePtMother/F");
        fTreeV0->Branch("fTreeVariableRapMother",&fTreeVariableRapMother,"fTreeVariableRapMother/F");
        fTreeV0->Branch("fTreeVariablePID",&fTreeVariablePID,"fTreeVariablePID/I");
        fTreeV0->Branch("fTreeVariablePIDPositive",&fTreeVariablePIDPositive,"fTreeVariablePIDPositive/I");
        fTreeV0->Branch("fTreeVariablePIDNegative",&fTreeVariablePIDNegative,"fTreeVariablePIDNegative/I");
        fTreeV0->Branch("fTreeVariablePIDMother",&fTreeVariablePIDMother,"fTreeVariablePIDMother/I");
        fTreeV0->Branch("fTreeVariablePrimaryStatus",&fTreeVariablePrimaryStatus,"fTreeVariablePrimaryStatus/I");
        fTreeV0->Branch("fTreeVariablePrimaryStatusMother",&fTreeVariablePrimaryStatusMother,"fTreeVariablePrimaryStatusMother/I");
        //------------------------------------------------
    }


    //------------------------------------------------
    // Particle Identification Setup
    //------------------------------------------------
    
    
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    if(fIsMC)
    {
        AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
        if(!mcH) return;
    }
    
    fPIDResponse = inputHandler->GetPIDResponse();
    inputHandler->SetNeedField();

    // Multiplicity
    if(! fESDtrackCuts ) {
        fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE,kFALSE);
        fESDtrackCuts->SetPtRange(0.15);
        fESDtrackCuts->SetEtaRange(-1.0, 1.0);
    }
    //Analysis Utils
    if(! fUtils ) {
        fUtils = new AliAnalysisUtils();
    }
    if(! fRand ){
        fRand = new TRandom3();
        // From TRandom3 reference:
        // if seed is 0 (default value) a TUUID is generated and
        // used to fill the first 8 integers of the seed array.
        fRand->SetSeed(0);
    }

    //------------------------------------------------
    // V0 Multiplicity Histograms
    //------------------------------------------------

    // Create histograms
    fListHist = new TList();
    fListHist->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner

    // event hists
    if(! fHistEventCounter ) {
        //Histogram Output: Event-by-Event
        fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",2,0,2);
        fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
        fHistEventCounter->GetXaxis()->SetBinLabel(2, "Selected");
        fListHist->Add(fHistEventCounter);
    }
    
    
    if(! fHistCentrality ) {
        //Histogram Output: Event-by-Event
        fHistCentrality = new TH1D( "fHistCentrality", "WARNING: no pileup rejection applied!;Centrality;Event Count",100,0,100);
        fListHist->Add(fHistCentrality);
    }
    
    

    const Int_t xNbins = 100;
    Double_t xBinEdge[xNbins+1];
    for(Int_t iBin = 0 ; iBin <= xNbins; iBin++){
        xBinEdge[iBin] = iBin - 0.5;
    }
    

    
    Double_t pidPtBins[20] = {1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9 , 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4 ,3.6, 3.7};
    fPtArray = new Double_t[20];
    fPtArray = pidPtBins;
    
    ///V0 Pt with Inv mass to check cuts
    if(! f3dHistInvMassVsPtVsCentLambda) {
        f3dHistInvMassVsPtVsCentLambda = new TH3F("f3dHistInvMassVsPtVsCentLambda","Inv Mass Lambda Vs Pt Vs Centrality",100,1.05,1.2,200,0,20,100,0,100);
        fListHist->Add(f3dHistInvMassVsPtVsCentLambda);
    }
    if(! f3dHistInvMassVsPtVsCentAntiLambda) {
        f3dHistInvMassVsPtVsCentAntiLambda = new TH3F("f3dHistInvMassVsPtVsCentAntiLambda","Inv Mass AntiLambda Vs Pt Vs Centrality",100,1.05,1.2,200,0,20,100,0,100);
        fListHist->Add(f3dHistInvMassVsPtVsCentAntiLambda);
    }
    
    
    ///VO pt
        if(! fhistLambdaRecPT) {
            fhistLambdaRecPT = new TH2F("fhistLambdaRecPT","fhistLambdaRecPT",xNbins, xBinEdge, fNptBins, fPtArray);
            fListHist->Add(fhistLambdaRecPT);
        }
    
    if(! fhistAntiLambdaRecPT) {
        fhistAntiLambdaRecPT = new TH2F("fhistAntiLambdaRecPT","fhistAntiLambdaRecPT",xNbins, xBinEdge, fNptBins, fPtArray);
        fListHist->Add(fhistAntiLambdaRecPT);
    }
    
    //Rec primaries & gen for L & A
    if(fIsMC){

        if(! fhistLambdaRecpri) {
            fhistLambdaRecpri = new TH2F("fhistLambdaRecpri","fhistLambdaRecpri",xNbins, xBinEdge, fNptBins, fPtArray);
            fListHist->Add(fhistLambdaRecpri);
        }
        if(! fhistAntiLambdaRecpri) {
            fhistAntiLambdaRecpri = new TH2F("fhistAntiLambdaRecpri","fhistAntiLambdaRecpri",xNbins, xBinEdge, fNptBins, fPtArray);
            fListHist->Add(fhistAntiLambdaRecpri);
        }
        
        //Generated 3D hists
        
        if(! fHistGeneratedPtVsYVsCentralityK0Short ) {
            //Histogram Output: Efficiency Denominator
            fHistGeneratedPtVsYVsCentralityK0Short = new TH3D( "fHistGeneratedPtVsYVsCentralityK0Short", ";pT;y;centrality",200,0,20,20,-1.0,1.0,100,0,100);
            fListHist->Add(fHistGeneratedPtVsYVsCentralityK0Short);
        }
        if(! fHistGeneratedPtVsYVsCentralityLambda ) {
            //Histogram Output: Efficiency Denominator
            fHistGeneratedPtVsYVsCentralityLambda = new TH3D( "fHistGeneratedPtVsYVsCentralityLambda", ";pT;y;centrality",200,0,20,20,-1.0,1.0,100,0,100);
            fListHist->Add(fHistGeneratedPtVsYVsCentralityLambda);
        }
        if(! fHistGeneratedPtVsYVsCentralityAntiLambda ) {
            //Histogram Output: Efficiency Denominator
            fHistGeneratedPtVsYVsCentralityAntiLambda = new TH3D( "fHistGeneratedPtVsYVsCentralityAntiLambda", ";pT;y;centrality",200,0,20,20,-1.0,1.0,100,0,100);
            fListHist->Add(fHistGeneratedPtVsYVsCentralityAntiLambda);
        }
        
        //Generated added 2D hists
        
        if(! fHistGeneratedCentralityVsLambda ) {
            fHistGeneratedCentralityVsLambda = new TH2F( "fHistGeneratedCentralityVsLambda", "fHistGeneratedCentralityVsLambda",xNbins, xBinEdge, fNptBins, fPtArray);
            fListHist->Add(fHistGeneratedCentralityVsLambda);
        }
        
        if(! fHistGeneratedCentralityVsAntiLambda ) {
            fHistGeneratedCentralityVsAntiLambda = new TH2F( "fHistGeneratedCentralityVsAntiLambda", "fHistGeneratedCentralityVsAntiLambda",xNbins, xBinEdge, fNptBins, fPtArray);
            fListHist->Add(fHistGeneratedCentralityVsAntiLambda);
        }
        
    
    }

    
    //THN SPARSE 
    
            const Int_t dim = 39; //1 centrality bin + (19 pt bins) * 2
            
            Int_t bin[dim]    = { 100,
                500, 500, 500,
                500, 500, 500, 500, 500, 500, 500, 500,
                200, 200, 200, 200, 200, 200, 200, 200,
                500, 500, 500,
                500, 500, 500, 500, 500, 500, 500, 500,
                200, 200, 200, 200, 200, 200, 200, 200 };
            
            Double_t min[dim] = { -0.5,
                -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
                -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
                -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
                -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
                -0.5, -0.5};
            
            Double_t max[dim] = { 99.5,
                499.5, 499.5, 499.5,
                499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5,
                199.5, 199.5, 199.5, 199.5, 199.5, 199.5, 199.5, 199.5,
                499.5, 499.5, 499.5,
                499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5,
                199.5, 199.5, 199.5, 199.5, 199.5, 199.5, 199.5, 199.5 };
            
            fPtBinNplusNminusCh = new THnSparseI("fPtBinNplusNminusCh","cent-nplus-nminus", dim, bin, min, max);
            fListHist->Add(fPtBinNplusNminusCh);
            
            if( fIsMC ){
                fPtBinNplusNminusChTruth = new THnSparseI("fPtBinNplusNminusChTruth","cent-nplus-nminus", dim, bin, min, max);
                fListHist->Add(fPtBinNplusNminusChTruth);
            }

    
    PostData(1, fListHist);
//     PostData(2, fTreeEvent);
//     PostData(3, fTreeV0);

}// end UserCreateOutputObjects


//________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskNetLambdaMCTrad::UserExec(Option_t *)
{
    const Int_t dim = fNptBins*2;
    Int_t ptCh[dim];
    Int_t ptChMC[dim];
    for(Int_t idx = 0; idx < dim; idx++){
        ptCh[idx] = 0.;
        ptChMC[idx] = 0;
    }
    

    AliESDEvent *lESDevent = 0x0;
    AliMCEvent  *lMCevent  = 0x0;
    AliStack    *lMCstack  = 0x0;

   
    lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
    if (!lESDevent) {
        AliWarning("ERROR: lESDevent not available \n");
        return;
    }

    //================================================================================================================
    // Monte Carlo-related information
    lMCevent = MCEvent();
    if (!lMCevent) {
        Printf("ERROR: Could not retrieve MC event \n");
        cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
        return;
    }
    lMCstack = lMCevent->Stack();
    if (!lMCstack) {
        Printf("ERROR: Could not retrieve MC stack \n");
        cout << "Name of the file with pb :" <<  fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
        return;
    }
    //=================================================================================================================
    
    //Get VZERO Information for multiplicity later
    AliVVZERO* esdV0 = lESDevent->GetVZEROData();
    if (!esdV0) {
        AliError("AliVVZERO not available");
        return;
    }

    Double_t lMagneticField = -10;
    lMagneticField = lESDevent->GetMagneticField( );


    fHistEventCounter->Fill(0.5);

    //------------------------------------------------
    // Primary Vertex Requirements Section:
    //------------------------------------------------

    //classical Proton-proton like selection
    const AliESDVertex *lPrimaryBestESDVtx     = lESDevent->GetPrimaryVertex();
    const AliESDVertex *lPrimaryTrackingESDVtx = lESDevent->GetPrimaryVertexTracks();
    const AliESDVertex *lPrimarySPDVtx         = lESDevent->GetPrimaryVertexSPD();

    Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
    lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );
    
    if(lBestPrimaryVtxPos[2] < -10. || lBestPrimaryVtxPos[2] > 10.) return;


    //------------------------------------------------
    // Multiplicity Information Acquistion
    //------------------------------------------------

    Float_t lPercentile = 500;
    Float_t lPercentileEmbeddedSelection = 500;
    Int_t lEvSelCode = 100;
    AliMultSelection *MultSelection = (AliMultSelection*) lESDevent -> FindListObject("MultSelection");
    if( !MultSelection) {
        //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
        AliWarning("AliMultSelection object not found!");
    } else {
        //V0M Multiplicity Percentile
        lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
        lPercentileEmbeddedSelection = MultSelection->GetMultiplicityPercentile("V0M", kTRUE );
        //Event Selection Code
        lEvSelCode = MultSelection->GetEvSelCode();
    }
    
    fMVPileupFlag = kFALSE;
    fMVPileupFlag = MultSelection->GetThisEventIsNotPileupMV();
    
    fCentrality = lPercentile;
    
    if( lEvSelCode != 0 ) {
        PostData(1, fListHist);
//         PostData(2, fTreeEvent);
//         PostData(3, fTreeV0);
        return;

    }

    
    if( fCentrality < 0 || fCentrality >=80 ) return;
    
    fHistEventCounter->Fill(1.5);
    fHistCentrality->Fill(fCentrality);
    
    //Event-level fill
    if ( fkSaveEventTree ) fTreeEvent->Fill() ;


    //-------------------------------------------------------------

    //----- Loop on Generated Particles ---------------------------
    Int_t    lThisPDG  = 0;
    Double_t lThisRap  = 0;
    Double_t lThisPt   = 0;
    
    if (fIsMC) {

    for (Int_t ilab = 0;  ilab < (lMCstack->GetNtrack()); ilab++)
    {   // This is the begining of the loop on tracks
        
        TParticle* lPart = 0x0;
        lPart = lMCstack->Particle( ilab );
        if(!lPart) {
            Printf("Generated loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", ilab );
            continue;
        }
        
        lThisPDG = lPart->GetPdgCode();
        
        //This if is necessary in some situations (rapidity calculation and PYTHIA junctions, etc)
        if ( (TMath::Abs(lThisPDG) == 3312) || (TMath::Abs(lThisPDG) == 3334) || (TMath::Abs(lThisPDG) == 3122) || lThisPDG == 310 )
        {
            lThisRap   = MyRapidity(lPart->Energy(),lPart->Pz());
            lThisPt    = lPart->Pt();
            
            //Use Physical Primaries only for filling These Histos
            if ( lMCstack->IsPhysicalPrimary(ilab)!=kTRUE ) continue;
            if (TMath::Abs(lThisRap) > 0.5) continue;
            
            Int_t iptbinMC = GetPtBin(lThisPt);
            if( iptbinMC < 0 || iptbinMC > fNptBins-1 ) continue;
            
            if( lThisPDG ==   310 ) {
                fHistGeneratedPtVsYVsCentralityK0Short       -> Fill (lThisPt, lThisRap, lPercentileEmbeddedSelection);
            }
            if( lThisPDG ==  3122 ) {
                fHistGeneratedPtVsYVsCentralityLambda       -> Fill (lThisPt, lThisRap, lPercentileEmbeddedSelection);
                fHistGeneratedCentralityVsLambda            -> Fill (fCentrality, lThisPt);
                ptChMC[iptbinMC] += 1;

            }
            if( lThisPDG == -3122 ) {
                fHistGeneratedPtVsYVsCentralityAntiLambda       -> Fill (lThisPt, lThisRap, lPercentileEmbeddedSelection);
                fHistGeneratedCentralityVsAntiLambda            -> Fill (fCentrality, lThisPt);
                ptChMC[iptbinMC+fNptBins] += 1;

            }
           
        }
    }//End of loop on tracks
        
        Double_t ptContainerMC[dim+1];
        ptContainerMC[0] = (Double_t) fCentrality;
        for(Int_t i = 1; i <= dim; i++){
            ptContainerMC[i] = ptChMC[i-1];
        }
        
        fPtBinNplusNminusChTruth->Fill(ptContainerMC);
    }
    
    
    //------------------------------------------------
    // Fill V0 Tree as needed
    //------------------------------------------------

    //Variable definition
    Int_t    lOnFlyStatus = 0;// nv0sOn = 0, nv0sOff = 0;
    Double_t lChi2V0 = 0;
    Double_t lDcaV0Daughters = 0, lDcaV0ToPrimVertex = 0;
    Double_t lDcaPosToPrimVertex = 0, lDcaNegToPrimVertex = 0;
    Double_t lV0CosineOfPointingAngle = 0;
    Double_t lV0Radius = 0, v0DecayLength = 0, lPt = 0;
    Double_t lRapK0Short = 0, lRapLambda = 0;
    Double_t lInvMassK0s = 0, lInvMassLambda = 0, lInvMassAntiLambda = 0;
    Double_t lAlphaV0 = 0, lPtArmV0 = 0;

    Double_t fMinV0Pt = 0.5;
    Double_t fMaxV0Pt = 10;

    Int_t nv0s = 0;
    nv0s = lESDevent->GetNumberOfV0s();

    for (Int_t iV0 = 0; iV0 < nv0s; iV0++) //extra-crazy test
    {   // This is the begining of the V0 loop
        AliESDv0 *v0 = ((AliESDEvent*)lESDevent)->GetV0(iV0);
        if (!v0) continue;

        Double_t tDecayVertexV0[3];
        v0->GetXYZ(tDecayVertexV0[0],tDecayVertexV0[1],tDecayVertexV0[2]);

        Double_t tV0mom[3];
        v0->GetPxPyPz( tV0mom[0],tV0mom[1],tV0mom[2] );
        Double_t lV0TotalMomentum = TMath::Sqrt(
                                        tV0mom[0]*tV0mom[0]+tV0mom[1]*tV0mom[1]+tV0mom[2]*tV0mom[2] );

        lV0Radius = TMath::Sqrt(tDecayVertexV0[0]*tDecayVertexV0[0]+tDecayVertexV0[1]*tDecayVertexV0[1]);
        
        v0DecayLength = TMath::Sqrt(TMath::Power(tDecayVertexV0[0] - lBestPrimaryVtxPos[0],2) +
                                    TMath::Power(tDecayVertexV0[1] - lBestPrimaryVtxPos[1],2) +
                                    TMath::Power(tDecayVertexV0[2] - lBestPrimaryVtxPos[2],2 ));

        lPt = v0->Pt();
        lRapK0Short = v0->RapK0Short();
        lRapLambda  = v0->RapLambda();
        if ((lPt<fMinV0Pt)||(fMaxV0Pt<lPt)) continue;

        UInt_t lKeyPos = (UInt_t)TMath::Abs(v0->GetPindex());
        UInt_t lKeyNeg = (UInt_t)TMath::Abs(v0->GetNindex());

        Double_t lMomPos[3];
        v0->GetPPxPyPz(lMomPos[0],lMomPos[1],lMomPos[2]);
        Double_t lMomNeg[3];
        v0->GetNPxPyPz(lMomNeg[0],lMomNeg[1],lMomNeg[2]);

        AliESDtrack *pTrack=((AliESDEvent*)lESDevent)->GetTrack(lKeyPos);
        AliESDtrack *nTrack=((AliESDEvent*)lESDevent)->GetTrack(lKeyNeg);
        if (!pTrack || !nTrack) {
            Printf("ERROR: Could not retreive one of the daughter track");
            continue;
        }
        fTreeVariablePosPIDForTracking = pTrack->GetPIDForTracking();
        fTreeVariableNegPIDForTracking = nTrack->GetPIDForTracking();
        
        const AliExternalTrackParam *innernegv0=nTrack->GetInnerParam();
        const AliExternalTrackParam *innerposv0=pTrack->GetInnerParam();
        Float_t lThisPosInnerP = -1;
        Float_t lThisNegInnerP = -1;
        if(innerposv0)  { lThisPosInnerP  = innerposv0 ->GetP(); }
        if(innernegv0)  { lThisNegInnerP  = innernegv0 ->GetP(); }
        Float_t lThisPosdEdx = pTrack -> GetTPCsignal();
        Float_t lThisNegdEdx = nTrack -> GetTPCsignal();
        
        fTreeVariablePosdEdx = lThisPosdEdx;
        fTreeVariableNegdEdx = lThisNegdEdx;
        
        fTreeVariablePosInnerP = lThisPosInnerP;
        fTreeVariableNegInnerP = lThisNegInnerP;
        
        fTreeVariableNegEta = nTrack->Eta();
        fTreeVariablePosEta = pTrack->Eta();

        if ( pTrack->GetSign() == nTrack->GetSign())
        {
            continue;
        }

        //________________________________________________________________________
        // Track quality cuts
        Float_t lPosTrackCrossedRows = pTrack->GetTPCClusterInfo(2,1);
        Float_t lNegTrackCrossedRows = nTrack->GetTPCClusterInfo(2,1);
        fTreeVariableLeastNbrCrossedRows = (Int_t) lPosTrackCrossedRows;
        if( lNegTrackCrossedRows < fTreeVariableLeastNbrCrossedRows )
            fTreeVariableLeastNbrCrossedRows = (Int_t) lNegTrackCrossedRows;

        // TPC refit condition (done during reconstruction for Offline but not for On-the-fly)
        if( !(pTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
        if( !(nTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;

        //Get status flags
        fTreeVariablePosTrackStatus = pTrack->GetStatus();
        fTreeVariableNegTrackStatus = nTrack->GetStatus();
        
        fTreeVariablePosDCAz = GetDCAz(pTrack);
        fTreeVariableNegDCAz = GetDCAz(nTrack);

        if ( ( ( pTrack->GetTPCClusterInfo(2,1) ) < 70 ) || ( ( nTrack->GetTPCClusterInfo(2,1) ) < 70 ) ) continue;

        //GetKinkIndex condition
        if( pTrack->GetKinkIndex(0)>0 || nTrack->GetKinkIndex(0)>0 ) continue;

        //Findable clusters > 0 condition
        if( pTrack->GetTPCNclsF()<=0 || nTrack->GetTPCNclsF()<=0 ) continue;

        //Compute ratio Crossed Rows / Findable clusters
        //Note: above test avoids division by zero!
        Float_t lPosTrackCrossedRowsOverFindable = lPosTrackCrossedRows / ((double)(pTrack->GetTPCNclsF()));
        Float_t lNegTrackCrossedRowsOverFindable = lNegTrackCrossedRows / ((double)(nTrack->GetTPCNclsF()));

        fTreeVariableLeastRatioCrossedRowsOverFindable = lPosTrackCrossedRowsOverFindable;
        if( lNegTrackCrossedRowsOverFindable < fTreeVariableLeastRatioCrossedRowsOverFindable )
            fTreeVariableLeastRatioCrossedRowsOverFindable = lNegTrackCrossedRowsOverFindable;

        //Lowest Cut Level for Ratio Crossed Rows / Findable = 0.8, set here
        if ( fTreeVariableLeastRatioCrossedRowsOverFindable < 0.8 ) continue;

        //Extra track quality: Chi2/cluster for cross-checks
        Float_t lBiggestChi2PerCluster = -1;
        
        Float_t lPosChi2PerCluster = 1000;
        Float_t lNegChi2PerCluster = 1000;
        
        if( pTrack->GetTPCNcls() > 0 ) lPosChi2PerCluster = pTrack->GetTPCchi2() / ((Float_t)pTrack->GetTPCNcls());
        if( nTrack->GetTPCNcls() > 0 ) lNegChi2PerCluster = nTrack->GetTPCchi2() / ((Float_t)nTrack->GetTPCNcls());
        
        if ( lPosChi2PerCluster  > lBiggestChi2PerCluster ) lBiggestChi2PerCluster = lPosChi2PerCluster;
        if ( lNegChi2PerCluster  > lBiggestChi2PerCluster ) lBiggestChi2PerCluster = lNegChi2PerCluster;
        
        fTreeVariableMaxChi2PerCluster = lBiggestChi2PerCluster;
        
        //Extra track quality: min track length
        Float_t lSmallestTrackLength = 1000;
        Float_t lPosTrackLength = -1;
        Float_t lNegTrackLength = -1;
        
        if (pTrack->GetInnerParam()) lPosTrackLength = pTrack->GetLengthInActiveZone(1, 2.0, 220.0, lESDevent->GetMagneticField());
        if (nTrack->GetInnerParam()) lNegTrackLength = nTrack->GetLengthInActiveZone(1, 2.0, 220.0, lESDevent->GetMagneticField());
        
        if ( lPosTrackLength  < lSmallestTrackLength ) lSmallestTrackLength = lPosTrackLength;
        if ( lNegTrackLength  < lSmallestTrackLength ) lSmallestTrackLength = lNegTrackLength;
        
        fTreeVariableMinTrackLength = lSmallestTrackLength;
        
        //End track Quality Cuts
        //__________________________________________________________________________________________________________________________

        lDcaPosToPrimVertex = TMath::Abs(pTrack->GetD(lBestPrimaryVtxPos[0],
                                         lBestPrimaryVtxPos[1],
                                         lMagneticField) );

        lDcaNegToPrimVertex = TMath::Abs(nTrack->GetD(lBestPrimaryVtxPos[0],
                                         lBestPrimaryVtxPos[1],
                                         lMagneticField) );

        lOnFlyStatus = v0->GetOnFlyStatus();
        lChi2V0 = v0->GetChi2V0();
        lDcaV0Daughters = v0->GetDcaV0Daughters();
        lDcaV0ToPrimVertex = v0->GetD(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]);
        lV0CosineOfPointingAngle = v0->GetV0CosineOfPointingAngle(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]);
        fTreeVariableV0CosineOfPointingAngle=lV0CosineOfPointingAngle;

        // Getting invariant mass infos directly from ESD
        v0->ChangeMassHypothesis(310);
        lInvMassK0s = v0->GetEffMass();
        v0->ChangeMassHypothesis(3122);
        lInvMassLambda = v0->GetEffMass();
        v0->ChangeMassHypothesis(-3122);
        lInvMassAntiLambda = v0->GetEffMass();
        lAlphaV0 = v0->AlphaV0();
        lPtArmV0 = v0->PtArmV0();
        
        fTreeVariableMVPileupFlag = fMVPileupFlag;

        fTreeVariablePt = v0->Pt();
        fTreeVariableChi2V0 = lChi2V0;
        fTreeVariableDcaV0ToPrimVertex = lDcaV0ToPrimVertex;
        fTreeVariableDcaV0Daughters = lDcaV0Daughters;
        fTreeVariableV0CosineOfPointingAngle = lV0CosineOfPointingAngle;
        fTreeVariableV0Radius = lV0Radius;
        fTreeVariableDcaPosToPrimVertex = lDcaPosToPrimVertex;
        fTreeVariableDcaNegToPrimVertex = lDcaNegToPrimVertex;
        fTreeVariableInvMassK0s = lInvMassK0s;
        fTreeVariableInvMassLambda = lInvMassLambda;
        fTreeVariableInvMassAntiLambda = lInvMassAntiLambda;
        fTreeVariableRapK0Short = lRapK0Short;
        fTreeVariableRapLambda = lRapLambda;
        fTreeVariableAlphaV0 = lAlphaV0;
        fTreeVariablePtArmV0 = lPtArmV0;

        //Official means of acquiring N-sigmas
        fTreeVariableNSigmasPosProton = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kProton );
        fTreeVariableNSigmasPosPion   = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kPion );
        fTreeVariableNSigmasNegProton = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kProton );
        fTreeVariableNSigmasNegPion   = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kPion );

        //This requires an Invariant Mass Hypothesis afterwards
        fTreeVariableDistOverTotMom = TMath::Sqrt(
                                          TMath::Power( tDecayVertexV0[0] - lBestPrimaryVtxPos[0] , 2) +
                                          TMath::Power( tDecayVertexV0[1] - lBestPrimaryVtxPos[1] , 2) +
                                          TMath::Power( tDecayVertexV0[2] - lBestPrimaryVtxPos[2] , 2)
                                      );
        fTreeVariableDistOverTotMom /= (lV0TotalMomentum+1e-10); //avoid division by zero, to be sure

        //Copy Multiplicity information
        fTreeVariableCentrality = fCentrality;

        //===============================================================================================
        // V0 Monte Carlo Association starts here
        //===============================================================================================
        
        //---> Set Everything to "I don't know" before starting
        //---> This is strictly necessary!
        
        fTreeVariablePIDPositive = 0;
        fTreeVariablePIDNegative = 0;
        
        fTreeVariablePtMother = -1;
        fTreeVariableRapMother = -100;
        fTreeVariablePtMC = -1;
        fTreeVariableRapMC = -100;
        
        fTreeVariablePID = -1;
        fTreeVariablePIDMother = -1;
        
        fTreeVariablePrimaryStatus = 0;
        fTreeVariablePrimaryStatusMother = 0;
        
        //fTreeVariablePosTransvMomentumMC = -1;
        //fTreeVariableNegTransvMomentumMC = -1;
        
        Int_t lblPosV0Dghter = (Int_t) TMath::Abs( pTrack->GetLabel() );
        Int_t lblNegV0Dghter = (Int_t) TMath::Abs( nTrack->GetLabel() );
        
        TParticle* mcPosV0Dghter = lMCstack->Particle( lblPosV0Dghter );
        if(!mcPosV0Dghter) continue;

        TParticle* mcNegV0Dghter = lMCstack->Particle( lblNegV0Dghter );
        if(!mcNegV0Dghter) continue;
        
        Int_t lPIDPositive = mcPosV0Dghter -> GetPdgCode();
        Int_t lPIDNegative = mcNegV0Dghter -> GetPdgCode();
        
        //fTreeVariablePosTransvMomentumMC = mcPosV0Dghter->Pt();
        //fTreeVariableNegTransvMomentumMC = mcNegV0Dghter->Pt();
        
        fTreeVariablePIDPositive = lPIDPositive;
        fTreeVariablePIDNegative = lPIDNegative;
        
        Int_t lblMotherPosV0Dghter = mcPosV0Dghter->GetFirstMother() ;
        Int_t lblMotherNegV0Dghter = mcNegV0Dghter->GetFirstMother();
        
        if( lblMotherPosV0Dghter == lblMotherNegV0Dghter && lblMotherPosV0Dghter > -1 ) {
            //either label is fine, they're equal at this stage
            TParticle* pThisV0 = lMCstack->Particle( lblMotherPosV0Dghter );
            if(!pThisV0) continue;

            //Set tree variables
            fTreeVariablePID   = pThisV0->GetPdgCode(); //PDG Code
            fTreeVariablePtMC  = pThisV0->Pt(); //Perfect Pt
            
            //Only Interested if it's a Lambda, AntiLambda or K0s
            //Avoid the Junction Bug! PYTHIA has particles with Px=Py=Pz=E=0 occasionally,
            //having particle code 88 (unrecognized by PDG), for documentation purposes.
            //Even ROOT's TParticle::Y() is not prepared to deal with that exception!
            //Note that TParticle::Pt() is immune (that would just return 0)...
            //Though granted that that should be extremely rare in this precise condition...
            if( TMath::Abs(fTreeVariablePID) == 3122 || fTreeVariablePID==310 ) {
                fTreeVariableRapMC = pThisV0->Y(); //Perfect Y
            }
            if( lMCstack->IsPhysicalPrimary       (lblMotherPosV0Dghter) ) fTreeVariablePrimaryStatus = 1; //Is Primary!
            if( lMCstack->IsSecondaryFromWeakDecay(lblMotherPosV0Dghter) ) fTreeVariablePrimaryStatus = 2; //Weak Decay!
            if( lMCstack->IsSecondaryFromMaterial (lblMotherPosV0Dghter) ) fTreeVariablePrimaryStatus = 3; //Material Int!
            
            //Now we try to acquire the V0 parent particle, if possible
            Int_t lblThisV0Parent = pThisV0->GetFirstMother();
            if ( lblThisV0Parent > -1 ) { //if it has a parent, get it and store specs
                TParticle* pThisV0Parent = lMCstack->Particle( lblThisV0Parent );
                fTreeVariablePIDMother   = pThisV0Parent->GetPdgCode(); //V0 Mother PDG
                fTreeVariablePtMother    = pThisV0Parent->Pt();         //V0 Mother Pt
                //NOTE: Fill only for charged xi
                if ( TMath::Abs(fTreeVariablePIDMother)==3312) fTreeVariableRapMother   = pThisV0Parent->Y();         //V0 Mother Pt
                //Primary Status for the V0 Mother particle
                if( lMCstack->IsPhysicalPrimary       (lblThisV0Parent) ) fTreeVariablePrimaryStatusMother = 1; //Is Primary!
                if( lMCstack->IsSecondaryFromWeakDecay(lblThisV0Parent) ) fTreeVariablePrimaryStatusMother = 2; //Weak Decay!
                if( lMCstack->IsSecondaryFromMaterial (lblThisV0Parent) ) fTreeVariablePrimaryStatusMother = 3; //Material Int!
            }
        }
        
        //------------------------------------------------
        // Fill Tree!
        //------------------------------------------------

        // The conditionals are meant to decrease excessive
        // memory usage!

        //First Selection: Reject OnFly
        if( lOnFlyStatus == 0 ) {
            //Second Selection: rough 20-sigma band, parametric.
            //K0Short: Enough to parametrize peak broadening with linear function.
            Double_t lUpperLimitK0Short = (5.63707e-01) + (1.14979e-02)*fTreeVariablePt;
            Double_t lLowerLimitK0Short = (4.30006e-01) - (1.10029e-02)*fTreeVariablePt;
            //Lambda: Linear (for higher pt) plus exponential (for low-pt broadening)
            //[0]+[1]*x+[2]*TMath::Exp(-[3]*x)
            Double_t lUpperLimitLambda = (1.13688e+00) + (5.27838e-03)*fTreeVariablePt + (8.42220e-02)*TMath::Exp(-(3.80595e+00)*fTreeVariablePt);
            Double_t lLowerLimitLambda = (1.09501e+00) - (5.23272e-03)*fTreeVariablePt - (7.52690e-02)*TMath::Exp(-(3.46339e+00)*fTreeVariablePt);
            //Do Selection
            if(
                //Case 1: Lambda Selection
                (fTreeVariableInvMassLambda    < lUpperLimitLambda  && fTreeVariableInvMassLambda     > lLowerLimitLambda &&
                 (!fkPreselectDedx || (fkPreselectDedx&&TMath::Abs(fTreeVariableNSigmasPosProton) < 3.0 && TMath::Abs(fTreeVariableNSigmasNegPion) < 3.0) ) &&
                 (!fkPreselectPID || fTreeVariablePID == 3122 )
                )
                ||
                //Case 2: AntiLambda Selection
                (fTreeVariableInvMassAntiLambda < lUpperLimitLambda  && fTreeVariableInvMassAntiLambda > lLowerLimitLambda &&
                 (!fkPreselectDedx || (fkPreselectDedx&&TMath::Abs(fTreeVariableNSigmasNegProton) < 3.0 && TMath::Abs(fTreeVariableNSigmasPosPion) < 3.0) ) &&
                 (!fkPreselectPID || fTreeVariablePID == -3122 )
                )
                ||
                //Case 3: K0Short Selection
                (fTreeVariableInvMassK0s        < lUpperLimitK0Short && fTreeVariableInvMassK0s        > lLowerLimitK0Short &&
                 (!fkPreselectDedx || (fkPreselectDedx&&TMath::Abs(fTreeVariableNSigmasNegPion)   < 3.0 && TMath::Abs(fTreeVariableNSigmasPosPion) < 3.0) ) &&
                 (!fkPreselectPID || fTreeVariablePID == 310 )
                ) ) {
                //Pre-selection in case this is AA...
                    
                    //Random denial
//                    Bool_t lKeepV0 = kTRUE;
//                    if(fkDownScaleV0 && ( fRand->Uniform() > fDownScaleFactorV0 )) lKeepV0 = kFALSE;
//
                if ( TMath::Abs(fTreeVariableNegEta)<0.8 && TMath::Abs(fTreeVariablePosEta)<0.8 && fkSaveV0Tree ) fTreeV0->Fill();
            }
        
        
            Int_t iptbin = GetPtBin(lPt);
            if( iptbin < 0 || iptbin > fNptBins-1 ) continue;
            
            
            // Topological cuts
            
        if(TMath::Abs(fTreeVariableNegEta)       <= 0.8               &&
           TMath::Abs(fTreeVariablePosEta)       <= 0.8               &&
           lV0Radius                 >= 5.0               &&
           lV0Radius                 <= 100.              &&
           lDcaV0Daughters           <= 0.8            &&
           lV0CosineOfPointingAngle    >= 0.999
           )
        {

            //Lambda
            if(
               lDcaNegToPrimVertex       >= 0.1     &&
               lDcaPosToPrimVertex       >= 0.2       &&
               fTreeVariablePID == 3122      &&
               TMath::Abs(fTreeVariableNSigmasNegPion)   <= 3. &&
               TMath::Abs(fTreeVariableNSigmasPosProton) <= 3.
               )
            {
                f3dHistInvMassVsPtVsCentLambda->Fill(lInvMassLambda,fTreeVariablePt,fCentrality);
                fhistLambdaRecPT->Fill(fCentrality,fTreeVariablePt);
                ptCh[iptbin] += 1;
                if(fIsMC)
                {
                   if( fTreeVariablePrimaryStatus == 1 )
                      fhistLambdaRecpri->Fill(fCentrality,fTreeVariablePtMC);
                }

            }
            
            //Anti-Lambda

            if(
               lDcaNegToPrimVertex       >= 0.2     &&
               lDcaPosToPrimVertex       >= 0.1       &&
               fTreeVariablePID == -3122      &&
               TMath::Abs(fTreeVariableNSigmasPosPion)   <= 3. &&
               TMath::Abs(fTreeVariableNSigmasNegProton) <= 3.
               )
            {
                f3dHistInvMassVsPtVsCentAntiLambda->Fill(lInvMassAntiLambda,fTreeVariablePt,fCentrality);
                fhistAntiLambdaRecPT->Fill(fCentrality,fTreeVariablePt);
                ptCh[iptbin+fNptBins] += 1;
                if(fIsMC){
                    if( fTreeVariablePrimaryStatus == 1 )
                    fhistAntiLambdaRecpri->Fill(fCentrality,fTreeVariablePtMC);
                }

            }
        }
    }
        
        
        
        
        //------------------------------------------------
        // Fill V0 tree over.
        //------------------------------------------------

        
    }// This is the end of the V0 loop

    Double_t ptContainer[dim+1];
    ptContainer[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++){
        ptContainer[i] = ptCh[i-1];
    }
    
    fPtBinNplusNminusCh->Fill(ptContainer);
    
    
    
    // Post output data.
    PostData(1, fListHist);
//     PostData(2, fTreeEvent);
//     PostData(3, fTreeV0);
}

//________________________________________________________________________
void AliAnalysisTaskNetLambdaMCTrad::Terminate(Option_t *)
{
    // Draw result to the screen
    // Called once at the end of the query

    TList *cRetrievedList = 0x0;
    cRetrievedList = (TList*)GetOutputData(1);
    if(!cRetrievedList) {
        Printf("ERROR - AliAnalysisTaskStrangenessVsMultiplicityMCRun2 : ouput data container list not available\n");
        return;
    }

    fHistEventCounter = dynamic_cast<TH1D*> (  cRetrievedList->FindObject("fHistEventCounter")  );
    if (!fHistEventCounter) {
        Printf("ERROR - AliAnalysisTaskStrangenessVsMultiplicityMCRun2 : fHistEventCounter not available");
        return;
    }

    TCanvas *canCheck = new TCanvas("AliAnalysisTaskStrangenessVsMultiplicityMCRun2","V0 Multiplicity",10,10,510,510);
    canCheck->cd(1)->SetLogy();

    fHistEventCounter->SetMarkerStyle(22);
    fHistEventCounter->DrawCopy("E");
}

//________________________________________________________________________
Double_t AliAnalysisTaskNetLambdaMCTrad::MyRapidity(Double_t rE, Double_t rPz) const
{
    // Local calculation for rapidity
    Double_t ReturnValue = -100;
    if( (rE-rPz+1.e-13) != 0 && (rE+rPz) != 0 ) {
        ReturnValue =  0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
    }
    return ReturnValue;
}

//---------------------------------------------------------------
Int_t AliAnalysisTaskNetLambdaMCTrad::GetPtBin(Double_t pt){
    
    Int_t bin = -1;
    
    Double_t pidPtBins[20] = {1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9 , 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4 ,3.6, 3.7};
        for(Int_t iBin = 0; iBin < fNptBins; iBin++){
            
            if( iBin == fNptBins-1){
                if( pt >= pidPtBins[iBin] && pt <= pidPtBins[iBin+1]){
                    bin = iBin;
                    break;
                }
            }
            else{
                if( pt >= pidPtBins[iBin] && pt < pidPtBins[iBin+1]){
                    bin = iBin;
                    break;
                    
                }
            }
        }//for
    
    return bin;
    
}
//___________________________________________________________


//________________________________________________________________________
Float_t AliAnalysisTaskNetLambdaMCTrad::GetDCAz(AliESDtrack *lTrack)
//Encapsulation of DCAz calculation
{
    Float_t b[2];
    Float_t bCov[3];
    lTrack->GetImpactParameters(b,bCov);
    if (bCov[0]<=0 || bCov[2]<=0) {
        AliDebug(1, "Estimated b resolution lower or equal to zero!");
        bCov[0]=0; bCov[2]=0;
    }
    Float_t dcaToVertexXY = b[0];
    Float_t dcaToVertexZ = b[1];
    
    return dcaToVertexZ;
}

//________________________________________________________________________
Float_t AliAnalysisTaskNetLambdaMCTrad::GetCosPA(AliESDtrack *lPosTrack, AliESDtrack *lNegTrack, AliESDEvent *lEvent)
//Encapsulation of CosPA calculation (warning: considers AliESDtrack clones)
{
    Float_t lCosPA = -1;
    AliESDtrack* lNegClone = (AliESDtrack*) lNegTrack->Clone("lNegClone"); //need clone, in order not to change track parameters
    AliESDtrack* lPosClone = (AliESDtrack*) lPosTrack->Clone("lPosClone"); //need clone, in order not to change track parameters
    
    //Get Magnetic field and primary vertex
    Double_t b=lEvent->GetMagneticField();
    const AliESDVertex *vtxT3D=lEvent->GetPrimaryVertex();
    Double_t xPrimaryVertex=vtxT3D->GetX();
    Double_t yPrimaryVertex=vtxT3D->GetY();
    Double_t zPrimaryVertex=vtxT3D->GetZ();
    
    //Get ExternalTrackParam
    AliExternalTrackParam nt(*lNegClone), pt(*lPosClone);
    
    //Find DCA
    Double_t xn, xp, dca=lNegClone->GetDCA(lPosClone,b,xn,xp);
    
    //Propagate to it
    nt.PropagateTo(xn,b); pt.PropagateTo(xp,b);
    
    //Create V0 object to do propagation
    AliESDv0 vertex(nt,1,pt,2); //Never mind indices, won't use
    
    //Get CosPA
    lCosPA=vertex.GetV0CosineOfPointingAngle(xPrimaryVertex,yPrimaryVertex,zPrimaryVertex);
    
    //Cleanup
    delete lNegClone;
    delete lPosClone;
    
    //Return value
    return lCosPA;
}


