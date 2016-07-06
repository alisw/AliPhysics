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
// Modified version of AliAnalysisTaskCheckCascade.cxx.
// This is a 'hybrid' output version, in that it uses a classic TTree
// ROOT object to store the candidates, plus a couple of histograms filled on
// a per-event basis for storing variables too numerous to put in a tree.
//
//
//  --- Algorithm Description
//   1. Loop over primaries in stack to acquire generated charged Xi
//   2. Loop over stack to find Cascades, fill TH3Fs "PrimRawPt"s for Efficiency
//   3. Perform Physics Selection
//   4. Perform Primary Vertex |z|<10cm selection
//   5. Perform Primary Vertex NoTPCOnly vertexing selection
//   6. Perform Pileup Rejection
//   7. Analysis Loops:
//    7a. Fill TH3Fs "PrimAnalysisPt" for control purposes only
//
//  Please Report Any Bugs!
//
//   --- David Dobrigkeit Chinellato
//        (david.chinellato@gmail.com)
//
// +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

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
#include "AliLog.h"

#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliLightV0vertexer.h"
#include "AliLightCascadeVertexer.h"

#include "AliCFContainer.h"
#include "AliMultiplicity.h"

#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliAODMCParticle.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDUtils.h"
#include "AliGenEventHeader.h"
#include "AliAnalysisUtils.h"

#include "AliAnalysisTaskExtractCascadePbPbRun2.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskExtractCascadePbPbRun2)

AliAnalysisTaskExtractCascadePbPbRun2::AliAnalysisTaskExtractCascadePbPbRun2()
: AliAnalysisTaskSE(), fListHist(0), fTreeCascade(0), fPIDResponse(0), fESDtrackCuts(0), fUtils(0),
fkRunVertexers ( kFALSE ),
fkSaveTree ( kTRUE ), 
fkSaveRawdEdxSignals ( kFALSE ),
fkSwitchCharges( kFALSE ),
fkSelectCentrality (kFALSE),
fCentSel_Low(0.0),
fCentSel_High(0.0),
fLowPtCutoff(0.0),
fCascadeMassWindow(0.075),
//------------------------------------------------
// Tree Variables
//------------------------------------------------

fTreeCascVarCharge(0),
fTreeCascVarMassAsXi(0),
fTreeCascVarMassAsOmega(0),
fTreeCascVarPt(0),
fTreeCascVarPtMC(0),
fTreeCascVarRapMC(0),
fTreeCascVarRapXi(0),
fTreeCascVarRapOmega(0),
fTreeCascVarNegEta(0),
fTreeCascVarPosEta(0),
fTreeCascVarBachEta(0),
fTreeCascVarDCACascDaughters(0),
fTreeCascVarDCABachToPrimVtx(0),
fTreeCascVarDCAV0Daughters(0),
fTreeCascVarDCAV0ToPrimVtx(0),
fTreeCascVarDCAPosToPrimVtx(0),
fTreeCascVarDCANegToPrimVtx(0),
fTreeCascVarCascCosPointingAngle(0),
fTreeCascVarCascRadius(0),
fTreeCascVarV0Mass(0),
fTreeCascVarV0CosPointingAngle(0),
fTreeCascVarV0CosPointingAngleSpecial(0),
fTreeCascVarV0Radius(0),
fTreeCascVarLeastNbrClusters(0),
fTreeCascVarCentrality(0), 
fTreeCascVarDistOverTotMom(0),
fTreeCascVarNegNSigmaPion(0),
fTreeCascVarNegNSigmaProton(0),
fTreeCascVarPosNSigmaPion(0),
fTreeCascVarPosNSigmaProton(0),
fTreeCascVarBachNSigmaPion(0),
fTreeCascVarBachNSigmaKaon(0),
fTreeCascVarRunNumber(0),
fTreeCascVarNegInDistortedRegion(0),
fTreeCascVarPosInDistortedRegion(0),
fTreeCascVarBachInDistortedRegion(0),
fTreeCascVarNegInnerP(0),
fTreeCascVarPosInnerP(0),
fTreeCascVarBachInnerP(0),
fTreeCascVarNegdEdx(0),
fTreeCascVarPosdEdx(0),
fTreeCascVarBachdEdx(0),

//------------------------------------------------
// HISTOGRAMS
// --- Filled on an Event-by-event basis
//------------------------------------------------
fHistEventCounter(0), 
fHistCentrality(0),
fHistdEdx(0), 
fHistdEdxPionsFromLambda(0), 
fHistdEdxProtonsFromLambda(0), 
fHistdEdxPionsFromK0s(0)
{
    // Dummy Constructor
}

AliAnalysisTaskExtractCascadePbPbRun2::AliAnalysisTaskExtractCascadePbPbRun2(const char *name, TString lExtraOptions)
: AliAnalysisTaskSE(name), fListHist(0), fTreeCascade(0), fPIDResponse(0), fESDtrackCuts(0), fUtils(0),
fkRunVertexers ( kFALSE ),
fkSaveTree ( kTRUE ), 
fkSaveRawdEdxSignals ( kFALSE ),
fkSwitchCharges( kFALSE ),
fkSelectCentrality (kFALSE),
fCentSel_Low(0.0),
fCentSel_High(0.0),
fLowPtCutoff(0.0),
fCascadeMassWindow(0.075),
//------------------------------------------------
// Tree Variables
//------------------------------------------------

fTreeCascVarCharge(0),
fTreeCascVarMassAsXi(0),
fTreeCascVarMassAsOmega(0),
fTreeCascVarPt(0),
fTreeCascVarPtMC(0),
fTreeCascVarRapMC(0),
fTreeCascVarRapXi(0),
fTreeCascVarRapOmega(0),
fTreeCascVarNegEta(0),
fTreeCascVarPosEta(0),
fTreeCascVarBachEta(0),
fTreeCascVarDCACascDaughters(0),
fTreeCascVarDCABachToPrimVtx(0),
fTreeCascVarDCAV0Daughters(0),
fTreeCascVarDCAV0ToPrimVtx(0),
fTreeCascVarDCAPosToPrimVtx(0),
fTreeCascVarDCANegToPrimVtx(0),
fTreeCascVarCascCosPointingAngle(0),
fTreeCascVarCascRadius(0),
fTreeCascVarV0Mass(0),
fTreeCascVarV0CosPointingAngle(0),
fTreeCascVarV0CosPointingAngleSpecial(0),
fTreeCascVarV0Radius(0),
fTreeCascVarLeastNbrClusters(0),
fTreeCascVarCentrality(0),
fTreeCascVarDistOverTotMom(0),
fTreeCascVarNegNSigmaPion(0),
fTreeCascVarNegNSigmaProton(0),
fTreeCascVarPosNSigmaPion(0),
fTreeCascVarPosNSigmaProton(0),
fTreeCascVarBachNSigmaPion(0),
fTreeCascVarBachNSigmaKaon(0),
fTreeCascVarRunNumber(0),
fTreeCascVarNegInDistortedRegion(0),
fTreeCascVarPosInDistortedRegion(0),
fTreeCascVarBachInDistortedRegion(0),
fTreeCascVarNegInnerP(0),
fTreeCascVarPosInnerP(0),
fTreeCascVarBachInnerP(0),
fTreeCascVarNegdEdx(0),
fTreeCascVarPosdEdx(0),
fTreeCascVarBachdEdx(0),

//------------------------------------------------
// HISTOGRAMS
// --- Filled on an Event-by-event basis
//------------------------------------------------
fHistEventCounter(0), 
fHistCentrality(0),
fHistdEdx(0), 
fHistdEdxPionsFromLambda(0), 
fHistdEdxProtonsFromLambda(0), 
fHistdEdxPionsFromK0s(0)

{
    // Constructor
    
    //Set Variables for re-running the cascade vertexers (as done for MS paper)
    if ( lExtraOptions.Contains("dEdx") )   fkSaveRawdEdxSignals = kTRUE;  
    if ( lExtraOptions.Contains("NoTree") ) fkSaveTree = kFALSE;  
  
    // New Loose : 1st step for the 7 TeV pp analysis
    
    fV0VertexerSels[0] =  33.  ;  // max allowed chi2
    fV0VertexerSels[1] =   0.02;  // min allowed impact parameter for the 1st daughter (LHC09a4 : 0.05)
    fV0VertexerSels[2] =   0.02;  // min allowed impact parameter for the 2nd daughter (LHC09a4 : 0.05)
    fV0VertexerSels[3] =   2.0 ;  // max allowed DCA between the daughter tracks       (LHC09a4 : 0.5)
    fV0VertexerSels[4] =   0.95;  // min allowed cosine of V0's pointing angle         (LHC09a4 : 0.99)
    fV0VertexerSels[5] =   1.0 ;  // min radius of the fiducial volume                 (LHC09a4 : 0.2)
    fV0VertexerSels[6] = 200.  ;  // max radius of the fiducial volume                 (LHC09a4 : 100.0)
    
    fCascadeVertexerSels[0] =  33.   ;  // max allowed chi2 (same as PDC07)
    fCascadeVertexerSels[1] =   0.05 ;  // min allowed V0 impact parameter                    (PDC07 : 0.05   / LHC09a4 : 0.025 )
    fCascadeVertexerSels[2] =   0.010;  // "window" around the Lambda mass                    (PDC07 : 0.008  / LHC09a4 : 0.010 )
    fCascadeVertexerSels[3] =   0.03 ;  // min allowed bachelor's impact parameter            (PDC07 : 0.035  / LHC09a4 : 0.025 )
    fCascadeVertexerSels[4] =   2.0  ;  // max allowed DCA between the V0 and the bachelor    (PDC07 : 0.1    / LHC09a4 : 0.2   )
    fCascadeVertexerSels[5] =   0.95 ;  // min allowed cosine of the cascade pointing angle   (PDC07 : 0.9985 / LHC09a4 : 0.998 )
    fCascadeVertexerSels[6] =   0.4  ;  // min radius of the fiducial volume                  (PDC07 : 0.9    / LHC09a4 : 0.2   )
    fCascadeVertexerSels[7] = 100.   ;  // max radius of the fiducial volume                  (PDC07 : 100    / LHC09a4 : 100   )
    
    // Output slot #0 writes into a TList container (Cascade)
    DefineOutput(1, TList::Class());
    DefineOutput(2, TTree::Class());
}


AliAnalysisTaskExtractCascadePbPbRun2::~AliAnalysisTaskExtractCascadePbPbRun2()
{
    //------------------------------------------------
    // DESTRUCTOR
    //------------------------------------------------
    
    if (fListHist){
        delete fListHist;
        fListHist = 0x0;
    }
    if (fTreeCascade){
        delete fTreeCascade;
        fTreeCascade = 0x0;
    }
    //cleanup esd track cuts object too...
    if (fESDtrackCuts){
        delete fESDtrackCuts;
        fESDtrackCuts = 0x0;
    }
    if (fUtils){
        delete fUtils;
        fUtils = 0x0;
    }
    
}

//________________________________________________________________________
void AliAnalysisTaskExtractCascadePbPbRun2::UserCreateOutputObjects()
{
    OpenFile(2);
    // Called once
    
    //------------------------------------------------
    
    if ( fkSaveTree ){ 
    fTreeCascade = new TTree("fTreeCascade","CascadeCandidates");
    
    //------------------------------------------------
    // fTreeCascade Branch definitions - Cascade Tree
    //------------------------------------------------
    
    //------------------------------------------------
    // fTreeCascade Branch definitions
    //------------------------------------------------
    
    //-----------BASIC-INFO---------------------------
    fTreeCascade->Branch("fTreeCascVarCharge",&fTreeCascVarCharge,"fTreeCascVarCharge/I");
    fTreeCascade->Branch("fTreeCascVarMassAsXi",&fTreeCascVarMassAsXi,"fTreeCascVarMassAsXi/F");
    fTreeCascade->Branch("fTreeCascVarMassAsOmega",&fTreeCascVarMassAsOmega,"fTreeCascVarMassAsOmega/F");
    fTreeCascade->Branch("fTreeCascVarPt",&fTreeCascVarPt,"fTreeCascVarPt/F");
    fTreeCascade->Branch("fTreeCascVarRapXi",&fTreeCascVarRapXi,"fTreeCascVarRapXi/F");
    fTreeCascade->Branch("fTreeCascVarRapOmega",&fTreeCascVarRapOmega,"fTreeCascVarRapOmega/F");
    fTreeCascade->Branch("fTreeCascVarNegEta",&fTreeCascVarNegEta,"fTreeCascVarNegEta/F");
    fTreeCascade->Branch("fTreeCascVarPosEta",&fTreeCascVarPosEta,"fTreeCascVarPosEta/F");
    fTreeCascade->Branch("fTreeCascVarBachEta",&fTreeCascVarBachEta,"fTreeCascVarBachEta/F");
    //-----------INFO-FOR-CUTS------------------------
    fTreeCascade->Branch("fTreeCascVarDCACascDaughters",&fTreeCascVarDCACascDaughters,"fTreeCascVarDCACascDaughters/F");
    fTreeCascade->Branch("fTreeCascVarDCABachToPrimVtx",&fTreeCascVarDCABachToPrimVtx,"fTreeCascVarDCABachToPrimVtx/F");
    fTreeCascade->Branch("fTreeCascVarDCAV0Daughters",&fTreeCascVarDCAV0Daughters,"fTreeCascVarDCAV0Daughters/F");
    fTreeCascade->Branch("fTreeCascVarDCAV0ToPrimVtx",&fTreeCascVarDCAV0ToPrimVtx,"fTreeCascVarDCAV0ToPrimVtx/F");
    fTreeCascade->Branch("fTreeCascVarDCAPosToPrimVtx",&fTreeCascVarDCAPosToPrimVtx,"fTreeCascVarDCAPosToPrimVtx/F");
    fTreeCascade->Branch("fTreeCascVarDCANegToPrimVtx",&fTreeCascVarDCANegToPrimVtx,"fTreeCascVarDCANegToPrimVtx/F");
    fTreeCascade->Branch("fTreeCascVarCascCosPointingAngle",&fTreeCascVarCascCosPointingAngle,"fTreeCascVarCascCosPointingAngle/F");
    fTreeCascade->Branch("fTreeCascVarCascRadius",&fTreeCascVarCascRadius,"fTreeCascVarCascRadius/F");
    fTreeCascade->Branch("fTreeCascVarV0Mass",&fTreeCascVarV0Mass,"fTreeCascVarV0Mass/F");
    fTreeCascade->Branch("fTreeCascVarV0CosPointingAngle",&fTreeCascVarV0CosPointingAngle,"fTreeCascVarV0CosPointingAngle/F");
    fTreeCascade->Branch("fTreeCascVarV0CosPointingAngleSpecial",&fTreeCascVarV0CosPointingAngleSpecial,"fTreeCascVarV0CosPointingAngleSpecial/F");
    fTreeCascade->Branch("fTreeCascVarV0Radius",&fTreeCascVarV0Radius,"fTreeCascVarV0Radius/F");
    fTreeCascade->Branch("fTreeCascVarLeastNbrClusters",&fTreeCascVarLeastNbrClusters,"fTreeCascVarLeastNbrClusters/I");
    //-----------MULTIPLICITY-INFO--------------------
    fTreeCascade->Branch("fTreeCascVarCentrality",&fTreeCascVarCentrality,"fTreeCascVarCentrality/F");
    //-----------DECAY-LENGTH-INFO--------------------
    fTreeCascade->Branch("fTreeCascVarDistOverTotMom",&fTreeCascVarDistOverTotMom,"fTreeCascVarDistOverTotMom/F");
    //------------------------------------------------
    fTreeCascade->Branch("fTreeCascVarNegNSigmaPion",&fTreeCascVarNegNSigmaPion,"fTreeCascVarNegNSigmaPion/F");
    fTreeCascade->Branch("fTreeCascVarNegNSigmaProton",&fTreeCascVarNegNSigmaProton,"fTreeCascVarNegNSigmaProton/F");
    fTreeCascade->Branch("fTreeCascVarPosNSigmaPion",&fTreeCascVarPosNSigmaPion,"fTreeCascVarPosNSigmaPion/F");
    fTreeCascade->Branch("fTreeCascVarPosNSigmaProton",&fTreeCascVarPosNSigmaProton,"fTreeCascVarPosNSigmaProton/F");
    fTreeCascade->Branch("fTreeCascVarBachNSigmaPion",&fTreeCascVarBachNSigmaPion,"fTreeCascVarBachNSigmaPion/F");
    fTreeCascade->Branch("fTreeCascVarBachNSigmaKaon",&fTreeCascVarBachNSigmaKaon,"fTreeCascVarBachNSigmaKaon/F");
    //-----------RUN-NUMBER---------------------------
    fTreeCascade->Branch("fTreeCascVarRunNumber",&fTreeCascVarRunNumber,"fTreeCascVarRunNumber/I");   
    //-----------DISTORTED-TPC-REGIONS----------------
    fTreeCascade->Branch("fTreeCascVarNegInDistortedRegion",&fTreeCascVarNegInDistortedRegion,"fTreeCascVarNegInDistortedRegion/O");   
    fTreeCascade->Branch("fTreeCascVarPosInDistortedRegion",&fTreeCascVarPosInDistortedRegion,"fTreeCascVarPosInDistortedRegion/O");   
    fTreeCascade->Branch("fTreeCascVarBachInDistortedRegion",&fTreeCascVarBachInDistortedRegion,"fTreeCascVarBachInDistortedRegion/O");   
    
    //-----------TPC-DEDX-INFO------------------------
    if ( fkSaveRawdEdxSignals ){ 
      fTreeCascade->Branch("fTreeCascVarPosInnerP",&fTreeCascVarPosInnerP,"fTreeCascVarPosInnerP/F");  
      fTreeCascade->Branch("fTreeCascVarNegInnerP",&fTreeCascVarNegInnerP,"fTreeCascVarNegInnerP/F");  
      fTreeCascade->Branch("fTreeCascVarBachInnerP",&fTreeCascVarBachInnerP,"fTreeCascVarBachInnerP/F");  
      fTreeCascade->Branch("fTreeCascVarPosdEdx",&fTreeCascVarPosdEdx,"fTreeCascVarPosdEdx/F");  
      fTreeCascade->Branch("fTreeCascVarNegdEdx",&fTreeCascVarNegdEdx,"fTreeCascVarNegdEdx/F");  
      fTreeCascade->Branch("fTreeCascVarBachdEdx",&fTreeCascVarBachdEdx,"fTreeCascVarBachdEdx/F");  
    }
    }
    
    //------------------------------------------------
    // Particle Identification Setup
    //------------------------------------------------
    
    AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    
    // Multiplicity
    
    if(! fESDtrackCuts ){
        fESDtrackCuts = new AliESDtrackCuts();
	
	//Set of silly track selections 
	fESDtrackCuts->SetAcceptKinkDaughters(kFALSE);
	fESDtrackCuts->SetMinNClustersTPC(80);
	fESDtrackCuts->SetMaxChi2PerClusterTPC(4);
	fESDtrackCuts->SetMaxDCAToVertexXY(.5);
	fESDtrackCuts->SetMaxDCAToVertexZ(.5);
	fESDtrackCuts->SetRequireTPCRefit(kTRUE);
    }
    if(! fUtils ){
        fUtils = new AliAnalysisUtils();
    }
    
    //------------------------------------------------
    // V0 Multiplicity Histograms
    //------------------------------------------------
    
    // Create histograms
    OpenFile(1);
    fListHist = new TList();
    fListHist->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner
    
    if(! fHistCentrality) {
        fHistCentrality = new TH1F("fHistCentrality", "Centrality Distribution;Centrality;Events", 250, 0, 250);
        fListHist->Add(fHistCentrality);
    }

    if(! fHistEventCounter) {
        fHistEventCounter = new TH1F("fHistEventCounter", "Event Counters;Selection Step;Events", 2, 0, 2);
	fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
	fHistEventCounter->GetXaxis()->SetBinLabel(2, "Selected");
        fListHist->Add(fHistEventCounter);
    }

    if (fkSaveRawdEdxSignals) { 
      if(! fHistdEdx) {
	  fHistdEdx = new TH2F("fHistdEdx", ";Momentum at Inner Wall;TPC Signal", 400,0,4,500,0,500);
	  fListHist->Add(fHistdEdx);
      }      
       if(! fHistdEdxPionsFromLambda) {
	  fHistdEdxPionsFromLambda = new TH2F("fHistdEdxPionsFromLambda", ";Momentum at Inner Wall;TPC Signal", 400,0,4,500,0,500);
	  fListHist->Add(fHistdEdxPionsFromLambda);
      }           
       if(! fHistdEdxProtonsFromLambda) {
	  fHistdEdxProtonsFromLambda = new TH2F("fHistdEdxProtonsFromLambda", ";Momentum at Inner Wall;TPC Signal", 400,0,4,500,0,500);
	  fListHist->Add(fHistdEdxProtonsFromLambda);
      }           
      if(! fHistdEdxPionsFromK0s) {
	  fHistdEdxPionsFromK0s = new TH2F("fHistdEdxPionsFromK0s", ";Momentum at Inner Wall;TPC Signal", 400,0,4,500,0,500);
	  fListHist->Add(fHistdEdxPionsFromK0s);
      } 
    }
    
    
    //List of Histograms: Normal
    PostData(1, fListHist);
    
    //TTree Object: Saved to base directory. Should cache to disk while saving.
    //(Important to avoid excessive memory usage, particularly when merging)
    PostData(2, fTreeCascade);
    
}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskExtractCascadePbPbRun2::UserExec(Option_t *)
{
    // Main loop
    // Called for each event
    
    AliESDEvent *lESDevent = 0x0;
    
    Int_t    lNumberOfV0s                      = -1;
    Double_t lTrkgPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
    Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
    Double_t lMagneticField                 = -10.;
    
    // Connect to the InputEvent
    // After these lines, we should have an ESD/AOD event + the number of V0s in it.
    
    // Appropriate for ESD analysis!
    
    lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
    if (!lESDevent) {
        AliWarning("ERROR: lESDevent not available \n");
        return;
    }

    //--- Acquisition of exact event ID
    fTreeCascVarRunNumber = lESDevent->GetRunNumber();

    //------------------------------------------------
    // Multiplicity Information Acquistion
    // New Multiplicity Selection Framework (Run 2)
    //------------------------------------------------

    Float_t lPercentile = 300; //not acquired
    Int_t lEvSelCode = 300;
    AliMultSelection *MultSelection = 0x0;
    MultSelection = (AliMultSelection * ) lESDevent->FindListObject("MultSelection");
    MultSelection = (AliMultSelection * ) lESDevent->FindListObject("MultSelection");
    if( !MultSelection) {
        //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
        AliWarning("AliMultSelection object not found!");
    } else {
        //V0M Multiplicity Percentile
        lPercentile = MultSelection->GetMultiplicityPercentile("V0M");
        //Event Selection Code
        lEvSelCode = MultSelection->GetEvSelCode();
    }

    fHistEventCounter -> Fill (0.5);

    //------------------------------------------------
    // Reject events which don't pass standard
    // AliMultSelection event selection criteria
    //------------------------------------------------
    if (lEvSelCode>1) {
        PostData(1, fListHist);
        PostData(2, fTreeCascade);
        return;
    }
    
    if( fkSelectCentrality ){
        if( lPercentile < fCentSel_Low || lPercentile >= fCentSel_High ){
            //Event is outside desired centrality centrality in V0M!
            PostData(1, fListHist);
            PostData(2, fTreeCascade);
            return;
        }
    }
    
    //Set variable for filling tree afterwards!
    fTreeCascVarCentrality = lPercentile;
    fHistCentrality->Fill ( lPercentile );
    
    //------------------------------------------------
    // Physics Selection
    //------------------------------------------------

    /* OBSOLETED, included in event selections above! 
    UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    Bool_t isSelected = 0;
    isSelected = (maskIsSelected & AliVEvent::kMB) == AliVEvent::kMB;
    
    //pA triggering: CINT7
    if( fkSwitchINT7 ) isSelected = (maskIsSelected & AliVEvent::kINT7) == AliVEvent::kINT7;
    
    //Standard Min-Bias Selection
    if ( ! isSelected ) {
        PostData(1, fListHist);
        PostData(2, fTreeCascade);
        return;
    }
    */ 

    //------------------------------------------------
    // Rerun cascade vertexer in its light version
    //------------------------------------------------
    
    if( fkRunVertexers && fkSaveTree ){
        lESDevent->ResetCascades();
        lESDevent->ResetV0s();
        
        AliLightV0vertexer lV0vtxer;
        AliLightCascadeVertexer lCascVtxer;
        
        //Do wrong charge combination (experimental!)
        lCascVtxer.SetSwitchCharges( fkSwitchCharges );
        
        lV0vtxer.SetDefaultCuts(fV0VertexerSels);
        lCascVtxer.SetDefaultCuts(fCascadeVertexerSels);
        
        lV0vtxer.Tracks2V0vertices(lESDevent);
        lCascVtxer.V0sTracks2CascadeVertices(lESDevent);
    }

    
    //------------------------------------------------
    // Getting: Primary Vertex + MagField Info
    //------------------------------------------------
    
    const AliESDVertex *lPrimaryTrackingESDVtx = lESDevent->GetPrimaryVertexTracks();
    // get the vtx stored in ESD found with tracks
    lPrimaryTrackingESDVtx->GetXYZ( lTrkgPrimaryVtxPos );
    
    const AliESDVertex *lPrimaryBestESDVtx = lESDevent->GetPrimaryVertex();
    // get the best primary vertex available for the event
    // As done in AliCascadeVertexer, we keep the one which is the best one available.
    // between : Tracking vertex > SPD vertex > TPC vertex > default SPD vertex
    // This one will be used for next calculations (DCA essentially)
    lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );
    
    Double_t lPrimaryVtxPosition[3];
    const AliVVertex *primaryVtx = lESDevent->GetPrimaryVertex();
    lPrimaryVtxPosition[0] = primaryVtx->GetX();
    lPrimaryVtxPosition[1] = primaryVtx->GetY();
    lPrimaryVtxPosition[2] = primaryVtx->GetZ();

    //------------------------------------------------
    // Getting information for playing around: 
    //  Loop over all tracks 
    //------------------------------------------------
    
    if (fkSaveRawdEdxSignals) { 
      //All-Track Inclusive (just in case) 
      for (Int_t i=0;i<lESDevent->GetNumberOfTracks();++i) {
	AliESDtrack *trackOne = 0x0;
	trackOne = lESDevent->GetTrack(i);
	if ( !fESDtrackCuts->AcceptTrack(trackOne) ) continue; 
	const AliExternalTrackParam *trackOneParams=trackOne->GetInnerParam();
	Float_t lThisInnerP = -1; 
	if(trackOneParams)  { lThisInnerP  = trackOneParams ->GetP(); }
	Float_t lThisTPCSignal = trackOne->GetTPCsignal();
	fHistdEdx -> Fill( lThisInnerP , lThisTPCSignal ); 	
      }
      
      //Clean Sample of V0s 
      Int_t nv0s = lESDevent->GetNumberOfV0s();
      for (Int_t iV0 = 0; iV0 < nv0s; iV0++) //extra-crazy test
      {   // This is the begining of the V0 loop
        AliESDv0 *v0 = ((AliESDEvent*)lESDevent)->GetV0(iV0);
        if (!v0) continue;
	
	//Get Daughters
	
        UInt_t lKeyPos = (UInt_t)TMath::Abs(v0->GetPindex());
        UInt_t lKeyNeg = (UInt_t)TMath::Abs(v0->GetNindex());

        AliESDtrack *pTrack=((AliESDEvent*)lESDevent)->GetTrack(lKeyPos);
        AliESDtrack *nTrack=((AliESDEvent*)lESDevent)->GetTrack(lKeyNeg);
        if (!pTrack || !nTrack) {
            Printf("ERROR: Could not retreive one of the daughter track");
            continue;
        }
	//Daughter Eta for Eta selection, afterwards
        if ( TMath::Abs( nTrack->Eta() ) > 0.9 || TMath::Abs( pTrack->Eta() ) > 0.9 ) continue;

        // Filter like-sign V0 (next: add counter and distribution)
        if ( pTrack->GetSign() == nTrack->GetSign()) {
            continue;
        }
	//________________________________________________________________________
        // Track quality cuts
        Float_t lPosTrackCrossedRows = pTrack->GetTPCClusterInfo(2,1);
        Float_t lNegTrackCrossedRows = nTrack->GetTPCClusterInfo(2,1);
	if ( lPosTrackCrossedRows < 70 || lNegTrackCrossedRows < 70 ) continue; 

        // TPC refit condition (done during reconstruction for Offline but not for On-the-fly)
        if( !(pTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
        if( !(nTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
	
	//GetKinkIndex condition
        if( pTrack->GetKinkIndex(0)>0 || nTrack->GetKinkIndex(0)>0 ) continue;

        //Findable clusters > 0 condition
        if( pTrack->GetTPCNclsF()<=0 || nTrack->GetTPCNclsF()<=0 ) continue;
	
	Float_t lV0CosineOfPointingAngle = v0->GetV0CosineOfPointingAngle(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1],lPrimaryVtxPosition[2]);
	
	if ( lV0CosineOfPointingAngle < 0.999 ) continue; //who cares? not me!
		
	Int_t lOnFlyStatus = v0->GetOnFlyStatus();
	if ( lOnFlyStatus != 0 ) continue; 
	
        Float_t lDcaV0Daughters = v0->GetDcaV0Daughters();
	if ( lDcaV0Daughters > 0.3 ) continue; //still not me!
		
        // Getting invariant mass infos directly from ESD
        v0->ChangeMassHypothesis(310);
        Float_t lInvMassK0s = v0->GetEffMass();
        v0->ChangeMassHypothesis(3122);
        Float_t lInvMassLambda = v0->GetEffMass();
        v0->ChangeMassHypothesis(-3122);
        Float_t lInvMassAntiLambda = v0->GetEffMass();
	
	const AliExternalTrackParam *innernegv0=nTrack->GetInnerParam();
	const AliExternalTrackParam *innerposv0=pTrack->GetInnerParam();
	Float_t lThisPosInnerP = -1; 
	Float_t lThisNegInnerP = -1; 
	if(innerposv0)  { lThisPosInnerP  = innerposv0 ->GetP(); }
	if(innernegv0)  { lThisNegInnerP  = innernegv0 ->GetP(); }
	Float_t lThisPosdEdx = pTrack -> GetTPCsignal(); 
	Float_t lThisNegdEdx = nTrack -> GetTPCsignal(); 
	
	if ( TMath::Abs( lInvMassLambda - 1.1157 ) < 0.0025 && TMath::Abs( lInvMassK0s - 0.497 ) > 0.005 && TMath::Abs( lInvMassAntiLambda - 1.1157 ) > 0.003 ){ 
	  fHistdEdxPionsFromLambda -> Fill ( lThisNegInnerP, lThisNegdEdx );
	  fHistdEdxProtonsFromLambda -> Fill ( lThisPosInnerP, lThisPosdEdx );
	}
	if ( TMath::Abs( lInvMassAntiLambda - 1.1157 ) < 0.0025 && TMath::Abs( lInvMassK0s - 0.497 ) > 0.005 && TMath::Abs( lInvMassLambda - 1.1157 ) > 0.003 ){ 
	  fHistdEdxPionsFromLambda -> Fill ( lThisPosInnerP, lThisPosdEdx );
	  fHistdEdxProtonsFromLambda -> Fill ( lThisNegInnerP, lThisNegdEdx );
	}	
	if ( TMath::Abs( lInvMassK0s - 0.497 ) < 0.0025 ){ 
	  fHistdEdxPionsFromK0s -> Fill( lThisPosInnerP , lThisPosdEdx ); 
	  fHistdEdxPionsFromK0s -> Fill( lThisNegInnerP , lThisNegdEdx ); 
	}
      }
    } 
    
    //------------------------------------------------
    // MAIN CASCADE LOOP STARTS HERE
    //------------------------------------------------
    // Code Credit: Antonin Maire (thanks^100)
    // ---> This is an adaptation

    Long_t ncascades = 0;
    ncascades = lESDevent->GetNumberOfCascades();

    for (Int_t iXi = 0; iXi < ncascades; iXi++){
        //------------------------------------------------
        // Initializations
        //------------------------------------------------
        //Double_t lTrkgPrimaryVtxRadius3D = -500.0;
        //Double_t lBestPrimaryVtxRadius3D = -500.0;
                
        // - 1st part of initialisation : variables needed to store AliESDCascade data members
        Double_t lEffMassXi      = 0. ;
        //Double_t lChi2Xi         = -1. ;
        Double_t lDcaXiDaughters = -1. ;
        Double_t lXiCosineOfPointingAngle = -1. ;
        Double_t lPosXi[3] = { -1000.0, -1000.0, -1000.0 };
        Double_t lXiRadius = -1000. ;
        
        // - 2nd part of initialisation : Nbr of clusters within TPC for the 3 daughter cascade tracks
        Int_t    lPosTPCClusters    = -1; // For ESD only ...//FIXME : wait for availability in AOD
        Int_t    lNegTPCClusters    = -1; // For ESD only ...
        Int_t    lBachTPCClusters   = -1; // For ESD only ...
        
        // - 3rd part of initialisation : about V0 part in cascades
        Double_t lInvMassLambdaAsCascDghter = 0.;
        //Double_t lV0Chi2Xi         = -1. ;
        Double_t lDcaV0DaughtersXi = -1.;
		
        Double_t lDcaBachToPrimVertexXi = -1., lDcaV0ToPrimVertexXi = -1.;
        Double_t lDcaPosToPrimVertexXi  = -1.;
        Double_t lDcaNegToPrimVertexXi  = -1.;
        Double_t lV0CosineOfPointingAngleXi = -1. ;
        Double_t lV0CosineOfPointingAngleXiSpecial = -1. ;
        Double_t lPosV0Xi[3] = { -1000. , -1000., -1000. }; // Position of VO coming from cascade
        Double_t lV0RadiusXi = -1000.0;
        Double_t lV0quality  = 0.;
        
        // - 4th part of initialisation : Effective masses
        Double_t lInvMassXiMinus    = 0.;
        Double_t lInvMassXiPlus     = 0.;
        Double_t lInvMassOmegaMinus = 0.;
        Double_t lInvMassOmegaPlus  = 0.;
        
        // - 6th part of initialisation : extra info for QA
        Double_t lXiMomX       = 0. , lXiMomY = 0., lXiMomZ = 0.;
        Double_t lXiTransvMom  = 0. ;
        Double_t lXiTransvMomMC= 0. ;
        Double_t lXiTotMom     = 0. ;
		
        Double_t lBachMomX       = 0., lBachMomY  = 0., lBachMomZ   = 0.;
        //Double_t lBachTransvMom  = 0.;
        //Double_t lBachTotMom     = 0.;
        
        fTreeCascVarNegNSigmaPion   = -100;
        fTreeCascVarNegNSigmaProton = -100;
        fTreeCascVarPosNSigmaPion   = -100;
        fTreeCascVarPosNSigmaProton = -100;
        fTreeCascVarBachNSigmaPion  = -100;
        fTreeCascVarBachNSigmaKaon  = -100;
        
        Short_t  lChargeXi = -2;
        //Double_t lV0toXiCosineOfPointingAngle = 0. ;
        
        Double_t lRapXi   = -20.0, lRapOmega = -20.0, lRapMC = -20.0; //  lEta = -20.0, lTheta = 360., lPhi = 720. ;
        //Double_t lAlphaXi = -200., lPtArmXi  = -200.0;
	    
        // -------------------------------------
        // II.ESD - Calculation Part dedicated to Xi vertices (ESD)
        
        AliESDcascade *xi = lESDevent->GetCascade(iXi);
        if (!xi) continue;
        
        
        // - II.Step 1 : around primary vertex
        //-------------
        //lTrkgPrimaryVtxRadius3D = TMath::Sqrt(  lTrkgPrimaryVtxPos[0] * lTrkgPrimaryVtxPos[0] +
        //                                        lTrkgPrimaryVtxPos[1] * lTrkgPrimaryVtxPos[1] +
        //                                        lTrkgPrimaryVtxPos[2] * lTrkgPrimaryVtxPos[2] );
        
        //lBestPrimaryVtxRadius3D = TMath::Sqrt(  lBestPrimaryVtxPos[0] * lBestPrimaryVtxPos[0] +
        //                                        lBestPrimaryVtxPos[1] * lBestPrimaryVtxPos[1] +
        //                                        lBestPrimaryVtxPos[2] * lBestPrimaryVtxPos[2] );
        
		// - II.Step 2 : Assigning the necessary variables for specific AliESDcascade data members (ESD)
		//-------------
        lV0quality = 0.;
        xi->ChangeMassHypothesis(lV0quality , 3312); // default working hypothesis : cascade = Xi- decay
        
        lEffMassXi  			= xi->GetEffMassXi();
        //lChi2Xi 			    = xi->GetChi2Xi();
        lDcaXiDaughters 	= xi->GetDcaXiDaughters();
        lXiCosineOfPointingAngle 	            = xi->GetCascadeCosineOfPointingAngle( lBestPrimaryVtxPos[0],
                                                                                      lBestPrimaryVtxPos[1],
                                                                                      lBestPrimaryVtxPos[2] );
        // Take care : the best available vertex should be used (like in AliCascadeVertexer)
        
        xi->GetXYZcascade( lPosXi[0],  lPosXi[1], lPosXi[2] );
        lXiRadius			= TMath::Sqrt( lPosXi[0]*lPosXi[0]  +  lPosXi[1]*lPosXi[1] );
        
		// - II.Step 3 : around the tracks : Bach + V0 (ESD)
		// ~ Necessary variables for ESDcascade data members coming from the ESDv0 part (inheritance)
		//-------------
		
        UInt_t lIdxPosXi 	= (UInt_t) TMath::Abs( xi->GetPindex() );
        UInt_t lIdxNegXi 	= (UInt_t) TMath::Abs( xi->GetNindex() );
        UInt_t lBachIdx 	= (UInt_t) TMath::Abs( xi->GetBindex() );
        // Care track label can be negative in MC production (linked with the track quality)
        // However = normally, not the case for track index ...
        
        // FIXME : rejection of a double use of a daughter track (nothing but just a crosscheck of what is done in the cascade vertexer)
        if(lBachIdx == lIdxNegXi) {
            AliWarning("Pb / Idx(Bach. track) = Idx(Neg. track) ... continue!"); continue;
        }
        if(lBachIdx == lIdxPosXi) {
            AliWarning("Pb / Idx(Bach. track) = Idx(Pos. track) ... continue!"); continue;
        }
        
        AliESDtrack *pTrackXi		= lESDevent->GetTrack( lIdxPosXi );
        AliESDtrack *nTrackXi		= lESDevent->GetTrack( lIdxNegXi );
        AliESDtrack *bachTrackXi	= lESDevent->GetTrack( lBachIdx );
        
        if (!pTrackXi || !nTrackXi || !bachTrackXi ) {
            AliWarning("ERROR: Could not retrieve one of the 3 ESD daughter tracks of the cascade ...");
            continue;
        }
        
        fTreeCascVarPosEta = pTrackXi->Eta();
        fTreeCascVarNegEta = nTrackXi->Eta();
        fTreeCascVarBachEta = bachTrackXi->Eta();
        
        fTreeCascVarPosInDistortedRegion = AliESDtrackCuts::IsTrackInDistortedTpcRegion( pTrackXi );
        fTreeCascVarNegInDistortedRegion = AliESDtrackCuts::IsTrackInDistortedTpcRegion( nTrackXi );
        fTreeCascVarBachInDistortedRegion = AliESDtrackCuts::IsTrackInDistortedTpcRegion( bachTrackXi );

	//Acquisition of TPC raw signal information 
	if ( fkSaveRawdEdxSignals ) {
	  //Step 1: Acquire TPC inner wall total momentum 
	  const AliExternalTrackParam *innerneg=nTrackXi->GetInnerParam();
	  const AliExternalTrackParam *innerpos=pTrackXi->GetInnerParam();
	  const AliExternalTrackParam *innerbach=bachTrackXi->GetInnerParam();
	  fTreeCascVarPosInnerP = -1; 
	  fTreeCascVarNegInnerP = -1; 
	  fTreeCascVarBachInnerP = -1; 
	  if(innerpos)  { fTreeCascVarPosInnerP  = innerpos ->GetP(); }
	  if(innerneg)  { fTreeCascVarNegInnerP  = innerneg ->GetP(); }
	  if(innerbach) { fTreeCascVarBachInnerP = innerbach->GetP(); }
	  
	  //Step 2: Acquire TPC Signals
	  fTreeCascVarPosdEdx = pTrackXi->GetTPCsignal();
	  fTreeCascVarNegdEdx = nTrackXi->GetTPCsignal();
	  fTreeCascVarBachdEdx = bachTrackXi->GetTPCsignal();
	}
	
	
        Double_t lBMom[3], lNMom[3], lPMom[3];
        xi->GetBPxPyPz( lBMom[0], lBMom[1], lBMom[2] );
        xi->GetPPxPyPz( lPMom[0], lPMom[1], lPMom[2] );
        xi->GetNPxPyPz( lNMom[0], lNMom[1], lNMom[2] );
        
        //------------------------------------------------
        // TPC dEdx information
        //------------------------------------------------
        fTreeCascVarNegNSigmaPion   = fPIDResponse->NumberOfSigmasTPC( nTrackXi, AliPID::kPion   );
        fTreeCascVarNegNSigmaProton = fPIDResponse->NumberOfSigmasTPC( nTrackXi, AliPID::kProton );
        fTreeCascVarPosNSigmaPion   = fPIDResponse->NumberOfSigmasTPC( pTrackXi, AliPID::kPion );
        fTreeCascVarPosNSigmaProton = fPIDResponse->NumberOfSigmasTPC( pTrackXi, AliPID::kProton );
        fTreeCascVarBachNSigmaPion  = fPIDResponse->NumberOfSigmasTPC( bachTrackXi, AliPID::kPion );
        fTreeCascVarBachNSigmaKaon  = fPIDResponse->NumberOfSigmasTPC( bachTrackXi, AliPID::kKaon );
        
        //------------------------------------------------
        // TPC Number of clusters info
        // --- modified to save the smallest number
        // --- of TPC clusters for the 3 tracks
        //------------------------------------------------
        
        lPosTPCClusters   = pTrackXi->GetTPCNcls();
        lNegTPCClusters   = nTrackXi->GetTPCNcls();
        lBachTPCClusters  = bachTrackXi->GetTPCNcls();
                
        // 1 - Poor quality related to TPCrefit
        ULong_t pStatus    = pTrackXi->GetStatus();
        ULong_t nStatus    = nTrackXi->GetStatus();
        ULong_t bachStatus = bachTrackXi->GetStatus();
                
        if ((pStatus&AliESDtrack::kTPCrefit)    == 0) { /*AliWarning("Pb / V0 Pos. track has no TPCrefit ... continue!"); */continue; }
        if ((nStatus&AliESDtrack::kTPCrefit)    == 0) { /*AliWarning("Pb / V0 Neg. track has no TPCrefit ... continue!"); */continue; }
        if ((bachStatus&AliESDtrack::kTPCrefit) == 0) { /*AliWarning("Pb / Bach.   track has no TPCrefit ... continue!"); */continue; }
                
        // 2 - Poor quality related to TPC clusters: lowest cut of 70 clusters
        if(lPosTPCClusters  < 70) { /*AliWarning("Pb / V0 Pos. track has less than 70 TPC clusters ... continue!"); */continue; }
        if(lNegTPCClusters  < 70) { /*AliWarning("Pb / V0 Neg. track has less than 70 TPC clusters ... continue!"); */continue; }
        if(lBachTPCClusters < 70) { /*AliWarning("Pb / Bach.   track has less than 70 TPC clusters ... continue!"); */continue; }
        Int_t leastnumberofclusters = 1000;
        if( lPosTPCClusters < leastnumberofclusters ) leastnumberofclusters = lPosTPCClusters;
        if( lNegTPCClusters < leastnumberofclusters ) leastnumberofclusters = lNegTPCClusters;
        if( lBachTPCClusters < leastnumberofclusters ) leastnumberofclusters = lBachTPCClusters;
        
        lInvMassLambdaAsCascDghter	= xi->GetEffMass();
        // This value shouldn't change, whatever the working hyp. is : Xi-, Xi+, Omega-, Omega+
        lDcaV0DaughtersXi 		= xi->GetDcaV0Daughters();
        //lV0Chi2Xi 			= xi->GetChi2V0();
        
        lV0CosineOfPointingAngleXi 	= xi->GetV0CosineOfPointingAngle( lBestPrimaryVtxPos[0],
                                                                     lBestPrimaryVtxPos[1],
                                                                     lBestPrimaryVtxPos[2] );
        //Modification: V0 CosPA wrt to Cascade decay vertex
        lV0CosineOfPointingAngleXiSpecial 	= xi->GetV0CosineOfPointingAngle( lPosXi[0],
                                                                             lPosXi[1],
                                                                             lPosXi[2] );
        
        lDcaV0ToPrimVertexXi 		= xi->GetD( lBestPrimaryVtxPos[0],
                                               lBestPrimaryVtxPos[1],
                                               lBestPrimaryVtxPos[2] );
		
        lDcaBachToPrimVertexXi = TMath::Abs( bachTrackXi->GetD(	lBestPrimaryVtxPos[0],
                                                               lBestPrimaryVtxPos[1],
                                                               lMagneticField  ) );
        // Note : AliExternalTrackParam::GetD returns an algebraic value ...
		
        xi->GetXYZ( lPosV0Xi[0],  lPosV0Xi[1], lPosV0Xi[2] );
        lV0RadiusXi		= TMath::Sqrt( lPosV0Xi[0]*lPosV0Xi[0]  +  lPosV0Xi[1]*lPosV0Xi[1] );
                
        lDcaPosToPrimVertexXi 	= TMath::Abs( pTrackXi	->GetD(	lBestPrimaryVtxPos[0],
                                                               lBestPrimaryVtxPos[1],
                                                               lMagneticField  )     );
        
        lDcaNegToPrimVertexXi 	= TMath::Abs( nTrackXi	->GetD(	lBestPrimaryVtxPos[0],
                                                               lBestPrimaryVtxPos[1],
                                                               lMagneticField  )     );
		
        // - II.Step 4 : around effective masses (ESD)
        // ~ change mass hypotheses to cover all the possibilities :  Xi-/+, Omega -/+
		
        if( bachTrackXi->Charge() < 0 )	{
            lV0quality = 0.;
            xi->ChangeMassHypothesis(lV0quality , 3312);
            // Calculate the effective mass of the Xi- candidate.
            // pdg code 3312 = Xi-
            lInvMassXiMinus = xi->GetEffMassXi();
            
            lV0quality = 0.;
            xi->ChangeMassHypothesis(lV0quality , 3334);
            // Calculate the effective mass of the Xi- candidate.
            // pdg code 3334 = Omega-
            lInvMassOmegaMinus = xi->GetEffMassXi();
            
            lV0quality = 0.;
            xi->ChangeMassHypothesis(lV0quality , 3312); 	// Back to default hyp.
        }// end if negative bachelor
        
        
        if( bachTrackXi->Charge() >  0 ){
            lV0quality = 0.;
            xi->ChangeMassHypothesis(lV0quality , -3312);
            // Calculate the effective mass of the Xi+ candidate.
            // pdg code -3312 = Xi+
            lInvMassXiPlus = xi->GetEffMassXi();
            
            lV0quality = 0.;
            xi->ChangeMassHypothesis(lV0quality , -3334);
            // Calculate the effective mass of the Xi+ candidate.
            // pdg code -3334  = Omega+
            lInvMassOmegaPlus = xi->GetEffMassXi();
            
            lV0quality = 0.;
            xi->ChangeMassHypothesis(lV0quality , -3312); 	// Back to "default" hyp.
        }// end if positive bachelor
        // - II.Step 6 : extra info for QA (ESD)
        // miscellaneous pieces of info that may help regarding data quality assessment.
        //-------------
        
        xi->GetPxPyPz( lXiMomX, lXiMomY, lXiMomZ );
        lXiTransvMom  	= TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY );
        lXiTotMom  	= TMath::Sqrt( lXiMomX*lXiMomX   + lXiMomY*lXiMomY   + lXiMomZ*lXiMomZ );
		
        xi->GetBPxPyPz(  lBachMomX,  lBachMomY,  lBachMomZ );
        //lBachTransvMom  = TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY );
        //lBachTotMom  	= TMath::Sqrt( lBachMomX*lBachMomX   + lBachMomY*lBachMomY  +  lBachMomZ*lBachMomZ  );
        
        lChargeXi = xi->Charge();
        
        //lV0toXiCosineOfPointingAngle = xi->GetV0CosineOfPointingAngle( lPosXi[0], lPosXi[1], lPosXi[2] );
        
        lRapXi    = xi->RapXi();
        lRapOmega = xi->RapOmega();
        //lEta      = xi->Eta();
        //lTheta    = xi->Theta() *180.0/TMath::Pi();
        //lPhi      = xi->Phi()   *180.0/TMath::Pi();
        //lAlphaXi  = xi->AlphaXi();
        //lPtArmXi  = xi->PtArmXi();
                
        //------------------------------------------------
        // Set Variables for adding to tree
        //------------------------------------------------
        
        /* 1*/		fTreeCascVarCharge	= lChargeXi;
        /* 2*/		if(lInvMassXiMinus!=0)    fTreeCascVarMassAsXi = lInvMassXiMinus;
        /* 2*/		if(lInvMassXiPlus!=0)     fTreeCascVarMassAsXi = lInvMassXiPlus;
        /* 3*/		if(lInvMassOmegaMinus!=0) fTreeCascVarMassAsOmega = lInvMassOmegaMinus;
        /* 3*/		if(lInvMassOmegaPlus!=0)  fTreeCascVarMassAsOmega = lInvMassOmegaPlus;
        /* 4*/		fTreeCascVarPt = lXiTransvMom;
        /* 4*/		fTreeCascVarPtMC = lXiTransvMomMC;
        /* 5*/		fTreeCascVarRapXi = lRapXi ;
        /* 5*/		fTreeCascVarRapMC = lRapMC ;
        /* 6*/		fTreeCascVarRapOmega = lRapOmega ;
        /* 7*/		fTreeCascVarDCACascDaughters = lDcaXiDaughters;
        /* 8*/		fTreeCascVarDCABachToPrimVtx = lDcaBachToPrimVertexXi;
        /* 9*/		fTreeCascVarDCAV0Daughters = lDcaV0DaughtersXi;
        /*10*/		fTreeCascVarDCAV0ToPrimVtx = lDcaV0ToPrimVertexXi;
        /*11*/		fTreeCascVarDCAPosToPrimVtx = lDcaPosToPrimVertexXi;
        /*12*/		fTreeCascVarDCANegToPrimVtx = lDcaNegToPrimVertexXi;
        /*13*/		fTreeCascVarCascCosPointingAngle = lXiCosineOfPointingAngle;
        /*14*/		fTreeCascVarCascRadius = lXiRadius;
        /*15*/		fTreeCascVarV0Mass = lInvMassLambdaAsCascDghter;
        /*16*/		fTreeCascVarV0CosPointingAngle = lV0CosineOfPointingAngleXi;
        /*16*/		fTreeCascVarV0CosPointingAngleSpecial = lV0CosineOfPointingAngleXiSpecial;
        /*17*/		fTreeCascVarV0Radius = lV0RadiusXi;
        /*20*/		fTreeCascVarLeastNbrClusters = leastnumberofclusters;
        /*21*/		fTreeCascVarCentrality = lPercentile; //multiplicity, whatever that may be
        /*23*/		fTreeCascVarDistOverTotMom = TMath::Sqrt(
                                                             TMath::Power( lPosXi[0] - lBestPrimaryVtxPos[0] , 2) +
                                                             TMath::Power( lPosXi[1] - lBestPrimaryVtxPos[1] , 2) +
                                                             TMath::Power( lPosXi[2] - lBestPrimaryVtxPos[2] , 2)
                                                             );
        /*23*/		fTreeCascVarDistOverTotMom /= (lXiTotMom+1e-13);
        
        //All vars not specified here: specified elsewhere!
        
        //------------------------------------------------
        // Fill Tree!
        //------------------------------------------------
        
        // The conditional is meant to decrease excessive
        // memory usage! Be careful when loosening the
        // cut!
        
        //Xi    Mass window: 150MeV wide
        //Omega mass window: 150MeV wide

        if( (fTreeCascVarMassAsXi<1.321+fCascadeMassWindow&&fTreeCascVarMassAsXi>1.321-fCascadeMassWindow) ||
                (fTreeCascVarMassAsOmega<1.672+fCascadeMassWindow&&fTreeCascVarMassAsOmega>1.672-fCascadeMassWindow) ) {
            //This cascade is useless until proven otherwise
            Bool_t lSaveThisCascade = kFALSE;
            
            //Extra selections in case this is supposed to be super-filtered
            //Inspired on tricks used for the V0 analysis in Pb-Pb
            if (TMath::Abs(fTreeCascVarNegEta) < 0.8 &&
                TMath::Abs(fTreeCascVarPosEta) < 0.8 &&
                TMath::Abs(fTreeCascVarBachEta) < 0.8 &&
                fTreeCascVarPt > fLowPtCutoff) { //beware ptMC and ptreco differences
                
                //Extra selections applied on a case-by-case basis:
                // (1) XiMinus
                if( fTreeCascVarCharge == -1 &&
                   TMath::Abs(fTreeCascVarMassAsXi-1.321)<fCascadeMassWindow &&
                   TMath::Abs(fTreeCascVarPosNSigmaProton) <= 4 &&
                   TMath::Abs(fTreeCascVarNegNSigmaPion  ) <= 4 &&
                   TMath::Abs(fTreeCascVarBachNSigmaPion ) <= 4 &&
                   TMath::Abs(fTreeCascVarRapXi          ) <= 0.5 ) {
                    lSaveThisCascade = kTRUE;
                }
                // (2) XiPlus
                if( fTreeCascVarCharge == +1 &&
                   TMath::Abs(fTreeCascVarMassAsXi-1.321)<fCascadeMassWindow &&
                   TMath::Abs(fTreeCascVarMassAsXi-1.321)<fCascadeMassWindow &&
                   TMath::Abs(fTreeCascVarPosNSigmaPion  ) <= 4 &&
                   TMath::Abs(fTreeCascVarNegNSigmaProton) <= 4 &&
                   TMath::Abs(fTreeCascVarBachNSigmaPion ) <= 4 &&
                   TMath::Abs(fTreeCascVarRapXi          ) <= 0.5 ) {
                    lSaveThisCascade = kTRUE;
                }
                // (3) OmegaMinus
                if( fTreeCascVarCharge == -1 &&
                   TMath::Abs(fTreeCascVarMassAsOmega-1.672)<fCascadeMassWindow &&
                   TMath::Abs(fTreeCascVarPosNSigmaProton) <= 4 &&
                   TMath::Abs(fTreeCascVarNegNSigmaPion  ) <= 4 &&
                   TMath::Abs(fTreeCascVarBachNSigmaKaon ) <= 4 &&
                   TMath::Abs(fTreeCascVarRapOmega       ) <= 0.5 ) {
                    lSaveThisCascade = kTRUE;
                }
                // (4) OmegaPlus
                if( fTreeCascVarCharge == +1 &&
                   TMath::Abs(fTreeCascVarMassAsOmega-1.672)<fCascadeMassWindow &&
                   TMath::Abs(fTreeCascVarPosNSigmaPion    ) <= 4 &&
                   TMath::Abs(fTreeCascVarNegNSigmaProton  ) <= 4 &&
                   TMath::Abs(fTreeCascVarBachNSigmaKaon   ) <= 4 &&
                   TMath::Abs(fTreeCascVarRapOmega         ) <= 0.5 ) {
                    lSaveThisCascade = kTRUE;
                }
            }
            if (lSaveThisCascade && fkSaveTree ) fTreeCascade -> Fill() ;
        }
        
        //------------------------------------------------
        // Fill tree over.
        //------------------------------------------------
    }// end of the Cascade loop (ESD or AOD)
    
    // Post output data.
    PostData(1, fListHist);
    PostData(2, fTreeCascade);
}

//________________________________________________________________________
void AliAnalysisTaskExtractCascadePbPbRun2::Terminate(Option_t *)
{
    // Draw result to the screen
    // Called once at the end of the query
    
    TList *cRetrievedList = 0x0;
    cRetrievedList = (TList*)GetOutputData(1);
    if(!cRetrievedList){
        Printf("ERROR - AliAnalysisTaskExtractCascadePbPbRun2 : ouput data container list not available\n");
        return;
    }
	
    fHistCentrality = dynamic_cast<TH1F*> (  cRetrievedList->FindObject("fHistCentrality")  );
    if (!fHistCentrality) {
        Printf("ERROR - AliAnalysisTaskExtractCascadePbPbRun2 : fHistCentrality not available");
        return;
    }
    
    TCanvas *canCheck = new TCanvas("AliAnalysisTaskExtractCascadePbPbRun2","V0 Multiplicity",10,10,510,510);
    canCheck->cd(1)->SetLogy();
    
    fHistCentrality->SetMarkerStyle(22);
    fHistCentrality->DrawCopy("E");
}

//----------------------------------------------------------------------------

Double_t AliAnalysisTaskExtractCascadePbPbRun2::MyRapidity(Double_t rE, Double_t rPz) const
{
    // Local calculation for rapidity
    Double_t ReturnValue = -100;
    if( (rE-rPz+1.e-13) != 0 && (rE+rPz) != 0 ){
        ReturnValue =  0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
    }
    return ReturnValue;
}
