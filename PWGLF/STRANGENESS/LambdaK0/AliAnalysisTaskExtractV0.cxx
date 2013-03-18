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
// --- Added bits of code for checking V0s 
//      (from AliAnalysisTaskCheckStrange.cxx)
//
//  --- Algorithm Description 
//   1. Perform Physics Selection
//   2. Perform Primary Vertex |z|<10cm selection
//   3. Perform Primary Vertex NoTPCOnly vertexing selection (>0 contrib.)
//   4. Perform Pileup Rejection
//   5. Analysis Loops: 
//    5a. Fill TTree object with V0 information, candidates
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
#include "THnSparse.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "AliLog.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"

#include "AliCFContainer.h"
#include "AliMultiplicity.h"

#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDUtils.h"
#include "AliESDHeader.h"

#include "AliAnalysisUtils.h"
#include "AliAnalysisTaskExtractV0.h"

//debugging purposes
#include "TObjectTable.h"

ClassImp(AliAnalysisTaskExtractV0)

AliAnalysisTaskExtractV0::AliAnalysisTaskExtractV0() 
  : AliAnalysisTaskSE(), fListHistV0(0), fTree(0), fPIDResponse(0),fESDtrackCuts(0), fUtils(0),
  fkIsNuclear     ( kFALSE ), 
  fkSwitchINT7    ( kFALSE ),
  fkUseOnTheFly   ( kFALSE ),
  fkTakeAllTracks ( kFALSE ),
  fCentralityEstimator("V0M"),
  fkLightWeight   ( kFALSE ),
  fkFastOnly      ( "" ),
  fkpAVertexSelection( kFALSE ),
  fkRunV0Vertexer ( kFALSE ),
  fkRejectPileup  ( kTRUE ),
  fkSpecialExecution( kFALSE ),
//------------------------------------------------
// Initialize 
	fTreeVariableChi2V0(0),
	fTreeVariableDcaV0Daughters(0),
	fTreeVariableDcaV0ToPrimVertex(0),
	fTreeVariableDcaPosToPrimVertex(0),
	fTreeVariableDcaNegToPrimVertex(0),
	fTreeVariableV0CosineOfPointingAngle(0),
	fTreeVariableV0Radius(0),
	fTreeVariablePt(0),
	fTreeVariableRapK0Short(0),
	fTreeVariableRapLambda(0),
	fTreeVariableInvMassK0s(0),
	fTreeVariableInvMassLambda(0),
	fTreeVariableInvMassAntiLambda(0),
	fTreeVariableAlphaV0(0),
	fTreeVariablePtArmV0(0),
	fTreeVariableNegTotMomentum(0),
	fTreeVariablePosTotMomentum(0),
	fTreeVariableNegdEdxSig(0),
	fTreeVariablePosdEdxSig(0),
	fTreeVariableNegEta(0),
	fTreeVariablePosEta(0),

	fTreeVariableNSigmasPosProton(0),
	fTreeVariableNSigmasPosPion(0),
	fTreeVariableNSigmasNegProton(0),
	fTreeVariableNSigmasNegPion(0),
	
	fTreeVariableDistOverTotMom(0),
	fTreeVariableLeastNbrCrossedRows(0),
	fTreeVariableLeastRatioCrossedRowsOverFindable(0),
	fTreeVariableMultiplicity(0),
	fTreeVariableMultiplicityV0A(0),
	fTreeVariableMultiplicityZNA(0),
	fTreeVariableMultiplicityTRK(0),
	fTreeVariableMultiplicitySPD(0),
  
  fTreeVariableRunNumber(0),
  fTreeVariableEventNumber(0),
  
  fTreeVariableV0x(0),
  fTreeVariableV0y(0),
  fTreeVariableV0z(0),

  fTreeVariableV0Px(0),
  fTreeVariableV0Py(0),
  fTreeVariableV0Pz(0),

  fTreeVariablePVx(0),
  fTreeVariablePVy(0),
  fTreeVariablePVz(0),

  fTreeVariableNegTrackStatus(0),
  fTreeVariablePosTrackStatus(0),

  fTreeVariableNegTPCSignal(0),
  fTreeVariablePosTPCSignal(0),
  fTreeVariableNegInnerP(0),
  fTreeVariablePosInnerP(0),

  fTreeVariableNegPx(0),
  fTreeVariableNegPy(0),
  fTreeVariableNegPz(0),
  fTreeVariablePosPx(0),
  fTreeVariablePosPy(0),
  fTreeVariablePosPz(0),


//------------------------------------------------
// HISTOGRAMS
// --- Filled on an Event-by-event basis
//------------------------------------------------
   fHistV0MultiplicityBeforeTrigSel(0),
   fHistV0MultiplicityForTrigEvt(0),
   fHistV0MultiplicityForSelEvt(0),
   fHistV0MultiplicityForSelEvtNoTPCOnly(0),
   fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup(0),
   fHistMultiplicityBeforeTrigSel(0),
   fHistMultiplicityForTrigEvt(0),
   fHistMultiplicity(0),
   fHistMultiplicityNoTPCOnly(0),
   fHistMultiplicityNoTPCOnlyNoPileup(0),

//V0A Centrality
fHistMultiplicityV0ABeforeTrigSel(0),
fHistMultiplicityV0AForTrigEvt(0),
fHistMultiplicityV0A(0),
fHistMultiplicityV0ANoTPCOnly(0),
fHistMultiplicityV0ANoTPCOnlyNoPileup(0),

//ZNA Centrality
fHistMultiplicityZNABeforeTrigSel(0),
fHistMultiplicityZNAForTrigEvt(0),
fHistMultiplicityZNA(0),
fHistMultiplicityZNANoTPCOnly(0),
fHistMultiplicityZNANoTPCOnlyNoPileup(0),

//TRK Centrality
fHistMultiplicityTRKBeforeTrigSel(0),
fHistMultiplicityTRKForTrigEvt(0),
fHistMultiplicityTRK(0),
fHistMultiplicityTRKNoTPCOnly(0),
fHistMultiplicityTRKNoTPCOnlyNoPileup(0),

//SPD Centrality
fHistMultiplicitySPDBeforeTrigSel(0),
fHistMultiplicitySPDForTrigEvt(0),
fHistMultiplicitySPD(0),
fHistMultiplicitySPDNoTPCOnly(0),
fHistMultiplicitySPDNoTPCOnlyNoPileup(0),

  //Raw Data for Vertex Z position estimator change
	f2dHistMultiplicityVsVertexZBeforeTrigSel(0),
	f2dHistMultiplicityVsVertexZForTrigEvt(0),
	f2dHistMultiplicityVsVertexZ(0),
	f2dHistMultiplicityVsVertexZNoTPCOnly(0),
	f2dHistMultiplicityVsVertexZNoTPCOnlyNoPileup(0),

   fHistPVx(0),
   fHistPVy(0),
   fHistPVz(0),
   fHistPVxAnalysis(0),
   fHistPVyAnalysis(0),
   fHistPVzAnalysis(0),
   fHistSwappedV0Counter(0)
{
  // Dummy Constructor
  for(Int_t iV0selIdx   = 0; iV0selIdx   < 7; iV0selIdx++   ) { fV0Sels          [iV0selIdx   ] = -1.; }
}

AliAnalysisTaskExtractV0::AliAnalysisTaskExtractV0(const char *name) 
  : AliAnalysisTaskSE(name), fListHistV0(0), fTree(0), fPIDResponse(0),fESDtrackCuts(0), fUtils(0),
  fkIsNuclear     ( kFALSE ), 
  fkSwitchINT7    ( kFALSE ),
  fkUseOnTheFly   ( kFALSE ),
  fkTakeAllTracks ( kFALSE ),
  fCentralityEstimator("V0M"),
  fkLightWeight   ( kFALSE ),
  fkFastOnly      ( "" ),
  fkpAVertexSelection( kFALSE ),
  fkRunV0Vertexer ( kFALSE ),
  fkRejectPileup  ( kTRUE ),
  fkSpecialExecution( kFALSE ),
//------------------------------------------------
// Initialize 
	fTreeVariableChi2V0(0),
	fTreeVariableDcaV0Daughters(0),
	fTreeVariableDcaV0ToPrimVertex(0),
	fTreeVariableDcaPosToPrimVertex(0),
	fTreeVariableDcaNegToPrimVertex(0),
	fTreeVariableV0CosineOfPointingAngle(0),
	fTreeVariableV0Radius(0),
	fTreeVariablePt(0),
	fTreeVariableRapK0Short(0),
	fTreeVariableRapLambda(0),
	fTreeVariableInvMassK0s(0),
	fTreeVariableInvMassLambda(0),
	fTreeVariableInvMassAntiLambda(0),
	fTreeVariableAlphaV0(0),
	fTreeVariablePtArmV0(0),
	fTreeVariableNegTotMomentum(0),
	fTreeVariablePosTotMomentum(0),
	fTreeVariableNegdEdxSig(0),
	fTreeVariablePosdEdxSig(0),
	fTreeVariableNegEta(0),
	fTreeVariablePosEta(0),

	fTreeVariableNSigmasPosProton(0),
	fTreeVariableNSigmasPosPion(0),
	fTreeVariableNSigmasNegProton(0),
	fTreeVariableNSigmasNegPion(0),
	
	fTreeVariableDistOverTotMom(0),
	fTreeVariableLeastNbrCrossedRows(0),
	fTreeVariableLeastRatioCrossedRowsOverFindable(0),
	fTreeVariableMultiplicity(0),
  fTreeVariableMultiplicityV0A(0),
  fTreeVariableMultiplicityZNA(0),
  fTreeVariableMultiplicityTRK(0),
  fTreeVariableMultiplicitySPD(0),

  fTreeVariableRunNumber(0),
  fTreeVariableEventNumber(0),
  
  fTreeVariableV0x(0),
  fTreeVariableV0y(0),
  fTreeVariableV0z(0),

  fTreeVariableV0Px(0),
  fTreeVariableV0Py(0),
  fTreeVariableV0Pz(0),

  fTreeVariablePVx(0),
  fTreeVariablePVy(0),
  fTreeVariablePVz(0),

  fTreeVariableNegTrackStatus(0),
  fTreeVariablePosTrackStatus(0),

  fTreeVariableNegTPCSignal(0),
  fTreeVariablePosTPCSignal(0),
  fTreeVariableNegInnerP(0),
  fTreeVariablePosInnerP(0),

  fTreeVariableNegPx(0),
  fTreeVariableNegPy(0),
  fTreeVariableNegPz(0),
  fTreeVariablePosPx(0),
  fTreeVariablePosPy(0),
  fTreeVariablePosPz(0),

//------------------------------------------------
// HISTOGRAMS
// --- Filled on an Event-by-event basis
//------------------------------------------------
   fHistV0MultiplicityBeforeTrigSel(0),
   fHistV0MultiplicityForTrigEvt(0),
   fHistV0MultiplicityForSelEvt(0),
   fHistV0MultiplicityForSelEvtNoTPCOnly(0),
   fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup(0),
   fHistMultiplicityBeforeTrigSel(0),
   fHistMultiplicityForTrigEvt(0),
   fHistMultiplicity(0),
   fHistMultiplicityNoTPCOnly(0),
   fHistMultiplicityNoTPCOnlyNoPileup(0),


//V0A Centrality
fHistMultiplicityV0ABeforeTrigSel(0),
fHistMultiplicityV0AForTrigEvt(0),
fHistMultiplicityV0A(0),
fHistMultiplicityV0ANoTPCOnly(0),
fHistMultiplicityV0ANoTPCOnlyNoPileup(0),

//ZNA Centrality
fHistMultiplicityZNABeforeTrigSel(0),
fHistMultiplicityZNAForTrigEvt(0),
fHistMultiplicityZNA(0),
fHistMultiplicityZNANoTPCOnly(0),
fHistMultiplicityZNANoTPCOnlyNoPileup(0),

//TRK Centrality
fHistMultiplicityTRKBeforeTrigSel(0),
fHistMultiplicityTRKForTrigEvt(0),
fHistMultiplicityTRK(0),
fHistMultiplicityTRKNoTPCOnly(0),
fHistMultiplicityTRKNoTPCOnlyNoPileup(0),

//SPD Centrality
fHistMultiplicitySPDBeforeTrigSel(0),
fHistMultiplicitySPDForTrigEvt(0),
fHistMultiplicitySPD(0),
fHistMultiplicitySPDNoTPCOnly(0),
fHistMultiplicitySPDNoTPCOnlyNoPileup(0),

  //Raw Data for Vertex Z position estimator change
	f2dHistMultiplicityVsVertexZBeforeTrigSel(0),
	f2dHistMultiplicityVsVertexZForTrigEvt(0),
	f2dHistMultiplicityVsVertexZ(0),
	f2dHistMultiplicityVsVertexZNoTPCOnly(0),
	f2dHistMultiplicityVsVertexZNoTPCOnlyNoPileup(0),

   fHistPVx(0),
   fHistPVy(0),
   fHistPVz(0),
   fHistPVxAnalysis(0),
   fHistPVyAnalysis(0),
   fHistPVzAnalysis(0),
   fHistSwappedV0Counter(0)
{
  // Constructor
  // Set Loose cuts or not here...
  // REALLY LOOSE? Be careful when attempting to run over PbPb if fkRunV0Vertexer is set! 
  fV0Sels[0] =  33.  ;  // max allowed chi2
  fV0Sels[1] =   0.02;  // min allowed impact parameter for the 1st daughter (LHC09a4 : 0.05)
  fV0Sels[2] =   0.02;  // min allowed impact parameter for the 2nd daughter (LHC09a4 : 0.05)
  fV0Sels[3] =   2.0 ;  // max allowed DCA between the daughter tracks       (LHC09a4 : 0.5)
  fV0Sels[4] =   0.95;  // min allowed cosine of V0's pointing angle         (LHC09a4 : 0.99)
  fV0Sels[5] =   0.5 ;  // min radius of the fiducial volume                 (LHC09a4 : 0.2)
  fV0Sels[6] = 200.  ;  // max radius of the fiducial volume                 (LHC09a4 : 100.0)
  
  // Output slot #0 writes into a TList container (Lambda Histos and fTree)
   DefineOutput(1, TList::Class());
   DefineOutput(2, TTree::Class());
}


AliAnalysisTaskExtractV0::~AliAnalysisTaskExtractV0()
{
//------------------------------------------------
// DESTRUCTOR
//------------------------------------------------

   if (fListHistV0){
      delete fListHistV0;
      fListHistV0 = 0x0;
   }
   if (fTree){
      delete fTree;
      fTree = 0x0;
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
void AliAnalysisTaskExtractV0::UserCreateOutputObjects()
{

   //Create File-resident Tree, please.
   OpenFile(2);
   // Called once
   fTree = new TTree("fTree","V0Candidates");

//------------------------------------------------
// fTree Branch definitions
//------------------------------------------------

//-----------BASIC-INFO---------------------------
/* 1*/  fTree->Branch("fTreeVariableChi2V0",&fTreeVariableChi2V0,"fTreeVariableChi2V0/F");
/* 2*/  fTree->Branch("fTreeVariableDcaV0Daughters",&fTreeVariableDcaV0Daughters,"fTreeVariableDcaV0Daughters/F");
/* 3*/	fTree->Branch("fTreeVariableDcaPosToPrimVertex",&fTreeVariableDcaPosToPrimVertex,"fTreeVariableDcaPosToPrimVertex/F");
/* 4*/	fTree->Branch("fTreeVariableDcaNegToPrimVertex",&fTreeVariableDcaNegToPrimVertex,"fTreeVariableDcaNegToPrimVertex/F");
/* 5*/	fTree->Branch("fTreeVariableV0Radius",&fTreeVariableV0Radius,"fTreeVariableV0Radius/F");
/* 6*/	fTree->Branch("fTreeVariablePt",&fTreeVariablePt,"fTreeVariablePt/F");
/* 7*/	fTree->Branch("fTreeVariableRapK0Short",&fTreeVariableRapK0Short,"fTreeVariableRapK0Short/F");
/* 8*/	fTree->Branch("fTreeVariableRapLambda",&fTreeVariableRapLambda,"fTreeVariableRapLambda/F");
/* 9*/	fTree->Branch("fTreeVariableInvMassK0s",&fTreeVariableInvMassK0s,"fTreeVariableInvMassK0s/F");
/*10*/	fTree->Branch("fTreeVariableInvMassLambda",&fTreeVariableInvMassLambda,"fTreeVariableInvMassLambda/F");
/*11*/	fTree->Branch("fTreeVariableInvMassAntiLambda",&fTreeVariableInvMassAntiLambda,"fTreeVariableInvMassAntiLambda/F");
/*12*/	fTree->Branch("fTreeVariableV0CosineOfPointingAngle",&fTreeVariableV0CosineOfPointingAngle,"fTreeVariableV0CosineOfPointingAngle/F");
/*13*/	fTree->Branch("fTreeVariableAlphaV0",&fTreeVariableAlphaV0,"fTreeVariableAlphaV0/F");
/*14*/	fTree->Branch("fTreeVariablePtArmV0",&fTreeVariablePtArmV0,"fTreeVariablePtArmV0/F");
/*15*/	fTree->Branch("fTreeVariableLeastNbrCrossedRows",&fTreeVariableLeastNbrCrossedRows,"fTreeVariableLeastNbrCrossedRows/I");
/*16*/	fTree->Branch("fTreeVariableLeastRatioCrossedRowsOverFindable",&fTreeVariableLeastRatioCrossedRowsOverFindable,"fTreeVariableLeastRatioCrossedRowsOverFindable/F");
//-----------MULTIPLICITY-INFO--------------------
/*17*/	fTree->Branch("fTreeVariableMultiplicity",&fTreeVariableMultiplicity,"fTreeVariableMultiplicity/I");
  /*17*/	fTree->Branch("fTreeVariableMultiplicityV0A",&fTreeVariableMultiplicityV0A,"fTreeVariableMultiplicityV0A/I");
  /*17*/	fTree->Branch("fTreeVariableMultiplicityZNA",&fTreeVariableMultiplicityZNA,"fTreeVariableMultiplicityZNA/I");
  /*17*/	fTree->Branch("fTreeVariableMultiplicityTRK",&fTreeVariableMultiplicityTRK,"fTreeVariableMultiplicityTRK/I");
  /*17*/	fTree->Branch("fTreeVariableMultiplicitySPD",&fTreeVariableMultiplicitySPD,"fTreeVariableMultiplicitySPD/I");
//------------------------------------------------
/*18*/	fTree->Branch("fTreeVariableDistOverTotMom",&fTreeVariableDistOverTotMom,"fTreeVariableDistOverTotMom/F");
/*19*/	fTree->Branch("fTreeVariableNSigmasPosProton",&fTreeVariableNSigmasPosProton,"fTreeVariableNSigmasPosProton/F");
/*20*/	fTree->Branch("fTreeVariableNSigmasPosPion",&fTreeVariableNSigmasPosPion,"fTreeVariableNSigmasPosPion/F");
/*21*/	fTree->Branch("fTreeVariableNSigmasNegProton",&fTreeVariableNSigmasNegProton,"fTreeVariableNSigmasNegProton/F");
/*22*/	fTree->Branch("fTreeVariableNSigmasNegPion",&fTreeVariableNSigmasNegPion,"fTreeVariableNSigmasNegPion/F");
/*23*/	fTree->Branch("fTreeVariableNegEta",&fTreeVariableNegEta,"fTreeVariableNegEta/F");
/*24*/	fTree->Branch("fTreeVariablePosEta",&fTreeVariablePosEta,"fTreeVariablePosEta/F");
/*25*/	fTree->Branch("fTreeVariableRunNumber",&fTreeVariableRunNumber,"fTreeVariableRunNumber/I");
/*26*/	fTree->Branch("fTreeVariableEventNumber",&fTreeVariableEventNumber,"fTreeVariableEventNumber/l");
  
  if( fkLightWeight == kFALSE ){
//-----------FOR CTAU DEBUGGING: Full Phase Space + Decay Position Info 
        fTree->Branch("fTreeVariablePVx",&fTreeVariablePVx,"fTreeVariablePVx/F");
        fTree->Branch("fTreeVariablePVy",&fTreeVariablePVy,"fTreeVariablePVy/F");
        fTree->Branch("fTreeVariablePVz",&fTreeVariablePVz,"fTreeVariablePVz/F");

        fTree->Branch("fTreeVariableV0x",&fTreeVariableV0x,"fTreeVariableV0x/F");
        fTree->Branch("fTreeVariableV0y",&fTreeVariableV0y,"fTreeVariableV0y/F");
        fTree->Branch("fTreeVariableV0z",&fTreeVariableV0z,"fTreeVariableV0z/F");

        fTree->Branch("fTreeVariableV0Px",&fTreeVariableV0Px,"fTreeVariableV0Px/F");
        fTree->Branch("fTreeVariableV0Py",&fTreeVariableV0Py,"fTreeVariableV0Py/F");
        fTree->Branch("fTreeVariableV0Pz",&fTreeVariableV0Pz,"fTreeVariableV0Pz/F");
  
        fTree->Branch("fTreeVariableNegTrackStatus",&fTreeVariableNegTrackStatus,"fTreeVariableNegTrackStatus/l");
        fTree->Branch("fTreeVariablePosTrackStatus",&fTreeVariablePosTrackStatus,"fTreeVariablePosTrackStatus/l");
  }
  if( fkSpecialExecution == kTRUE ){
    fTree->Branch("fTreeVariablePosTPCSignal",&fTreeVariablePosTPCSignal,"fTreeVariablePosTPCSignal/F");
    fTree->Branch("fTreeVariableNegTPCSignal",&fTreeVariableNegTPCSignal,"fTreeVariableNegTPCSignal/F");
    fTree->Branch("fTreeVariablePosInnerP",&fTreeVariablePosInnerP,"fTreeVariablePosInnerP/F");
    fTree->Branch("fTreeVariableNegInnerP",&fTreeVariableNegInnerP,"fTreeVariableNegInnerP/F");
    
    fTree->Branch("fTreeVariablePosPx",&fTreeVariablePosPx,"fTreeVariablePosPx/F");
    fTree->Branch("fTreeVariablePosPy",&fTreeVariablePosPy,"fTreeVariablePosPy/F");
    fTree->Branch("fTreeVariablePosPz",&fTreeVariablePosPz,"fTreeVariablePosPz/F");
    fTree->Branch("fTreeVariableNegPx",&fTreeVariableNegPx,"fTreeVariableNegPx/F");
    fTree->Branch("fTreeVariableNegPy",&fTreeVariableNegPy,"fTreeVariableNegPy/F");
    fTree->Branch("fTreeVariableNegPz",&fTreeVariableNegPz,"fTreeVariableNegPz/F");
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
  }
  if(! fUtils ){
    fUtils = new AliAnalysisUtils();
  }

//------------------------------------------------
// V0 Multiplicity Histograms
//------------------------------------------------

   // Create histograms
   //Create File-resident Tree, please.
   OpenFile(1);

   fListHistV0 = new TList();
   fListHistV0->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner

   if(! fHistV0MultiplicityBeforeTrigSel) {
      fHistV0MultiplicityBeforeTrigSel = new TH1F("fHistV0MultiplicityBeforeTrigSel", 
         "V0s per event (before Trig. Sel.);Nbr of V0s/Evt;Events", 
         25, 0, 25);
      fListHistV0->Add(fHistV0MultiplicityBeforeTrigSel);
   }
           
   if(! fHistV0MultiplicityForTrigEvt) {
      fHistV0MultiplicityForTrigEvt = new TH1F("fHistV0MultiplicityForTrigEvt", 
         "V0s per event (for triggered evt);Nbr of V0s/Evt;Events", 
         25, 0, 25);
      fListHistV0->Add(fHistV0MultiplicityForTrigEvt);
   }

   if(! fHistV0MultiplicityForSelEvt) {
      fHistV0MultiplicityForSelEvt = new TH1F("fHistV0MultiplicityForSelEvt", 
         "V0s per event;Nbr of V0s/Evt;Events", 
         25, 0, 25);
      fListHistV0->Add(fHistV0MultiplicityForSelEvt);
   }

   if(! fHistV0MultiplicityForSelEvtNoTPCOnly) {
      fHistV0MultiplicityForSelEvtNoTPCOnly = new TH1F("fHistV0MultiplicityForSelEvtNoTPCOnly", 
         "V0s per event;Nbr of V0s/Evt;Events", 
         25, 0, 25);
      fListHistV0->Add(fHistV0MultiplicityForSelEvtNoTPCOnly);
   }

   if(! fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup) {
      fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup = new TH1F("fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup", 
         "V0s per event;Nbr of V0s/Evt;Events", 
         25, 0, 25);
      fListHistV0->Add(fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup);
   }

//------------------------------------------------
// Track Multiplicity Histograms
//------------------------------------------------

  //Default V0M Centrality (if PbPb) 
   if(! fHistMultiplicityBeforeTrigSel) {
      fHistMultiplicityBeforeTrigSel = new TH1F("fHistMultiplicityBeforeTrigSel", 
         "Tracks per event;Nbr of Tracks;Events", 
         200, 0, 200); 		
      fListHistV0->Add(fHistMultiplicityBeforeTrigSel);
   }
   if(! fHistMultiplicityForTrigEvt) {
      fHistMultiplicityForTrigEvt = new TH1F("fHistMultiplicityForTrigEvt", 
         "Tracks per event;Nbr of Tracks;Events", 
         200, 0, 200); 		
      fListHistV0->Add(fHistMultiplicityForTrigEvt);
   }
   if(! fHistMultiplicity) {
      fHistMultiplicity = new TH1F("fHistMultiplicity", 
         "Tracks per event;Nbr of Tracks;Events", 
         200, 0, 200); 		
      fListHistV0->Add(fHistMultiplicity);
   }
   if(! fHistMultiplicityNoTPCOnly) {
      fHistMultiplicityNoTPCOnly = new TH1F("fHistMultiplicityNoTPCOnly", 
         "Tracks per event;Nbr of Tracks;Events", 
         200, 0, 200); 		
      fListHistV0->Add(fHistMultiplicityNoTPCOnly);
   }
   if(! fHistMultiplicityNoTPCOnlyNoPileup) {
      fHistMultiplicityNoTPCOnlyNoPileup = new TH1F("fHistMultiplicityNoTPCOnlyNoPileup", 
         "Tracks per event;Nbr of Tracks;Events", 
         200, 0, 200); 		
      fListHistV0->Add(fHistMultiplicityNoTPCOnlyNoPileup);
   }
  
  //V0A Centrality (if PbPb / pPb)
  if(! fHistMultiplicityV0ABeforeTrigSel) {
    fHistMultiplicityV0ABeforeTrigSel = new TH1F("fHistMultiplicityV0ABeforeTrigSel",
                                              "Centrality Distribution: V0A;V0A Centrality;Events",
                                              200, 0, 200);
    fListHistV0->Add(fHistMultiplicityV0ABeforeTrigSel);
  }
  if(! fHistMultiplicityV0AForTrigEvt) {
    fHistMultiplicityV0AForTrigEvt = new TH1F("fHistMultiplicityV0AForTrigEvt",
                                              "Centrality Distribution: V0A;V0A Centrality;Events",
                                           200, 0, 200);
    fListHistV0->Add(fHistMultiplicityV0AForTrigEvt);
  }
  if(! fHistMultiplicityV0A) {
    fHistMultiplicityV0A = new TH1F("fHistMultiplicityV0A",
                                              "Centrality Distribution: V0A;V0A Centrality;Events",
                                 200, 0, 200);
    fListHistV0->Add(fHistMultiplicityV0A);
  }
  if(! fHistMultiplicityV0ANoTPCOnly) {
    fHistMultiplicityV0ANoTPCOnly = new TH1F("fHistMultiplicityV0ANoTPCOnly",
                                              "Centrality Distribution: V0A;V0A Centrality;Events",
                                          200, 0, 200);
    fListHistV0->Add(fHistMultiplicityV0ANoTPCOnly);
  }
  if(! fHistMultiplicityV0ANoTPCOnlyNoPileup) {
    fHistMultiplicityV0ANoTPCOnlyNoPileup = new TH1F("fHistMultiplicityV0ANoTPCOnlyNoPileup",
                                              "Centrality Distribution: V0A;V0A Centrality;Events",
                                                  200, 0, 200); 		
    fListHistV0->Add(fHistMultiplicityV0ANoTPCOnlyNoPileup);
  }
  
  //ZNA Centrality (if PbPb / pPb)
  if(! fHistMultiplicityZNABeforeTrigSel) {
    fHistMultiplicityZNABeforeTrigSel = new TH1F("fHistMultiplicityZNABeforeTrigSel",
                                                 "Centrality Distribution: ZNA;ZNA Centrality;Events",
                                                 200, 0, 200);
    fListHistV0->Add(fHistMultiplicityZNABeforeTrigSel);
  }
  if(! fHistMultiplicityZNAForTrigEvt) {
    fHistMultiplicityZNAForTrigEvt = new TH1F("fHistMultiplicityZNAForTrigEvt",
                                              "Centrality Distribution: ZNA;ZNA Centrality;Events",
                                              200, 0, 200);
    fListHistV0->Add(fHistMultiplicityZNAForTrigEvt);
  }
  if(! fHistMultiplicityZNA) {
    fHistMultiplicityZNA = new TH1F("fHistMultiplicityZNA",
                                    "Centrality Distribution: ZNA;ZNA Centrality;Events",
                                    200, 0, 200);
    fListHistV0->Add(fHistMultiplicityZNA);
  }
  if(! fHistMultiplicityZNANoTPCOnly) {
    fHistMultiplicityZNANoTPCOnly = new TH1F("fHistMultiplicityZNANoTPCOnly",
                                             "Centrality Distribution: ZNA;ZNA Centrality;Events",
                                             200, 0, 200);
    fListHistV0->Add(fHistMultiplicityZNANoTPCOnly);
  }
  if(! fHistMultiplicityZNANoTPCOnlyNoPileup) {
    fHistMultiplicityZNANoTPCOnlyNoPileup = new TH1F("fHistMultiplicityZNANoTPCOnlyNoPileup",
                                                     "Centrality Distribution: ZNA;ZNA Centrality;Events",
                                                     200, 0, 200);
    fListHistV0->Add(fHistMultiplicityZNANoTPCOnlyNoPileup);
  }
  
  //TRK Centrality (if PbPb / pPb)
  if(! fHistMultiplicityTRKBeforeTrigSel) {
    fHistMultiplicityTRKBeforeTrigSel = new TH1F("fHistMultiplicityTRKBeforeTrigSel",
                                                 "Centrality Distribution: TRK;TRK Centrality;Events",
                                                 200, 0, 200);
    fListHistV0->Add(fHistMultiplicityTRKBeforeTrigSel);
  }
  if(! fHistMultiplicityTRKForTrigEvt) {
    fHistMultiplicityTRKForTrigEvt = new TH1F("fHistMultiplicityTRKForTrigEvt",
                                              "Centrality Distribution: TRK;TRK Centrality;Events",
                                              200, 0, 200);
    fListHistV0->Add(fHistMultiplicityTRKForTrigEvt);
  }
  if(! fHistMultiplicityTRK) {
    fHistMultiplicityTRK = new TH1F("fHistMultiplicityTRK",
                                    "Centrality Distribution: TRK;TRK Centrality;Events",
                                    200, 0, 200);
    fListHistV0->Add(fHistMultiplicityTRK);
  }
  if(! fHistMultiplicityTRKNoTPCOnly) {
    fHistMultiplicityTRKNoTPCOnly = new TH1F("fHistMultiplicityTRKNoTPCOnly",
                                             "Centrality Distribution: TRK;TRK Centrality;Events",
                                             200, 0, 200);
    fListHistV0->Add(fHistMultiplicityTRKNoTPCOnly);
  }
  if(! fHistMultiplicityTRKNoTPCOnlyNoPileup) {
    fHistMultiplicityTRKNoTPCOnlyNoPileup = new TH1F("fHistMultiplicityTRKNoTPCOnlyNoPileup",
                                                     "Centrality Distribution: TRK;TRK Centrality;Events",
                                                     200, 0, 200);
    fListHistV0->Add(fHistMultiplicityTRKNoTPCOnlyNoPileup);
  }
  
  //SPD Centrality (if PbPb / pPb)
  if(! fHistMultiplicitySPDBeforeTrigSel) {
    fHistMultiplicitySPDBeforeTrigSel = new TH1F("fHistMultiplicitySPDBeforeTrigSel",
                                                 "Centrality Distribution: SPD;SPD Centrality;Events",
                                                 200, 0, 200);
    fListHistV0->Add(fHistMultiplicitySPDBeforeTrigSel);
  }
  if(! fHistMultiplicitySPDForTrigEvt) {
    fHistMultiplicitySPDForTrigEvt = new TH1F("fHistMultiplicitySPDForTrigEvt",
                                              "Centrality Distribution: SPD;SPD Centrality;Events",
                                              200, 0, 200);
    fListHistV0->Add(fHistMultiplicitySPDForTrigEvt);
  }
  if(! fHistMultiplicitySPD) {
    fHistMultiplicitySPD = new TH1F("fHistMultiplicitySPD",
                                    "Centrality Distribution: SPD;SPD Centrality;Events",
                                    200, 0, 200);
    fListHistV0->Add(fHistMultiplicitySPD);
  }
  if(! fHistMultiplicitySPDNoTPCOnly) {
    fHistMultiplicitySPDNoTPCOnly = new TH1F("fHistMultiplicitySPDNoTPCOnly",
                                             "Centrality Distribution: SPD;SPD Centrality;Events",
                                             200, 0, 200);
    fListHistV0->Add(fHistMultiplicitySPDNoTPCOnly);
  }
  if(! fHistMultiplicitySPDNoTPCOnlyNoPileup) {
    fHistMultiplicitySPDNoTPCOnlyNoPileup = new TH1F("fHistMultiplicitySPDNoTPCOnlyNoPileup",
                                                     "Centrality Distribution: SPD;SPD Centrality;Events",
                                                     200, 0, 200);
    fListHistV0->Add(fHistMultiplicitySPDNoTPCOnlyNoPileup);
  }
  
  //Raw Data for Vertex Z position estimator change
	//TH2F    *f2dHistMultiplicityVsVertexZBeforeTrigSel; 	        //! multiplicity distribution    
	//TH2F    *f2dHistMultiplicityVsVertexZForTrigEvt;  		        //! multiplicity distribution
	//TH2F    *f2dHistMultiplicityVsVertexZ;     					        //! multiplicity distribution
	//TH2F    *f2dHistMultiplicityVsVertexZNoTPCOnly;			        //! multiplicity distribution
	//TH2F    *f2dHistMultiplicityVsVertexZNoTPCOnlyNoPileup;			//! multiplicity distribution

   if(! f2dHistMultiplicityVsVertexZBeforeTrigSel) {
      f2dHistMultiplicityVsVertexZBeforeTrigSel = new TH2F("f2dHistMultiplicityVsVertexZBeforeTrigSel", 
         "Tracks per event", 200, 0, 200,400, -20, 20); 		
      fListHistV0->Add(f2dHistMultiplicityVsVertexZBeforeTrigSel);
   }
   if(! f2dHistMultiplicityVsVertexZForTrigEvt) {
      f2dHistMultiplicityVsVertexZForTrigEvt = new TH2F("f2dHistMultiplicityVsVertexZForTrigEvt", 
         "Tracks per event", 200, 0, 200, 400, -20, 20); 		
      fListHistV0->Add(f2dHistMultiplicityVsVertexZForTrigEvt);
   }
   if(! f2dHistMultiplicityVsVertexZ) {
      f2dHistMultiplicityVsVertexZ = new TH2F("f2dHistMultiplicityVsVertexZ", 
         "Tracks per event", 200, 0, 200, 400, -20, 20); 		
      fListHistV0->Add(f2dHistMultiplicityVsVertexZ);
   }
   if(! f2dHistMultiplicityVsVertexZNoTPCOnly) {
      f2dHistMultiplicityVsVertexZNoTPCOnly = new TH2F("f2dHistMultiplicityVsVertexZNoTPCOnly", 
         "Tracks per event", 200, 0, 200, 400, -20, 20); 		
      fListHistV0->Add(f2dHistMultiplicityVsVertexZNoTPCOnly);
   }
   if(! f2dHistMultiplicityVsVertexZNoTPCOnlyNoPileup) {
      f2dHistMultiplicityVsVertexZNoTPCOnlyNoPileup = new TH2F("f2dHistMultiplicityVsVertexZNoTPCOnlyNoPileup", 
         "Tracks per event", 200, 0, 200, 400, -20, 20); 		
      fListHistV0->Add(f2dHistMultiplicityVsVertexZNoTPCOnlyNoPileup);
   }

//-----------------------------------------------------

   if(! fHistPVx) {
      fHistPVx = new TH1F("fHistPVx", 
         "PV x position;Nbr of Evts;x", 
         2000, -0.5, 0.5); 		
      fListHistV0->Add(fHistPVx);
   }
   if(! fHistPVy) {
      fHistPVy = new TH1F("fHistPVy", 
         "PV y position;Nbr of Evts;y", 
         2000, -0.5, 0.5); 		
      fListHistV0->Add(fHistPVy);
   }
   if(! fHistPVz) {
      fHistPVz = new TH1F("fHistPVz", 
         "PV z position;Nbr of Evts;z", 
         400, -20, 20); 		
      fListHistV0->Add(fHistPVz);
   }
   if(! fHistPVxAnalysis) {
      fHistPVxAnalysis = new TH1F("fHistPVxAnalysis", 
         "PV x position;Nbr of Evts;x", 
         2000, -0.5, 0.5);		
      fListHistV0->Add(fHistPVxAnalysis);
   }
   if(! fHistPVyAnalysis) {
      fHistPVyAnalysis = new TH1F("fHistPVyAnalysis", 
         "PV y position;Nbr of Evts;y", 
         2000, -0.5, 0.5);		
      fListHistV0->Add(fHistPVyAnalysis);
   }
   if(! fHistPVzAnalysis) {
      fHistPVzAnalysis = new TH1F("fHistPVzAnalysis", 
         "PV z position;Nbr of Evts;z", 
         400, -20, 20); 		
      fListHistV0->Add(fHistPVzAnalysis);
   }
   if(! fHistSwappedV0Counter) {
      fHistSwappedV0Counter = new TH1F("fHistSwappedV0Counter", 
         "Swap or not histo;Swapped (1) or not (0); count", 
         2, 0, 2); 		
      fListHistV0->Add(fHistSwappedV0Counter);
   }
   //Regular output: Histograms
   PostData(1, fListHistV0);
   //TTree Object: Saved to base directory. Should cache to disk while saving. 
   //(Important to avoid excessive memory usage, particularly when merging)
   PostData(2, fTree);

}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskExtractV0::UserExec(Option_t *) 
{
   // Main loop
   // Called for each event
   //gObjectTable->Print();
   AliESDEvent *lESDevent = 0x0;

   //AliAODEvent *lAODevent = 0x0;
   Int_t    nV0s                        = -1;

   Double_t lTrkgPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
   Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
   Double_t lMagneticField                 = -10.;
	
   // Connect to the InputEvent	
   // After these lines, we should have an ESD/AOD event + the number of cascades in it.

   // Appropriate for ESD analysis! 
		
   lESDevent = dynamic_cast<AliESDEvent*>( InputEvent() );
   if (!lESDevent) {
      AliWarning("ERROR: lESDevent not available \n");
      return;
   }
  
  //------------------------------------------------
  // Rerun V0 vertexer, if asked for
  // --- WARNING: Be careful when using in PbPb
  //------------------------------------------------
  if( fkRunV0Vertexer ){
    lESDevent->ResetV0s();
    AliV0vertexer lV0vtxer;
    lV0vtxer.SetDefaultCuts(fV0Sels);
    lV0vtxer.Tracks2V0vertices(lESDevent);
  }
  
   fTreeVariableRunNumber = lESDevent->GetRunNumber();
   fTreeVariableEventNumber =  
    ( ( ((ULong64_t)lESDevent->GetPeriodNumber() ) << 36 ) |
      ( ((ULong64_t)lESDevent->GetOrbitNumber () ) << 12 ) |
        ((ULong64_t)lESDevent->GetBunchCrossNumber() )  );

   //REVISED multiplicity estimator after 'multiplicity day' (2011)
   Int_t lMultiplicity = -100;
   Int_t lMultiplicityV0A = -100;
   Int_t lMultiplicityZNA = -100;
   Int_t lMultiplicityTRK = -100;
   Int_t lMultiplicitySPD = -100;  

   if(fkIsNuclear == kFALSE) lMultiplicity = fESDtrackCuts->GetReferenceMultiplicity(lESDevent, AliESDtrackCuts::kTrackletsITSTPC,0.5);

   //---> If this is a nuclear collision, then go nuclear on "multiplicity" variable...
   //---> Warning: Experimental
   if(fkIsNuclear == kTRUE){ 
      AliCentrality* centrality;
      centrality = lESDevent->GetCentrality();
      lMultiplicity = ( ( Int_t ) ( centrality->GetCentralityPercentile(   fCentralityEstimator.Data() ) ) );
      lMultiplicityV0A = ( ( Int_t ) ( centrality->GetCentralityPercentile(   "V0A" ) ) );
      lMultiplicityZNA = ( ( Int_t ) ( centrality->GetCentralityPercentile(   "ZNA" ) ) );
      lMultiplicityTRK = ( ( Int_t ) ( centrality->GetCentralityPercentile(   "TRK" ) ) );
      lMultiplicitySPD = ( ( Int_t ) ( centrality->GetCentralityPercentile(   "CL1" ) ) );
      if (centrality->GetQuality()>1) {
        PostData(1, fListHistV0);
        PostData(2, fTree);
        return;
      }
   }

   //Set variable for filling tree afterwards!
   //---> pp case......: GetReferenceMultiplicity
   //---> Pb-Pb case...: Centrality by V0M
  fTreeVariableMultiplicity = lMultiplicity;
  fTreeVariableMultiplicityV0A = lMultiplicityV0A;
  fTreeVariableMultiplicityZNA = lMultiplicityZNA;
  fTreeVariableMultiplicityTRK = lMultiplicityTRK;
  fTreeVariableMultiplicitySPD = lMultiplicitySPD;
  
  fHistMultiplicityBeforeTrigSel->Fill ( lMultiplicity );
  fHistMultiplicityV0ABeforeTrigSel->Fill ( lMultiplicityV0A );
  fHistMultiplicityZNABeforeTrigSel->Fill ( lMultiplicityZNA );
  fHistMultiplicityTRKBeforeTrigSel->Fill ( lMultiplicityTRK );
  fHistMultiplicitySPDBeforeTrigSel->Fill ( lMultiplicitySPD );
  
  fHistV0MultiplicityBeforeTrigSel->Fill ( lESDevent->GetNumberOfV0s() );

//------------------------------------------------
// Physics Selection
//------------------------------------------------
  
  // new method
  UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  Bool_t isSelected = 0;
  Bool_t isSelectedExtra = kTRUE; //extra sel, default YES
  isSelected = (maskIsSelected & AliVEvent::kMB) == AliVEvent::kMB;
  
  //pA triggering: CINT7
  if( fkSwitchINT7 ) isSelected = (maskIsSelected & AliVEvent::kINT7) == AliVEvent::kINT7;
  
  //Extra selection applies if with/without SDD is to be dealth with
  if( fkFastOnly == "kFastOnly"){
    //If not kFastOnly, isSelectedExtra will be kFALSE; procedure will reject it
    isSelectedExtra = (maskIsSelected & AliVEvent::kFastOnly) == AliVEvent::kFastOnly;
  }
  if( fkFastOnly == "NotkFastOnly"){
    //If not kFastOnly, isSelectedExtra will be kTRUE; procedure will accept it
    isSelectedExtra = !( (maskIsSelected & AliVEvent::kFastOnly) == AliVEvent::kFastOnly );
  }
  
  //Standard Min-Bias Selection
  if ( !isSelected ) {
    PostData(1, fListHistV0);
    PostData(2, fTree);
    return;
  }
  //Check if goes through extra selections
  //isSelectedExtra will be true in case -> fkFastOnly==""
  //isSelectedExtra will be true in case -> fkFastOnly=="kFastOnly"    && bit kFastOnly ON
  //isSelectedExtra will be true in case -> fkFastOnly=="NotkFastOnly" && bit kFastOnly OFF
  if ( !isSelectedExtra ) {
    PostData(1, fListHistV0);
    PostData(2, fTree);
    return;
  }
  
  
//------------------------------------------------
// After Trigger Selection
//------------------------------------------------

	nV0s = lESDevent->GetNumberOfV0s();

  fHistV0MultiplicityForTrigEvt     ->Fill( nV0s       );
  fHistMultiplicityForTrigEvt       ->Fill( lMultiplicity  );
  fHistMultiplicityV0AForTrigEvt       ->Fill( lMultiplicityV0A  );
  fHistMultiplicityZNAForTrigEvt       ->Fill( lMultiplicityZNA  );
  fHistMultiplicityTRKForTrigEvt       ->Fill( lMultiplicityTRK  );
  fHistMultiplicitySPDForTrigEvt       ->Fill( lMultiplicitySPD  );
  
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

   Double_t tPrimaryVtxPosition[3];
   const AliVVertex *primaryVtx = lESDevent->GetPrimaryVertex();
   tPrimaryVtxPosition[0] = primaryVtx->GetX();
   tPrimaryVtxPosition[1] = primaryVtx->GetY();
   tPrimaryVtxPosition[2] = primaryVtx->GetZ();
   fHistPVx->Fill( tPrimaryVtxPosition[0] );
   fHistPVy->Fill( tPrimaryVtxPosition[1] );
   fHistPVz->Fill( tPrimaryVtxPosition[2] );

   f2dHistMultiplicityVsVertexZForTrigEvt->Fill( lMultiplicity, tPrimaryVtxPosition[2] );

//------------------------------------------------
// Primary Vertex Z position: SKIP
//------------------------------------------------

  if(TMath::Abs(lBestPrimaryVtxPos[2]) > 10.0 && fkpAVertexSelection==kFALSE) {
    AliWarning("Pb / | Z position of Best Prim Vtx | > 10.0 cm ... return !");
    PostData(1, fListHistV0);
    PostData(2, fTree);
    return;
  }
  if(fkpAVertexSelection==kTRUE && fUtils->IsFirstEventInChunk(lESDevent)) {
    AliWarning("Pb / | This is the first event in the chunk!");
    PostData(1, fListHistV0);
    PostData(2, fTree);
    return;
  }
  if(fkpAVertexSelection==kTRUE && !fUtils->IsVertexSelected2013pA(lESDevent)) {
    AliWarning("Pb / | Vertex not selected by 2013 pA criteria!");
    PostData(1, fListHistV0);
    PostData(2, fTree);
    return;
  }
  
  f2dHistMultiplicityVsVertexZ->Fill( lMultiplicity, tPrimaryVtxPosition[2] );
  
  lMagneticField = lESDevent->GetMagneticField( );
  fHistV0MultiplicityForSelEvt ->Fill( nV0s );
  fHistMultiplicity->Fill(lMultiplicity);
  fHistMultiplicityV0A->Fill(lMultiplicityV0A);
  fHistMultiplicityZNA->Fill(lMultiplicityZNA);
  fHistMultiplicityTRK->Fill(lMultiplicityTRK);
  fHistMultiplicitySPD->Fill(lMultiplicitySPD);

//------------------------------------------------
// Only look at events with well-established PV
//------------------------------------------------
	
   const AliESDVertex *lPrimaryTrackingESDVtxCheck = lESDevent->GetPrimaryVertexTracks();
   const AliESDVertex *lPrimarySPDVtx = lESDevent->GetPrimaryVertexSPD();
   if (!lPrimarySPDVtx->GetStatus() && !lPrimaryTrackingESDVtxCheck->GetStatus() && fkpAVertexSelection==kFALSE ){
      AliWarning("Pb / No SPD prim. vertex nor prim. Tracking vertex ... return !");
      PostData(1, fListHistV0);
      PostData(2, fTree);
      return;
   }


   f2dHistMultiplicityVsVertexZNoTPCOnly->Fill( lMultiplicity, tPrimaryVtxPosition[2] );
  fHistV0MultiplicityForSelEvtNoTPCOnly ->Fill( nV0s );
  fHistMultiplicityNoTPCOnly->Fill(lMultiplicity);
  fHistMultiplicityV0ANoTPCOnly->Fill(lMultiplicityV0A);
  fHistMultiplicityZNANoTPCOnly->Fill(lMultiplicityZNA);
  fHistMultiplicityTRKNoTPCOnly->Fill(lMultiplicityTRK);
  fHistMultiplicitySPDNoTPCOnly->Fill(lMultiplicitySPD);

//------------------------------------------------
// Pileup Rejection
//------------------------------------------------

   // FIXME : quality selection regarding pile-up rejection 
   if(lESDevent->IsPileupFromSPD() && !fkIsNuclear && fkRejectPileup){// minContributors=3, minZdist=0.8, nSigmaZdist=3., nSigmaDiamXY=2., nSigmaDiamZ=5.  -> see http://alisoft.cern.ch/viewvc/trunk/STEER/AliESDEvent.h?root=AliRoot&r1=41914&r2=42199&pathrev=42199
      PostData(1, fListHistV0);
      PostData(2, fTree); 
      return;
   }

  f2dHistMultiplicityVsVertexZNoTPCOnlyNoPileup->Fill( lMultiplicity, tPrimaryVtxPosition[2] );
  fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup ->Fill( nV0s );
  fHistMultiplicityNoTPCOnlyNoPileup->Fill(lMultiplicity);
  fHistMultiplicityV0ANoTPCOnlyNoPileup->Fill(lMultiplicityV0A);
  fHistMultiplicityZNANoTPCOnlyNoPileup->Fill(lMultiplicityZNA);
  fHistMultiplicityTRKNoTPCOnlyNoPileup->Fill(lMultiplicityTRK);
  fHistMultiplicitySPDNoTPCOnlyNoPileup->Fill(lMultiplicitySPD);
  
//------------------------------------------------
// MAIN LAMBDA LOOP STARTS HERE
//------------------------------------------------

   if( ! (lESDevent->GetPrimaryVertex()->IsFromVertexerZ() )	 ){ 
      fHistPVxAnalysis->Fill( tPrimaryVtxPosition[0] );
      fHistPVyAnalysis->Fill( tPrimaryVtxPosition[1] );
      fHistPVzAnalysis->Fill( tPrimaryVtxPosition[2] );
   }

  fTreeVariablePVx = tPrimaryVtxPosition[0];
  fTreeVariablePVy = tPrimaryVtxPosition[1];
  fTreeVariablePVz = tPrimaryVtxPosition[2];

//Variable definition
   Int_t    lOnFlyStatus = 0;// nv0sOn = 0, nv0sOff = 0;
   Double_t lChi2V0 = 0;
   Double_t lDcaV0Daughters = 0, lDcaV0ToPrimVertex = 0;
   Double_t lDcaPosToPrimVertex = 0, lDcaNegToPrimVertex = 0;
   Double_t lV0CosineOfPointingAngle = 0;
   Double_t lV0Radius = 0, lPt = 0;
   Double_t lRapK0Short = 0, lRapLambda = 0;
   Double_t lInvMassK0s = 0, lInvMassLambda = 0, lInvMassAntiLambda = 0;
   Double_t lAlphaV0 = 0, lPtArmV0 = 0;

   Double_t fMinV0Pt = 0; 
   Double_t fMaxV0Pt = 100; 

   Int_t nv0s = 0;
   nv0s = lESDevent->GetNumberOfV0s();

   //for (Int_t iV0 = 0; iV0 < nv0s; iV0++) 
   for (Int_t iV0 = 0; iV0 < nv0s; iV0++) //extra-crazy test
   {// This is the begining of the V0 loop
      AliESDv0 *v0 = ((AliESDEvent*)lESDevent)->GetV0(iV0);
      if (!v0) continue;

      //---> Fix On-the-Fly candidates, count how many swapped
      if( v0->GetParamN()->Charge() > 0 && v0->GetParamP()->Charge() < 0 ){
        fHistSwappedV0Counter -> Fill( 1 );
      }else{
        fHistSwappedV0Counter -> Fill( 0 ); 
      }
      if ( fkUseOnTheFly ) CheckChargeV0(v0); 

      Double_t tDecayVertexV0[3]; v0->GetXYZ(tDecayVertexV0[0],tDecayVertexV0[1],tDecayVertexV0[2]); 

      Double_t tV0mom[3];
      v0->GetPxPyPz( tV0mom[0],tV0mom[1],tV0mom[2] ); 
      Double_t lV0TotalMomentum = TMath::Sqrt(
      tV0mom[0]*tV0mom[0]+tV0mom[1]*tV0mom[1]+tV0mom[2]*tV0mom[2] );

      lV0Radius = TMath::Sqrt(tDecayVertexV0[0]*tDecayVertexV0[0]+tDecayVertexV0[1]*tDecayVertexV0[1]);

      //Set Variables for later filling
      fTreeVariableV0x = tDecayVertexV0[0];
      fTreeVariableV0y = tDecayVertexV0[1];
      fTreeVariableV0z = tDecayVertexV0[2];

      //Set Variables for later filling
      fTreeVariableV0Px = tV0mom[0];
      fTreeVariableV0Py = tV0mom[1];
      fTreeVariableV0Pz = tV0mom[2];

      lPt = v0->Pt();
      lRapK0Short = v0->RapK0Short();
      lRapLambda  = v0->RapLambda();
      if ((lPt<fMinV0Pt)||(fMaxV0Pt<lPt)) continue;

      UInt_t lKeyPos = (UInt_t)TMath::Abs(v0->GetPindex());
      UInt_t lKeyNeg = (UInt_t)TMath::Abs(v0->GetNindex());

      Double_t lMomPos[3]; v0->GetPPxPyPz(lMomPos[0],lMomPos[1],lMomPos[2]);
      Double_t lMomNeg[3]; v0->GetNPxPyPz(lMomNeg[0],lMomNeg[1],lMomNeg[2]);

      AliESDtrack *pTrack=((AliESDEvent*)lESDevent)->GetTrack(lKeyPos);
      AliESDtrack *nTrack=((AliESDEvent*)lESDevent)->GetTrack(lKeyNeg);
      if (!pTrack || !nTrack) {
         Printf("ERROR: Could not retreive one of the daughter track");
         continue;
      }

      //Daughter Eta for Eta selection, afterwards
      fTreeVariableNegEta = nTrack->Eta();
      fTreeVariablePosEta = pTrack->Eta();
     
     if( fkSpecialExecution ){
       fTreeVariableNegPx = lMomNeg[0]; fTreeVariableNegPy = lMomNeg[1]; fTreeVariableNegPz = lMomNeg[2];
       fTreeVariablePosPx = lMomPos[0]; fTreeVariablePosPy = lMomPos[1]; fTreeVariablePosPz = lMomPos[2];
       //Need to acquire info to do dEdx cut afterwards...
       fTreeVariablePosTPCSignal = pTrack->GetTPCsignal();
       fTreeVariableNegTPCSignal = nTrack->GetTPCsignal();
       //Get Inner momentum if possible
       fTreeVariableNegInnerP = nTrack->GetP();
       fTreeVariablePosInnerP = pTrack->GetP();
       const AliExternalTrackParam *innerpos=pTrack->GetInnerParam();
       const AliExternalTrackParam *innerneg=nTrack->GetInnerParam();
       if(innerpos) { fTreeVariablePosInnerP = innerpos->GetP(); }
       if(innerneg) { fTreeVariableNegInnerP = innerneg->GetP(); }
     }

      // Filter like-sign V0 (next: add counter and distribution)
      if ( pTrack->GetSign() == nTrack->GetSign()){
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
     
      if ( ( ( ( pTrack->GetTPCClusterInfo(2,1) ) < 70 ) || ( ( nTrack->GetTPCClusterInfo(2,1) ) < 70 ) )&&(fkTakeAllTracks==kFALSE) ) continue;
	
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
      if ( (fTreeVariableLeastRatioCrossedRowsOverFindable < 0.8)&&(fkTakeAllTracks==kFALSE) ) continue;

      //End track Quality Cuts
      //________________________________________________________________________

      lDcaPosToPrimVertex = TMath::Abs(pTrack->GetD(tPrimaryVtxPosition[0],
							tPrimaryVtxPosition[1],
							lMagneticField) );

      lDcaNegToPrimVertex = TMath::Abs(nTrack->GetD(tPrimaryVtxPosition[0],
							tPrimaryVtxPosition[1],
							lMagneticField) );

      lOnFlyStatus = v0->GetOnFlyStatus();
      lChi2V0 = v0->GetChi2V0();
      lDcaV0Daughters = v0->GetDcaV0Daughters();
      lDcaV0ToPrimVertex = v0->GetD(tPrimaryVtxPosition[0],tPrimaryVtxPosition[1],tPrimaryVtxPosition[2]);
      lV0CosineOfPointingAngle = v0->GetV0CosineOfPointingAngle(tPrimaryVtxPosition[0],tPrimaryVtxPosition[1],tPrimaryVtxPosition[2]);
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

//------------------------------------------------
// Fill Tree! 
//------------------------------------------------
     
     // The conditionals are meant to decrease excessive
     // memory usage!
     
     //First Selection: Reject OnFly
     if( (lOnFlyStatus == 0 && fkUseOnTheFly == kFALSE) || (lOnFlyStatus != 0 && fkUseOnTheFly == kTRUE ) ){
       //Second Selection: rough 20-sigma band, parametric.
       //K0Short: Enough to parametrize peak broadening with linear function.
       Double_t lUpperLimitK0Short = (5.63707e-01) + (1.14979e-02)*fTreeVariablePt;
       Double_t lLowerLimitK0Short = (4.30006e-01) - (1.10029e-02)*fTreeVariablePt;
       //Lambda: Linear (for higher pt) plus exponential (for low-pt broadening)
       //[0]+[1]*x+[2]*TMath::Exp(-[3]*x)
       Double_t lUpperLimitLambda = (1.13688e+00) + (5.27838e-03)*fTreeVariablePt + (8.42220e-02)*TMath::Exp(-(3.80595e+00)*fTreeVariablePt);
       Double_t lLowerLimitLambda = (1.09501e+00) - (5.23272e-03)*fTreeVariablePt - (7.52690e-02)*TMath::Exp(-(3.46339e+00)*fTreeVariablePt);
       //Do Selection
       if( (fTreeVariableInvMassLambda     < lUpperLimitLambda  && fTreeVariableInvMassLambda     > lLowerLimitLambda     ) ||
          (fTreeVariableInvMassAntiLambda < lUpperLimitLambda  && fTreeVariableInvMassAntiLambda > lLowerLimitLambda     ) ||
          (fTreeVariableInvMassK0s        < lUpperLimitK0Short && fTreeVariableInvMassK0s        > lLowerLimitK0Short    ) ){
         //Pre-selection in case this is AA...
         if( fkIsNuclear == kFALSE ) fTree->Fill();
         if( fkIsNuclear == kTRUE){
           //If this is a nuclear collision___________________
           // ... pre-filter with TPC, daughter eta selection
           if( (fTreeVariableInvMassLambda     < lUpperLimitLambda  && fTreeVariableInvMassLambda     > lLowerLimitLambda
                && TMath::Abs(fTreeVariableNSigmasPosProton) < 6.0 && TMath::Abs(fTreeVariableNSigmasNegPion) < 6.0 ) ||
              (fTreeVariableInvMassAntiLambda < lUpperLimitLambda  && fTreeVariableInvMassAntiLambda > lLowerLimitLambda
               && TMath::Abs(fTreeVariableNSigmasNegProton) < 6.0 && TMath::Abs(fTreeVariableNSigmasPosPion) < 6.0 ) ||
              (fTreeVariableInvMassK0s        < lUpperLimitK0Short && fTreeVariableInvMassK0s        > lLowerLimitK0Short
               && TMath::Abs(fTreeVariableNSigmasNegPion)   < 6.0 && TMath::Abs(fTreeVariableNSigmasPosPion) < 6.0 ) ){
                //insane test
                if ( TMath::Abs(fTreeVariableNegEta)<0.8 && TMath::Abs(fTreeVariablePosEta)<0.8 && !fkSpecialExecution) fTree->Fill();
              }
         }//end nuclear_____________________________________
       }
     }
     
     if(lOnFlyStatus == 0 && fkSpecialExecution){
       if(
          (fTreeVariableNSigmasPosProton > 6 && TMath::Abs(fTreeVariableNSigmasNegPion)< 6) ||
          (fTreeVariableNSigmasNegProton > 6 && TMath::Abs(fTreeVariableNSigmasPosPion)< 6) ){
         fTree->Fill();
       }
     }

//------------------------------------------------
// Fill tree over.
//------------------------------------------------

   }// This is the end of the V0 loop
  
  // Post output data.
    PostData(1, fListHistV0);
    PostData(2, fTree);


}

//________________________________________________________________________
void AliAnalysisTaskExtractV0::Terminate(Option_t *)
{
   // Draw result to the screen
   // Called once at the end of the query
   // This will draw the V0 candidate multiplicity, whose 
   // number of entries corresponds to the number of triggered events. 
   TList *cRetrievedList = 0x0;
   cRetrievedList = (TList*)GetOutputData(1);
   if(!cRetrievedList){
      Printf("ERROR - AliAnalysisTaskExtractV0 : ouput data container list not available\n");
      return;
   }		
   fHistV0MultiplicityForTrigEvt = dynamic_cast<TH1F*> (  cRetrievedList->FindObject("fHistV0MultiplicityForTrigEvt")  );
   if (!fHistV0MultiplicityForTrigEvt) {
      Printf("ERROR - AliAnalysisTaskExtractV0 : fHistV0MultiplicityForTrigEvt not available");
      return;
   }
   TCanvas *canCheck = new TCanvas("AliAnalysisTaskExtractV0","V0 Multiplicity",10,10,510,510);
   canCheck->cd(1)->SetLogy();
   fHistV0MultiplicityForTrigEvt->SetMarkerStyle(22);
   fHistV0MultiplicityForTrigEvt->DrawCopy("E");
}

//________________________________________________________________________
void AliAnalysisTaskExtractV0::CheckChargeV0(AliESDv0 *v0)
{
   // This function checks charge of negative and positive daughter tracks. 
   // If incorrectly defined (onfly vertexer), swaps out. 
   if( v0->GetParamN()->Charge() > 0 && v0->GetParamP()->Charge() < 0 ){
      //V0 daughter track swapping is required! Note: everything is swapped here... P->N, N->P
      Long_t lCorrectNidx = v0->GetPindex();
      Long_t lCorrectPidx = v0->GetNindex();
      Double32_t	lCorrectNmom[3];
      Double32_t	lCorrectPmom[3];
      v0->GetPPxPyPz( lCorrectNmom[0], lCorrectNmom[1], lCorrectNmom[2] );
      v0->GetNPxPyPz( lCorrectPmom[0], lCorrectPmom[1], lCorrectPmom[2] );

      AliExternalTrackParam	lCorrectParamN(
        v0->GetParamP()->GetX() , 
        v0->GetParamP()->GetAlpha() , 
        v0->GetParamP()->GetParameter() , 
        v0->GetParamP()->GetCovariance() 
      );
      AliExternalTrackParam	lCorrectParamP(
        v0->GetParamN()->GetX() , 
        v0->GetParamN()->GetAlpha() , 
        v0->GetParamN()->GetParameter() , 
        v0->GetParamN()->GetCovariance() 
      );
      lCorrectParamN.SetMostProbablePt( v0->GetParamP()->GetMostProbablePt() );
      lCorrectParamP.SetMostProbablePt( v0->GetParamN()->GetMostProbablePt() );

      //Get Variables___________________________________________________
      Double_t lDcaV0Daughters = v0 -> GetDcaV0Daughters();
      Double_t lCosPALocal     = v0 -> GetV0CosineOfPointingAngle(); 
      Bool_t lOnFlyStatusLocal = v0 -> GetOnFlyStatus();

      //Create Replacement Object_______________________________________
      AliESDv0 *v0correct = new AliESDv0(lCorrectParamN,lCorrectNidx,lCorrectParamP,lCorrectPidx);
      v0correct->SetDcaV0Daughters          ( lDcaV0Daughters   );
      v0correct->SetV0CosineOfPointingAngle ( lCosPALocal       );
      v0correct->ChangeMassHypothesis       ( kK0Short          );
      v0correct->SetOnFlyStatus             ( lOnFlyStatusLocal );

      //Reverse Cluster info..._________________________________________
      v0correct->SetClusters( v0->GetClusters( 1 ), v0->GetClusters ( 0 ) );

      *v0 = *v0correct;
      //Proper cleanup..._______________________________________________
      v0correct->Delete();
      v0correct = 0x0;

      //Just another cross-check and output_____________________________
      if( v0->GetParamN()->Charge() > 0 && v0->GetParamP()->Charge() < 0 ) {
        AliWarning("Found Swapped Charges, tried to correct but something FAILED!");
      }else{
        //AliWarning("Found Swapped Charges and fixed.");
      }
      //________________________________________________________________
   }else{
      //Don't touch it! ---
      //Printf("Ah, nice. Charges are already ordered...");
   }
   return;
} 
