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
#include "AliMultSelection.h"

#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDUtils.h"
#include "AliESDHeader.h"
#include "AliAODTrack.h"
#include "AliAnalysisTaskExtractV0AODRun2.h"

//debugging purposes
#include "TObjectTable.h"

ClassImp(AliAnalysisTaskExtractV0AODRun2)

AliAnalysisTaskExtractV0AODRun2::AliAnalysisTaskExtractV0AODRun2() 
  : AliAnalysisTaskSE(), fListHistV0(0), fTree(0), fPIDResponse(0),
   fkIsNuclear   ( kFALSE ), 
   fkLowEnergyPP ( kFALSE ),
   fkUseOnTheFly ( kFALSE ),
   fTriggerMask  ( "kMB"  ),
    fkLowE ( kFALSE ),
    fkSaveAllInvMasses (kFALSE),
    fkPreSelect (kFALSE),
//------------------------------------------------
// Variables for tree 
//------------------------------------------------

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
    fTreeVariableAlphaV0BoostAsK0(0),
	 fTreeVariablePtArmV0(0),
	 fTreeVariableNegEta(0),
	 fTreeVariablePosEta(0),

	 fTreeVariableNSigmasPosProton(0),
	 fTreeVariableNSigmasPosPion(0),
	 fTreeVariableNSigmasNegProton(0),
	 fTreeVariableNSigmasNegPion(0),

    fTreeVariableNegTPCSignal(0),
    fTreeVariablePosTPCSignal(0),
    fTreeVariableNegInnerP(0),
    fTreeVariablePosInnerP(0),

	 fTreeVariableDistOverTotMom(0),
   fTreeVariableLeastNbrCrossedRows(0),
	 fTreeVariableLeastRatioCrossedRowsOverFindable(0),
	 fTreeVariableCentrality(0),

   fTreeVariableRunNumber(0),
   fTreeVariablePeriodNumber(0),
    fTreeVariableOrbitNumber(0),
    fTreeVariableBunchCrossNumber(0),

      fTreeVariableV0Px(0),
      fTreeVariableV0Py(0),
      fTreeVariableV0Pz(0),

    fTreeVariableDecayMomLambda(0),
    fTreeVariableDecayMomAntiLambda(0),
    fTreeVariableDecayMomK0Short(0),

    fTreeVariablePrimVX(0),
    fTreeVariablePrimVY(0),
    fTreeVariablePrimVZ(0),

    fTreeVariableCowboy(0),


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
    fHistCentralityBeforeTrigSel(0),
    fHistCentralityForTrigEvt(0),
    fHistCentrality(0),
    fHistCentralityNoTPCOnly(0),
    fHistCentralityNoTPCOnlyNoPileup(0),
   fHistPVx(0),
   fHistPVy(0),
   fHistPVz(0),
   fHistPVxAnalysis(0),
   fHistPVyAnalysis(0),
   fHistPVzAnalysis(0),
   fHistSwappedV0Counter(0),

    fHistRunNumber(0),
    fHistPeriodNumber(0),
    fHistOrbitNumber(0),
    fHistBunchCrossNumber(0),

f2dHistTPCSignal(0),
f2dHistTPCSignalPionsFromK0(0),
f2dHistTPCSignalPionsFromLambda(0),
f2dHistTPCSignalProtonsFromLambda(0),

    f3dHistInvMassVsPtVsCentLambda(0),
    f3dHistInvMassVsPtVsCentAntiLambda(0),
    f3dHistInvMassVsPtVsCentK0Short(0),

    f3dHistInvMassVsPtVsCentLambdaCowboy(0),
    f3dHistInvMassVsPtVsCentAntiLambdaCowboy(0),
    f3dHistInvMassVsPtVsCentK0ShortCowboy(0),

    f3dHistInvMassVsPtVsCentLambdaSailor(0),
    f3dHistInvMassVsPtVsCentAntiLambdaSailor(0),
    f3dHistInvMassVsPtVsCentK0ShortSailor(0)
{
  // Dummy Constructor
}

AliAnalysisTaskExtractV0AODRun2::AliAnalysisTaskExtractV0AODRun2(const char *name) 
  : AliAnalysisTaskSE(name), fListHistV0(0), fTree(0), fPIDResponse(0),
   fkIsNuclear   ( kFALSE ), 
   fkLowEnergyPP ( kFALSE ),
   fkUseOnTheFly ( kFALSE ),
   fTriggerMask  ( "kMB"  ),
    fkLowE ( kFALSE ),
    fkSaveAllInvMasses (kFALSE),
    fkPreSelect (kFALSE),
//------------------------------------------------
// Variables for tree 
//------------------------------------------------

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
    fTreeVariableAlphaV0BoostAsK0(0),
	 fTreeVariablePtArmV0(0),
	 fTreeVariableNegEta(0),
	 fTreeVariablePosEta(0),

	 fTreeVariableNSigmasPosProton(0),
	 fTreeVariableNSigmasPosPion(0),
	 fTreeVariableNSigmasNegProton(0),
	 fTreeVariableNSigmasNegPion(0),

    fTreeVariableNegTPCSignal(0),
    fTreeVariablePosTPCSignal(0),
    fTreeVariableNegInnerP(0),
    fTreeVariablePosInnerP(0),

	 fTreeVariableDistOverTotMom(0),
   fTreeVariableLeastNbrCrossedRows(0),
	 fTreeVariableLeastRatioCrossedRowsOverFindable(0),
	 fTreeVariableCentrality(0),

   fTreeVariableRunNumber(0),
   fTreeVariablePeriodNumber(0),
    fTreeVariableOrbitNumber(0),
    fTreeVariableBunchCrossNumber(0),

      fTreeVariableV0Px(0),
      fTreeVariableV0Py(0),
      fTreeVariableV0Pz(0),

    fTreeVariableDecayMomLambda(0),
    fTreeVariableDecayMomAntiLambda(0),
    fTreeVariableDecayMomK0Short(0),

    fTreeVariablePrimVX(0),
    fTreeVariablePrimVY(0),
    fTreeVariablePrimVZ(0),

    fTreeVariableCowboy(0),


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
    fHistCentralityBeforeTrigSel(0),
    fHistCentralityForTrigEvt(0),
    fHistCentrality(0),
    fHistCentralityNoTPCOnly(0),
    fHistCentralityNoTPCOnlyNoPileup(0),
   fHistPVx(0),
   fHistPVy(0),
   fHistPVz(0),
   fHistPVxAnalysis(0),
   fHistPVyAnalysis(0),
   fHistPVzAnalysis(0),
   fHistSwappedV0Counter(0),

    fHistRunNumber(0),
    fHistPeriodNumber(0),
    fHistOrbitNumber(0),
    fHistBunchCrossNumber(0),

f2dHistTPCSignal(0),
f2dHistTPCSignalPionsFromK0(0),
f2dHistTPCSignalPionsFromLambda(0),
f2dHistTPCSignalProtonsFromLambda(0),

    f3dHistInvMassVsPtVsCentLambda(0),
    f3dHistInvMassVsPtVsCentAntiLambda(0),
    f3dHistInvMassVsPtVsCentK0Short(0),

    f3dHistInvMassVsPtVsCentLambdaCowboy(0),
    f3dHistInvMassVsPtVsCentAntiLambdaCowboy(0),
    f3dHistInvMassVsPtVsCentK0ShortCowboy(0),

    f3dHistInvMassVsPtVsCentLambdaSailor(0),
    f3dHistInvMassVsPtVsCentAntiLambdaSailor(0),
    f3dHistInvMassVsPtVsCentK0ShortSailor(0)
{
  // Constructor
  // Output slot #0 writes into a TList container (Lambda Histos and fTree)
   DefineOutput(1, TList::Class());
   DefineOutput(2, TTree::Class());
}


AliAnalysisTaskExtractV0AODRun2::~AliAnalysisTaskExtractV0AODRun2()
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
}



//________________________________________________________________________
void AliAnalysisTaskExtractV0AODRun2::UserCreateOutputObjects()
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
    fTree->Branch("fTreeVariableAlphaV0BoostAsK0",&fTreeVariableAlphaV0BoostAsK0,"fTreeVariableAlphaV0BoostAsK0/F");
/*14*/	fTree->Branch("fTreeVariablePtArmV0",&fTreeVariablePtArmV0,"fTreeVariablePtArmV0/F");
/*15*/	fTree->Branch("fTreeVariableLeastNbrCrossedRows",&fTreeVariableLeastNbrCrossedRows,"fTreeVariableLeastNbrCrossedRows/I");
/*16*/	fTree->Branch("fTreeVariableLeastRatioCrossedRowsOverFindable",&fTreeVariableLeastRatioCrossedRowsOverFindable,"fTreeVariableLeastRatioCrossedRowsOverFindable/F");
//-----------CENTRALITY-INFO--------------------
/*17*/	fTree->Branch("fTreeVariableCentrality",&fTreeVariableCentrality,"fTreeVariableCentrality/I");
//------------------------------------------------
/*18*/	fTree->Branch("fTreeVariableDistOverTotMom",&fTreeVariableDistOverTotMom,"fTreeVariableDistOverTotMom/F");
    //------------------------------------------------
/*19*/	fTree->Branch("fTreeVariableNSigmasPosProton",&fTreeVariableNSigmasPosProton,"fTreeVariableNSigmasPosProton/F");
/*20*/	fTree->Branch("fTreeVariableNSigmasPosPion",&fTreeVariableNSigmasPosPion,"fTreeVariableNSigmasPosPion/F");
/*21*/	fTree->Branch("fTreeVariableNSigmasNegProton",&fTreeVariableNSigmasNegProton,"fTreeVariableNSigmasNegProton/F");
/*22*/	fTree->Branch("fTreeVariableNSigmasNegPion",&fTreeVariableNSigmasNegPion,"fTreeVariableNSigmasNegPion/F");
    //------------------------------------------------
/*23*/	fTree->Branch("fTreeVariableNegEta",&fTreeVariableNegEta,"fTreeVariableNegEta/F");
/*24*/	fTree->Branch("fTreeVariablePosEta",&fTreeVariablePosEta,"fTreeVariablePosEta/F");
    //------------------------------------------------
/*25*/	fTree->Branch("fTreeVariableRunNumber",&fTreeVariableRunNumber,"fTreeVariableRunNumber/I");
/*26*/	fTree->Branch("fTreeVariablePeriodNumber",&fTreeVariablePeriodNumber,"fTreeVariablePeriodNumber/I");
    fTree->Branch("fTreeVariableOrbitNumber",&fTreeVariableOrbitNumber,"fTreeVariableOrbitNumber/I");
    fTree->Branch("fTreeVariableBunchCrossNumber",&fTreeVariableBunchCrossNumber,"fTreeVariableBunchCrossNumber/I");

    //------------------------------------------------
    fTree->Branch("fTreeVariableNegTPCSignal",&fTreeVariableNegTPCSignal,"fTreeVariableNegTPCSignal/F");
    fTree->Branch("fTreeVariablePosTPCSignal",&fTreeVariablePosTPCSignal,"fTreeVariablePosTPCSignal/F");
    fTree->Branch("fTreeVariableNegInnerP",&fTreeVariableNegInnerP,"fTreeVariableNegInnerP/F");
    fTree->Branch("fTreeVariablePosInnerP",&fTreeVariablePosInnerP,"fTreeVariablePosInnerP/F");
    
        fTree->Branch("fTreeVariableV0Px",&fTreeVariableV0Px,"fTreeVariableV0Px/F");
        fTree->Branch("fTreeVariableV0Py",&fTreeVariableV0Py,"fTreeVariableV0Py/F");
        fTree->Branch("fTreeVariableV0Pz",&fTreeVariableV0Pz,"fTreeVariableV0Pz/F");
    
    fTree->Branch("fTreeVariableDecayMomLambda",&fTreeVariableDecayMomLambda,"fTreeVariableDecayMomLambda/F");
    fTree->Branch("fTreeVariableDecayMomAntiLambda",&fTreeVariableDecayMomAntiLambda,"fTreeVariableDecayMomAntiLambda/F");
    fTree->Branch("fTreeVariableDecayMomK0Short",&fTreeVariableDecayMomK0Short,"fTreeVariableDecayMomK0Short/F");

    fTree->Branch("fTreeVariablePrimVX",&fTreeVariablePrimVX,"fTreeVariablePrimVX/F");
    fTree->Branch("fTreeVariablePrimVY",&fTreeVariablePrimVY,"fTreeVariablePrimVY/F");
    fTree->Branch("fTreeVariablePrimVZ",&fTreeVariablePrimVZ,"fTreeVariablePrimVZ/F");

    fTree->Branch("fTreeVariableCowboy",&fTreeVariableCowboy,"fTreeVariableCowboy/O");
    

//------------------------------------------------
// Particle Identification Setup
//------------------------------------------------

   AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
   AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
   fPIDResponse = inputHandler->GetPIDResponse();

//Skipped. Will use Local setup.

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
         2500, 0, 2500);
      fListHistV0->Add(fHistV0MultiplicityBeforeTrigSel);
   }
           
   if(! fHistV0MultiplicityForTrigEvt) {
      fHistV0MultiplicityForTrigEvt = new TH1F("fHistV0MultiplicityForTrigEvt", 
         "V0s per event (for triggered evt);Nbr of V0s/Evt;Events", 
         2500, 0, 2500);
      fListHistV0->Add(fHistV0MultiplicityForTrigEvt);
   }

   if(! fHistV0MultiplicityForSelEvt) {
      fHistV0MultiplicityForSelEvt = new TH1F("fHistV0MultiplicityForSelEvt", 
         "V0s per event;Nbr of V0s/Evt;Events", 
         2500, 0, 2500);
      fListHistV0->Add(fHistV0MultiplicityForSelEvt);
   }

   if(! fHistV0MultiplicityForSelEvtNoTPCOnly) {
      fHistV0MultiplicityForSelEvtNoTPCOnly = new TH1F("fHistV0MultiplicityForSelEvtNoTPCOnly", 
         "V0s per event;Nbr of V0s/Evt;Events", 
         2500, 0, 2500);
      fListHistV0->Add(fHistV0MultiplicityForSelEvtNoTPCOnly);
   }

   if(! fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup) {
      fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup = new TH1F("fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup", 
         "V0s per event;Nbr of V0s/Evt;Events", 
         2500, 0, 2500);
      fListHistV0->Add(fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup);
   }

//------------------------------------------------
// Track Multiplicity Histograms
//------------------------------------------------

   if(! fHistMultiplicityBeforeTrigSel) {
      fHistMultiplicityBeforeTrigSel = new TH1F("fHistMultiplicityBeforeTrigSel", 
         "Tracks per event;Nbr of Tracks;Events", 
         20000, 0, 20000);
      fListHistV0->Add(fHistMultiplicityBeforeTrigSel);
   }
   if(! fHistMultiplicityForTrigEvt) {
      fHistMultiplicityForTrigEvt = new TH1F("fHistMultiplicityForTrigEvt", 
         "Tracks per event;Nbr of Tracks;Events", 
         20000, 0, 20000);
      fListHistV0->Add(fHistMultiplicityForTrigEvt);
   }
   if(! fHistMultiplicity) {
      fHistMultiplicity = new TH1F("fHistMultiplicity", 
         "Tracks per event;Nbr of Tracks;Events", 
         20000, 0, 20000);
      fListHistV0->Add(fHistMultiplicity);
   }
   if(! fHistMultiplicityNoTPCOnly) {
      fHistMultiplicityNoTPCOnly = new TH1F("fHistMultiplicityNoTPCOnly", 
         "Tracks per event;Nbr of Tracks;Events", 
         20000, 0, 20000);
      fListHistV0->Add(fHistMultiplicityNoTPCOnly);
   }
   if(! fHistMultiplicityNoTPCOnlyNoPileup) {
      fHistMultiplicityNoTPCOnlyNoPileup = new TH1F("fHistMultiplicityNoTPCOnlyNoPileup", 
         "Tracks per event;Nbr of Tracks;Events", 
         20000, 0, 20000);
      fListHistV0->Add(fHistMultiplicityNoTPCOnlyNoPileup);
   }

//------------------------------------------------
// Centrality Histograms
//------------------------------------------------
    
    if(! fHistCentralityBeforeTrigSel) {
        fHistCentralityBeforeTrigSel = new TH1F("fHistCentralityBeforeTrigSel",
                                                  "Centrality;Centrality;Events",
                                                  100, 0, 100);
        fListHistV0->Add(fHistCentralityBeforeTrigSel);
    }
    if(! fHistCentralityForTrigEvt) {
        fHistCentralityForTrigEvt = new TH1F("fHistCentralityForTrigEvt",
                                               "Centrality;Centrality;Events",
                                               100, 0, 100);
        fListHistV0->Add(fHistCentralityForTrigEvt);
    }
    if(! fHistCentrality) {
        fHistCentrality = new TH1F("fHistCentrality",
                                     "Centrality;Centrality;Events",
                                     100, 0, 100);
        fListHistV0->Add(fHistCentrality);
    }
    if(! fHistCentralityNoTPCOnly) {
        fHistCentralityNoTPCOnly = new TH1F("fHistCentralityNoTPCOnly",
                                              "Centrality;Centrality;Events",
                                              100, 0, 100);
        fListHistV0->Add(fHistCentralityNoTPCOnly);
    }
    if(! fHistCentralityNoTPCOnlyNoPileup) {
        fHistCentralityNoTPCOnlyNoPileup = new TH1F("fHistCentralityNoTPCOnlyNoPileup", 
                                                      "Centrality;Centrality;Events", 
                                                      100, 0, 100);
        fListHistV0->Add(fHistCentralityNoTPCOnlyNoPileup);
    }

    
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
    
//----------------------------------------
// event id histograms
//----------------------------------------
    
    if(! fHistRunNumber) {
        fHistRunNumber = new TH1D("fHistRunNumber",
                                  "Run Number",
                                  2000, 224000, 226000);
        fListHistV0->Add(fHistRunNumber);
    }
    if(! fHistPeriodNumber) {
        fHistPeriodNumber = new TH1D("fHistPeriodNumber",
                                     "Period Number",
                                     20, 0, 20);
        fListHistV0->Add(fHistPeriodNumber);
    }
    if(! fHistOrbitNumber) {
        fHistOrbitNumber = new TH1D("fHistOrbitNumber",
                                    "Orbit Number",
                                    1000, 0., 3.e7);
        fListHistV0->Add(fHistOrbitNumber);
    }
    if(! fHistBunchCrossNumber) {
        fHistBunchCrossNumber = new TH1D("fHistBunchCrossNumber",
                                         "Bunch Cross Number",
                                         4000, 0, 4000);
        fListHistV0->Add(fHistBunchCrossNumber);
    }
    
//----------------------------------------
// TPC signal histograms
//----------------------------------------
    
    const Double_t lPowerLowX = -1.5;
    const Double_t lPowerHighX = 1.5;
    const Int_t lBinNumbX = 1000*(lPowerHighX-lPowerLowX);
    Double_t lBinLimitsX[lBinNumbX+1];
    for(Int_t i = 0; i < lBinNumbX+1; i++) {
        lBinLimitsX[i] = TMath::Power(10, lPowerLowX + (Double_t)i/lBinNumbX*(lPowerHighX-lPowerLowX));
    }
    
    const Double_t lPowerLowY = 1.;
    const Double_t lPowerHighY = 4.;
    const Int_t lBinNumbY = 1000*(lPowerHighY-lPowerLowY);
    Double_t lBinLimitsY[lBinNumbY+1];
    for(Int_t i = 0; i < lBinNumbY+1; i++) {
        lBinLimitsY[i] = TMath::Power(10, lPowerLowY + (Double_t)i/lBinNumbY*(lPowerHighY-lPowerLowY));
    }
    
    if(! f2dHistTPCSignal) {
        f2dHistTPCSignal = new TH2D("f2dHistTPCSignal","TPC Signal;TPC Inner Wall Momentum;dE/dx",lBinNumbX,lBinLimitsX,lBinNumbY,lBinLimitsY);
        fListHistV0->Add(f2dHistTPCSignal);
    }
    if(! f2dHistTPCSignalPionsFromK0) {
        f2dHistTPCSignalPionsFromK0 = new TH2D("f2dHistTPCSignalPionsFromK0","TPC Signal of Pions from K0;TPC Inner Wall Momentum;dE/dx",lBinNumbX,lBinLimitsX,lBinNumbY,lBinLimitsY);
        fListHistV0->Add(f2dHistTPCSignalPionsFromK0);
    }
    if(! f2dHistTPCSignalPionsFromLambda) {
        f2dHistTPCSignalPionsFromLambda = new TH2D("f2dHistTPCSignalPionsFromLambda","TPC Signal of Pions from Lambda;TPC Inner Wall Momentum;dE/dx",lBinNumbX,lBinLimitsX,lBinNumbY,lBinLimitsY);
        fListHistV0->Add(f2dHistTPCSignalPionsFromLambda);
    }
    if(! f2dHistTPCSignalProtonsFromLambda) {
        f2dHistTPCSignalProtonsFromLambda = new TH2D("f2dHistTPCSignalProtonsFromLambda","TPC Signal of Protons from Lambda;TPC Inner Wall Momentum;dE/dx",lBinNumbX,lBinLimitsX,lBinNumbY,lBinLimitsY);
        fListHistV0->Add(f2dHistTPCSignalProtonsFromLambda);
    }
    
//----------------------------------------
// 3d histograms
//----------------------------------------
    
    const Int_t mBinNumb = 3000;
    Double_t mBinLimits[mBinNumb+1];
    for(Int_t i = 0; i < mBinNumb+1; i++)
    {
        //cast is required so that dividing returns decimals
        mBinLimits[i] = (Double_t)i*3/mBinNumb;
    }
    Double_t ptbinlimits[] = {0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4,4.4,4.9,5.5,6.0,7,8,10};
    Long_t ptbinnumb = sizeof(ptbinlimits)/sizeof(Double_t) - 1;
    Double_t centBinLimits[] = {0.,5.,10.,20.,40.,60.,80.,90.};
    Int_t nCentBin = sizeof(centBinLimits)/sizeof(Double_t) - 1;
    
    if(! f3dHistInvMassVsPtVsCentLambda) {
        f3dHistInvMassVsPtVsCentLambda = new TH3F("f3dHistInvMassVsPtVsCentLambda","Inv Mass Lambda Vs Pt Vs Centrality",mBinNumb,mBinLimits,ptbinnumb,ptbinlimits,nCentBin,centBinLimits);
        fListHistV0->Add(f3dHistInvMassVsPtVsCentLambda);
    }
    if(! f3dHistInvMassVsPtVsCentAntiLambda) {
        f3dHistInvMassVsPtVsCentAntiLambda = new TH3F("f3dHistInvMassVsPtVsCentAntiLambda","Inv Mass AntiLambda Vs Pt Vs Centrality",mBinNumb,mBinLimits,ptbinnumb,ptbinlimits,nCentBin,centBinLimits);
        fListHistV0->Add(f3dHistInvMassVsPtVsCentAntiLambda);
    }
    if(! f3dHistInvMassVsPtVsCentK0Short) {
        f3dHistInvMassVsPtVsCentK0Short = new TH3F("f3dHistInvMassVsPtVsCentK0Short","Inv Mass K0Short Vs Pt Vs Centrality",mBinNumb,mBinLimits,ptbinnumb,ptbinlimits,nCentBin,centBinLimits);
        fListHistV0->Add(f3dHistInvMassVsPtVsCentK0Short);
    }

    //histograms for cowboys
    if(! f3dHistInvMassVsPtVsCentLambdaCowboy) {
        f3dHistInvMassVsPtVsCentLambdaCowboy = new TH3F("f3dHistInvMassVsPtVsCentLambdaCowboy","Inv Mass Lambda Cowboy Vs Pt Vs Centrality",mBinNumb,mBinLimits,ptbinnumb,ptbinlimits,nCentBin,centBinLimits);
        fListHistV0->Add(f3dHistInvMassVsPtVsCentLambdaCowboy);
    }
    if(! f3dHistInvMassVsPtVsCentAntiLambdaCowboy) {
        f3dHistInvMassVsPtVsCentAntiLambdaCowboy = new TH3F("f3dHistInvMassVsPtVsCentAntiLambdaCowboy","Inv Mass AntiLambda Cowboy Vs Pt Vs Centrality",mBinNumb,mBinLimits,ptbinnumb,ptbinlimits,nCentBin,centBinLimits);
        fListHistV0->Add(f3dHistInvMassVsPtVsCentAntiLambdaCowboy);
    }
    if(! f3dHistInvMassVsPtVsCentK0ShortCowboy) {
        f3dHistInvMassVsPtVsCentK0ShortCowboy = new TH3F("f3dHistInvMassVsPtVsCentK0ShortCowboy","Inv Mass K0Short Cowboy Vs Pt Vs Centrality",mBinNumb,mBinLimits,ptbinnumb,ptbinlimits,nCentBin,centBinLimits);
        fListHistV0->Add(f3dHistInvMassVsPtVsCentK0ShortCowboy);
    }

    //histograms for sailors
    if(! f3dHistInvMassVsPtVsCentLambdaSailor) {
        f3dHistInvMassVsPtVsCentLambdaSailor = new TH3F("f3dHistInvMassVsPtVsCentLambdaSailor","Inv Mass Lambda Sailor Vs Pt Vs Centrality",mBinNumb,mBinLimits,ptbinnumb,ptbinlimits,nCentBin,centBinLimits);
        fListHistV0->Add(f3dHistInvMassVsPtVsCentLambdaSailor);
    }
    if(! f3dHistInvMassVsPtVsCentAntiLambdaSailor) {
        f3dHistInvMassVsPtVsCentAntiLambdaSailor = new TH3F("f3dHistInvMassVsPtVsCentAntiLambdaSailor","Inv Mass AntiLambda Sailor Vs Pt Vs Centrality",mBinNumb,mBinLimits,ptbinnumb,ptbinlimits,nCentBin,centBinLimits);
        fListHistV0->Add(f3dHistInvMassVsPtVsCentAntiLambdaSailor);
    }
    if(! f3dHistInvMassVsPtVsCentK0ShortSailor) {
        f3dHistInvMassVsPtVsCentK0ShortSailor = new TH3F("f3dHistInvMassVsPtVsCentK0ShortSailor","Inv Mass K0Short Sailor Vs Pt Vs Centrality",mBinNumb,mBinLimits,ptbinnumb,ptbinlimits,nCentBin,centBinLimits);
        fListHistV0->Add(f3dHistInvMassVsPtVsCentK0ShortSailor);
    }
    
    
   //Regular output: Histograms
   PostData(1, fListHistV0);
   //TTree Object: Saved to base directory. Should cache to disk while saving. 
   //(Important to avoid excessive memory usage, particularly when merging)
   PostData(2, fTree);

}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskExtractV0AODRun2::UserExec(Option_t *) 
{

   // Main loop
   // Called for each event
   //gObjectTable->Print();
   AliAODEvent *lAODevent = 0x0;

   //AliAODEvent *lAODevent = 0x0;
   Int_t    nV0s                        = -1;

   Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
   Double_t lMagneticField                 = -10.;
	
   // Connect to the InputEvent	
   // After these lines, we should have an ESD/AOD event + the number of cascades in it.

   // Appropriate for ESD analysis! 
		
   lAODevent = dynamic_cast<AliAODEvent*>( InputEvent() );
   if (!lAODevent) {
      AliWarning("ERROR: lAODevent not available \n");
      return;
   }
   fTreeVariableRunNumber = lAODevent->GetRunNumber();
   //fTreeVariableEventNumber = ( ( ((ULong64_t)lAODevent->GetPeriodNumber() ) << 36 ) | ( ((ULong64_t)lAODevent->GetOrbitNumber () ) << 12 ) | ((ULong64_t)lAODevent->GetBunchCrossNumber() )  );
    
    
    fTreeVariablePeriodNumber = lAODevent->GetPeriodNumber();
    fTreeVariableOrbitNumber = lAODevent->GetOrbitNumber();
    fTreeVariableBunchCrossNumber = lAODevent->GetBunchCrossNumber();
    //ULong64_t lGlobalOrbitNumber = (((((ULong64_t)fTreeVariableRunNumber-224000) << 4) | (ULong64_t)lPeriodNumber) << 25) | (ULong64_t)lOrbitNumber;
    
    fHistRunNumber->Fill(fTreeVariableRunNumber);
    fHistPeriodNumber->Fill(fTreeVariablePeriodNumber);
    fHistOrbitNumber->Fill(fTreeVariableOrbitNumber);
    fHistBunchCrossNumber->Fill(fTreeVariableBunchCrossNumber);
    
    /*
    cout << "       Run Number: " << fTreeVariableRunNumber << endl;
    cout << "       Period Number: " << lPeriodNumber << endl;
    cout << "       Orbit Number: " << lOrbitNumber << endl;
    cout << "       Bunch Cross Number: " << lBunchCrossNumber << endl;
    cout << "       Event Type: " << lAODevent->GetEventType() << endl;
    */

   //REVISED multiplicity estimator after 'multiplicity day' (2011)
   Double_t lCentrality = -100;

   if(fkIsNuclear == kFALSE) lCentrality = -1;

   //---> If this is a nuclear collision, then go nuclear on "multiplicity" variable...
   //---> Warning: Experimental
   if(fkIsNuclear == kTRUE){
       if(!fkLowE) {
           AliMultSelection* lCentSelect = 0x0;
           lCentSelect = (AliMultSelection*)lAODevent->FindListObject("MultSelection");
           if(!lCentSelect) {
               AliWarning("AliMultSelection object not found!");
           }
           lCentrality = lCentSelect->GetMultiplicityPercentile("V0M");
       } else {
           AliCentrality* centrality;
           centrality = lAODevent->GetCentrality();
           lCentrality = centrality->GetCentralityPercentile("V0M");
           if (centrality->GetQuality()>1) {
               PostData(1, fListHistV0);
               PostData(2, fTree);
               return;
           }
       }
   }

   //Set variable for filling tree afterwards!
   //---> pp case......: GetReferenceMultiplicity
   //---> Pb-Pb case...: Centrality by V0M
   fTreeVariableCentrality = (Int_t)lCentrality;
    Int_t lMultiplicity = lAODevent->GetNumberOfTracks();
    fHistMultiplicityBeforeTrigSel->Fill ( lMultiplicity );
   fHistCentralityBeforeTrigSel->Fill ( lCentrality );
   fHistV0MultiplicityBeforeTrigSel->Fill ( lAODevent->GetNumberOfV0s() );

//------------------------------------------------
// Physics Selection
//------------------------------------------------
    
// new method        
//   UInt_t maskIsSelected = 1;
   UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
   Bool_t isSelected = 0;
   //kMB: default selection, also if fTriggerMask is something not understood...
   isSelected = (maskIsSelected & AliVEvent::kMB) == AliVEvent::kMB;
    
   if( fTriggerMask == "kINT7" )
     isSelected = (maskIsSelected & AliVEvent::kINT7) == AliVEvent::kINT7;
   if( fTriggerMask == "kINT8" )
     isSelected = (maskIsSelected & AliVEvent::kINT8) == AliVEvent::kINT8;
   if( fTriggerMask == "kAnyINT" )
     isSelected = (maskIsSelected & AliVEvent::kAnyINT) == AliVEvent::kAnyINT;

   //pp at 2.76TeV: special case, ignore FastOnly
   if ( (fkLowEnergyPP == kTRUE) && (maskIsSelected& AliVEvent::kFastOnly) ){
      PostData(1, fListHistV0);
      PostData(2, fTree);
      return;
   } 
   //Standard Min-Bias Selection
   if ( ! isSelected ) { 
      PostData(1, fListHistV0);
      PostData(2, fTree);
      return;
   }

//------------------------------------------------
// After Trigger Selection
//------------------------------------------------

	nV0s = lAODevent->GetNumberOfV0s();

  fHistV0MultiplicityForTrigEvt     ->Fill( nV0s       );
  fHistMultiplicityForTrigEvt       ->Fill( lMultiplicity  );
    fHistCentralityForTrigEvt       ->Fill( lCentrality  );

//------------------------------------------------
// Getting: Primary Vertex + MagField Info
//------------------------------------------------
        
	const AliAODVertex *lPrimaryBestAODVtx = lAODevent->GetPrimaryVertex();	
  // get the best primary vertex available for the event
	// As done in AliCascadeVertexer, we keep the one which is the best one available.
	// between : Tracking vertex > SPD vertex > TPC vertex > default SPD vertex
	// This one will be used for next calculations (DCA essentially)
   lPrimaryBestAODVtx->GetXYZ( lBestPrimaryVtxPos );

   Double_t tPrimaryVtxPosition[3];
   const AliVVertex *primaryVtx = lAODevent->GetPrimaryVertex();
   tPrimaryVtxPosition[0] = primaryVtx->GetX();
   tPrimaryVtxPosition[1] = primaryVtx->GetY();
   tPrimaryVtxPosition[2] = primaryVtx->GetZ();
   fHistPVx->Fill( tPrimaryVtxPosition[0] );
   fHistPVy->Fill( tPrimaryVtxPosition[1] );
   fHistPVz->Fill( tPrimaryVtxPosition[2] );
    
    lMagneticField = lAODevent->GetMagneticField( );

//------------------------------------------------
// Primary Vertex Z position: SKIP
//------------------------------------------------

   if(TMath::Abs(lBestPrimaryVtxPos[2]) > 10.0 ) { 
      AliWarning("Pb / | Z position of Best Prim Vtx | > 10.0 cm ... return !"); 
      PostData(1, fListHistV0);
      PostData(2, fTree);
      return; 
   }

//------------------------------------------------
// Centrality cut
//------------------------------------------------
    
    if((lCentrality < 0.) || (lCentrality > 90.)) {
        AliWarning("Pb / Centrality not in (0,90) ... return !");
        PostData(1, fListHistV0);
        PostData(2, fTree);
        return; 
    }
    
   fHistV0MultiplicityForSelEvt ->Fill( nV0s );
   fHistMultiplicity->Fill(lMultiplicity);
    fHistCentrality->Fill(lCentrality);

//------------------------------------------------
// Only look at events with well-established PV
//------------------------------------------------
	
   const AliAODVertex *lPrimaryTrackingAODVtxCheck = lAODevent->GetPrimaryVertexTPC();
   const AliAODVertex *lPrimarySPDVtx = lAODevent->GetPrimaryVertexSPD();
   if (!lPrimarySPDVtx && !lPrimaryTrackingAODVtxCheck ){
      AliWarning("Pb / No SPD prim. vertex nor prim. Tracking vertex ... return !");
      PostData(1, fListHistV0);
      PostData(2, fTree);
      return;
   }

   fHistV0MultiplicityForSelEvtNoTPCOnly ->Fill( nV0s );
   fHistMultiplicityNoTPCOnly->Fill(lMultiplicity);
    fHistCentralityNoTPCOnly->Fill(lCentrality);

//------------------------------------------------
// Pileup Rejection
//------------------------------------------------

   // FIXME : quality selection regarding pile-up rejection 
   if(lAODevent->IsPileupFromSPD() && !fkIsNuclear){// minContributors=3, minZdist=0.8, nSigmaZdist=3., nSigmaDiamXY=2., nSigmaDiamZ=5.  -> see http://alisoft.cern.ch/viewvc/trunk/STEER/AliESDEvent.h?root=AliRoot&r1=41914&r2=42199&pathrev=42199
      PostData(1, fListHistV0);
      PostData(2, fTree); 
      return;
   }
   fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup ->Fill( nV0s );
   fHistMultiplicityNoTPCOnlyNoPileup->Fill(lMultiplicity);
    fHistCentralityNoTPCOnlyNoPileup->Fill(lCentrality);
  
//------------------------------------------------
// MAIN LAMBDA LOOP STARTS HERE
//------------------------------------------------

   //if( ! (lAODevent->GetPrimaryVertex()->IsFromVertexerZ() )	 ){ 
      fHistPVxAnalysis->Fill( tPrimaryVtxPosition[0] );
      fHistPVyAnalysis->Fill( tPrimaryVtxPosition[1] );
      fHistPVzAnalysis->Fill( tPrimaryVtxPosition[2] );
   //}
    
    fTreeVariablePrimVX = tPrimaryVtxPosition[0];
    fTreeVariablePrimVY = tPrimaryVtxPosition[1];
    fTreeVariablePrimVZ = tPrimaryVtxPosition[2];

//Variable definition
   Int_t    lOnFlyStatus = 0;// nv0sOn = 0, nv0sOff = 0;
   Double_t lChi2V0 = 0;
   Double_t lDcaV0Daughters = 0, lDcaV0ToPrimVertex = 0;
   Double_t lDcaPosToPrimVertex = 0, lDcaNegToPrimVertex = 0;
   Double_t lV0CosineOfPointingAngle = 0;
   Double_t lV0Radius = 0, lPt = 0;
   Double_t lRapK0Short = 0, lRapLambda = 0;
   Double_t lInvMassK0s = 0, lInvMassLambda = 0, lInvMassAntiLambda = 0;
   Double_t lAlphaV0 = 0, lAlphaV0BoostAsK0 = 0, lPtArmV0 = 0;

   Double_t fMinV0Pt = 0; 
   Double_t fMaxV0Pt = 100; 

    Int_t nTracks = lAODevent->GetNumberOfTracks();
   Int_t nv0s = 0;
   nv0s = lAODevent->GetNumberOfV0s();
    
    for (Int_t i = 0; i < nTracks; i++)
    {//begin track selection
        AliAODTrack* tr = (AliAODTrack*)lAODevent->GetTrack(i);
        
        // TPC refit
        if (!tr->IsOn(AliAODTrack::kTPCrefit)) continue;
        // Minimum number of clusters
        Float_t nCrossedRowsTPC = tr->GetTPCClusterInfo(2,1);
        if (nCrossedRowsTPC < 70) continue;
        Int_t findable = tr->GetTPCNclsF();
        if (findable <= 0) continue;
        if (nCrossedRowsTPC/findable < 0.8) continue;
        if (tr->GetKinkIndex(0) > 0) continue;
        
        Double_t lTPCSignal = tr->GetTPCsignal();
        Double_t lTPCMomentum = tr->GetTPCmomentum();
        f2dHistTPCSignal->Fill(lTPCMomentum,lTPCSignal);
    }//end track selection

   for (Int_t iV0 = 0; iV0 < nv0s; iV0++) //extra-crazy test
   {// This is the begining of the V0 loop
      AliAODv0 *v0 = lAODevent->GetV0(iV0);
      if (!v0) continue;

      //Obsolete at AOD level... 
      //---> Fix On-the-Fly candidates, count how many swapped
      //if( v0->GetParamN()->Charge() > 0 && v0->GetParamP()->Charge() < 0 ){
      //  fHistSwappedV0Counter -> Fill( 1 );
      //}else{
      //  fHistSwappedV0Counter -> Fill( 0 ); 
      //}
      //if ( fkUseOnTheFly ) CheckChargeV0(v0); 
      
      Double_t tDecayVertexV0[3]; v0->GetXYZ(tDecayVertexV0); 
      //Double_t tV0mom[3];
	//nefunguje v AOD, iba v ESD:
//      v0->GetPxPyPz( tV0mom ); 
	//namiesto toho:
	//tV0mom[0] = v0->Px();
	//tV0mom[1] = v0->Py();
	//tV0mom[2] = v0->Pz();
      //Double_t lV0TotalMomentum = TMath::Sqrt(tV0mom[0]*tV0mom[0]+tV0mom[1]*tV0mom[1]+tV0mom[2]*tV0mom[2]);
       Double_t lV0TotalMomentum = v0->P();

      lV0Radius = TMath::Sqrt(tDecayVertexV0[0]*tDecayVertexV0[0]+tDecayVertexV0[1]*tDecayVertexV0[1]);

        //Set Variables for later filling
        fTreeVariableV0Px = v0->Px();
        fTreeVariableV0Py = v0->Py();
        fTreeVariableV0Pz = v0->Pz();

      lPt = v0->Pt();
      lRapK0Short = v0->RapK0Short();
      lRapLambda  = v0->RapLambda();
      if ((lPt<fMinV0Pt)||(fMaxV0Pt<lPt)) continue;

      //UInt_t lKeyPos = (UInt_t)TMath::Abs(v0->GetPosID());
      //UInt_t lKeyNeg = (UInt_t)TMath::Abs(v0->GetPosID());

      Double_t lMomPos[3]; //v0->GetPPxPyPz(lMomPos[0],lMomPos[1],lMomPos[2]);
      Double_t lMomNeg[3]; //v0->GetNPxPyPz(lMomNeg[0],lMomNeg[1],lMomNeg[2]);
      lMomPos[0] = v0->MomPosX();
      lMomPos[1] = v0->MomPosY();
      lMomPos[2] = v0->MomPosZ();
      lMomNeg[0] = v0->MomNegX();
      lMomNeg[1] = v0->MomNegY();
      lMomNeg[2] = v0->MomNegZ();

      AliAODTrack *pTrack=(AliAODTrack *)v0->GetDaughter(0); //0->Positive Daughter
      AliAODTrack *nTrack=(AliAODTrack *)v0->GetDaughter(1); //1->Negative Daughter
      if (!pTrack || !nTrack) {
         Printf("ERROR: Could not retreive one of the daughter track");
         continue;
      }

      //Daughter Eta for Eta selection, afterwards
      fTreeVariableNegEta = nTrack->Eta();
      fTreeVariablePosEta = pTrack->Eta();

      // Filter like-sign V0 (next: add counter and distribution)
      if ( pTrack->Charge() == nTrack->Charge()){
         continue;
      } 

      //Quick test this far! 
      

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

      if ( fTreeVariableLeastNbrCrossedRows < 70 ) continue;
	
      //GetKinkIndex condition
      //if( pTrack->GetKinkIndex(0)>0 || nTrack->GetKinkIndex(0)>0 ) continue;

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

      //End track Quality Cuts
      //________________________________________________________________________

      
      //ESD code: obsolete 
      //lDcaPosToPrimVertex = TMath::Abs(pTrack->GetD(tPrimaryVtxPosition[0],
			//				tPrimaryVtxPosition[1],
			//				lMagneticField) );

      //lDcaNegToPrimVertex = TMath::Abs(nTrack->GetD(tPrimaryVtxPosition[0],
			//				tPrimaryVtxPosition[1],
			//				lMagneticField) );

      lDcaPosToPrimVertex = v0->DcaPosToPrimVertex();
      lDcaNegToPrimVertex = v0->DcaNegToPrimVertex();

        
      lOnFlyStatus = v0->GetOnFlyStatus();
      lChi2V0 = v0->Chi2V0();
      lDcaV0Daughters = v0->DcaV0Daughters();
      lDcaV0ToPrimVertex = v0->DcaV0ToPrimVertex();
      lV0CosineOfPointingAngle = v0->CosPointingAngle(tPrimaryVtxPosition);

      // Getting invariant mass infos directly from ESD
      lInvMassK0s        = v0->MassK0Short();
      lInvMassLambda     = v0->MassLambda();
      lInvMassAntiLambda = v0->MassAntiLambda();
      lAlphaV0 = v0->AlphaV0();
      lPtArmV0 = v0->PtArmV0();
       
       
       // get daughter track momentum in v0 rest frame
       const Double_t kMassLambda = 1.115683;
       const Double_t kMassProton = 0.938272;
       const Double_t kMassK0Short = 0.497614;
       const Double_t kMassPion = 0.139570;
       Double_t lLongMomPos = lV0TotalMomentum*(1+lAlphaV0)/2;
       Double_t lEV0Lambda = TMath::Sqrt(lV0TotalMomentum*lV0TotalMomentum + kMassLambda*kMassLambda);
       Double_t lEV0K0Short = TMath::Sqrt(lV0TotalMomentum*lV0TotalMomentum + kMassK0Short*kMassK0Short);
       Double_t lPPos = pTrack->P();
       Double_t lEPosProton = TMath::Sqrt(lPPos*lPPos + kMassProton*kMassProton);
       Double_t lEPosPion = TMath::Sqrt(lPPos*lPPos + kMassPion*kMassPion);
       Double_t fTreeVariableDecayMomLambda = TMath::Sqrt(TMath::Power(lEV0Lambda/kMassLambda*lLongMomPos - lV0TotalMomentum/kMassLambda*lEPosProton,2) + lPtArmV0*lPtArmV0);
       Double_t fTreeVariableDecayMomAntiLambda = TMath::Sqrt(TMath::Power(lEV0Lambda/kMassLambda*lLongMomPos - lV0TotalMomentum/kMassLambda*lEPosPion,2) + lPtArmV0*lPtArmV0);
       Double_t fTreeVariableDecayMomK0Short = TMath::Sqrt(TMath::Power(lEV0K0Short/kMassK0Short*lLongMomPos - lV0TotalMomentum/kMassK0Short*lEPosPion,2) + lPtArmV0*lPtArmV0);
       
       //get ArmPod alpha of V0 boosted as K0
       Double_t lPNeg = nTrack->P();
       Double_t lENegPion = TMath::Sqrt(lPNeg*lPNeg + kMassPion*kMassPion);
       lAlphaV0BoostAsK0 = (lAlphaV0*lV0TotalMomentum + lEPosPion - lENegPion)/(lV0TotalMomentum + lEV0K0Short);


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
       fTreeVariableAlphaV0BoostAsK0 = lAlphaV0BoostAsK0;
      fTreeVariablePtArmV0 = lPtArmV0;

      //Official means of acquiring N-sigmas 
      fTreeVariableNSigmasPosProton = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kProton );
      fTreeVariableNSigmasPosPion   = fPIDResponse->NumberOfSigmasTPC( pTrack, AliPID::kPion );
      fTreeVariableNSigmasNegProton = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kProton );
      fTreeVariableNSigmasNegPion   = fPIDResponse->NumberOfSigmasTPC( nTrack, AliPID::kPion );
       
       //TPC Signal Information
       fTreeVariableNegTPCSignal = nTrack->GetTPCsignal();
       fTreeVariablePosTPCSignal = pTrack->GetTPCsignal();
       fTreeVariableNegInnerP = nTrack->GetTPCmomentum();
       fTreeVariablePosInnerP = pTrack->GetTPCmomentum();

//This requires an Invariant Mass Hypothesis afterwards
      fTreeVariableDistOverTotMom = TMath::Sqrt(
						TMath::Power( tDecayVertexV0[0] - lBestPrimaryVtxPos[0] , 2) +
						TMath::Power( tDecayVertexV0[1] - lBestPrimaryVtxPos[1] , 2) +
						TMath::Power( tDecayVertexV0[2] - lBestPrimaryVtxPos[2] , 2)
					);
      fTreeVariableDistOverTotMom /= (lV0TotalMomentum+1e-10); //avoid division by zero, to be sure
       
       fTreeVariableCowboy = kFALSE;
       if ( lMagneticField * ( pTrack->Px()*nTrack->Py() - pTrack->Py()*nTrack->Px() ) < 0 ) fTreeVariableCowboy = kTRUE;

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
          Bool_t lSel = kFALSE;
          //Lambda
          if(fkSaveAllInvMasses || (!fkSaveAllInvMasses && (fTreeVariableInvMassLambda < lUpperLimitLambda && fTreeVariableInvMassLambda > lLowerLimitLambda))) {
              if(!fkPreSelect) lSel = kTRUE;
              else {
                  if(TMath::Abs(fTreeVariableNSigmasPosProton) < 10.0 &&
                     TMath::Abs(fTreeVariableNSigmasNegPion) < 10.0
                     )
                      lSel = kTRUE;
              }
          }
          //AntiLambda
          if(fkSaveAllInvMasses || (!fkSaveAllInvMasses && (fTreeVariableInvMassAntiLambda < lUpperLimitLambda && fTreeVariableInvMassAntiLambda > lLowerLimitLambda))) {
              if(!fkPreSelect) lSel = kTRUE;
              else {
                  if(TMath::Abs(fTreeVariableNSigmasNegProton) < 10.0 &&
                     TMath::Abs(fTreeVariableNSigmasPosPion) < 10.0
                     )
                      lSel = kTRUE;
              }
          }
          //K0Short
          if(fkSaveAllInvMasses || (!fkSaveAllInvMasses && (fTreeVariableInvMassK0s < lUpperLimitK0Short && fTreeVariableInvMassK0s > lLowerLimitK0Short))) {
              if(!fkPreSelect) lSel = kTRUE;
              else {
                  if(TMath::Abs(fTreeVariableNSigmasNegPion) < 10.0 &&
                     TMath::Abs(fTreeVariableNSigmasPosPion) < 10.0
                     )
                      lSel = kTRUE;
              }
          }
          if(lSel)
          {
              if(!fkPreSelect) fTree->Fill();
              else {
                  if(TMath::Abs(fTreeVariableNegEta)<0.8 &&
                     TMath::Abs(fTreeVariablePosEta)<0.8 &&
                     lPtArmV0 > 0.03 &&
                     lV0Radius > 1.
                      )
                      fTree->Fill();
              }
          }
          
          //Default cuts
          if(TMath::Abs(fTreeVariableNegEta)       <= 0.8               &&
             TMath::Abs(fTreeVariablePosEta)       <= 0.8               &&
             lV0Radius                 >= 5.0               &&
             lV0Radius                 <= 100.              &&
             lDcaNegToPrimVertex       >= 0.1                &&
             lDcaPosToPrimVertex       >= 0.1                &&
             lDcaV0Daughters           <= 1.0            &&
             lV0CosineOfPointingAngle    >= 0.998
             ) {
              Double_t kMassLambda = 1.115683;
              Double_t kMassK0Short = 0.4976;
              
              //Lambda
              if(TMath::Abs(lRapLambda)<0.5 &&
                 kMassLambda*fTreeVariableDistOverTotMom <= 3*7.89     &&
                 //TMath::Abs(lInvMassK0s - kMassK0Short) > 0.010 &&
                 TMath::Abs(fTreeVariableNSigmasNegPion)   <= 6. &&
                 TMath::Abs(fTreeVariableNSigmasPosProton) <= 6.
                 )
              {
                  f3dHistInvMassVsPtVsCentLambda->Fill(lInvMassLambda,lPt,lCentrality);
                  if(fTreeVariableCowboy) f3dHistInvMassVsPtVsCentLambdaCowboy->Fill(lInvMassLambda,lPt,lCentrality);
                  else f3dHistInvMassVsPtVsCentLambdaSailor->Fill(lInvMassLambda,lPt,lCentrality);
                  
                  if (TMath::Abs(lInvMassK0s-kMassK0Short) > 0.0025 &&
                      TMath::Abs(lInvMassLambda-kMassLambda) < 0.001 &&
                      TMath::Abs(lInvMassAntiLambda-kMassLambda) > 0.002 &&
                      lDcaV0Daughters                       < 0.3 &&
                      lV0CosineOfPointingAngle              > 0.999 &&
                      pTrack->Pt()                          > 0.15 &&
                      nTrack->Pt()                          > 0.15
                      ) {
                      f2dHistTPCSignalProtonsFromLambda->Fill(fTreeVariablePosInnerP,fTreeVariablePosTPCSignal);
                      f2dHistTPCSignalPionsFromLambda->Fill(fTreeVariableNegInnerP,fTreeVariableNegTPCSignal);
                  }
              }
              //AntiLambda
              if(TMath::Abs(lRapLambda)<0.5 &&
                 kMassLambda*fTreeVariableDistOverTotMom <= 3*7.89     &&
                 //TMath::Abs(lInvMassK0s - kMassK0Short) > 0.010 &&
                 TMath::Abs(fTreeVariableNSigmasPosPion)   <= 6. &&
                 TMath::Abs(fTreeVariableNSigmasNegProton) <= 6.
                 )
              {
                  f3dHistInvMassVsPtVsCentAntiLambda->Fill(lInvMassAntiLambda,lPt,lCentrality);
                  if(fTreeVariableCowboy) f3dHistInvMassVsPtVsCentAntiLambdaCowboy->Fill(lInvMassAntiLambda,lPt,lCentrality);
                  else f3dHistInvMassVsPtVsCentAntiLambdaSailor->Fill(lInvMassAntiLambda,lPt,lCentrality);
                  
                  if (TMath::Abs(lInvMassK0s-kMassK0Short) > 0.0025 &&
                      TMath::Abs(lInvMassLambda-kMassLambda) > 0.002 &&
                      TMath::Abs(lInvMassAntiLambda-kMassLambda) < 0.001 &&
                      lDcaV0Daughters                       < 0.3 &&
                      lV0CosineOfPointingAngle              > 0.999 &&
                      pTrack->Pt()                          > 0.15 &&
                      nTrack->Pt()                          > 0.15
                      ) {
                      f2dHistTPCSignalPionsFromLambda->Fill(fTreeVariablePosInnerP,fTreeVariablePosTPCSignal);
                      f2dHistTPCSignalProtonsFromLambda->Fill(fTreeVariableNegInnerP,fTreeVariableNegTPCSignal);
                  }
              }
              //K0Short
              if(TMath::Abs(lRapK0Short)<0.5 &&
                 kMassK0Short*fTreeVariableDistOverTotMom <= 3*2.68     &&
                 //TMath::Abs(lInvMassLambda - kMassLambda) > 0.005 &&
                 //TMath::Abs(lInvMassAntiLambda - kMassLambda) > 0.005 &&
                 ( lPtArmV0*5>TMath::Abs(lAlphaV0) ) &&
                 TMath::Abs(fTreeVariableNSigmasPosPion)   <= 6. &&
                 TMath::Abs(fTreeVariableNSigmasNegPion)   <= 6.
                 )
              {
                  f3dHistInvMassVsPtVsCentK0Short->Fill(lInvMassK0s,lPt,lCentrality);
                  if(fTreeVariableCowboy) f3dHistInvMassVsPtVsCentK0ShortCowboy->Fill(lInvMassK0s,lPt,lCentrality);
                  else f3dHistInvMassVsPtVsCentK0ShortSailor->Fill(lInvMassK0s,lPt,lCentrality);
                  
                  if (TMath::Abs(lInvMassK0s-kMassK0Short) < 0.0015 &&
                      TMath::Abs(lInvMassLambda-kMassLambda) > 0.002 &&
                      TMath::Abs(lInvMassAntiLambda-kMassLambda) > 0.002 &&
                      lDcaV0Daughters                       < 0.3 &&
                      lV0CosineOfPointingAngle              > 0.999 &&
                      pTrack->Pt()                          > 0.15 &&
                      nTrack->Pt()                          > 0.15
                      ) {
                      f2dHistTPCSignalPionsFromK0->Fill(fTreeVariablePosInnerP,fTreeVariablePosTPCSignal);
                      f2dHistTPCSignalPionsFromK0->Fill(fTreeVariableNegInnerP,fTreeVariableNegTPCSignal);
                  }
              }
          } // end if cuts
          
      } // end if onfly selection

//------------------------------------------------
// Fill tree over.
//------------------------------------------------

   }// This is the end of the V0 loop

  // Post output data.
    PostData(1, fListHistV0);
    PostData(2, fTree);


}

//________________________________________________________________________
void AliAnalysisTaskExtractV0AODRun2::Terminate(Option_t *)
{
   // Draw result to the screen
   // Called once at the end of the query
   // This will draw the V0 candidate multiplicity, whose
   // number of entries corresponds to the number of triggered events. 
   TList *cRetrievedList = 0x0;
   cRetrievedList = (TList*)GetOutputData(1);
   if(!cRetrievedList){
      Printf("ERROR - AliAnalysisTaskExtractV0AODRun2 : ouput data container list not available\n");
      return;
   }		
   fHistV0MultiplicityForTrigEvt = dynamic_cast<TH1F*> (  cRetrievedList->FindObject("fHistV0MultiplicityForTrigEvt")  );
   if (!fHistV0MultiplicityForTrigEvt) {
      Printf("ERROR - AliAnalysisTaskExtractV0AODRun2 : fHistV0MultiplicityForTrigEvt not available");
      return;
   }
   TCanvas *canCheck = new TCanvas("AliAnalysisTaskExtractV0AODRun2","V0 Multiplicity",10,10,510,510);
   canCheck->cd(1)->SetLogy();
   fHistV0MultiplicityForTrigEvt->SetMarkerStyle(22);
   fHistV0MultiplicityForTrigEvt->DrawCopy("E");
}

//________________________________________________________________________
void AliAnalysisTaskExtractV0AODRun2::CheckChargeV0(AliESDv0 *v0)
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
