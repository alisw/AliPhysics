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
// --- Adapted to look for lambdas as well, using code from 
//        AliAnalysisTaskCheckPerformanceStrange.cxx
//
//  --- Algorithm Description 
//   1. Loop over primaries in stack to acquire generated charged Xi
//   2. Loop over stack to find V0s, fill TH3Fs "PrimRawPt"s for Efficiency
//   3. Perform Physics Selection
//   4. Perform Primary Vertex |z|<10cm selection
//   5. Perform Primary Vertex NoTPCOnly vertexing selection (>0 contrib.)
//   6. Perform Pileup Rejection
//   7. Analysis Loops: 
//    7a. Fill TH3Fs "PrimAnalysisPt" for control purposes only
//    7b. Fill TTree object with V0 information, candidates
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
#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"
#include "AliESDpid.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include "AliCFContainer.h"
#include "AliMultiplicity.h"
#include "AliAODMCParticle.h"
#include "AliESDcascade.h"
#include "AliAODcascade.h"
#include "AliESDUtils.h"

#include "AliAnalysisTaskExtractPerformanceV0.h"

ClassImp(AliAnalysisTaskExtractPerformanceV0)

AliAnalysisTaskExtractPerformanceV0::AliAnalysisTaskExtractPerformanceV0() 
  : AliAnalysisTaskSE(), fListHistV0(0), fTree(0), fPIDResponse(0),
   fkIsNuclear   ( kFALSE ), 
   fkLowEnergyPP ( kFALSE ),
   fkUseOnTheFly ( kFALSE ),

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

//------------------------------------------------
// PARTICLE HISTOGRAMS
// --- Filled on a Particle-by-Particle basis
//------------------------------------------------
   f3dHistPrimAnalysisPtVsYVsMultLambda(0),
   f3dHistPrimAnalysisPtVsYVsMultAntiLambda(0),
   f3dHistPrimAnalysisPtVsYVsMultK0Short(0),
   f3dHistPrimRawPtVsYVsMultLambda(0),
   f3dHistPrimRawPtVsYVsMultAntiLambda(0),
   f3dHistPrimRawPtVsYVsMultK0Short(0),
   f3dHistPrimRawPtVsYVsDecayLengthLambda(0),
   f3dHistPrimRawPtVsYVsDecayLengthAntiLambda(0),
   f3dHistPrimRawPtVsYVsDecayLengthK0Short(0),
   f3dHistGenPtVsYVsMultXiMinus(0),
   f3dHistGenPtVsYVsMultXiPlus(0),
   fHistPVx(0),
   fHistPVy(0),
   fHistPVz(0),
   fHistPVxAnalysis(0),
   fHistPVyAnalysis(0),
   fHistPVzAnalysis(0),
   fHistPVxAnalysisHasHighPtLambda(0),
   fHistPVyAnalysisHasHighPtLambda(0),
   fHistPVzAnalysisHasHighPtLambda(0)
{
  // Dummy Constructor
}

AliAnalysisTaskExtractPerformanceV0::AliAnalysisTaskExtractPerformanceV0(const char *name) 
  : AliAnalysisTaskSE(name), fListHistV0(0), fTree(0), fPIDResponse(0),
   fkIsNuclear   ( kFALSE ), 
   fkLowEnergyPP ( kFALSE ),
   fkUseOnTheFly ( kFALSE ),
     
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


//------------------------------------------------
// PARTICLE HISTOGRAMS
// --- Filled on a Particle-by-Particle basis
//------------------------------------------------
   f3dHistPrimAnalysisPtVsYVsMultLambda(0),
   f3dHistPrimAnalysisPtVsYVsMultAntiLambda(0),
   f3dHistPrimAnalysisPtVsYVsMultK0Short(0),
   f3dHistPrimRawPtVsYVsMultLambda(0),
   f3dHistPrimRawPtVsYVsMultAntiLambda(0),
   f3dHistPrimRawPtVsYVsMultK0Short(0),
   f3dHistPrimRawPtVsYVsDecayLengthLambda(0),
   f3dHistPrimRawPtVsYVsDecayLengthAntiLambda(0),
   f3dHistPrimRawPtVsYVsDecayLengthK0Short(0),
   f3dHistGenPtVsYVsMultXiMinus(0),
   f3dHistGenPtVsYVsMultXiPlus(0),
   fHistPVx(0),
   fHistPVy(0),
   fHistPVz(0),
   fHistPVxAnalysis(0),
   fHistPVyAnalysis(0),
   fHistPVzAnalysis(0),
   fHistPVxAnalysisHasHighPtLambda(0),
   fHistPVyAnalysisHasHighPtLambda(0),
   fHistPVzAnalysisHasHighPtLambda(0)
{
   // Constructor
   // Output slot #0 writes into a TList container (Cascade)
   DefineOutput(1, TList::Class());
   DefineOutput(2, TTree::Class());
}


AliAnalysisTaskExtractPerformanceV0::~AliAnalysisTaskExtractPerformanceV0()
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
void AliAnalysisTaskExtractPerformanceV0::UserCreateOutputObjects()
{

   OpenFile(2);	
   // Called once
   fTree = new TTree("fTree","V0Candidates");

//------------------------------------------------
// fTree Branch definitions
//------------------------------------------------

//-----------BASIC-INFO---------------------------
/* 1*/   fTree->Branch("fTreeVariablePrimaryStatus",&fTreeVariablePrimaryStatus,"fTreeVariablePrimaryStatus/I");	
/* 1*/   fTree->Branch("fTreeVariablePrimaryStatusMother",&fTreeVariablePrimaryStatusMother,"fTreeVariablePrimaryStatusMother/I");	
/* 2*/   fTree->Branch("fTreeVariableChi2V0",&fTreeVariableChi2V0,"Chi2V0/F");
/* 3*/   fTree->Branch("fTreeVariableDcaV0Daughters",&fTreeVariableDcaV0Daughters,"fTreeVariableDcaV0Daughters/F");
/* 4*/   fTree->Branch("fTreeVariableDcaPosToPrimVertex",&fTreeVariableDcaPosToPrimVertex,"fTreeVariableDcaPosToPrimVertex/F");
/* 5*/   fTree->Branch("fTreeVariableDcaNegToPrimVertex",&fTreeVariableDcaNegToPrimVertex,"fTreeVariableDcaNegToPrimVertex/F");
/* 6*/   fTree->Branch("fTreeVariableV0Radius",&fTreeVariableV0Radius,"fTreeVariableV0Radius/F");
/* 7*/   fTree->Branch("fTreeVariablePt",&fTreeVariablePt,"fTreeVariablePt/F");
/* 7*/   fTree->Branch("fTreeVariablePtMC",&fTreeVariablePtMC,"fTreeVariablePtMC/F");
/* 8*/   fTree->Branch("fTreeVariableRapK0Short",&fTreeVariableRapK0Short,"fTreeVariableRapK0Short/F");
/* 9*/   fTree->Branch("fTreeVariableRapLambda",&fTreeVariableRapLambda,"fTreeVariableRapLambda/F");
/*10*/   fTree->Branch("fTreeVariableRapMC",&fTreeVariableRapMC,"fTreeVariableRapMC/F");
/*11*/   fTree->Branch("fTreeVariableInvMassK0s",&fTreeVariableInvMassK0s,"fTreeVariableInvMassK0s/F");
/*12*/   fTree->Branch("fTreeVariableInvMassLambda",&fTreeVariableInvMassLambda,"fTreeVariableInvMassLambda/F");
/*13*/   fTree->Branch("fTreeVariableInvMassAntiLambda",&fTreeVariableInvMassAntiLambda,"fTreeVariableInvMassAntiLambda/F");
/*14*/   fTree->Branch("fTreeVariableAlphaV0",&fTreeVariableAlphaV0,"fTreeVariableAlphaV0/F");
/*15*/   fTree->Branch("fTreeVariablePtArmV0",&fTreeVariablePtArmV0,"fTreeVariablePtArmV0/F");
/*16*/   fTree->Branch("fTreeVariableNegTransvMomentum",&fTreeVariableNegTransvMomentum,"fTreeVariableNegTransvMomentum/F");
/*17*/   fTree->Branch("fTreeVariablePosTransvMomentum",&fTreeVariablePosTransvMomentum,"fTreeVariablePosTransvMomentum/F");
/*18*/   fTree->Branch("fTreeVariableNegTransvMomentumMC",&fTreeVariableNegTransvMomentumMC,"fTreeVariableNegTransvMomentumMC/F");
/*19*/   fTree->Branch("fTreeVariablePosTransvMomentumMC",&fTreeVariablePosTransvMomentumMC,"fTreeVariablePosTransvMomentumMC/F");
/*20*/   fTree->Branch("fTreeVariableLeastNbrCrossedRows",&fTreeVariableLeastNbrCrossedRows,"fTreeVariableLeastNbrCrossedRows/I");
/*21*/   fTree->Branch("fTreeVariableLeastRatioCrossedRowsOverFindable",&fTreeVariableLeastRatioCrossedRowsOverFindable,"fTreeVariableLeastRatioCrossedRowsOverFindable/F");
/*22*/   fTree->Branch("fTreeVariablePID",&fTreeVariablePID,"fTreeVariablePID/I");
/*23*/   fTree->Branch("fTreeVariablePIDPositive",&fTreeVariablePIDPositive,"fTreeVariablePIDPositive/I");
/*24*/   fTree->Branch("fTreeVariablePIDNegative",&fTreeVariablePIDNegative,"fTreeVariablePIDNegative/I");
/*25*/   fTree->Branch("fTreeVariablePIDMother",&fTreeVariablePIDMother,"fTreeVariablePIDMother/I");
/*26*/   fTree->Branch("fTreeVariablePtXiMother",&fTreeVariablePtMother,"fTreeVariablePtMother/F");
/*27*/   fTree->Branch("fTreeVariableV0CosineOfPointingAngle",&fTreeVariableV0CosineOfPointingAngle,"fTreeVariableV0CosineOfPointingAngle/F");
//-----------MULTIPLICITY-INFO--------------------
/*28*/   fTree->Branch("fTreeVariableMultiplicity",&fTreeVariableMultiplicity,"fTreeVariableMultiplicity/I");
//------------------------------------------------
/*29*/   fTree->Branch("fTreeVariableDistOverTotMom",&fTreeVariableDistOverTotMom,"fTreeVariableDistOverTotMom/F");
/*30*/   fTree->Branch("fTreeVariableNSigmasPosProton",&fTreeVariableNSigmasPosProton,"fTreeVariableNSigmasPosProton/F");
/*31*/   fTree->Branch("fTreeVariableNSigmasPosPion",&fTreeVariableNSigmasPosPion,"fTreeVariableNSigmasPosPion/F");
/*32*/   fTree->Branch("fTreeVariableNSigmasNegProton",&fTreeVariableNSigmasNegProton,"fTreeVariableNSigmasNegProton/F");
/*33*/   fTree->Branch("fTreeVariableNSigmasNegPion",&fTreeVariableNSigmasNegPion,"fTreeVariableNSigmasNegPion/F");
//------------------------------------------------
/*34*/   fTree->Branch("fTreeVariableNegEta",&fTreeVariableNegEta,"fTreeVariableNegEta/F");
/*35*/   fTree->Branch("fTreeVariablePosEta",&fTreeVariablePosEta,"fTreeVariablePosEta/F");
/*36*/   fTree->Branch("fTreeVariableV0CreationRadius",&fTreeVariableV0CreationRadius,"fTreeVariableV0CreationRadius/F");
/*37*/   fTree->Branch("fTreeVariableIndexStatus",&fTreeVariableIndexStatus,"fTreeVariableIndexStatus/I");
/*38*/   fTree->Branch("fTreeVariableIndexStatusMother",&fTreeVariableIndexStatusMother,"fTreeVariableIndexStatusMother/I");

//------------------------------------------------
// Particle Identification Setup
//------------------------------------------------

   AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
   AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
   fPIDResponse = inputHandler->GetPIDResponse();

//------------------------------------------------
// V0 Multiplicity Histograms
//------------------------------------------------

   // Create histograms
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

//------------------------------------------------
// Generated Particle Histograms
//------------------------------------------------

   Int_t lCustomNBins = 200; 
   Double_t lCustomPtUpperLimit = 20; 
   Int_t lCustomNBinsMultiplicity = 100;

//----------------------------------
// Raw Generated (Pre-physics-selection)
//----------------------------------

//--- 3D Histo (Pt, Y, Multiplicity)  

   if(! f3dHistPrimRawPtVsYVsMultLambda) {
      f3dHistPrimRawPtVsYVsMultLambda = new TH3F( "f3dHistPrimRawPtVsYVsMultLambda", "Pt_{lambda} Vs Y_{#Lambda} Vs Multiplicity; Pt_{lambda} (GeV/c); Y_{#Lambda} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimRawPtVsYVsMultLambda);
   }
   if(! f3dHistPrimRawPtVsYVsMultAntiLambda) {
      f3dHistPrimRawPtVsYVsMultAntiLambda = new TH3F( "f3dHistPrimRawPtVsYVsMultAntiLambda", "Pt_{antilambda} Vs Y_{#Lambda} Vs Multiplicity; Pt_{antilambda} (GeV/c); Y_{#Lambda} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimRawPtVsYVsMultAntiLambda);
   }
   if(! f3dHistPrimRawPtVsYVsMultK0Short) {
      f3dHistPrimRawPtVsYVsMultK0Short = new TH3F( "f3dHistPrimRawPtVsYVsMultK0Short", "Pt_{K0S} Vs Y_{K0S} Vs Multiplicity; Pt_{K0S} (GeV/c); Y_{K0S} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimRawPtVsYVsMultK0Short);
   }

//--- 3D Histo (Pt, Y, Proper Decay Length)

   if(! f3dHistPrimRawPtVsYVsDecayLengthLambda) {
      f3dHistPrimRawPtVsYVsDecayLengthLambda = new TH3F( "f3dHistPrimRawPtVsYVsDecayLengthLambda", "Pt_{lambda} Vs Y_{#Lambda} Vs DecayLength; Pt_{lambda} (GeV/c); Y_{#Lambda} ; DecayLength", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,200,0,50);
      fListHistV0->Add(f3dHistPrimRawPtVsYVsDecayLengthLambda);
   }
   if(! f3dHistPrimRawPtVsYVsDecayLengthAntiLambda) {
      f3dHistPrimRawPtVsYVsDecayLengthAntiLambda = new TH3F( "f3dHistPrimRawPtVsYVsDecayLengthAntiLambda", "Pt_{antilambda} Vs Y_{#Lambda} Vs DecayLength; Pt_{antilambda} (GeV/c); Y_{#Lambda} ; DecayLength", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,200,0,50);
      fListHistV0->Add(f3dHistPrimRawPtVsYVsDecayLengthAntiLambda);
   }
   if(! f3dHistPrimRawPtVsYVsDecayLengthK0Short) {
      f3dHistPrimRawPtVsYVsDecayLengthK0Short = new TH3F( "f3dHistPrimRawPtVsYVsDecayLengthK0Short", "Pt_{K0S} Vs Y_{K0S} Vs DecayLength; Pt_{K0S} (GeV/c); Y_{K0S} ; DecayLength", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,200,0,50);
      fListHistV0->Add(f3dHistPrimRawPtVsYVsDecayLengthK0Short);
   }

//--- 3D Histo (Pt, Y, Multiplicity) for generated charged Xi (feeddown)

   if(! f3dHistGenPtVsYVsMultXiMinus) {
      f3dHistGenPtVsYVsMultXiMinus = new TH3F( "f3dHistGenPtVsYVsMultXiMinus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistGenPtVsYVsMultXiMinus);
   }
   if(! f3dHistGenPtVsYVsMultXiPlus) {
      f3dHistGenPtVsYVsMultXiPlus = new TH3F( "f3dHistGenPtVsYVsMultXiPlus", "Pt_{#Xi} Vs Y_{#Xi} Vs Multiplicity; Pt_{cascade} (GeV/c); Y_{#Xi} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistGenPtVsYVsMultXiPlus);
   }

//----------------------------------
// Histos at analysis level 
//----------------------------------

   if(! f3dHistPrimAnalysisPtVsYVsMultLambda) {
      f3dHistPrimAnalysisPtVsYVsMultLambda = new TH3F( "f3dHistPrimAnalysisPtVsYVsMultLambda", "Pt_{lambda} Vs Y_{#Lambda} Vs Multiplicity; Pt_{lambda} (GeV/c); Y_{#Lambda} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimAnalysisPtVsYVsMultLambda);
   }
   if(! f3dHistPrimAnalysisPtVsYVsMultAntiLambda) {
      f3dHistPrimAnalysisPtVsYVsMultAntiLambda = new TH3F( "f3dHistPrimAnalysisPtVsYVsMultAntiLambda", "Pt_{antilambda} Vs Y_{#Lambda} Vs Multiplicity; Pt_{antilambda} (GeV/c); Y_{#Lambda} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimAnalysisPtVsYVsMultAntiLambda);
   }
   if(! f3dHistPrimAnalysisPtVsYVsMultK0Short) {
      f3dHistPrimAnalysisPtVsYVsMultK0Short = new TH3F( "f3dHistPrimAnalysisPtVsYVsMultK0Short", "Pt_{K0S} Vs Y_{K0S} Vs Multiplicity; Pt_{K0S} (GeV/c); Y_{K0S} ; Mult", lCustomNBins, 0., lCustomPtUpperLimit, 48, -1.2,1.2,lCustomNBinsMultiplicity,0,lCustomNBinsMultiplicity);
      fListHistV0->Add(f3dHistPrimAnalysisPtVsYVsMultK0Short);
   }

//----------------------------------
// Primary Vertex Position Histos
//----------------------------------

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

   if(! fHistPVxAnalysisHasHighPtLambda) {
         fHistPVxAnalysisHasHighPtLambda = new TH1F("fHistPVxAnalysisHasHighPtLambda", 
            "PV x position;Nbr of Evts;x", 
            2000, -0.5, 0.5);       
      fListHistV0->Add(fHistPVxAnalysisHasHighPtLambda);
   }
   if(! fHistPVyAnalysisHasHighPtLambda) {
         fHistPVyAnalysisHasHighPtLambda = new TH1F("fHistPVyAnalysisHasHighPtLambda", 
            "PV y position;Nbr of Evts;y", 
            2000, -0.5, 0.5);       
      fListHistV0->Add(fHistPVyAnalysisHasHighPtLambda);
   }
   if(! fHistPVzAnalysisHasHighPtLambda) {
         fHistPVzAnalysisHasHighPtLambda = new TH1F("fHistPVzAnalysisHasHighPtLambda", 
            "PV z position;Nbr of Evts;z", 
            400, -20, 20);       
      fListHistV0->Add(fHistPVzAnalysisHasHighPtLambda);
   }
   //List of Histograms: Normal
   PostData(1, fListHistV0);

   //TTree Object: Saved to base directory. Should cache to disk while saving. 
   //(Important to avoid excessive memory usage, particularly when merging)
   PostData(2, fTree);

}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskExtractPerformanceV0::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

   AliESDEvent *lESDevent = 0x0;
   AliMCEvent  *lMCevent  = 0x0; 
   AliStack    *lMCstack  = 0x0; 

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
        
//------------------------------------------------
// Multiplicity Information Acquistion
//------------------------------------------------

  //REVISED multiplicity estimator after 'multiplicity day' (2011)
   Int_t lMultiplicity = AliESDtrackCuts::GetReferenceMultiplicity( lESDevent );

   //---> If this is a nuclear collision, then go nuclear on "multiplicity" variable...
   //---> Warning: Experimental
   if(fkIsNuclear == kTRUE){ 
      AliCentrality* centrality;
      centrality = lESDevent->GetCentrality();
      lMultiplicity = ( ( Int_t ) ( centrality->GetCentralityPercentile( "V0M" ) ) );
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

   fHistV0MultiplicityBeforeTrigSel->Fill ( lESDevent->GetNumberOfV0s() );
   fHistMultiplicityBeforeTrigSel->Fill ( lMultiplicity );
        
//------------------------------------------------
// MC Information Acquistion
//------------------------------------------------

   Int_t iNumberOfPrimaries = -1;
   iNumberOfPrimaries = lMCstack->GetNprimary();
   if(iNumberOfPrimaries < 1) return; 
   Bool_t lHasHighPtLambda = kFALSE;

//------------------------------------------------
// Variable Definition
//------------------------------------------------

   Int_t lNbMCPrimary        = 0;

   Int_t lPdgcodeCurrentPart = 0;
   Double_t lRapCurrentPart  = 0;
   Double_t lPtCurrentPart   = 0;
  
   //Int_t lComeFromSigma      = 0;

   // current mc particle 's mother
   //Int_t iCurrentMother  = 0;
   lNbMCPrimary = lMCstack->GetNprimary();

//------------------------------------------------
// Pre-Physics Selection
//------------------------------------------------

//----- Loop on primary Xi, Omega --------------------------------------------------------------
   for (Int_t iCurrentLabelStack = 0; iCurrentLabelStack < lNbMCPrimary; iCurrentLabelStack++) 
   {// This is the begining of the loop on primaries
      
      TParticle* lCurrentParticlePrimary = 0x0; 
      lCurrentParticlePrimary = lMCstack->Particle( iCurrentLabelStack );
      if(!lCurrentParticlePrimary){
         Printf("Cascade loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iCurrentLabelStack );
         continue;
      }
      if ( TMath::Abs(lCurrentParticlePrimary->GetPdgCode()) == 3312 ){ 
         Double_t lRapXiMCPrimary = 0.5*TMath::Log((lCurrentParticlePrimary->Energy() + lCurrentParticlePrimary->Pz()) / (lCurrentParticlePrimary->Energy() - lCurrentParticlePrimary->Pz() +1.e-13));

         //=================================================================================
         // Xi Histograms for Feeddown - Primary Charged Xis
         if( lCurrentParticlePrimary->GetPdgCode() == 3312 ){ 
            lPtCurrentPart    = lCurrentParticlePrimary->Pt();
            f3dHistGenPtVsYVsMultXiMinus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicity);
         }
         if( lCurrentParticlePrimary->GetPdgCode() == -3312 ){ 
            lPtCurrentPart    = lCurrentParticlePrimary->Pt();
            f3dHistGenPtVsYVsMultXiPlus->Fill(lPtCurrentPart, lRapXiMCPrimary, lMultiplicity);
         }
      } 
   }
//----- End Loop on primary Xi, Omega ----------------------------------------------------------

//----- Loop on Lambda, K0Short ----------------------------------------------------------------
   for (Int_t iCurrentLabelStack = 0;  iCurrentLabelStack < (lMCstack->GetNtrack()); iCurrentLabelStack++) 
   {// This is the begining of the loop on tracks
      
      TParticle* lCurrentParticleForLambdaCheck = 0x0; 
      lCurrentParticleForLambdaCheck = lMCstack->Particle( iCurrentLabelStack );
      if(!lCurrentParticleForLambdaCheck){
         Printf("V0s loop %d - MC TParticle pointer to current stack particle = 0x0 ! Skip ...\n", iCurrentLabelStack );
         continue;
      }

      //=================================================================================
      //Single-Strange checks
      // Keep only K0s, Lambda and AntiLambda:
      lPdgcodeCurrentPart = lCurrentParticleForLambdaCheck->GetPdgCode();	      

      if ( (lCurrentParticleForLambdaCheck->GetPdgCode() == 310   ) ||
           (lCurrentParticleForLambdaCheck->GetPdgCode() == 3122  ) ||
           (lCurrentParticleForLambdaCheck->GetPdgCode() == -3122 ) )
	   {
         lRapCurrentPart   = MyRapidity(lCurrentParticleForLambdaCheck->Energy(),lCurrentParticleForLambdaCheck->Pz());
         lPtCurrentPart    = lCurrentParticleForLambdaCheck->Pt();

         //Use Physical Primaries only for filling PrimRaw Histograms!
         if ( lMCstack->IsPhysicalPrimary(iCurrentLabelStack)!=kTRUE ) continue;

         if( lPdgcodeCurrentPart == 3122 ){
            f3dHistPrimRawPtVsYVsMultLambda->Fill(lPtCurrentPart, lRapCurrentPart, lMultiplicity);
            if( TMath::Abs( lCurrentParticleForLambdaCheck->Eta() )<1.2 && lPtCurrentPart>2 ){
               lHasHighPtLambda = kTRUE; //Keep track of events with Lambda within |eta|<1.2 and pt>2
            }
         }
         if( lPdgcodeCurrentPart == -3122 ){
            f3dHistPrimRawPtVsYVsMultAntiLambda->Fill(lPtCurrentPart, lRapCurrentPart, lMultiplicity);
         }
         if( lPdgcodeCurrentPart == 310 ){
            f3dHistPrimRawPtVsYVsMultK0Short->Fill(lPtCurrentPart, lRapCurrentPart, lMultiplicity);
         }
         //Decay Length Acquisition=====================================================
         Double_t decaylength = -1; 
         Double_t lV0Mass = -1; 
          
         if( !(lCurrentParticleForLambdaCheck->GetDaughter(0) < 0) ) {
            TParticle* lDght0ofV0 = lMCstack->Particle(  lCurrentParticleForLambdaCheck->GetDaughter(0) ); //get first daughter
            if(lDght0ofV0){ // skip if not defined. 
               decaylength = TMath::Sqrt(
				        TMath::Power( lCurrentParticleForLambdaCheck->Vx() - lDght0ofV0->Vx() , 2) + 
				        TMath::Power( lCurrentParticleForLambdaCheck->Vy() - lDght0ofV0->Vy() , 2) + 
				        TMath::Power( lCurrentParticleForLambdaCheck->Vz() - lDght0ofV0->Vz() , 2)
               );
               //Need to correct for relativitity! Involves multiplying by mass and dividing by momentum. 
               if(TMath::Abs( lPdgcodeCurrentPart ) == 3122 ) { lV0Mass = 1.115683; }
               if(TMath::Abs( lPdgcodeCurrentPart ) == 310 ) { lV0Mass = 0.497614; }
               decaylength = ( lV0Mass * decaylength ) / ( lCurrentParticleForLambdaCheck->P() + 1e-10 );
            }
         }
         if( lPdgcodeCurrentPart == 3122) f3dHistPrimRawPtVsYVsDecayLengthLambda ->Fill( lPtCurrentPart, lRapCurrentPart , decaylength ); 
         if( lPdgcodeCurrentPart == -3122) f3dHistPrimRawPtVsYVsDecayLengthAntiLambda ->Fill( lPtCurrentPart, lRapCurrentPart , decaylength ); 
         if( lPdgcodeCurrentPart == 310) f3dHistPrimRawPtVsYVsDecayLengthK0Short ->Fill( lPtCurrentPart, lRapCurrentPart , decaylength ); 
      }
   }//End of loop on tracks
//----- End Loop on Lambda, K0Short ------------------------------------------------------------

// ---> Set Variables to Zero again
// ---> Variable Definition

   lPdgcodeCurrentPart = 0;
   lRapCurrentPart  = 0;
   lPtCurrentPart   = 0;

//------------------------------------------------
// Physics Selection
//------------------------------------------------

   UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
   Bool_t isSelected = 0;
   isSelected = (maskIsSelected & AliVEvent::kMB) == AliVEvent::kMB;

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

   lNumberOfV0s          = lESDevent->GetNumberOfV0s();
  
   //Set variable for filling tree afterwards!
   fHistV0MultiplicityForTrigEvt->Fill(lNumberOfV0s);
   fHistMultiplicityForTrigEvt->Fill ( lMultiplicity );

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
   fHistPVx->Fill( lPrimaryVtxPosition[0] );
   fHistPVy->Fill( lPrimaryVtxPosition[1] );
   fHistPVz->Fill( lPrimaryVtxPosition[2] );

//------------------------------------------------
// Primary Vertex Z position: SKIP
//------------------------------------------------

   if(TMath::Abs(lBestPrimaryVtxPos[2]) > 10.0 ) { 
      AliWarning("Pb / | Z position of Best Prim Vtx | > 10.0 cm ... return !"); 
      PostData(1, fListHistV0);
      PostData(2, fTree);
      return; 
   }

   lMagneticField = lESDevent->GetMagneticField( );
   fHistV0MultiplicityForSelEvt ->Fill( lNumberOfV0s );
   fHistMultiplicity->Fill(lMultiplicity);

//------------------------------------------------
// SKIP: Events with well-established PVtx
//------------------------------------------------
	
   const AliESDVertex *lPrimaryTrackingESDVtxCheck = lESDevent->GetPrimaryVertexTracks();
   const AliESDVertex *lPrimarySPDVtx = lESDevent->GetPrimaryVertexSPD();
   if (!lPrimarySPDVtx->GetStatus() && !lPrimaryTrackingESDVtxCheck->GetStatus() ){
      AliWarning("Pb / No SPD prim. vertex nor prim. Tracking vertex ... return !");
      PostData(1, fListHistV0);
      PostData(2, fTree);
      return;
   }
   fHistV0MultiplicityForSelEvtNoTPCOnly ->Fill( lNumberOfV0s );
   fHistMultiplicityNoTPCOnly->Fill(lMultiplicity);

//------------------------------------------------
// Pileup Rejection Studies
//------------------------------------------------

   // FIXME : quality selection regarding pile-up rejection 
   if(lESDevent->IsPileupFromSPD() ){// minContributors=3, minZdist=0.8, nSigmaZdist=3., nSigmaDiamXY=2., nSigmaDiamZ=5.  -> see http://alisoft.cern.ch/viewvc/trunk/STEER/AliESDEvent.h?root=AliRoot&r1=41914&r2=42199&pathrev=42199
      AliWarning("Pb / This is tagged as Pileup from SPD... return !");
      PostData(1, fListHistV0);
      PostData(2, fTree);
      return;
   }
   fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup ->Fill( lNumberOfV0s );
   fHistMultiplicityNoTPCOnlyNoPileup->Fill(lMultiplicity);

   //Do control histograms without the IsFromVertexerZ events, but consider them in analysis...
   if( ! (lESDevent->GetPrimaryVertex()->IsFromVertexerZ() )	 ){ 
      fHistPVxAnalysis->Fill( lPrimaryVtxPosition[0] );
      fHistPVyAnalysis->Fill( lPrimaryVtxPosition[1] );
      fHistPVzAnalysis->Fill( lPrimaryVtxPosition[2] );
      if ( lHasHighPtLambda == kTRUE ){ 
         fHistPVxAnalysisHasHighPtLambda->Fill( lPrimaryVtxPosition[0] );
         fHistPVyAnalysisHasHighPtLambda->Fill( lPrimaryVtxPosition[1] );
         fHistPVzAnalysisHasHighPtLambda->Fill( lPrimaryVtxPosition[2] );
      }
   }

//------------------------------------------------
// stack loop starts here
//------------------------------------------------

//---> Loop over ALL PARTICLES
 
   for (Int_t iMc = 0; iMc < (lMCstack->GetNtrack()); iMc++) {  
      TParticle *p0 = lMCstack->Particle(iMc); 
      if (!p0) {
         //Printf("ERROR: particle with label %d not found in lMCstack (mc loop)", iMc);
         continue;
      }
      lPdgcodeCurrentPart = p0->GetPdgCode();

      // Keep only K0s, Lambda and AntiLambda:
      if ( (lPdgcodeCurrentPart != 310 ) && (lPdgcodeCurrentPart != 3122 ) && (lPdgcodeCurrentPart != -3122 ) ) continue;
	
      lRapCurrentPart   = MyRapidity(p0->Energy(),p0->Pz());
      lPtCurrentPart    = p0->Pt();

        //Use Physical Primaries only for filling PrimRaw Histograms!
      if ( lMCstack->IsPhysicalPrimary(iMc)!=kTRUE ) continue;

      if( lPdgcodeCurrentPart == 3122 ){
         f3dHistPrimAnalysisPtVsYVsMultLambda->Fill(lPtCurrentPart, lRapCurrentPart, lMultiplicity);
      }
      if( lPdgcodeCurrentPart == -3122 ){
         f3dHistPrimAnalysisPtVsYVsMultAntiLambda->Fill(lPtCurrentPart, lRapCurrentPart, lMultiplicity);
      }
      if( lPdgcodeCurrentPart == 310 ){
         f3dHistPrimAnalysisPtVsYVsMultK0Short->Fill(lPtCurrentPart, lRapCurrentPart, lMultiplicity);
      }
   }

//------------------------------------------------
// MAIN LAMBDA LOOP STARTS HERE
//------------------------------------------------

   //Variable definition
   Int_t    lOnFlyStatus = 0;
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
   
   for (Int_t iV0 = 0; iV0 < nv0s; iV0++) 
	{// This is the begining of the V0 loop
      AliESDv0 *v0 = ((AliESDEvent*)lESDevent)->GetV0(iV0);
      if (!v0) continue;

      Double_t tV0mom[3];
      v0->GetPxPyPz( tV0mom[0],tV0mom[1],tV0mom[2] ); 
      Double_t lV0TotalMomentum = TMath::Sqrt(
         tV0mom[0]*tV0mom[0]+tV0mom[1]*tV0mom[1]+tV0mom[2]*tV0mom[2] );

      Double_t tDecayVertexV0[3]; v0->GetXYZ(tDecayVertexV0[0],tDecayVertexV0[1],tDecayVertexV0[2]); 
      lV0Radius = TMath::Sqrt(tDecayVertexV0[0]*tDecayVertexV0[0]+tDecayVertexV0[1]*tDecayVertexV0[1]);
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

      fTreeVariableNegEta = nTrack->Eta();
      fTreeVariablePosEta = pTrack->Eta();

      // Filter like-sign V0 (next: add counter and distribution)
      if ( pTrack->GetSign() == nTrack->GetSign()){
         continue;
      } 

      //________________________________________________________________________
      // Track quality cuts 
      Float_t lPosTrackCrossedRows = pTrack->GetTPCClusterInfo(2,1);
      Float_t lNegTrackCrossedRows = nTrack->GetTPCClusterInfo(2,1);
      fTreeVariableLeastNbrCrossedRows = lPosTrackCrossedRows;
      if( lNegTrackCrossedRows < fTreeVariableLeastNbrCrossedRows )
         fTreeVariableLeastNbrCrossedRows = lNegTrackCrossedRows;

      // TPC refit condition (done during reconstruction for Offline but not for On-the-fly)
      if( !(pTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;      
      if( !(nTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;

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

      //End track Quality Cuts
      //________________________________________________________________________

      lDcaPosToPrimVertex = TMath::Abs(pTrack->GetD(lPrimaryVtxPosition[0],
							lPrimaryVtxPosition[1],
							lMagneticField) );

      lDcaNegToPrimVertex = TMath::Abs(nTrack->GetD(lPrimaryVtxPosition[0],
							lPrimaryVtxPosition[1],
							lMagneticField) );

      lOnFlyStatus = v0->GetOnFlyStatus();
      lChi2V0 = v0->GetChi2V0();
      lDcaV0Daughters = v0->GetDcaV0Daughters();
      lDcaV0ToPrimVertex = v0->GetD(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1],lPrimaryVtxPosition[2]);
      lV0CosineOfPointingAngle = v0->GetV0CosineOfPointingAngle(lPrimaryVtxPosition[0],lPrimaryVtxPosition[1],lPrimaryVtxPosition[2]);
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

      //fTreeVariableOnFlyStatus = lOnFlyStatus;
      //fHistV0OnFlyStatus->Fill(lOnFlyStatus);

//===============================================
// Monte Carlo Association starts here
//===============================================

      //---> Set Everything to "I don't know" before starting

      fTreeVariablePIDPositive = 0;
      fTreeVariablePIDNegative = 0;

      fTreeVariableIndexStatus = 0;
      fTreeVariableIndexStatusMother = 0;

      fTreeVariablePtMother = -1;
      fTreeVariablePtMC = -1;
      fTreeVariableRapMC = -100;

      fTreeVariablePID = -1; 
      fTreeVariablePIDMother = -1;

      fTreeVariablePrimaryStatus = 0; 
      fTreeVariablePrimaryStatusMother = 0; 
      fTreeVariableV0CreationRadius = -1;

      Int_t lblPosV0Dghter = (Int_t) TMath::Abs( pTrack->GetLabel() );  
      Int_t lblNegV0Dghter = (Int_t) TMath::Abs( nTrack->GetLabel() );
		
      TParticle* mcPosV0Dghter = lMCstack->Particle( lblPosV0Dghter );
      TParticle* mcNegV0Dghter = lMCstack->Particle( lblNegV0Dghter );
	    
      fTreeVariablePosTransvMomentumMC = mcPosV0Dghter->Pt();
      fTreeVariableNegTransvMomentumMC = mcNegV0Dghter->Pt();

      Int_t lPIDPositive = mcPosV0Dghter -> GetPdgCode();
      Int_t lPIDNegative = mcNegV0Dghter -> GetPdgCode();

      fTreeVariablePIDPositive = lPIDPositive;
      fTreeVariablePIDNegative = lPIDNegative;

      Int_t lblMotherPosV0Dghter = mcPosV0Dghter->GetFirstMother() ; 
      Int_t lblMotherNegV0Dghter = mcNegV0Dghter->GetFirstMother();

      if( lblMotherPosV0Dghter == lblMotherNegV0Dghter && lblMotherPosV0Dghter > -1 ){
         //either label is fine, they're equal at this stage
         TParticle* pThisV0 = lMCstack->Particle( lblMotherPosV0Dghter ); 
         //Set tree variables
         fTreeVariablePID   = pThisV0->GetPdgCode(); //PDG Code
         fTreeVariablePtMC  = pThisV0->Pt(); //Perfect Pt
         //Only Interested if it's a Lambda, AntiLambda or K0s 
         //Avoid the Junction Bug! PYTHIA has particles with Px=Py=Pz=E=0 occasionally, 
         //having particle code 88 (unrecognized by PDG), for documentation purposes.
         //Even ROOT's TParticle::Y() is not prepared to deal with that exception!
         //Note that TParticle::Pt() is immune (that would just return 0)...
         //Though granted that that should be extremely rare in this precise condition...
         if( TMath::Abs(fTreeVariablePID) == 3122 || fTreeVariablePID==310 ){
            fTreeVariableRapMC = pThisV0->Y(); //Perfect Y
         }
         fTreeVariableV0CreationRadius = pThisV0->R(); // Creation Radius
         if( lblMotherPosV0Dghter  < lNbMCPrimary ) fTreeVariableIndexStatus = 1; //looks primary
         if( lblMotherPosV0Dghter >= lNbMCPrimary ) fTreeVariableIndexStatus = 2; //looks secondary
         if( lMCstack->IsPhysicalPrimary       (lblMotherPosV0Dghter) ) fTreeVariablePrimaryStatus = 1; //Is Primary!
         if( lMCstack->IsSecondaryFromWeakDecay(lblMotherPosV0Dghter) ) fTreeVariablePrimaryStatus = 2; //Weak Decay!
         if( lMCstack->IsSecondaryFromMaterial (lblMotherPosV0Dghter) ) fTreeVariablePrimaryStatus = 3; //Material Int!
         
         //Now we try to acquire the V0 parent particle, if possible
         Int_t lblThisV0Parent = pThisV0->GetFirstMother();
         if ( lblThisV0Parent > -1 ){ //if it has a parent, get it and store specs
            TParticle* pThisV0Parent = lMCstack->Particle( lblThisV0Parent );
            fTreeVariablePIDMother   = pThisV0Parent->GetPdgCode(); //V0 Mother PDG
            fTreeVariablePtMother    = pThisV0Parent->Pt();         //V0 Mother Pt
            //Primary Status for the V0 Mother particle 
            if( lblThisV0Parent  < lNbMCPrimary ) fTreeVariableIndexStatusMother = 1; //looks primary
            if( lblThisV0Parent >= lNbMCPrimary ) fTreeVariableIndexStatusMother = 2; //looks secondary
            if( lMCstack->IsPhysicalPrimary       (lblThisV0Parent) ) fTreeVariablePrimaryStatusMother = 1; //Is Primary!
            if( lMCstack->IsSecondaryFromWeakDecay(lblThisV0Parent) ) fTreeVariablePrimaryStatusMother = 2; //Weak Decay!
            if( lMCstack->IsSecondaryFromMaterial (lblThisV0Parent) ) fTreeVariablePrimaryStatusMother = 3; //Material Int!
         }
      }

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

//tDecayVertexV0[0],tDecayVertexV0[1],tDecayVertexV0[2]
      fTreeVariableDistOverTotMom = TMath::Sqrt(
						TMath::Power( tDecayVertexV0[0] - lBestPrimaryVtxPos[0] , 2) +
						TMath::Power( tDecayVertexV0[1] - lBestPrimaryVtxPos[1] , 2) +
						TMath::Power( tDecayVertexV0[2] - lBestPrimaryVtxPos[2] , 2)
					);
      fTreeVariableDistOverTotMom /= (lV0TotalMomentum + 1e-10); //avoid division by zero, to be sure

      Double_t lMomentumPosTemp[3];
      pTrack->GetPxPyPz(lMomentumPosTemp);
      Double_t lPtPosTemporary = sqrt(pow(lMomentumPosTemp[0],2) + pow(lMomentumPosTemp[1],2));

      Double_t lMomentumNegTemp[3];
      nTrack->GetPxPyPz(lMomentumNegTemp);
      Double_t lPtNegTemporary = sqrt(pow(lMomentumNegTemp[0],2) + pow(lMomentumNegTemp[1],2));

      fTreeVariablePosTransvMomentum = lPtPosTemporary;
      fTreeVariableNegTransvMomentum = lPtNegTemporary;


//------------------------------------------------
// Fill Tree! 
//------------------------------------------------

      // The conditionals are meant to decrease excessive
      // memory usage! 

      //Modified version: Keep only OnFlyStatus == 0
      //Keep only if included in a parametric InvMass Region 20 sigmas away from peak

      //First Selection: Reject OnFly
      if( (lOnFlyStatus == 0 && fkUseOnTheFly == kFALSE) || (lOnFlyStatus != 1 && fkUseOnTheFly == kTRUE ) ){
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
void AliAnalysisTaskExtractPerformanceV0::Terminate(Option_t *)
{
   // Draw result to the screen
   // Called once at the end of the query

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

//----------------------------------------------------------------------------

Double_t AliAnalysisTaskExtractPerformanceV0::MyRapidity(Double_t rE, Double_t rPz) const
{
   // Local calculation for rapidity
   return 0.5*TMath::Log((rE+rPz)/(rE-rPz+1.e-13));
} 
//----------------------------------------------------------------------------

