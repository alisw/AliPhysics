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
#include "AliAnalysisTaskExtractV0pPb.h"

//debugging purposes
#include "TObjectTable.h"

ClassImp(AliAnalysisTaskExtractV0pPb)

AliAnalysisTaskExtractV0pPb::AliAnalysisTaskExtractV0pPb() 
  : AliAnalysisTaskSE(), fListHistV0(0), fTree(0), fTreeEvents(0), fPIDResponse(0),
  fTPCdEdxSelection ( kTRUE ),
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
	fTreeVariableCentrality(0),

//Event Tree
  fTreeEventsCentrality(0), 

//------------------------------------------------
// HISTOGRAMS
// --- Filled on an Event-by-event basis
//------------------------------------------------
  //V0A Centrality
  fHistCentralityProcessed(0),
  fHistCentralityTrigEvt(0),
  fHistCentralityHasVtx(0),
  fHistCentralityVtxZ(0)
{
  // Dummy Constructor
}

AliAnalysisTaskExtractV0pPb::AliAnalysisTaskExtractV0pPb(const char *name) 
  : AliAnalysisTaskSE(name), fListHistV0(0), fTree(0), fTreeEvents(0), fPIDResponse(0),
  fTPCdEdxSelection ( kTRUE ),
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
	fTreeVariableCentrality(0),

//Event Tree
  fTreeEventsCentrality(0), 

//------------------------------------------------
// HISTOGRAMS
// --- Filled on an Event-by-event basis
//------------------------------------------------
  //V0A Centrality
  fHistCentralityProcessed(0),
  fHistCentralityTrigEvt(0),
  fHistCentralityHasVtx(0),
  fHistCentralityVtxZ(0)
{
  // Output slot #0 writes into a TList container (Lambda Histos and fTree)
   DefineOutput(1, TList::Class());
   DefineOutput(2, TTree::Class());
   DefineOutput(3, TTree::Class());
}


AliAnalysisTaskExtractV0pPb::~AliAnalysisTaskExtractV0pPb()
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
   if (fTreeEvents){
      delete fTreeEvents;
      fTree = 0x0;
   }
}



//________________________________________________________________________
void AliAnalysisTaskExtractV0pPb::UserCreateOutputObjects()
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
/*17*/	fTree->Branch("fTreeVariableCentrality",&fTreeVariableCentrality,"fTreeVariableCentrality/F");
//------------------------------------------------
/*18*/	fTree->Branch("fTreeVariableDistOverTotMom",&fTreeVariableDistOverTotMom,"fTreeVariableDistOverTotMom/F");
/*19*/	fTree->Branch("fTreeVariableNSigmasPosProton",&fTreeVariableNSigmasPosProton,"fTreeVariableNSigmasPosProton/F");
/*20*/	fTree->Branch("fTreeVariableNSigmasPosPion",&fTreeVariableNSigmasPosPion,"fTreeVariableNSigmasPosPion/F");
/*21*/	fTree->Branch("fTreeVariableNSigmasNegProton",&fTreeVariableNSigmasNegProton,"fTreeVariableNSigmasNegProton/F");
/*22*/	fTree->Branch("fTreeVariableNSigmasNegPion",&fTreeVariableNSigmasNegPion,"fTreeVariableNSigmasNegPion/F");
/*23*/	fTree->Branch("fTreeVariableNegEta",&fTreeVariableNegEta,"fTreeVariableNegEta/F");
/*24*/	fTree->Branch("fTreeVariablePosEta",&fTreeVariablePosEta,"fTreeVariablePosEta/F");
  
//------------------------------------------------
// fTreeEvents Branch definitions
//------------------------------------------------

   fTreeEvents = new TTree("fTreeEvents","Events");
/*01*/	fTree->Branch("fTreeEventsCentrality",&fTreeEventsCentrality,"fTreeEventsCentrality/F");

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
   //Create File-resident Tree, please.
   OpenFile(1);

   fListHistV0 = new TList();
   fListHistV0->SetOwner();  // See http://root.cern.ch/root/html/TCollection.html#TCollection:SetOwner

//------------------------------------------------
// V0A Centrality Histograms
//------------------------------------------------

  //Default V0A Centrality (if PbPb) 
   if(! fHistCentralityProcessed) {
      fHistCentralityProcessed = new TH1F("fHistCentralityProcessed", 
         "All processed Events;V0A Centrality;Events", 
         200, 0, 100); 		
      fListHistV0->Add(fHistCentralityProcessed);
   }
   if(! fHistCentralityTrigEvt) {
      fHistCentralityTrigEvt = new TH1F("fHistCentralityTrigEvt", 
         "PS selected Events;V0A Centrality;Events", 
         200, 0, 100); 		
      fListHistV0->Add(fHistCentralityTrigEvt);
   }
   if(! fHistCentralityHasVtx) {
      fHistCentralityHasVtx = new TH1F("fHistCentralityHasVtx", 
         "Events having Vertex;V0A Centrality;Events", 
         200, 0, 100); 		
      fListHistV0->Add(fHistCentralityHasVtx);
   }
   if(! fHistCentralityVtxZ) {
      fHistCentralityVtxZ = new TH1F("fHistCentralityVtxZ", 
         "Vertex |z|<10cm;V0A Centrality;Events", 
         200, 0, 100); 		
      fListHistV0->Add(fHistCentralityVtxZ);
   }
 
   //Regular output: Histograms
   PostData(1, fListHistV0);
   //TTree Object: Saved to base directory. Should cache to disk while saving. 
   //(Important to avoid excessive memory usage, particularly when merging)
   PostData(2, fTree);
   PostData(3, fTreeEvents);

}// end UserCreateOutputObjects


//________________________________________________________________________
void AliAnalysisTaskExtractV0pPb::UserExec(Option_t *) 
{
   // Main loop
   // Called for each event
   //gObjectTable->Print();
   AliESDEvent *lESDevent = 0x0;

   //AliAODEvent *lAODevent = 0x0;

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
  
   //REVISED multiplicity estimator after 'multiplicity day' (2011)
   Float_t lMultiplicity = -100;
   fTreeVariableCentrality = -100;

   //---> If this is a nuclear collision, then go nuclear on "multiplicity" variable...
   //---> Warning: Experimental
   AliCentrality* centrality;
   centrality = lESDevent->GetCentrality();
   if (centrality->GetQuality()>1) {
     PostData(1, fListHistV0);
     PostData(2, fTree);
     PostData(3, fTreeEvents);
     return;
   }
   lMultiplicity = centrality->GetCentralityPercentile( "V0A" );


   //Set variable for filling tree afterwards!
   //---> pp case......: GetReferenceMultiplicity
   //---> Pb-Pb case...: Centrality by V0M
  fTreeVariableCentrality = lMultiplicity;
  fTreeEventsCentrality = lMultiplicity;
 
  fHistCentralityProcessed->Fill ( fTreeVariableCentrality );


//------------------------------------------------
// Physics Selection
//------------------------------------------------
  
  // new method
  UInt_t maskIsSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  Bool_t isSelected = 0;
  isSelected = (maskIsSelected & AliVEvent::kINT7) == AliVEvent::kINT7;
  
  //Trigger Selection for kINT7
  if ( !isSelected ) {
    PostData(1, fListHistV0);
    PostData(2, fTree);
    PostData(3, fTreeEvents);
    return;
  }
  
//------------------------------------------------
// After Trigger Selection
//------------------------------------------------

  fHistCentralityTrigEvt -> Fill( fTreeVariableCentrality );
  
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

  //------------------------------------------------
  // Primary Vertex Requirements Section:
  //  ---> pp and PbPb: Only requires |z|<10cm
  //  ---> pPb: all requirements checked at this stage
  //------------------------------------------------
  
  //Roberto's PV selection criteria, implemented 17th April 2013
  
  /* vertex selection */
  Bool_t fHasVertex = kFALSE;
  
  const AliESDVertex *vertex = lESDevent->GetPrimaryVertexTracks();
  if (vertex->GetNContributors() < 1) {
    vertex = lESDevent->GetPrimaryVertexSPD();
    if (vertex->GetNContributors() < 1) fHasVertex = kFALSE;
    else fHasVertex = kTRUE;
    TString vtxTyp = vertex->GetTitle();
    Double_t cov[6]={0};
    vertex->GetCovarianceMatrix(cov);
    Double_t zRes = TMath::Sqrt(cov[5]);
    if (vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) fHasVertex = kFALSE;
  }
  else fHasVertex = kTRUE;
  
  //Is First event in chunk rejection: Still present!
  if(fHasVertex == kFALSE) {
    AliWarning("Pb / | PV does not satisfy selection criteria!");
    PostData(1, fListHistV0);
    PostData(2, fTree);
    PostData(3, fTreeEvents);
    return;
  }

  fHistCentralityHasVtx -> Fill ( fTreeVariableCentrality );
  
  //17 April Fix: Always do primary vertex Z selection, after pA vertex selection from Roberto
  if(TMath::Abs(lBestPrimaryVtxPos[2]) > 10.0) {
    AliWarning("Pb / | pPb case | Z position of Best Prim Vtx | > 10.0 cm ... return !");
    PostData(1, fListHistV0);
    PostData(2, fTree);
    PostData(3, fTreeEvents);
    return;
  }
  
  fHistCentralityVtxZ -> Fill ( fTreeVariableCentrality );  
  lMagneticField = lESDevent->GetMagneticField( );

//------------------------------------------------
// MAIN LAMBDA LOOP STARTS HERE
//------------------------------------------------

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

      Double_t tDecayVertexV0[3]; v0->GetXYZ(tDecayVertexV0[0],tDecayVertexV0[1],tDecayVertexV0[2]); 

      Double_t tV0mom[3];
      v0->GetPxPyPz( tV0mom[0],tV0mom[1],tV0mom[2] ); 
      Double_t lV0TotalMomentum = TMath::Sqrt(
      tV0mom[0]*tV0mom[0]+tV0mom[1]*tV0mom[1]+tV0mom[2]*tV0mom[2] );

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

      //Daughter Eta for Eta selection, afterwards
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
      fTreeVariableLeastNbrCrossedRows = (Int_t) lPosTrackCrossedRows;
      if( lNegTrackCrossedRows < fTreeVariableLeastNbrCrossedRows )
         fTreeVariableLeastNbrCrossedRows = (Int_t) lNegTrackCrossedRows;

      // TPC refit condition (done during reconstruction for Offline but not for On-the-fly)
      if( !(pTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
      if( !(nTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
     
      //Get status flags
      //fTreeVariablePosTrackStatus = pTrack->GetStatus();
      //fTreeVariableNegTrackStatus = nTrack->GetStatus();
     
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
      //if ( fTreeVariableLeastRatioCrossedRowsOverFindable < 0.8 ) continue;

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
     if( lOnFlyStatus == 0 ){
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
           //If this is a nuclear collision___________________
           // ... pre-filter with TPC, daughter eta selection
           if( (fTreeVariableInvMassLambda     < lUpperLimitLambda  && fTreeVariableInvMassLambda     > lLowerLimitLambda
                && ( fTPCdEdxSelection==kFALSE ||( fTPCdEdxSelection==kTRUE && TMath::Abs(fTreeVariableNSigmasPosProton) < 6.0 && TMath::Abs(fTreeVariableNSigmasNegPion) < 6.0 ) ) ) ||
              (fTreeVariableInvMassAntiLambda < lUpperLimitLambda  && fTreeVariableInvMassAntiLambda > lLowerLimitLambda
                && ( fTPCdEdxSelection==kFALSE ||( fTPCdEdxSelection==kTRUE && TMath::Abs(fTreeVariableNSigmasNegProton) < 6.0 && TMath::Abs(fTreeVariableNSigmasPosPion) < 6.0 ) ) ) ||
              (fTreeVariableInvMassK0s        < lUpperLimitK0Short && fTreeVariableInvMassK0s        > lLowerLimitK0Short
                && ( fTPCdEdxSelection==kFALSE ||( fTPCdEdxSelection==kTRUE && TMath::Abs(fTreeVariableNSigmasNegPion)   < 6.0 && TMath::Abs(fTreeVariableNSigmasPosPion) < 6.0 ) ) ) ){
                //insane test
                if ( TMath::Abs(fTreeVariableNegEta)<0.8 && TMath::Abs(fTreeVariablePosEta)<0.8) fTree->Fill();
              }
         }//end nuclear_____________________________________
       }

//------------------------------------------------
// Fill tree over.
//------------------------------------------------

   }// This is the end of the V0 loop
  
  // Post output data.
  PostData(1, fListHistV0);
  PostData(2, fTree);
  PostData(3, fTreeEvents);
}

//________________________________________________________________________
void AliAnalysisTaskExtractV0pPb::Terminate(Option_t *)
{
   // Draw result to the screen
   // Called once at the end of the query
   // This will draw the V0 candidate multiplicity, whose 
   // number of entries corresponds to the number of triggered events. 
   TList *cRetrievedList = 0x0;
   cRetrievedList = (TList*)GetOutputData(1);
   if(!cRetrievedList){
      Printf("ERROR - AliAnalysisTaskExtractV0pPb : ouput data container list not available\n");
      return;
   }		
   fHistCentralityProcessed = dynamic_cast<TH1F*> (  cRetrievedList->FindObject("fHistCentralityProcessed")  );
   if (!fHistCentralityProcessed) {
      Printf("ERROR - AliAnalysisTaskExtractV0pPb : fHistV0MultiplicityForTrigEvt not available");
      return;
   }
   TCanvas *canCheck = new TCanvas("AliAnalysisTaskExtractV0pPb","V0 Multiplicity",10,10,510,510);
   canCheck->cd(1)->SetLogy();
   fHistCentralityProcessed->SetMarkerStyle(22);
   fHistCentralityProcessed->DrawCopy("E");
}
