/**************************************************************************
 * Author: Michal Meres.                                               *
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

//-----------------------------------------------------------------
//                 AliAnalysisTaskProtonAbsorbtion class
//            This task is for absorption efficiency of protons and antiprotons from ESD only
//-----------------------------------------------------------------
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"

#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TDatabasePDG.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"

#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliAODEvent.h"

#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "TParticle.h"
#include "TMCProcess.h"

#include "AliV0vertexer.h"
#include "AliCascadeVertexer.h"

#include "AliESDv0.h"
#include "AliESDcascade.h"
#include "AliTPCPIDResponse.h"
#include "AliPhysicsSelection.h"
#include "AliPIDResponse.h"
#include "AliMultiplicity.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include "AliProdInfo.h"
#include "AliCFManager.h"
#include "AliCFContainer.h"
#include "AliCFTrackKineCuts.h"
#include "AliCFParticleGenCuts.h"
#include "AliCFAcceptanceCuts.h"

#include "AliAnalysisTaskProtonAbsorbtion.h"

ClassImp(AliAnalysisTaskProtonAbsorbtion)

AliAnalysisTaskProtonAbsorbtion::AliAnalysisTaskProtonAbsorbtion()
  : AliAnalysisTaskSE(), fAnalysisType("ESD"), fCollidingSystems(0), 
    fUsePhysicsSelection(1), fMaxPrimaryVtxPosZ(100.) ,fPhysicsSelection(0),fMultiplicityMode(0), fListHist(0),
     fCFManagerProtonsPrim(0), fCFManagerAntiProtonsPrim(0),
  fCFManagerProtonsPrimMulti(0), fCFManagerAntiProtonsPrimMulti(0),
fMinTPCClusters(0),fMinITSClusters(0), fMaxChi2PerTPCCluster(0.), fMaxDCAXY(0),
    fPIDMode(kSigma), fNBoundP(0), fNSigma1(0), fNSigma2(0), fNRatio1(0), fNRatio2(0),
nbinsPtPrim(6), fLowPtPrim(0.45), fHighPtPrim(1.05), fPIDResponse(0), fOADBPath(0), fList(0), containerProtonsPrim(0), containerAntiProtonsPrim(0), containerProtonsPrimMulti(0), containerAntiProtonsPrimMulti(0), nbinsYPrim(5), fLowYPrim(-0.5), fHighYPrim(0.5), nbinsYMulti(100), fLowYMulti(0), fHighYMulti(100), mintrackrefsTPC(1), MAXCent(100), MINCent(0), fDebugMode(kFALSE), fMaxDCAZ(1.), fMaxChi2PerITSCluster(36), fMaxDCAXYFlag(kTRUE), fMaxDCAZFlag(kFALSE)
{
  // Dummy constructor
}

//________________________________________________________________________
AliAnalysisTaskProtonAbsorbtion::AliAnalysisTaskProtonAbsorbtion(const char *name) 
  : AliAnalysisTaskSE(name), fAnalysisType("ESD"), fCollidingSystems(0), 
    fUsePhysicsSelection(1), fMaxPrimaryVtxPosZ(100.) ,fPhysicsSelection(0),fMultiplicityMode(0), fListHist(0),
     fCFManagerProtonsPrim(0), fCFManagerAntiProtonsPrim(0),
  fCFManagerProtonsPrimMulti(0), fCFManagerAntiProtonsPrimMulti(0),
 fMinTPCClusters(0),fMinITSClusters(0), fMaxChi2PerTPCCluster(0.), fMaxDCAXY(0),
    fPIDMode(kSigma), fNBoundP(0), fNSigma1(0), fNSigma2(0), fNRatio1(0), fNRatio2(0),
nbinsPtPrim(6), fLowPtPrim(0.45), fHighPtPrim(1.05), fPIDResponse(0), fOADBPath(0), fList(0), containerProtonsPrim(0), containerAntiProtonsPrim(0), containerProtonsPrimMulti(0), containerAntiProtonsPrimMulti(0), nbinsYPrim(5), fLowYPrim(-0.5), fHighYPrim(0.5), nbinsYMulti(100), fLowYMulti(0), fHighYMulti(100), mintrackrefsTPC(1), MAXCent(100), MINCent(0), fDebugMode(kFALSE), fMaxDCAZ(1.), fMaxChi2PerITSCluster(36), fMaxDCAXYFlag(kTRUE), fMaxDCAZFlag(kFALSE)
{
  // Constructor
  // Define output slots only here
   fPhysicsSelection = new AliPhysicsSelection();
  AliMCEventHandler *lmcEvtHandler  = dynamic_cast<AliMCEventHandler*>((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
  if(lmcEvtHandler) fPhysicsSelection->SetAnalyzeMC();

  AliInputEventHandler* handler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (handler) {
    handler->SetEventSelection(fPhysicsSelection);
    AliInfo("Physics Event Selection enabled.");
  } else {
    AliError("No input event handler connected to analysis manager. No Physics Event Selection.");
  }
  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class());
}
//________________________________________________________________________
AliAnalysisTaskProtonAbsorbtion::~AliAnalysisTaskProtonAbsorbtion(){

}
//________________________________________________________________________
void AliAnalysisTaskProtonAbsorbtion::UserCreateOutputObjects()
{
  //input hander
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=dynamic_cast<AliInputEventHandler*>(man->GetInputEventHandler());
  if (!inputHandler) AliFatal("Input handler needed");

  //pid response object
  AliMCEventHandler *lmcEvtHandler  = dynamic_cast<AliMCEventHandler*>((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
  if(lmcEvtHandler) inputHandler->CreatePIDResponse(kTRUE);
  else 
  inputHandler->CreatePIDResponse(kFALSE);
  fPIDResponse=inputHandler->GetPIDResponse();
  if (!fPIDResponse) AliFatal("PIDResponse object was not created");

  fPIDResponse->SetOADBPath(AliAnalysisManager::GetOADBPath());
  if (!fOADBPath.IsNull()) fPIDResponse->SetOADBPath(fOADBPath.Data());

  // Create histograms
  // Called once
 	fList = new TList();
	fList->SetOwner();

//  mintrackrefsTPC = 1;

//  nbinsYMulti = 100;
//  fLowYMulti = 0;
//  fHighYMulti = 100;

//  nbinsYPrim = 5;
//  fLowYPrim = -0.5;
//  fHighYPrim = 0.5;
//  nbinsPtPrim = 6;
//  fLowPtPrim = 0.45;
//  fHighPtPrim = 1.05;

Int_t iBinPrim[2];
  iBinPrim[0] = nbinsYPrim;
  iBinPrim[1] = nbinsPtPrim;
  Double_t *binLimYPrim = new Double_t[nbinsYPrim+1];
  Double_t *binLimPtPrim = new Double_t[nbinsPtPrim+1];
  //values for bin lower bounds
  for(Int_t i = 0; i <= nbinsYPrim; i++) 
    binLimYPrim[i]=(Double_t)fLowYPrim  + (fHighYPrim - fLowYPrim)  /nbinsYPrim*(Double_t)i;
  for(Int_t i = 0; i <= nbinsPtPrim; i++) 
    binLimPtPrim[i]=(Double_t)fLowPtPrim  + (fHighPtPrim - fLowPtPrim)  /nbinsPtPrim*(Double_t)i;

  Int_t iBinPrimMulti[2];
  iBinPrimMulti[0] = nbinsYMulti;
  iBinPrimMulti[1] = nbinsPtPrim;
  Double_t *binLimYMulti = new Double_t[nbinsYMulti+1];
  //values for bin lower bounds
  for(Int_t i = 0; i <= nbinsYMulti; i++) 
    binLimYMulti[i]=(Double_t)fLowYMulti  + (fHighYMulti - fLowYMulti)  /nbinsYMulti*(Double_t)i;

  //Proton container
  containerProtonsPrim = new AliCFContainer("containerProtonsPrim","container for protons",kNSteps,2,iBinPrim);
  containerProtonsPrim->SetBinLimits(0,binLimYPrim); //rapidity or eta
  containerProtonsPrim->SetBinLimits(1,binLimPtPrim); //pT
  fList->Add(containerProtonsPrim);
  //Anti-proton container
  containerAntiProtonsPrim = new AliCFContainer("containerAntiProtonsPrim","container for antiprotons",kNSteps,2,iBinPrim);
  containerAntiProtonsPrim->SetBinLimits(0,binLimYPrim); //rapidity or eta
  containerAntiProtonsPrim->SetBinLimits(1,binLimPtPrim); //pT
  fList->Add(containerAntiProtonsPrim);
  //Proton container - multi
  containerProtonsPrimMulti = new AliCFContainer("containerProtonsPrimMulti","container for protons",kNSteps,2,iBinPrimMulti);
  containerProtonsPrimMulti->SetBinLimits(0,binLimYMulti); //rapidity or eta
  containerProtonsPrimMulti->SetBinLimits(1,binLimPtPrim); //pT
  fList->Add(containerProtonsPrimMulti);
  //Anti-proton container - multi
  containerAntiProtonsPrimMulti = new AliCFContainer("containerAntiProtonsPrimMulti","container for antiprotons",kNSteps,2,iBinPrimMulti);
  containerAntiProtonsPrimMulti->SetBinLimits(0,binLimYMulti); //rapidity or eta
  containerAntiProtonsPrimMulti->SetBinLimits(1,binLimPtPrim); //pT
  fList->Add(containerAntiProtonsPrimMulti);




PostData(1, fList);
}

void AliAnalysisTaskProtonAbsorbtion::UserExec(Option_t *){

//Setting up the criteria for the generated particles
  AliCFTrackKineCuts *mcKineCutsPrim = new AliCFTrackKineCuts("mcKineCutsPrim","MC-level kinematic cuts");
  mcKineCutsPrim->SetPtRange(fLowPtPrim,fHighPtPrim);
  mcKineCutsPrim->SetRapidityRange(fLowYPrim,fHighYPrim);

  AliCFParticleGenCuts* mcGenCutsProtons = new AliCFParticleGenCuts("mcGenCutsProtons","MC particle generation cuts");
  mcGenCutsProtons->SetRequireIsPrimary();
  mcGenCutsProtons->SetRequirePdgCode(2212);

  AliCFParticleGenCuts* mcGenCutsAntiProtons = new AliCFParticleGenCuts("mcGenCutsAntiProtons","MC particle generation cuts");
  mcGenCutsAntiProtons->SetRequireIsPrimary();
  mcGenCutsAntiProtons->SetRequirePdgCode(-2212);

  TObjArray* mcListProtonsPrim = new TObjArray(0);
  mcListProtonsPrim->AddLast(mcKineCutsPrim);
  mcListProtonsPrim->AddLast(mcGenCutsProtons);
  TObjArray* mcListAntiProtonsPrim = new TObjArray(0);
  mcListAntiProtonsPrim->AddLast(mcKineCutsPrim);
  mcListAntiProtonsPrim->AddLast(mcGenCutsAntiProtons);

  //____________________________________________//
  //Setting up the criteria for the reconstructible particles
  //____________________________________________//
  AliCFAcceptanceCuts *mcAccCuts = new AliCFAcceptanceCuts("mcAccCuts","Acceptance cuts");
  mcAccCuts->SetMinNHitTPC(mintrackrefsTPC);
  TObjArray* accList = new TObjArray(0);
  accList->AddLast(mcAccCuts);

  //____________________________________________//
  //Setting up the criteria for the reconstructed tracks
  //____________________________________________//
  AliCFTrackKineCuts *recKineCutsPrim = new AliCFTrackKineCuts("recKineCutsPrim","rec-level kine cuts");
  recKineCutsPrim->SetPtRange(fLowPtPrim,fHighPtPrim);
  recKineCutsPrim->SetRapidityRange(fLowYPrim,fHighYPrim);
  recKineCutsPrim->SetChargeRec(1.0);

  AliCFTrackKineCuts *recKineCutsAntiPrim = new AliCFTrackKineCuts("recKineCutsAntiPrim","rec-level kine cuts");
  recKineCutsAntiPrim->SetPtRange(fLowPtPrim,fHighPtPrim);
  recKineCutsAntiPrim->SetRapidityRange(fLowYPrim,fHighYPrim);
  recKineCutsAntiPrim->SetChargeRec(-1.0);

  TObjArray* recListProtonsPrim = new TObjArray(0);
  recListProtonsPrim->AddLast(recKineCutsPrim);
  recListProtonsPrim->AddLast(mcGenCutsProtons);
  TObjArray* recListAntiProtonsPrim = new TObjArray(0);
  recListAntiProtonsPrim->AddLast(recKineCutsAntiPrim);
  recListAntiProtonsPrim->AddLast(mcGenCutsAntiProtons);

 //=========================================================//
  //CF manager - Protons
  fCFManagerProtonsPrim = new AliCFManager();
  fCFManagerProtonsPrim->SetParticleContainer(containerProtonsPrim);
  fCFManagerProtonsPrim->SetParticleCutsList(AliCFManager::kPartGenCuts,mcListProtonsPrim);
  fCFManagerProtonsPrim->SetParticleCutsList(AliCFManager::kPartAccCuts,accList);
  fCFManagerProtonsPrim->SetParticleCutsList(AliCFManager::kPartRecCuts,recListProtonsPrim);

  //CF manager - AntiProtons
  fCFManagerAntiProtonsPrim = new AliCFManager();
  fCFManagerAntiProtonsPrim->SetParticleContainer(containerAntiProtonsPrim);
  fCFManagerAntiProtonsPrim->SetParticleCutsList(AliCFManager::kPartGenCuts,mcListAntiProtonsPrim);
  fCFManagerAntiProtonsPrim->SetParticleCutsList(AliCFManager::kPartAccCuts,accList);
  fCFManagerAntiProtonsPrim->SetParticleCutsList(AliCFManager::kPartRecCuts,recListAntiProtonsPrim);

 //CF manager - Protons Multi
  fCFManagerProtonsPrimMulti = new AliCFManager();
  fCFManagerProtonsPrimMulti->SetParticleContainer(containerProtonsPrimMulti);
  fCFManagerProtonsPrimMulti->SetParticleCutsList(AliCFManager::kPartGenCuts,mcListProtonsPrim);
  fCFManagerProtonsPrimMulti->SetParticleCutsList(AliCFManager::kPartAccCuts,accList);
  fCFManagerProtonsPrimMulti->SetParticleCutsList(AliCFManager::kPartRecCuts,recListProtonsPrim);

  //CF manager - AntiProtons Multi
  fCFManagerAntiProtonsPrimMulti = new AliCFManager();
  fCFManagerAntiProtonsPrimMulti->SetParticleContainer(containerAntiProtonsPrimMulti);
  fCFManagerAntiProtonsPrimMulti->SetParticleCutsList(AliCFManager::kPartGenCuts,mcListAntiProtonsPrim);
  fCFManagerAntiProtonsPrimMulti->SetParticleCutsList(AliCFManager::kPartAccCuts,accList);
  fCFManagerAntiProtonsPrimMulti->SetParticleCutsList(AliCFManager::kPartRecCuts,recListAntiProtonsPrim);

// Main loop
  // Called for each event


// Create MC part and stack
           AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
           if (!eventHandler) {
	     Printf("ERROR: Could not retrieve MC event handler");
	     return;
           }
           AliMCEvent* mcEvent = eventHandler->MCEvent();
           if (!mcEvent) {
	     Printf("ERROR: Could not retrieve MC event");
	     return;
           }
	  AliStack* mcStack = mcEvent->Stack(); 

//--------------------------------------------------------------------------

// Create ESD part and stack
  AliVEvent* lEvent = InputEvent();
  if (!lEvent) {
    Printf("ERROR: Event not available");
    return;
  }

  AliESDEvent* lESDEvent = (AliESDEvent*)lEvent;



 AliCentrality *centrality = lESDEvent->GetCentrality();
  AliESDtrackCuts *fTrackCuts = new AliESDtrackCuts();
  if(!fTrackCuts){
	AliError("Can't get track cut object");
	}
Int_t nTracklets = 0;
nTracklets = centrality->GetCentralityPercentile("V0M");

  if(centrality->GetCentralityPercentile("V0M") > MAXCent || centrality->GetCentralityPercentile("V0M") < MINCent) return;


  Bool_t isSelected = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  if ((fUsePhysicsSelection)&&(!isSelected)) return;


// Create vertex
const AliVVertex *primaryVtx = lEvent->GetPrimaryVertex();
const AliESDVertex *trkVtx = lESDEvent->GetPrimaryVertex();
  if (!trkVtx || trkVtx->GetNContributors()<=0) return;
  TString vtxTtl = trkVtx->GetTitle();
  if (!vtxTtl.Contains("VertexerTracks")) return;
  Float_t zvtx = trkVtx->GetZ();
const AliESDVertex* spdVtx = lESDEvent->GetPrimaryVertexSPD();
  if (spdVtx->GetNContributors()<=0) return;
  TString vtxTyp = spdVtx->GetTitle();
  Double_t cov[6]={0};
  spdVtx->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if (vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) return;
  if (TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5) return;
//--------------------------------------------------------------------------
if (TMath::Abs(zvtx)>fMaxPrimaryVtxPosZ) return;

  fCFManagerProtonsPrim->SetMCEventInfo(mcEvent);
  fCFManagerAntiProtonsPrim->SetMCEventInfo(mcEvent);

  fCFManagerProtonsPrimMulti->SetMCEventInfo(mcEvent);
  fCFManagerAntiProtonsPrimMulti->SetMCEventInfo(mcEvent);

  //Dummy container
  Double_t containerInput[2];
  Double_t containerInputMulti[2];

//loop on the MC event
  for (Int_t ipart = 0; ipart < mcEvent->GetNumberOfTracks(); ipart++) { 
  AliMCParticle *mcPart  = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(ipart));
  if(!mcPart) continue;

Double_t vz = mcPart->Zv();
    if (TMath::Abs(vz) > 50.) continue;//exclude particles generated out of the acceptance
    if(TMath::Abs(mcPart->Eta()) > 1.0) continue;

    containerInput[0] = (Float_t)mcPart->Y();
    containerInput[1] = (Float_t)mcPart->Pt();
    containerInputMulti[0] = (Float_t)nTracklets;
    containerInputMulti[1] = (Float_t)mcPart->Pt();


    //Protons
    //check the MC-level cuts
    if (fCFManagerProtonsPrim->CheckParticleCuts(AliCFManager::kPartGenCuts,mcPart)) {
 	fCFManagerProtonsPrim->GetParticleContainer()->Fill(containerInput,kStepGenerated);
	fCFManagerProtonsPrimMulti->GetParticleContainer()->Fill(containerInputMulti,kStepGenerated);
        //check the Acceptance-level cuts
        if (fCFManagerProtonsPrim->CheckParticleCuts(AliCFManager::kPartAccCuts,mcPart)){ 
		fCFManagerProtonsPrim->GetParticleContainer()->Fill(containerInput,kStepReconstructible);
		fCFManagerProtonsPrimMulti->GetParticleContainer()->Fill(containerInputMulti,kStepReconstructible);
	}//MC acceptance cuts - protons
    }//MC level cuts - protons

    //Antiprotons
    //check the MC-level cuts
    if (fCFManagerAntiProtonsPrim->CheckParticleCuts(AliCFManager::kPartGenCuts,mcPart)) {
	fCFManagerAntiProtonsPrim->GetParticleContainer()->Fill(containerInput,kStepGenerated);
	fCFManagerAntiProtonsPrimMulti->GetParticleContainer()->Fill(containerInputMulti,kStepGenerated);
        //check the Acceptance-level cuts
        if (fCFManagerAntiProtonsPrim->CheckParticleCuts(AliCFManager::kPartAccCuts,mcPart)){
		fCFManagerAntiProtonsPrim->GetParticleContainer()->Fill(containerInput,kStepReconstructible);
		fCFManagerAntiProtonsPrimMulti->GetParticleContainer()->Fill(containerInputMulti,kStepReconstructible);
	}//MC acceptance cuts - antiprotons
    }//MC level cuts - antiprotons

  }//loop over MC particles

  //__________________________________________________________//
  //ESD track loop
 for (Int_t iTrack = 0; iTrack < lESDEvent->GetNumberOfTracks(); iTrack++) {
	   AliESDtrack *track = dynamic_cast<AliESDtrack *>(lESDEvent->GetTrack(iTrack));
    if(!track) continue;

    Int_t label = track->GetLabel();
    if (label < 0) continue;
    if(label > mcEvent->GetNumberOfTracks()) continue;

    AliMCParticle *mcPart  = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(label));
    if(!mcPart) continue;

    Double_t vz = mcPart->Zv();
    if (TMath::Abs(vz) > 50.) continue;//exclude particles generated out of the acceptance
    if(TMath::Abs(mcPart->Eta()) > 1.0) continue;

 AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
      if(!tpcTrack) continue;

    containerInput[0] = (Float_t)mcPart->Y();
    containerInput[1] = (Float_t)mcPart->Pt();
    containerInputMulti[0] = (Float_t)nTracklets;
    containerInputMulti[1] = (Float_t)mcPart->Pt();

    //Protons
    // check if this track was part of the signal
    if (fCFManagerProtonsPrim->CheckParticleCuts(AliCFManager::kPartGenCuts,mcPart)) {
	fCFManagerProtonsPrim->GetParticleContainer()->Fill(containerInput,kStepReconstructed);
	fCFManagerProtonsPrimMulti->GetParticleContainer()->Fill(containerInputMulti,kStepReconstructed);
	if(IsAccepted(track)){
		fCFManagerProtonsPrim->GetParticleContainer()->Fill(containerInput,kStepSurvived);
		fCFManagerProtonsPrimMulti->GetParticleContainer()->Fill(containerInputMulti,kStepSurvived);
	}//cut survived tracks   
    }//MC primaries - protons
    
    //Antiprotons    
    // check if this track was part of the signal
    if (fCFManagerAntiProtonsPrim->CheckParticleCuts(AliCFManager::kPartGenCuts,mcPart)) {
	fCFManagerAntiProtonsPrim->GetParticleContainer()->Fill(containerInput,kStepReconstructed);
	fCFManagerAntiProtonsPrimMulti->GetParticleContainer()->Fill(containerInputMulti,kStepReconstructed);
	if(IsAccepted(track)){
		fCFManagerAntiProtonsPrim->GetParticleContainer()->Fill(containerInput,kStepSurvived);
		fCFManagerAntiProtonsPrimMulti->GetParticleContainer()->Fill(containerInputMulti,kStepSurvived);
	}//cut survived tracks   
    }//MC primaries - antiprotons

}//track loop

}

//________________________________________________________________________
Bool_t AliAnalysisTaskProtonAbsorbtion::IsProton(AliESDtrack *track) {
  //Function that checks if a track is a proton
  
   Double_t gPt = 0.0, gP = 0.0, gEta = 0.0;
  
  //Ratio of the measured over the theoretical dE/dx a la STAR
  if(fPIDMode == kRatio) {
    AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
    if(tpcTrack) {
      gPt = tpcTrack->Pt();
      gP = track->GetInnerParam()->P();
      gEta = tpcTrack->Eta();
    }
    
    Double_t normalizeddEdx = -10.;
    if((track->GetTPCsignal() > 0.0) && (fPIDResponse->GetTPCResponse().GetExpectedSignal(gP,AliPID::kProton) > 0.0))
      normalizeddEdx = TMath::Log(track->GetTPCsignal()/fPIDResponse->GetTPCResponse().GetExpectedSignal(gP,AliPID::kProton));

    if (gP <= fNBoundP) if(normalizeddEdx >= fNRatio1) return kTRUE;
    if (gP >  fNBoundP) if(normalizeddEdx >= fNRatio2) return kTRUE;
  }//kRatio PID mode

  //Definition of an N-sigma area around the dE/dx vs P band
 else if(fPIDMode == kSigma) {
   
    Double_t nsigma = 100.0;
    
    Double_t mom = track->GetP();
    const AliExternalTrackParam *in = track->GetInnerParam();
    if (in)
      mom = in->GetP();

    nsigma = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton));
  
      if(nsigma <= fNSigma1) return kTRUE;
    //if (mom <= fNBoundP) if(nsigma <= fNSigma1) return kTRUE;
    //if (mom >  fNBoundP) if(nsigma <= fNSigma2) return kTRUE;
 }//kSigma PID method 

  return kFALSE;
}

Bool_t AliAnalysisTaskProtonAbsorbtion::IsAccepted(AliESDtrack* track) {
  // Checks if the track is excluded from the cuts
  //Int_t  fIdxInt[200];
  //Int_t nClustersITS = track->GetITSclusters(fIdxInt);
  //Int_t nClustersTPC = track->GetTPCclusters(fIdxInt);
  Int_t nClustersITS = track->GetITSclusters(0x0);
  Int_t nClustersTPC = track->GetTPCclusters(0x0);

  Float_t chi2PerClusterITS = -1;
  if (nClustersITS!=0)
    chi2PerClusterITS = track->GetITSchi2()/Float_t(nClustersITS);
  Float_t chi2PerClusterTPC = -1;
  if (nClustersTPC!=0)
    chi2PerClusterTPC = track->GetTPCchi2()/Float_t(nClustersTPC);

 // Double_t extCov[15];
 // track->GetExternalCovariance(extCov);

if(nClustersITS < fMinITSClusters) {
	if(fDebugMode) Printf("IsAccepted: Track rejected because it has %d ITS points (min. requested: %d)",nClustersITS,fMinITSClusters);
      return kFALSE;
    }
if(chi2PerClusterITS > fMaxChi2PerITSCluster) {
	if(fDebugMode) Printf("IsAccepted: Track rejected because it has a chi2 per ITS cluster %lf (max. requested: %lf)",chi2PerClusterITS,fMaxChi2PerITSCluster);
      return kFALSE;
    }
if(nClustersTPC < fMinTPCClusters) {
	if(fDebugMode) Printf("IsAccepted: Track rejected because it has %d TPC clusters (min. requested: %d)",nClustersTPC,fMinTPCClusters);
      return kFALSE;
    }
if(chi2PerClusterTPC > fMaxChi2PerTPCCluster) {
	if(fDebugMode) Printf("IsAccepted: Track rejected because it has a chi2 per TPC cluster %lf (max. requested: %lf)",chi2PerClusterTPC,fMaxChi2PerTPCCluster);
      return kFALSE; 
    }
if(track->GetTPCsignalN() < fMinTPCClusters) {
	if(fDebugMode) Printf("IsAccepted: Track rejected because it has %d TPC points for the calculation of the energy loss (min. requested: %d)",track->GetTPCsignalN(),fMinTPCClusters);
      return kFALSE;
    }
return kTRUE;
}

Bool_t AliAnalysisTaskProtonAbsorbtion::IsPrimary(AliESDEvent *esd,
					const AliESDVertex *vertex, 
					AliESDtrack* track) {
  // Checks if the track is a primary-like candidate
  Double_t gPt = 0.0, gPx = 0.0, gPy = 0.0, gPz = 0.0;
  Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};  //The impact parameters and their covariance.
  Double_t dca3D = 0.0;
  Float_t dcaXY = 0.0, dcaZ = 0.0;

    AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
    AliExternalTrackParam cParam;
    if(!tpcTrack) {
      gPt = 0.0; gPx = 0.0; gPy = 0.0; gPz = 0.0;
      dca[0] = -100.; dca[1] = -100.; dca3D = -100.;
      cov[0] = -100.; cov[1] = -100.; cov[2] = -100.;
    }
    else {
      gPt = tpcTrack->Pt();
      gPx = tpcTrack->Px();
      gPy = tpcTrack->Py();
      gPz = tpcTrack->Pz();

      track->RelateToVertex(vertex,
			    esd->GetMagneticField(),
			    100.,&cParam);
      track->GetImpactParameters(dcaXY,dcaZ);
      dca[0] = dcaXY; dca[1] = dcaZ;
    }
  //standalone TPC or hybrid TPC approaches
  dca3D = TMath::Sqrt(TMath::Power(dca[0],2) +
		      TMath::Power(dca[1],2));

   if(fMaxDCAXYFlag) { 
    if(TMath::Abs(dca[0]) > fMaxDCAXY) {
	if(fDebugMode)
	Printf("IsPrimary: Track rejected because it has a value of dca(xy) of %lf (max. requested: %lf)",TMath::Abs(dca[0]),fMaxDCAXY);
      return kFALSE;
    }
  }

   if(fMaxDCAZFlag) { 
    if(TMath::Abs(dca[1]) > fMaxDCAZ) {
	if(fDebugMode)
	Printf("IsPrimary: Track rejected because it has a value of dca(z) of %lf (max. requested: %lf)",TMath::Abs(dca[1]),fMaxDCAZ);
      return kFALSE;
    }
  }

return kTRUE;
  }

//____________________________________________________________________//
Double_t AliAnalysisTaskProtonAbsorbtion::Rapidity(Double_t gPx, 
					      Double_t gPy, 
					      Double_t gPz, Int_t fType) const {
  //returns the rapidity of the proton - to be removed
  Double_t fMass = 9.38270000000000048e-01;
  if (fType == 0) fMass = TDatabasePDG::Instance()->GetParticle(2212)->Mass();
  if (fType == 1) fMass = TDatabasePDG::Instance()->GetParticle(321)->Mass();
  
  Double_t gP = TMath::Sqrt(TMath::Power(gPx,2) + 
                           TMath::Power(gPy,2) + 
			   TMath::Power(gPz,2));
  Double_t energy = TMath::Sqrt(gP*gP + fMass*fMass);
  Double_t y = -999;
  if(energy != gPz) 
    y = 0.5*TMath::Log((energy + gPz)/(energy - gPz));

  return y;
}
