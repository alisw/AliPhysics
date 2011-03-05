#include <Riostream.h>
#include <TFile.h>
#include <TSystem.h>
#include <TF1.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TH1I.h>
#include <TParticle.h>
#include <string.h>

#include <AliExternalTrackParam.h>
#include <AliAODEvent.h>
#include <AliESDEvent.h>
//#include <AliLog.h>
#include <AliPID.h>
#include <AliStack.h>
#include <AliCFContainer.h>
#include <AliCFEffGrid.h>
#include <AliCFDataGrid.h>
#include <AliCFManager.h>
#include <AliCFTrackKineCuts.h>
#include <AliCFParticleGenCuts.h>
#include <AliCFAcceptanceCuts.h>
#include <AliMCEvent.h>
//#include <AliESDVertex.h>
class AliLog;
class AliESDVertex;

#include "AliProtonSpectraCorrection.h"
#include "AliProtonAnalysisBase.h"

ClassImp(AliProtonSpectraCorrection)

//____________________________________________________________________//
AliProtonSpectraCorrection::AliProtonSpectraCorrection() : 
  TObject(), fProtonAnalysisBase(0),
  fNBinsY(0), fMinY(0), fMaxY(0),
  fNBinsPt(0), fMinPt(0), fMaxPt(0),
  fCFManagerProtons(0), fCFManagerAntiProtons(0) {
  //fProtonContainer(0), fAntiProtonContainer(0) {
  //Default constructor
}

//____________________________________________________________________//
AliProtonSpectraCorrection::~AliProtonSpectraCorrection() {
  //Default destructor
  if(fCFManagerProtons) delete fCFManagerProtons;
  if(fCFManagerAntiProtons) delete fCFManagerAntiProtons;
  if(fProtonAnalysisBase) delete fProtonAnalysisBase;
  //if(fProtonContainer) delete fProtonContainer;
  //if(fAntiProtonContainer) delete fAntiProtonContainer;
}

//____________________________________________________________________//
void AliProtonSpectraCorrection::InitAnalysisHistograms(Int_t nbinsY, 
							   Float_t fLowY, 
							   Float_t fHighY, 
							   Int_t nbinsPt, 
							   Float_t fLowPt, 
							   Float_t fHighPt) {
  //Initializes the histograms
  fNBinsY = nbinsY;
  fMinY = fLowY;
  fMaxY = fHighY;
  fNBinsPt = nbinsPt;
  fMinPt = fLowPt;
  fMaxPt = fHighPt;
  const Int_t    mintrackrefsTPC = 1;

  //=========================================================//
  //setting up the containers
  Int_t iBin[2];
  iBin[0] = nbinsY;
  iBin[1] = nbinsPt;
  Double_t *binLimY = new Double_t[nbinsY+1];
  Double_t *binLimPt = new Double_t[nbinsPt+1];
  //values for bin lower bounds
  for(Int_t i = 0; i <= nbinsY; i++) 
    binLimY[i]=(Double_t)fLowY  + (fHighY - fLowY)  /nbinsY*(Double_t)i;
  for(Int_t i = 0; i <= nbinsPt; i++) 
    binLimPt[i]=(Double_t)fLowPt  + (fHighPt - fLowPt)  /nbinsPt*(Double_t)i;

  //Proton container
  AliCFContainer *containerProtons = new AliCFContainer("containerProtons",
							"container for protons",
							kNSteps,2,iBin);
  containerProtons->SetBinLimits(0,binLimY); //rapidity or eta
  containerProtons->SetBinLimits(1,binLimPt); //pT

  //Anti-proton container
  AliCFContainer *containerAntiProtons = new AliCFContainer("containerAntiProtons",
							    "container for antiprotons",
							    kNSteps,2,iBin);
  containerAntiProtons->SetBinLimits(0,binLimY); //rapidity or eta
  containerAntiProtons->SetBinLimits(1,binLimPt); //pT
  
  //=========================================================//
  //Setting up the criteria for the generated particles
  Char_t titleCF[256] = "mcKineCutsProtons";
  Char_t nameCF[256] = "MC-level kinematic cuts";
  AliCFTrackKineCuts *mcKineCutsProtons = new AliCFTrackKineCuts(titleCF,nameCF);
  mcKineCutsProtons->SetPtRange(fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    mcKineCutsProtons->SetEtaRange(fMinY,fMaxY);
  else
    mcKineCutsProtons->SetRapidityRange(fMinY,fMaxY);
  mcKineCutsProtons->SetChargeMC(1.0);

  strncpy(titleCF,"mcKineCutsAntiProtons",22);
  AliCFTrackKineCuts *mcKineCutsAntiProtons = new AliCFTrackKineCuts(titleCF,nameCF);
  mcKineCutsAntiProtons->SetPtRange(fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    mcKineCutsAntiProtons->SetEtaRange(fMinY,fMaxY);
  else
    mcKineCutsAntiProtons->SetRapidityRange(fMinY,fMaxY);
  mcKineCutsAntiProtons->SetChargeMC(-1.0);

  strncpy(titleCF,"mcGenCutsProtons",17);
  strncpy(nameCF,"MC particle generation cuts",28);
  AliCFParticleGenCuts* mcGenCutsProtons = new AliCFParticleGenCuts(titleCF,nameCF);
  mcGenCutsProtons->SetRequireIsPrimary();
  mcGenCutsProtons->SetRequirePdgCode(2212);

  strncpy(titleCF,"mcGenCutsAntiProtons",21);
  AliCFParticleGenCuts* mcGenCutsAntiProtons = new AliCFParticleGenCuts(titleCF,nameCF);
  mcGenCutsAntiProtons->SetRequireIsPrimary();
  mcGenCutsAntiProtons->SetRequirePdgCode(-2212);

  TObjArray* mcListProtons = new TObjArray(0);
  mcListProtons->AddLast(mcKineCutsProtons);
  mcListProtons->AddLast(mcGenCutsProtons);
  TObjArray* mcListAntiProtons = new TObjArray(0);
  mcListAntiProtons->AddLast(mcKineCutsAntiProtons);
  mcListAntiProtons->AddLast(mcGenCutsAntiProtons);

  //Setting up the criteria for the reconstructible particles
  strncpy(titleCF,"mcAccCuts",10);
  strncpy(nameCF,"Acceptance cuts",16);
  AliCFAcceptanceCuts *mcAccCuts = new AliCFAcceptanceCuts(titleCF,nameCF);
  mcAccCuts->SetMinNHitTPC(mintrackrefsTPC);
  TObjArray* accList = new TObjArray(0);
  accList->AddLast(mcAccCuts);

  //____________________________________________//
  //Setting up the criteria for the reconstructed tracks
  //____________________________________________//
  strncpy(titleCF,"recKineCutsProtons",19);
  strncpy(nameCF,"rec-level kine cuts",20);
  AliCFTrackKineCuts *recKineCutsProtons = new AliCFTrackKineCuts(titleCF,nameCF);								  
  recKineCutsProtons->SetPtRange(fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    recKineCutsProtons->SetEtaRange(fMinY,fMaxY);
  else
    recKineCutsProtons->SetRapidityRange(fMinY,fMaxY);
  recKineCutsProtons->SetChargeRec(1.0);

  //____________________________________________//
  strncpy(titleCF,"recKineCutsAntiProtons",23);
  AliCFTrackKineCuts *recKineCutsAntiProtons = new AliCFTrackKineCuts(titleCF,nameCF);
  recKineCutsAntiProtons->SetPtRange(fMinPt,fMaxPt);
  if(fProtonAnalysisBase->GetEtaMode()) 
    recKineCutsAntiProtons->SetEtaRange(fMinY,fMinY);
  else
    recKineCutsAntiProtons->SetRapidityRange(fMinY,fMaxY);
  recKineCutsAntiProtons->SetChargeRec(-1.0);

  //____________________________________________//
  TObjArray* recListProtons = new TObjArray(0);
  recListProtons->AddLast(recKineCutsProtons);
  recListProtons->AddLast(mcGenCutsProtons);

  TObjArray* recListAntiProtons = new TObjArray(0);
  recListAntiProtons->AddLast(recKineCutsAntiProtons);
  recListAntiProtons->AddLast(mcGenCutsAntiProtons);

  //=========================================================//
  //CF manager - Protons
  fCFManagerProtons = new AliCFManager();
  fCFManagerProtons->SetParticleContainer(containerProtons);
  fCFManagerProtons->SetParticleCutsList(AliCFManager::kPartGenCuts,mcListProtons);
  fCFManagerProtons->SetParticleCutsList(AliCFManager::kPartAccCuts,accList);
  fCFManagerProtons->SetParticleCutsList(AliCFManager::kPartRecCuts,recListProtons);

  //CF manager - Protons
  fCFManagerAntiProtons = new AliCFManager();
  fCFManagerAntiProtons->SetParticleContainer(containerAntiProtons);
  fCFManagerAntiProtons->SetParticleCutsList(AliCFManager::kPartGenCuts,mcListAntiProtons);
  fCFManagerAntiProtons->SetParticleCutsList(AliCFManager::kPartAccCuts,accList);
  fCFManagerAntiProtons->SetParticleCutsList(AliCFManager::kPartRecCuts,recListAntiProtons);
}

//_________________________________________________________________________//
void AliProtonSpectraCorrection::FillCorrectionMaps(AliESDEvent *esd, 
						    const AliESDVertex *vertex,
						    AliMCEvent *mcEvent) {
  //Function to fill the correction containers
  fCFManagerProtons->SetMCEventInfo(mcEvent);
  fCFManagerAntiProtons->SetMCEventInfo(mcEvent);
 
  //Dummy container
  Double_t containerInput[2];
  //__________________________________________________________//
  //loop on the MC event
  for (Int_t ipart = 0; ipart < mcEvent->GetNumberOfTracks(); ipart++) { 
    AliMCParticle *mcPart  = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(ipart));

    //Protons
    //check the MC-level cuts
    if (mcPart && fCFManagerProtons->CheckParticleCuts(AliCFManager::kPartGenCuts,
					     mcPart)) {
      containerInput[0] = (Float_t)mcPart->Eta();
      containerInput[1] = (Float_t)mcPart->Pt();
      //fill the container for Gen-level selection
      fCFManagerProtons->GetParticleContainer()->Fill(containerInput,
						      kStepGenerated);
      //check the Acceptance-level cuts
      if (fCFManagerProtons->CheckParticleCuts(AliCFManager::kPartAccCuts,
					       mcPart)) {
	//fill the container for Acceptance-level selection
	fCFManagerProtons->GetParticleContainer()->Fill(containerInput,
							kStepReconstructible);
      }//acceptance cuts - protons
    }//MC level cuts - protons

    //Antiprotons
    //check the MC-level cuts
    if (fCFManagerAntiProtons->CheckParticleCuts(AliCFManager::kPartGenCuts,
						 mcPart)) {
      containerInput[0] = (Float_t)mcPart->Eta();
      containerInput[1] = (Float_t)mcPart->Pt();
      //fill the container for Gen-level selection
      fCFManagerAntiProtons->GetParticleContainer()->Fill(containerInput,
							  kStepGenerated);
      //check the Acceptance-level cuts
      if (fCFManagerAntiProtons->CheckParticleCuts(AliCFManager::kPartAccCuts,
						   mcPart)) {
	//fill the container for Acceptance-level selection
	fCFManagerAntiProtons->GetParticleContainer()->Fill(containerInput,
							    kStepReconstructible);
      }//acceptance cuts - antiprotons
    }//MC level cuts - antiprotons
  }//loop over MC particles
  
  
  //__________________________________________________________//
  //ESD track loop
  for (Int_t iTrack = 0; iTrack < esd->GetNumberOfTracks(); iTrack++) {
    AliESDtrack *track = dynamic_cast<AliESDtrack *>(esd->GetTrack(iTrack));
    if(!track) continue;

    // is track associated to particle ?
    Int_t label = track->GetLabel();
    if (label<0) continue;
    AliMCParticle *mcPart  = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(label));

    if((fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kTPC)||(fProtonAnalysisBase->GetAnalysisMode()==AliProtonAnalysisBase::kHybrid)) {
      AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
      if(!tpcTrack) continue;
    }//Hybrid or standalone TPC
    
    //Protons
    // check if this track was part of the signal - primaries
    if (fCFManagerProtons->CheckParticleCuts(AliCFManager::kPartGenCuts,
					     mcPart)) {
      //fill the container - reconstructed protons
      containerInput[0] = mcPart->Eta();
      containerInput[1] = mcPart->Pt();
      fCFManagerProtons->GetParticleContainer()->Fill(containerInput,
						      kStepReconstructed);
      //fill the container - identified protons
      if(fProtonAnalysisBase->IsProton(track)) {
	containerInput[0] = mcPart->Eta();
	containerInput[1] = mcPart->Pt();
	fCFManagerProtons->GetParticleContainer()->Fill(containerInput,
							kStepIdentified);

	//fill the container - survived protons
	//track cuts
	if(!fProtonAnalysisBase->IsPrimary(esd,vertex,track)) continue;
	if(!fProtonAnalysisBase->IsAccepted(track)) continue;
	//track outside the analyzed y-Pt
	if(!fProtonAnalysisBase->IsInPhaseSpace(track)) continue; 
	containerInput[0] = mcPart->Eta();
	containerInput[1] = mcPart->Pt();
	fCFManagerProtons->GetParticleContainer()->Fill(containerInput,
							kStepSelected);
      }
    }//MC primaries - protons
    
    //Antiprotons    
    // check if this track was part of the signal - primaries
    if (fCFManagerAntiProtons->CheckParticleCuts(AliCFManager::kPartGenCuts,
						 mcPart)) {
      //fill the container - reconstructed antiprotons
      containerInput[0] = mcPart->Eta();
      containerInput[1] = mcPart->Pt();
      fCFManagerAntiProtons->GetParticleContainer()->Fill(containerInput,
							  kStepReconstructed);

      //fill the container - identified antiprotons
      if(fProtonAnalysisBase->IsProton(track)) {
	containerInput[0] = mcPart->Eta();
	containerInput[1] = mcPart->Pt();
	fCFManagerAntiProtons->GetParticleContainer()->Fill(containerInput,
							    kStepIdentified);

	//fill the container - survived antiprotons
	//track cuts
	if(!fProtonAnalysisBase->IsPrimary(esd,vertex,track)) continue;
	if(!fProtonAnalysisBase->IsAccepted(track)) continue;
	//track outside the analyzed y-Pt
	if(!fProtonAnalysisBase->IsInPhaseSpace(track)) continue; 
	containerInput[0] = mcPart->Eta();
	containerInput[1] = mcPart->Pt();
	fCFManagerAntiProtons->GetParticleContainer()->Fill(containerInput,
							    kStepSelected);
      }
    }//MC primaries - antiprotons
  }//track loop
  
  //if(fProtonAnalysisBase->GetDebugMode())
  //Printf("Initial number of tracks: %d | Identified (anti)protons: %d - %d | Survived (anti)protons: %d - %d",nTracks,nIdentifiedProtons,nIdentifiedAntiProtons,nSurvivedProtons,nSurvivedAntiProtons);
}

//_________________________________________________________________________//
void AliProtonSpectraCorrection::FillCorrectionMaps(AliAODEvent *fAOD) {
  // Analysis from AOD: to be implemented...
  // in the near future.
  fAOD->Print();
}
