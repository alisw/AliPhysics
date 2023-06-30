/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Baldo Sahlmueller, Friederike Bock                             *
 * Version 1.0                                                            *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// Photon specific study by Hikari Murakami                               
// Purpose: (1)calculate Nch based on V0 multiplicity (ALICE acceptance)
//          (2)define multiplicity class
//          (3)get photon distribution based on the multiplicity class

//////////////////////////////////////////////////////////////////
//----------------------------------------------------------------
// Class used to multiplicity dependent photon analysis 
//----------------------------------------------------------------
//////////////////////////////////////////////////////////////////
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TPDGCode.h"
#include "TProfile.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliAnalysisTaskGammaPythia.h"
#include "AliVParticle.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliEventplane.h"
#include "AliInputEventHandler.h"
#include <algorithm>
#include <array>
#include <vector>
#include <map>

ClassImp(AliAnalysisTaskGammaPythia)

//________________________________________________________________________
AliAnalysisTaskGammaPythia::AliAnalysisTaskGammaPythia(): AliAnalysisTaskSE(),
  fOutputContainer(nullptr),
  fHistNEvents(nullptr),
  fHistNEventsHM(nullptr),
  fHistXSection(nullptr),
  fHistPtHard(nullptr),
  fHistParticlePDG(nullptr),
  fHistMotherParticlePDG(nullptr),
  fHistGrandMotherParticlePDG(nullptr),
  fHistParticlevsMother(nullptr),
  fHistPtGamma(nullptr),
  fHistPtYGamma(nullptr),
  fHistPtMultGamma(nullptr),
  fHistMult(nullptr),
  fHistV0Mult(nullptr),
  fHistV0MultHM(nullptr),
  fIsMC(1),
  fMaxY(2),
  fNch_min(0),
  fNch_max(30),
  fMaxpT(100),
  fDoMultStudies(0),
  fNTracks(0),
  fNTracksInV0Acc(0),
  fNTracksInV0AccHM(0),
  fIsEvtINELgtZERO(0)
{

}

//________________________________________________________________________
AliAnalysisTaskGammaPythia::AliAnalysisTaskGammaPythia(const char *name):
  AliAnalysisTaskSE(name),
  fOutputContainer(nullptr),
  fHistNEvents(nullptr),
  fHistNEventsHM(nullptr),
  fHistXSection(nullptr),
  fHistPtHard(nullptr),
  fHistParticlePDG(nullptr),
  fHistMotherParticlePDG(nullptr),
  fHistGrandMotherParticlePDG(nullptr),
  fHistParticlevsMother(nullptr),
  fHistPtGamma(nullptr),
  fHistPtYGamma(nullptr),
  fHistPtMultGamma(nullptr),
  fHistMult(nullptr),
  fHistV0Mult(nullptr),
  fHistV0MultHM(nullptr),
  fIsMC(1),
  fMaxY(2),
  fNch_min(0),
  fNch_max(30),
  fMaxpT(100),
  fDoMultStudies(0),
  fNTracks(0),
  fNTracksInV0Acc(0),
  fNTracksInV0AccHM(0),
  fIsEvtINELgtZERO(0)
{
  // Define output slots here
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskGammaPythia::~AliAnalysisTaskGammaPythia()
{

}

//________________________________________________________________________
void AliAnalysisTaskGammaPythia::UserCreateOutputObjects(){

  // Create histograms
  if(fOutputContainer != nullptr){
    delete fOutputContainer;
    fOutputContainer          = nullptr;
  }
  if(fOutputContainer == nullptr){
    fOutputContainer          = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }

  fHistNEvents                		= new TH1F("NEvents", "", 3, -0.5, 2.5);
  fHistNEvents->Sumw2();
  fOutputContainer->Add(fHistNEvents);

  fHistNEventsHM                		= new TH1F("NEventsHM", "", 3, -0.5, 2.5);
  fHistNEventsHM->Sumw2();
  fOutputContainer->Add(fHistNEventsHM);
 
  fHistXSection               		= new TH1D("XSection", "", 1000000, 0, 1e4);  
  fHistXSection->Sumw2();
  fOutputContainer->Add(fHistXSection);

  fHistPtHard                 		= new TH1F("PtHard", "", 400, 0, 200);
  fHistPtHard->Sumw2();
  fOutputContainer->Add(fHistPtHard);

  fHistParticlePDG                 = new TH1I("ParticlePDG", "", 5000, 0, 5000);
  fOutputContainer->Add(fHistParticlePDG);

  fHistMotherParticlePDG                 = new TH1I("MotherParticlePDG", "", 5000, 0, 5000);
  fOutputContainer->Add(fHistMotherParticlePDG);

  fHistGrandMotherParticlePDG                 = new TH1I("GrandMotherParticlePDG", "", 5000, 0, 5000);
  fOutputContainer->Add(fHistGrandMotherParticlePDG);

  fHistParticlevsMother                 = new TH2I("ParticlevsMother", "", 5000, 0, 5000, 5000, 0, 5000);
  fOutputContainer->Add(fHistParticlevsMother);

  fHistPtGamma                		= new TH1F("Pt_Gamma","", fMaxpT*10, 0, fMaxpT);
  fHistPtGamma->Sumw2();
  fOutputContainer->Add(fHistPtGamma);

  fHistPtYGamma                		= new TH2F("Pt_Y_Gamma","", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYGamma->Sumw2();
  fOutputContainer->Add(fHistPtYGamma);

  fHistPtMultGamma                		= new TH2F("Pt_Mult_Gamma","", fMaxpT*10, 0, fMaxpT, 1000, -0.5, 1000 - 0.5);
  fHistPtMultGamma->Sumw2();
  fOutputContainer->Add(fHistPtMultGamma);

  fHistMult = new TH1D("Multiplicity", "", 1000, -0.5, 1000 - 0.5);
  fHistMult->Sumw2();
  fOutputContainer->Add(fHistMult);

  fHistV0Mult = new TH1D("V0Multiplicity", "", 1000, -0.5, 1000 - 0.5);
  fHistV0Mult->Sumw2();
  fOutputContainer->Add(fHistV0Mult);

  fHistV0MultHM = new TH1D("V0MultiplicityHM", "", 1000, -0.5, 1000 - 0.5);
  fHistV0MultHM->Sumw2();
  fOutputContainer->Add(fHistV0MultHM);

  PostData(1, fOutputContainer);
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaPythia::UserExec(Option_t *)
{

  fInputEvent = InputEvent();

  fMCEvent = MCEvent();
  if(fMCEvent == nullptr) fIsMC = 0;
  if (fIsMC==0){
    //    printf("UserExec()   fMCEvent=NULL \n");	
    return;
  }			
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();
  
  if (TMath::Abs(mcProdVtxZ) < 10 ){
    fHistNEvents->Fill(0);
  } else {
    fHistNEvents->Fill(1);
  }
  
  AliGenEventHeader* mcEH = fMCEvent->GenEventHeader();
  AliGenPythiaEventHeader *pyH  = dynamic_cast<AliGenPythiaEventHeader*>(mcEH);
  AliGenHijingEventHeader *hiH  = 0;
  AliGenDPMjetEventHeader *dpmH = 0;
  
  // it can be only one save some casts
  // assuming PYTHIA and HIJING are the most likely ones...
  if(!pyH){
    hiH = dynamic_cast<AliGenHijingEventHeader*>(mcEH);
    if(!hiH){
      dpmH = dynamic_cast<AliGenDPMjetEventHeader*>(mcEH);
    }
  }

  // fetch the trials on a event by event basis, not from pyxsec.root otherwise
  // we will get a problem when running on proof since Notify may be called
  // more than once per file
  // consider storing this information in the AOD output via AliAODHandler
  Float_t ntrials = 0;
  if (!pyH || !hiH || dpmH) {
    AliGenCocktailEventHeader *ccEH = dynamic_cast<AliGenCocktailEventHeader *>(mcEH);
    if (ccEH) {
      TList *genHeaders = ccEH->GetHeaders();
      for (int imch=0; imch<genHeaders->GetEntries(); imch++) {
        if(!pyH)pyH = dynamic_cast<AliGenPythiaEventHeader*>(genHeaders->At(imch));
        if(!hiH)hiH = dynamic_cast<AliGenHijingEventHeader*>(genHeaders->At(imch));
        if(!dpmH)dpmH = dynamic_cast<AliGenDPMjetEventHeader*>(genHeaders->At(imch));
      }
    }
  }
  
  // take the trials from the p+p event
  if(hiH)ntrials = hiH->Trials();
  if(dpmH)ntrials = dpmH->Trials();
  if(pyH)ntrials = pyH->Trials();
  if(ntrials)fHistNEvents->Fill(2,ntrials);
  
  Double_t xSection = 0;
  Double_t ptHard = 0;
  if (pyH) xSection = pyH->GetXsection();
  if (pyH) ptHard = pyH->GetPtHard();
  if (xSection) fHistXSection->Fill(xSection);
  if (ptHard) fHistPtHard->Fill(ptHard);
  
  ProcessMultiplicity();

  if((fNch_min <= fNTracksInV0Acc) && (fNTracksInV0Acc < fNch_max) ){
    fHistNEventsHM->Fill(0);
    fHistV0MultHM->Fill(fNTracksInV0Acc);
    // Loop over all primary MC particle
    for(Long_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {
      // fill primary histograms
      AliVParticle* particle     = nullptr;
      particle                    = (AliVParticle *)fMCEvent->GetTrack(i);
      //    printf("%d",particle->PdgCode());
      fHistParticlePDG->Fill(TMath::Abs(particle->PdgCode()));
      if (!particle) continue;
      Bool_t hasMother            = kFALSE;
      Bool_t particleIsPrimary    = kTRUE;

      if (particle->GetMother()>-1){
	hasMother         = kTRUE;
	particleIsPrimary = kFALSE;
      }
      AliVParticle*  motherParticle   = nullptr;
      if (hasMother){
	//      printf("%d",motherParticle->PdgCode());
	motherParticle  = (AliVParticle*)fMCEvent->GetTrack(particle->GetMother());
	fHistMotherParticlePDG->Fill(TMath::Abs(motherParticle->PdgCode()));
      }
      if (motherParticle) hasMother   = kTRUE;
      else                hasMother   = kFALSE;

      Bool_t motherIsPrimary                                = kFALSE;
      if(hasMother){
	if(motherParticle->GetMother()>-1) motherIsPrimary = kFALSE;
	else                                motherIsPrimary = kTRUE;
      }

      AliVParticle* grandMotherParticle  = nullptr;
      Bool_t hasGrandMother          = kFALSE;
      if (hasMother && !motherIsPrimary) {
	grandMotherParticle           = (AliVParticle*)fMCEvent->GetTrack(motherParticle->GetMother());
	hasGrandMother               = kTRUE;
	fHistGrandMotherParticlePDG->Fill(TMath::Abs(grandMotherParticle->PdgCode()));
      }

      Bool_t grandMotherIsPrimary                                     = kFALSE;
      if (hasGrandMother) {
	if(grandMotherParticle->GetMother()>-1)  grandMotherIsPrimary = kFALSE;
	else                                     grandMotherIsPrimary = kTRUE;
      }
      
      if (!(TMath::Abs(particle->E()-particle->Pz())>0.)) continue;
      Double_t yPre = (particle->E()+particle->Pz())/(particle->E()-particle->Pz());
      if( yPre <= 0 ) continue;
      
      Double_t y = 0.5*TMath::Log(yPre);
      
      Int_t PdgAnalyzedParticle = 22;
      // gamma from source
      if(TMath::Abs(particle->PdgCode())==PdgAnalyzedParticle){
	if(hasMother==kTRUE){
	  fHistParticlevsMother->Fill(TMath::Abs(motherParticle->PdgCode()),TMath::Abs(particle->PdgCode()));
	  if(hasGrandMother==kFALSE){
	    if (TMath::Abs(y) > fMaxY) continue;
	    fHistPtGamma->Fill(particle->Pt());
	    fHistPtYGamma->Fill(particle->Pt(), particle->Y());
	    fHistPtMultGamma->Fill(particle->Pt(), fNTracksInV0Acc);
	  }
	}
      }

    }//End of loop over all primary MC particle
  }

  PostData(1, fOutputContainer);
}


void AliAnalysisTaskGammaPythia::ProcessMultiplicity()
{
  // set number of tracks to 0
  fNTracks = 0;
  // set number of tracks in V0 acceptance to 0
  fNTracksInV0Acc = 0;
  // set INEL>0 to false
  fIsEvtINELgtZERO = false;

  // Loop over all primary MC particle
  for(Long_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {
    AliVParticle* particle     = nullptr;
    particle                    = (AliVParticle *)fMCEvent->GetTrack(i);
    if (!particle) continue;
    fNTracks++;
    // selected charged primary particles in V0 acceptance
    if(particle->IsPhysicalPrimary() && particle->Charge() != 0){
      if(IsInV0Acceptance(particle)){
        fNTracksInV0Acc++;
      }
      // check if event is INEL>0
      if(!fIsEvtINELgtZERO){
        if(std::abs(particle->Eta()) < 1){
          fIsEvtINELgtZERO = true;
        }
      }
    }
  }

  fHistMult->Fill(fNTracks);
  fHistV0Mult->Fill(fNTracksInV0Acc);

}


//________________________________________________________________________
bool AliAnalysisTaskGammaPythia::IsInV0Acceptance(AliVParticle* part) const {
  const Double_t kBoundaryEtaMinV0A = 2.8;
  const Double_t kBoundaryEtaMaxV0A = 5.1;
  const Double_t kBoundaryEtaMinV0C = -3.7;
  const Double_t kBoundaryEtaMaxV0C = -1.7;
  if (part->Eta() < kBoundaryEtaMaxV0A && part->Eta() > kBoundaryEtaMinV0A) return true;
  if (part->Eta() < kBoundaryEtaMaxV0C && part->Eta() > kBoundaryEtaMinV0C) return true;

  //  if ((part->Eta() < kBoundaryEtaMaxV0A && part->Eta() > kBoundaryEtaMinV0A) || ( (part->Eta() < kBoundaryEtaMaxV0C && part->Eta() > kBoundaryEtaMinV0C))) return true;

  return false;
}
//________________________________________________________________________
void AliAnalysisTaskGammaPythia::Terminate(const Option_t *)
{

  //fOutputContainer->Print(); // Will crash on GRID
}

