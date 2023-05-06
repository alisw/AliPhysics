/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Baldo Sahlmueller, Friederike Bock                             *
 * Version 1.0                                                            *
 * Photon specific study                                                  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

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
  fHistXSection(nullptr),
  fHistPtHard(nullptr),
  fHistPtGamma(nullptr),
  fHistPtYGamma(nullptr),
  fHistPtGammaClassI(nullptr),
  fHistPtGammaClassII(nullptr),
  fHistPtGammaClassIII(nullptr),
  fHistPtYGammaClassI(nullptr),
  fHistPtYGammaClassII(nullptr),
  fHistPtYGammaClassIII(nullptr),
  fHistV0Mult(nullptr),
  fIsMC(1),
  fMaxpT(100),
  fDoMultStudies(0),
  fNTracksInV0Acc(0),
  fIsEvtINELgtZERO(0)
{

}

//________________________________________________________________________
AliAnalysisTaskGammaPythia::AliAnalysisTaskGammaPythia(const char *name):
  AliAnalysisTaskSE(name),
  fOutputContainer(nullptr),
  fHistNEvents(nullptr),
  fHistXSection(nullptr),
  fHistPtHard(nullptr),
  fHistPtGamma(nullptr),
  fHistPtYGamma(nullptr),
  fHistPtGammaClassI(nullptr),
  fHistPtGammaClassII(nullptr),
  fHistPtGammaClassIII(nullptr),
  fHistPtYGammaClassI(nullptr),
  fHistPtYGammaClassII(nullptr),
  fHistPtYGammaClassIII(nullptr),
  fHistV0Mult(nullptr),
  fIsMC(1),
  fMaxpT(100),
  fDoMultStudies(0),
  fNTracksInV0Acc(0),
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

  fHistXSection               		= new TH1D("XSection", "", 1000000, 0, 1e4);  
  fHistXSection->Sumw2();
  fOutputContainer->Add(fHistXSection);

  fHistPtHard                 		= new TH1F("PtHard", "", 400, 0, 200);
  fHistPtHard->Sumw2();
  fOutputContainer->Add(fHistPtHard);

  fHistPtGamma                		= new TH1F("Pt_Gamma","", fMaxpT*10, 0, fMaxpT);
  fHistPtGamma->Sumw2();
  fOutputContainer->Add(fHistPtGamma);

  fHistPtGammaClassI                        = new TH1F("Pt_GammaClassI","", fMaxpT*10, 0, fMaxpT);
  fHistPtGammaClassI->Sumw2();
  fOutputContainer->Add(fHistPtGammaClassI);

  fHistPtGammaClassII                        = new TH1F("Pt_GammaClassII","", fMaxpT*10, 0, fMaxpT);
  fHistPtGammaClassII->Sumw2();
  fOutputContainer->Add(fHistPtGammaClassII);

  fHistPtGammaClassIII                        = new TH1F("Pt_GammaClassIII","", fMaxpT*10, 0, fMaxpT);
  fHistPtGammaClassIII->Sumw2();
  fOutputContainer->Add(fHistPtGammaClassIII);

  fHistPtYGamma                		= new TH2F("Pt_Y_Gamma","", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYGamma->Sumw2();
  fOutputContainer->Add(fHistPtYGamma);

  fHistPtYGammaClassI                	= new TH2F("Pt_Y_GammaClassI","", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYGammaClassI->Sumw2();
  fOutputContainer->Add(fHistPtYGammaClassI);

  fHistPtYGammaClassII                	= new TH2F("Pt_Y_GammaClassII","", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYGammaClassII->Sumw2();
  fOutputContainer->Add(fHistPtYGammaClassII);

  fHistPtYGammaClassIII                	= new TH2F("Pt_Y_GammaClassIII","", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYGammaClassIII->Sumw2();
  fOutputContainer->Add(fHistPtYGammaClassIII);

  fHistV0Mult = new TH1D("V0Multiplicity", "", 1000, -0.5, 1000 - 0.5);
  fHistV0Mult->Sumw2();
  fOutputContainer->Add(fHistV0Mult);

  PostData(1, fOutputContainer);
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaPythia::UserExec(Option_t *)
{

  fInputEvent = InputEvent();
  //   cout << "I found an Event" << endl;

  fMCEvent = MCEvent();
  if(fMCEvent == nullptr) fIsMC = 0;
  if (fIsMC==0) return;
  //   cout << "I found an MC header" << endl;

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
  ProcessMCParticles();

  PostData(1, fOutputContainer);
}

void AliAnalysisTaskGammaPythia::ProcessMultiplicity()
{
  // set number of tracks in V0 acceptance to 0
  fNTracksInV0Acc = 0;
  // set INEL>0 to false
  fIsEvtINELgtZERO = false;
  // Loop over all primary MC particle
  for(Long_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++){
    AliVParticle* particle = nullptr;
    particle               = (AliVParticle *)fMCEvent->GetTrack(i);
    if (!particle) continue;
    // selected charged primary particles in V0 acceptance
    if(particle->IsPhysicalPrimary() && particle->Charge() != 0){
      if(IsInV0Acceptance(particle)){
	fNTracksInV0Acc++;
      }
    }
  }
  fHistV0Mult->Fill(fNTracksInV0Acc);
}

//________________________________________________________________________
void AliAnalysisTaskGammaPythia::ProcessMCParticles()
{
  
  // Loop over all primary MC particle
  for(Long_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {
    // fill primary histograms
    AliVParticle* particle     = nullptr;
    particle                    = (AliVParticle *)fMCEvent->GetTrack(i);
    if (!particle) continue;
    Bool_t hasMother            = kFALSE;
    //     cout << i << "\t"<< particle->GetMother() << endl;
    if (particle->GetMother()>-1)
      hasMother                 = kTRUE;
    AliVParticle* motherParticle   = nullptr;
    if( hasMother )
      motherParticle            = (AliVParticle *)fMCEvent->GetTrack(particle->GetMother());
    if (motherParticle)
      hasMother                 = kTRUE;
    else
      hasMother                 = kFALSE;

    //    const std::array<int, 1> kAcceptPdgCodes = {kPdgGamma};
    //    if(std::find(kAcceptPdgCodes.begin(), kAcceptPdgCodes.end(), TMath::Abs(particle->PdgCode())) ==  kAcceptPdgCodes.end()) continue;  // species not supported

    if (!(TMath::Abs(particle->E()-particle->Pz())>0.)) continue;
    Double_t yPre = (particle->E()+particle->Pz())/(particle->E()-particle->Pz());
    //     cout << i << "\t"<< particle->GetPdgCode() << "\t"<< particle->Pz() << "\t" << particle->Energy()<< "\t" << particle->Energy()-particle->Pz() << "\t"<< yPre << endl;
    if( yPre <= 0 ) continue;

    Double_t y = 0.5*TMath::Log(yPre);

    if (y > 1.000) continue;
    if(particle->PdgCode()==22){

      //whatever you require hasMother or not, both histogram are filled  (2023/04/13 v5)
      // if (!hasMother){//Primary
      // }else{//Secondary
      // 	fHistPtGammaHM->Fill(particle->Pt());	
      // 	fHistPtYGammaHM->Fill(particle->Pt(), particle->Y());
      // }
      if(std::abs(particle->Y()) <= 0.8){
	//MB case 
	fHistPtGamma->Fill(particle->Pt());	
	fHistPtYGamma->Fill(particle->Pt(), particle->Y());		  
	if(fNTracksInV0Acc < 7.0){
	  //Low multiplicity < 7.0 
	  fHistPtGammaClassI->Fill(particle->Pt());	
	  fHistPtYGammaClassI->Fill(particle->Pt(), particle->Y());		   
	}else if( (fNTracksInV0Acc < 30) && (7.0 <= fNTracksInV0Acc) ){
	  //7.0 <= Nch < 30.0
	  fHistPtGammaClassII->Fill(particle->Pt());	
	  fHistPtYGammaClassII->Fill(particle->Pt(), particle->Y());		   
	}else if( 30.0 <= fNTracksInV0Acc ){
	  //high multiplicity 30 <= Nch
	  fHistPtGammaClassIII->Fill(particle->Pt());	
	  fHistPtYGammaClassIII->Fill(particle->Pt(), particle->Y());		   
	}
      }	
    }//if particle is photon 
  }//Number of Tracks
}//ProcessMCParticles
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
