/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Baldo Sahlmueller, Friederike Bock                             *
 * Modified by Hikari Murakami                                            *
 * Version 1.0                                                            *
 *                                                                        *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//----------------------------------------------------------------
// Class used to obtain Particle ratio
//----------------------------------------------------------------
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
#include "AliAnalysisTaskLMeePureMC.h"
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

ClassImp(AliAnalysisTaskLMeePureMC)

//________________________________________________________________________
AliAnalysisTaskLMeePureMC::AliAnalysisTaskLMeePureMC(): AliAnalysisTaskSE(),
  fOutputContainer(nullptr),
  fHistNEvents(nullptr),
  fHistXSection(nullptr),
  fHistPtHard(nullptr),
  fHistPtYPi0(nullptr),
  fHistPtYPiPl(nullptr),
  fHistPtYPiMi(nullptr),
  fHistPtYEta(nullptr),
  fHistPtYEtaPrime(nullptr),
  fHistPtYOmega(nullptr),
  fHistPtYRho0(nullptr),
  fHistPtYRhoPl(nullptr),
  fHistPtYRhoMi(nullptr),
  fHistPtYPhi(nullptr),
  fHistPtYJPsi(nullptr),
  fHistPtYSigma0(nullptr),
  fHistPtYK0s(nullptr),
  fHistPtYK0l(nullptr),
  fHistPtYK0star(nullptr),
  fHistPtYDeltaPlPl(nullptr),
  fHistPtYDeltaPl(nullptr),
  fHistPtYDeltaMi(nullptr),
  fHistPtYDelta0(nullptr),
  fHistPtYLambda(nullptr),
  fHistPtYKPl(nullptr),
  fHistPtYKMi(nullptr),
  fHistPtYPi0FromEta(nullptr),
  fHistPtYPi0FromLambda(nullptr),
  fHistPtYPi0FromK(nullptr),
  fHistPtYPiPlFromK(nullptr),
  fHistPtYPiPlFromEta(nullptr),
  fHistPtYPiMiFromK(nullptr),
  fHistPtYPiMiFromEta(nullptr),
  fIsMC(1)
{

}

//________________________________________________________________________
AliAnalysisTaskLMeePureMC::AliAnalysisTaskLMeePureMC(const char *name):
  AliAnalysisTaskSE(name),
  fOutputContainer(nullptr),
  fHistNEvents(nullptr),
  fHistXSection(nullptr),
  fHistPtHard(nullptr),
  fHistPtYPi0(nullptr),
  fHistPtYPiPl(nullptr),
  fHistPtYPiMi(nullptr),
  fHistPtYEta(nullptr),
  fHistPtYEtaPrime(nullptr),
  fHistPtYOmega(nullptr),
  fHistPtYRho0(nullptr),
  fHistPtYRhoPl(nullptr),
  fHistPtYRhoMi(nullptr),
  fHistPtYPhi(nullptr),
  fHistPtYJPsi(nullptr),
  fHistPtYSigma0(nullptr),
  fHistPtYK0s(nullptr),
  fHistPtYK0l(nullptr),
  fHistPtYK0star(nullptr),
  fHistPtYDeltaPlPl(nullptr),
  fHistPtYDeltaPl(nullptr),
  fHistPtYDeltaMi(nullptr),
  fHistPtYDelta0(nullptr),
  fHistPtYLambda(nullptr),
  fHistPtYKPl(nullptr),
  fHistPtYKMi(nullptr),
  fHistPtYPi0FromEta(nullptr),
  fHistPtYPi0FromLambda(nullptr),
  fHistPtYPi0FromK(nullptr),
  fHistPtYPiPlFromK(nullptr),
  fHistPtYPiPlFromEta(nullptr),
  fHistPtYPiMiFromK(nullptr),
  fHistPtYPiMiFromEta(nullptr),
  fIsMC(1)
{
  // Define output slots here
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskLMeePureMC::~AliAnalysisTaskLMeePureMC()
{

}

//________________________________________________________________________
void AliAnalysisTaskLMeePureMC::UserCreateOutputObjects(){

  // Create histograms
  if(fOutputContainer != nullptr){
    delete fOutputContainer;
    fOutputContainer          = nullptr;
  }
  if(fOutputContainer == nullptr){
    fOutputContainer          = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }

  fHistNEvents          = new TH1F("NEvents", "",               3, -0.5, 2.5);
  fHistXSection         = new TH1D("XSection", "",           1000000, 0, 1e4);
  fHistPtHard           = new TH1F("PtHard", "",                 400, 0, 200);

  //0 < pt < 20 2000 bin --> per 10MeV/c
  fHistPtYPi0           = new TH2F("Pt_Y_Pi0","",               2000, 0, 20, 200, -1.0, 1.0);  
  fHistPtYPiPl     	= new TH2F("Pt_Y_PiPl","",              2000, 0, 20, 200, -1.0, 1.0);
  fHistPtYPiMi         	= new TH2F("Pt_Y_PiMi","",              2000, 0, 20, 200, -1.0, 1.0);
  fHistPtYEta           = new TH2F("Pt_Y_Eta","",               2000, 0, 20, 200, -1.0, 1.0);
  fHistPtYEtaPrime     	= new TH2F("Pt_Y_EtaPrime","",          2000, 0, 20, 200, -1.0, 1.0);
  fHistPtYOmega        	= new TH2F("Pt_Y_Omega","",             2000, 0, 20, 200, -1.0, 1.0);
  fHistPtYRho0         	= new TH2F("Pt_Y_Rho0","",              2000, 0, 20, 200, -1.0, 1.0);
  fHistPtYRhoPl         = new TH2F("Pt_Y_RhoPl","",             2000, 0, 20, 200, -1.0, 1.0);
  fHistPtYRhoMi         = new TH2F("Pt_Y_RhoMi","",             2000, 0, 20, 200, -1.0, 1.0);
  fHistPtYPhi         	= new TH2F("Pt_Y_Phi","",               2000, 0, 20, 200, -1.0, 1.0); 
  fHistPtYJPsi        	= new TH2F("Pt_Y_JPsi","",              2000, 0, 20, 200, -1.0, 1.0);
  fHistPtYSigma0        = new TH2F("Pt_Y_Sigma0","",            2000, 0, 20, 200, -1.0, 1.0);
  fHistPtYK0s         	= new TH2F("Pt_Y_K0s","",               2000, 0, 20, 200, -1.0, 1.0); 
  fHistPtYK0l         	= new TH2F("Pt_Y_K0l","",               2000, 0, 20, 200, -1.0, 1.0); 
  fHistPtYK0star       	= new TH2F("Pt_Y_K0star","",            2000, 0, 20, 200, -1.0, 1.0);
  fHistPtYDeltaPlPl    	= new TH2F("Pt_Y_DeltaPlPl","",         2000, 0, 20, 200, -1.0, 1.0); 
  fHistPtYDeltaPl      	= new TH2F("Pt_Y_DeltaPl","",           2000, 0, 20, 200, -1.0, 1.0); 
  fHistPtYDeltaMi       = new TH2F("Pt_Y_DeltaMi","",           2000, 0, 20, 200, -1.0, 1.0); 
  fHistPtYDelta0      	= new TH2F("Pt_Y_Delta0","",            2000, 0, 20, 200, -1.0, 1.0);
  fHistPtYLambda      	= new TH2F("Pt_Y_Lambda","",            2000, 0, 20, 200, -1.0, 1.0);
  fHistPtYKPl         	= new TH2F("Pt_Y_KPl","",               2000, 0, 20, 200, -1.0, 1.0); 
  fHistPtYKMi           = new TH2F("Pt_Y_KMi","",               2000, 0, 20, 200, -1.0, 1.0); 
  fHistPtYPi0FromEta   	= new TH2F("Pt_Y_Pi0FromEta","",        2000, 0, 20, 200, -1.0, 1.0);
  fHistPtYPi0FromLambda = new TH2F("Pt_Y_Pi0FromLambda","",     2000, 0, 20, 200, -1.0, 1.0);
  fHistPtYPi0FromK     	= new TH2F("Pt_Y_Pi0FromK","",          2000, 0, 20, 200, -1.0, 1.0);
  fHistPtYPiPlFromK    	= new TH2F("Pt_Y_PiPlFromK","",         2000, 0, 20, 200, -1.0, 1.0);
  fHistPtYPiPlFromEta  	= new TH2F("Pt_Y_PiPlFromEta","",       2000, 0, 20, 200, -1.0, 1.0);
  fHistPtYPiMiFromK    	= new TH2F("Pt_Y_PiMiFromK","",         2000, 0, 20, 200, -1.0, 1.0);
  fHistPtYPiMiFromEta 	= new TH2F("Pt_Y_PiMiFromEta","",       2000, 0, 20, 200, -1.0, 1.0);

  fHistNEvents                  ->Sumw2();
  fHistXSection                 ->Sumw2();
  fHistPtHard                   ->Sumw2();
  fHistPtYPi0        		->Sumw2();
  fHistPtYPiPl        		->Sumw2();
  fHistPtYPiMi         		->Sumw2();
  fHistPtYEta         		->Sumw2();
  fHistPtYEtaPrime     		->Sumw2();
  fHistPtYOmega        		->Sumw2();
  fHistPtYRho0         		->Sumw2();
  fHistPtYRhoPl        		->Sumw2();
  fHistPtYRhoMi       		->Sumw2();
  fHistPtYPhi         		->Sumw2();
  fHistPtYJPsi        		->Sumw2();
  fHistPtYSigma0       		->Sumw2();
  fHistPtYK0s          		->Sumw2();
  fHistPtYK0l         		->Sumw2();
  fHistPtYK0star       		->Sumw2();
  fHistPtYDeltaPlPl    		->Sumw2();
  fHistPtYDeltaPl      		->Sumw2();
  fHistPtYDeltaMi     		->Sumw2();
  fHistPtYDelta0       		->Sumw2();
  fHistPtYLambda       		->Sumw2();
  fHistPtYKPl         		->Sumw2();
  fHistPtYKMi         		->Sumw2();
  fHistPtYPi0FromEta   		->Sumw2();
  fHistPtYPi0FromLambda		->Sumw2();
  fHistPtYPi0FromK        	->Sumw2();
  fHistPtYPiPlFromK             ->Sumw2();
  fHistPtYPiPlFromEta  		->Sumw2();
  fHistPtYPiMiFromK             ->Sumw2();
  fHistPtYPiMiFromEta  		->Sumw2();

  fOutputContainer->Add(fHistNEvents);
  fOutputContainer->Add(fHistXSection);
  fOutputContainer->Add(fHistPtHard);
  fOutputContainer->Add(fHistPtYPi0);
  fOutputContainer->Add(fHistPtYPiPl);
  fOutputContainer->Add(fHistPtYPiMi);
  fOutputContainer->Add(fHistPtYEta);  
  fOutputContainer->Add(fHistPtYEtaPrime);
  fOutputContainer->Add(fHistPtYOmega);
  fOutputContainer->Add(fHistPtYRho0);
  fOutputContainer->Add(fHistPtYRhoPl);
  fOutputContainer->Add(fHistPtYRhoMi);
  fOutputContainer->Add(fHistPtYPhi);
  fOutputContainer->Add(fHistPtYJPsi);
  fOutputContainer->Add(fHistPtYSigma0);
  fOutputContainer->Add(fHistPtYK0s);
  fOutputContainer->Add(fHistPtYK0l);
  fOutputContainer->Add(fHistPtYK0star);
  fOutputContainer->Add(fHistPtYDeltaPlPl);
  fOutputContainer->Add(fHistPtYDeltaPl);
  fOutputContainer->Add(fHistPtYDeltaMi);
  fOutputContainer->Add(fHistPtYDelta0);
  fOutputContainer->Add(fHistPtYLambda);
  fOutputContainer->Add(fHistPtYKPl);
  fOutputContainer->Add(fHistPtYKMi);
  fOutputContainer->Add(fHistPtYPi0FromEta);
  fOutputContainer->Add(fHistPtYPi0FromLambda);
  fOutputContainer->Add(fHistPtYPi0FromK);
  fOutputContainer->Add(fHistPtYPiPlFromK);
  fOutputContainer->Add(fHistPtYPiPlFromEta);
  fOutputContainer->Add(fHistPtYPiMiFromK);
  fOutputContainer->Add(fHistPtYPiMiFromEta);
  
  PostData(1, fOutputContainer);
}

//_____________________________________________________________________________
void AliAnalysisTaskLMeePureMC::UserExec(Option_t *)
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

  ProcessMCParticles();


  PostData(1, fOutputContainer);
}


//________________________________________________________________________
void AliAnalysisTaskLMeePureMC::ProcessMCParticles()
{
  // Loop over all primary MC particle
  for(Long_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {
    // fill primary histograms
    TParticle* particle         = nullptr;
    particle                    = (TParticle *)fMCEvent->Particle(i);
    if (!particle) continue;
    Bool_t hasMother            = kFALSE;
    //     cout << i << "\t"<< particle->GetMother(0) << endl;
    if (particle->GetMother(0)>-1)
      hasMother                 = kTRUE;
    TParticle* motherParticle   = nullptr;
    if( hasMother )
      motherParticle            = (TParticle *)fMCEvent->Particle(particle->GetMother(0));
    if (motherParticle)
      hasMother                 = kTRUE;
    else
      hasMother                 = kFALSE;

    const std::array<int, 19> kAcceptPdgCodes = {kPdgPi0, kPdgEta, kPdgEtaPrime, kPdgOmega, kPdgPiPlus, kPdgRho0, kPdgPhi, kPdgJPsi, kPdgSigma0, kPdgK0Short, kPdgDeltaPlus, kPdgDeltaPlusPlus, kPdgDeltaMinus, kPdgDelta0, kPdgRhoPlus, kPdgKStar, kPdgK0Long, kPdgLambda, kPdgKPlus};
    if(std::find(kAcceptPdgCodes.begin(), kAcceptPdgCodes.end(), TMath::Abs(particle->GetPdgCode())) ==  kAcceptPdgCodes.end()) continue;  // species not supported

    if (!(TMath::Abs(particle->Energy()-particle->Pz())>0.)) continue;
    Double_t yPre = (particle->Energy()+particle->Pz())/(particle->Energy()-particle->Pz());
    //     cout << i << "\t"<< particle->GetPdgCode() << "\t"<< particle->Pz() << "\t" << particle->Energy()<< "\t" << particle->Energy()-particle->Pz() << "\t"<< yPre << endl;
    if( yPre <= 0 ) continue;

    Double_t y = 0.5*TMath::Log(yPre);
    
    if (TMath::Abs(y) > 1.000) continue;
    switch(particle->GetPdgCode()){
    case kPdgPi0:
      fHistPtYPi0->Fill(particle->Pt(), particle->Y());
      if (hasMother){
        if (TMath::Abs(motherParticle->GetPdgCode()) == kPdgK0Short ||
            TMath::Abs(motherParticle->GetPdgCode()) == kPdgK0Long ||
            TMath::Abs(motherParticle->GetPdgCode()) == kPdgKPlus
	    )
          fHistPtYPi0FromK->Fill(particle->Pt(), particle->Y());
        if (TMath::Abs(motherParticle->GetPdgCode()) == kPdgLambda)
          fHistPtYPi0FromLambda->Fill(particle->Pt(), particle->Y());
        if (motherParticle->GetPdgCode() == kPdgEta)
          fHistPtYPi0FromEta->Fill(particle->Pt(), particle->Y());
      }
      break;
    case kPdgEta:
      fHistPtYEta->Fill(particle->Pt(), particle->Y());
      break;
    case kPdgEtaPrime:
      fHistPtYEtaPrime->Fill(particle->Pt(), particle->Y());
      break;
    case kPdgOmega:
      fHistPtYOmega->Fill(particle->Pt(), particle->Y());
      break;
    case kPdgPiPlus:
      fHistPtYPiPl->Fill(particle->Pt(), particle->Y());
      if (hasMother){
        if (TMath::Abs(motherParticle->GetPdgCode()) == kPdgK0Short ||
            TMath::Abs(motherParticle->GetPdgCode()) == kPdgK0Long ||
            TMath::Abs(motherParticle->GetPdgCode()) == kPdgKPlus
	    )
          fHistPtYPiPlFromK->Fill(particle->Pt(), particle->Y());
        if (motherParticle->GetPdgCode() == kPdgEta)
          fHistPtYPiPlFromEta->Fill(particle->Pt(), particle->Y());
      }
      break;
    case kPdgPiMinus:
      fHistPtYPiMi->Fill(particle->Pt(), particle->Y());
      if (hasMother){
        if (TMath::Abs(motherParticle->GetPdgCode()) == kPdgK0Short ||
            TMath::Abs(motherParticle->GetPdgCode()) == kPdgK0Long ||
            TMath::Abs(motherParticle->GetPdgCode()) == kPdgKPlus
	    )
          fHistPtYPiMiFromK->Fill(particle->Pt(), particle->Y());
        if (motherParticle->GetPdgCode() == kPdgEta)
          fHistPtYPiMiFromEta->Fill(particle->Pt(), particle->Y());
      }
      break;
    case kPdgRho0:
      fHistPtYRho0->Fill(particle->Pt(), particle->Y());
      break;
    case kPdgRhoPlus:
      fHistPtYRhoPl->Fill(particle->Pt(), particle->Y());
      break;
    case kPdgRhoMinus:
      fHistPtYRhoMi->Fill(particle->Pt(), particle->Y());
      break;
    case kPdgPhi:
      fHistPtYPhi->Fill(particle->Pt(), particle->Y());
      break;
    case kPdgJPsi:
      fHistPtYJPsi->Fill(particle->Pt(), particle->Y());
      break;
    case 3212:
      fHistPtYSigma0->Fill(particle->Pt(), particle->Y());
      break;
    case kPdgK0Short:
      fHistPtYK0s->Fill(particle->Pt(), particle->Y());
      break;
    case kPdgK0Long:
      fHistPtYK0l->Fill(particle->Pt(), particle->Y());
      break;
    case kPdgKStar:
      fHistPtYK0star->Fill(particle->Pt(), particle->Y());
      break;
    case kPdgDeltaPlusPlus:
      fHistPtYDeltaPlPl->Fill(particle->Pt(), particle->Y());
      break;
    case kPdgDeltaPlus:
      fHistPtYDeltaPl->Fill(particle->Pt(), particle->Y());
      break;
    case kPdgDeltaMinus:
      fHistPtYDeltaMi->Fill(particle->Pt(), particle->Y());
      break;
    case kPdgDelta0:
      fHistPtYDelta0->Fill(particle->Pt(), particle->Y());
      break;
    case kPdgLambda:
      fHistPtYLambda->Fill(particle->Pt(), particle->Y());
      break;
    case kPdgKPlus:
      fHistPtYKPl->Fill(particle->Pt(), particle->Y());
      break;
    case kPdgKMinus:
      fHistPtYKMi->Fill(particle->Pt(), particle->Y());
      break;
    }

  }// Loop over all primary MC particle
}

//________________________________________________________________________
void AliAnalysisTaskLMeePureMC::Terminate(const Option_t *)
{

  //fOutputContainer->Print(); // Will crash on GRID
}
