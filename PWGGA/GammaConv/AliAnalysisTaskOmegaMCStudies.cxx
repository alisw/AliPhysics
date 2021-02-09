/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                       *
 * Author: Baldo Sahlmueller, Friederike Bock                     *
 * Version 1.0                                 *
 *                                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its    *
 * documentation strictly for non-commercial purposes is hereby granted    *
 * without fee, provided that the above copyright notice appears in all    *
 * copies and that both the copyright notice and this permission notice    *
 * appear in the supporting documentation. The authors make no claims    *
 * about the suitability of this software for any purpose. It is      *
 * provided "as is" without express or implied warranty.               *
 **************************************************************************/

//////////////////////////////////////////////////////////////////
//----------------------------------------------------------------
// Class used to do MC studies on heavy mesons
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
#include "AliAnalysisTaskOmegaMCStudies.h"
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

ClassImp(AliAnalysisTaskOmegaMCStudies)

//________________________________________________________________________
AliAnalysisTaskOmegaMCStudies::AliAnalysisTaskOmegaMCStudies(): AliAnalysisTaskSE(),
  fOutputContainer(nullptr),
  fHistNEvents(nullptr),
  fHistXSection(nullptr),
  fHistPtHard(nullptr),
  fHistPtYOmega(nullptr),
  fHistPtYOmegaPiPiPi(nullptr),
  fHistPtYPi0(nullptr),
  fHistPtYEtaPrime(nullptr),
  fHistPtYEtaPrimeEtaPiPi(nullptr),
  fHistOmegaPtPi0Pt(nullptr),
  fHistEtaPrimePtEtaPt(nullptr),
  fIsMC(1),
  fMaxpT(20),
  fProcessCode(-1),
  fWeight(-1),
  fEventCounter(0),
  fNTotEvents(1),
  fNRejectEvents(10000)
{

}

//________________________________________________________________________
AliAnalysisTaskOmegaMCStudies::AliAnalysisTaskOmegaMCStudies(const char *name):
  AliAnalysisTaskSE(name),
  fOutputContainer(nullptr),
  fHistNEvents(nullptr),
  fHistXSection(nullptr),
  fHistPtHard(nullptr),
  fHistPtYOmega(nullptr),
  fHistPtYOmegaPiPiPi(nullptr),
  fHistPtYPi0(nullptr),
  fHistPtYEtaPrime(nullptr),
  fHistPtYEtaPrimeEtaPiPi(nullptr),
  fHistOmegaPtPi0Pt(nullptr),
  fHistEtaPrimePtEtaPt(nullptr),
  
  fIsMC(1),
  fMaxpT(20),
  fProcessCode(-1),
  fWeight(-1),
  fEventCounter(0),
  fNTotEvents(1),
  fNRejectEvents(10000)
{
  // Define output slots here
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskOmegaMCStudies::~AliAnalysisTaskOmegaMCStudies()
{

}

//________________________________________________________________________
void AliAnalysisTaskOmegaMCStudies::UserCreateOutputObjects(){

  // Create histograms
  if(fOutputContainer != nullptr){
    delete fOutputContainer;
    fOutputContainer          = nullptr;
  }
  if(fOutputContainer == nullptr){
    fOutputContainer          = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }

  fHistNEvents                		= new TH1F("NEvents", "NEvents", 3, -0.5, 2.5);
  fHistNEvents->Sumw2();
  fOutputContainer->Add(fHistNEvents);

  fHistXSection               		= new TH1D("XSection", "XSection", 1000, 0, 1e4);

  //   SetLogBinningXTH1(fHistXSection);
  fHistXSection->Sumw2();
  fOutputContainer->Add(fHistXSection);

  fHistPtHard                 		= new TH1F("PtHard", "PtHard", fMaxpT*50, 0, fMaxpT);
  fHistPtHard->Sumw2();
  fOutputContainer->Add(fHistPtHard);

// Spectra
  fHistPtYOmega                		= new TH2F("Pt_Y_Omega","Pt_Y_Omega", fMaxpT*50, 0, fMaxpT, fMaxpT*10, -1.0, 1.0);
  fHistPtYOmega->Sumw2();
  fOutputContainer->Add(fHistPtYOmega);

  fHistPtYOmegaPiPiPi                		= new TH2F("Pt_Y_OmegaPiPiPI","Pt_Y_OmegaPiPiPi", fMaxpT*50, 0, fMaxpT*10, fMaxpT, -1.0, 1.0);
  fHistPtYOmegaPiPiPi->Sumw2();
  fOutputContainer->Add(fHistPtYOmegaPiPiPi);

  fHistPtYPi0               		= new TH2F("Pt_Y_Pi0","Pt_Y_Pi0", fMaxpT*50, 0, fMaxpT, fMaxpT*10, -1.0, 1.0);
  fHistPtYPi0->Sumw2();
  fOutputContainer->Add(fHistPtYPi0);

  fHistPtYEtaPrime               		= new TH2F("Pt_Y_EtaPrime","Pt_Y_EtaPrime", fMaxpT*50, 0, fMaxpT, fMaxpT*10, -1.0, 1.0);
  fHistPtYEtaPrime->Sumw2();
  fOutputContainer->Add(fHistPtYEtaPrime);
  
  fHistPtYEtaPrimeEtaPiPi               		= new TH2F("Pt_Y_EtaPrimeEtaPiPi","Pt_Y_EtaPrimeEtaPiPi", fMaxpT*50, 0, fMaxpT, fMaxpT*10, -1.0, 1.0);
  fHistPtYEtaPrimeEtaPiPi->Sumw2();
  fOutputContainer->Add(fHistPtYEtaPrimeEtaPiPi);

  fHistOmegaPtPi0Pt                		= new TH2F("fHistOmegaPtPi0Pt","fHistOmegaPtPi0Pt", fMaxpT*50, 0, fMaxpT,  fMaxpT*50, 0, fMaxpT);
  fHistOmegaPtPi0Pt->Sumw2();
  fOutputContainer->Add(fHistOmegaPtPi0Pt);

  fHistEtaPrimePtEtaPt                		= new TH2F("fHistEtaPrimePtEtaPt","fHistEtaPrimePtEtaPt", fMaxpT*50, 0, fMaxpT,  fMaxpT*50, 0, fMaxpT);
  fHistEtaPrimePtEtaPt->Sumw2();
  fOutputContainer->Add(fHistEtaPrimePtEtaPt);
  
  PostData(1, fOutputContainer);
}

//_____________________________________________________________________________
void AliAnalysisTaskOmegaMCStudies::UserExec(Option_t *)
{
  fEventCounter++;

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

  if(pyH) {
    fProcessCode = pyH->ProcessType();
   // printf("x1 = %f x2 = %f Q = %f \n",fX1,fX2,fQFrac);
  }
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

  // calculate weight
  fWeight = xSection/(ntrials*(fNTotEvents));

  fEventCounter++;
  //printf("XSec = %f\t NTrials=%f \t  Weight = %f \n",xSection,ntrials,fWeight);
  ProcessMCParticles();


  PostData(1, fOutputContainer);
}


//________________________________________________________________________
void AliAnalysisTaskOmegaMCStudies::ProcessMCParticles()
{

  for(Long_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {

    
    // fill primary histograms
    TParticle* particle         = nullptr;
    particle                    = (TParticle *)fMCEvent->Particle(i);
    if (!particle) continue;

    switch(particle->GetPdgCode()){
    case kPdgOmega:
      fHistPtYOmega->Fill(particle->Pt(), particle->Y(),fWeight);
      if(IsPiPlusPiMinusPiZeroDecay(particle)){
        fHistPtYOmegaPiPiPi->Fill(particle->Pt(), particle->Y(),fWeight);
        TParticle* pi0 = (TParticle*) fMCEvent->Particle(ReturnPi0FromOmega(particle));
        fHistOmegaPtPi0Pt->Fill(particle->Pt(),pi0->Pt(),fWeight);
      }
      break; 
    case kPdgPi0:
      fHistPtYPi0->Fill(particle->Pt(), particle->Y(),fWeight);
      break; 
    case kPdgEtaPrime:
      fHistPtYEtaPrime->Fill(particle->Pt(), particle->Y(),fWeight);
      if(IsPiPlusPiMinusEtaDecay(particle)){
        fHistPtYEtaPrimeEtaPiPi->Fill(particle->Pt(), particle->Y(),fWeight);
        TParticle* eta = (TParticle*) fMCEvent->Particle(ReturnEtaFromEtaPrime(particle));
        fHistEtaPrimePtEtaPt->Fill(particle->Pt(),eta->Pt(),fWeight);
      }
      break; 
    }
  }
}



//________________________________________________________________________
// check if given particle decays to pi+pi-pi0
bool AliAnalysisTaskOmegaMCStudies::IsPiPlusPiMinusPiZeroDecay(TParticle* part) const{
   // check number of daughters
   Bool_t foundPi0     = kFALSE;
   Bool_t foundPiPlus  = kFALSE;
   Bool_t foundPiMinus = kFALSE;
   
   if(part->GetNDaughters()==3){
     for(Int_t i = part->GetFirstDaughter(); i <= part->GetLastDaughter();i++){
         TParticle* daughter = (TParticle*) fMCEvent->Particle(i);
         switch(daughter->GetPdgCode()){
           case kPdgPi0:     foundPi0      = kTRUE; break;
           case kPdgPiPlus:  foundPiPlus   = kTRUE; break;
           case kPdgPiMinus: foundPiMinus  = kTRUE; break;
         }  
     }
     if(foundPi0 && foundPiPlus && foundPiMinus){
       return true;
     } else{
       return false;
     }
   } else{
     return false;
   }
}

//________________________________________________________________________
// check if given particle decays to pi+pi-eta
bool AliAnalysisTaskOmegaMCStudies::IsPiPlusPiMinusEtaDecay(TParticle* part) const{
   // check number of daughters
   Bool_t foundEta     = kFALSE;
   Bool_t foundPiPlus  = kFALSE;
   Bool_t foundPiMinus = kFALSE;
   
   if(part->GetNDaughters()==3){
     for(Int_t i = part->GetFirstDaughter(); i <= part->GetLastDaughter();i++){
         TParticle* daughter = (TParticle*) fMCEvent->Particle(i);
         switch(daughter->GetPdgCode()){
           case kPdgEta:     foundEta      = kTRUE; break;
           case kPdgPiPlus:  foundPiPlus   = kTRUE; break;
           case kPdgPiMinus: foundPiMinus  = kTRUE; break;
         }  
     }
     if(foundEta && foundPiPlus && foundPiMinus){
       return true;
     } else{
       return false;
     }
   } else{
     return false;
   }
}

//________________________________________________________________________
// return stack position pi0 from omega->pi+pi-pi0 decay
Int_t AliAnalysisTaskOmegaMCStudies::ReturnPi0FromOmega(TParticle* part){
   // check number of daughters
   if(part->GetNDaughters()==3){
     for(Int_t i = part->GetFirstDaughter(); i <= part->GetLastDaughter();i++){
       TParticle* daughter = (TParticle*) fMCEvent->Particle(i);
       if(   (TMath::Abs(daughter->GetPdgCode())== kPdgPi0)){
           return i;
       } 

    }
    return -1;
   } else{
     return -1;
   }
}


//________________________________________________________________________
// return stack position of eta from eta'->pi+pi-eta decay
Int_t AliAnalysisTaskOmegaMCStudies::ReturnEtaFromEtaPrime(TParticle* part){
   // check number of daughters
   if(part->GetNDaughters()==3){
     for(Int_t i = part->GetFirstDaughter(); i <= part->GetLastDaughter();i++){
       TParticle* daughter = (TParticle*) fMCEvent->Particle(i);       
       if(   (TMath::Abs(daughter->GetPdgCode())== kPdgEta)){
           return i;
       } 

    }
    return -1;
   } else{
     return -1;
   }
}

//________________________________________________________________________
bool AliAnalysisTaskOmegaMCStudies::IsInPCMAcceptance(TParticle* part) const {
  const Double_t kBoundaryEta = 0.900001;
  if (//part->Pt() > 0.050 
  //&& 
  TMath::Abs(part->Eta()) < kBoundaryEta) return true;

  return false;
}

//________________________________________________________________________
bool AliAnalysisTaskOmegaMCStudies::IsInPHOSAcceptance(TParticle* part) const {
  const Double_t kBoundaryEtaMin = -0.13;
  const Double_t kBoundaryEtaMax = 0.13;
  const Double_t kBoundaryPhiMin = 4.54;
  const Double_t kBoundaryPhiMax = 5.59;
  //if (part->Pt() < 0.300) return false;
  if (part->Eta() > kBoundaryEtaMax || part->Eta() < kBoundaryEtaMin) return false;
  if (part->Phi() > kBoundaryPhiMax || part->Phi() < kBoundaryPhiMin) return false;
  return true;
}

//________________________________________________________________________
bool AliAnalysisTaskOmegaMCStudies::IsInEMCalAcceptance(TParticle* part) const {
  const Double_t kBoundaryEtaMin = -0.6687;
  const Double_t kBoundaryEtaMax = 0.66465;
  const Double_t kBoundaryPhiMin = 1.39626;
  const Double_t kBoundaryPhiMax = 3.15;
  //if (part->Pt() < 0.400) return false;
  if (part->Eta() > kBoundaryEtaMax || part->Eta() < kBoundaryEtaMin) return false;
  if (part->Phi() > kBoundaryPhiMax || part->Phi() < kBoundaryPhiMin) return false;
  return true;
}

//________________________________________________________________________
bool AliAnalysisTaskOmegaMCStudies::IsInFOCALAcceptance(TParticle* part) const {
  const Double_t kBoundaryEtaMin = 3.2;
  const Double_t kBoundaryEtaMax = 5.3;
  //if (part->Pt() < 0.400) return false;
  if (part->Eta() > kBoundaryEtaMax || part->Eta() < kBoundaryEtaMin) return false;
  return true;
}

//________________________________________________________________________
bool AliAnalysisTaskOmegaMCStudies::IsInLHCbAcceptance(TParticle* part) const {
  const Double_t kBoundaryEtaMin = 1.9;
  const Double_t kBoundaryEtaMax = 5.1;
  //if (part->Pt() < 0.400) return false;
  if (part->Eta() > kBoundaryEtaMax || part->Eta() < kBoundaryEtaMin) return false;
  return true;
}

//________________________________________________________________________
bool AliAnalysisTaskOmegaMCStudies::IsInMidAcceptance(TParticle* part) const {
  const Double_t kBoundaryEtaMin = -1;
  const Double_t kBoundaryEtaMax = 1.;
  //if (part->Pt() < 0.400) return false;
  if (part->Eta() > kBoundaryEtaMax || part->Eta() < kBoundaryEtaMin) return false;
  return true;
}
//________________________________________________________________________
void AliAnalysisTaskOmegaMCStudies::Terminate(const Option_t *)
{

  //fOutputContainer->Print(); // Will crash on GRID
}


//_________________________________________________________________________________
void AliAnalysisTaskOmegaMCStudies::SetLogBinningXTH1(TH1* histoRebin){
  TAxis *axisafter = histoRebin->GetXaxis();
  Int_t bins = axisafter->GetNbins();
  Double_t from = axisafter->GetXmin();
  Double_t to = axisafter->GetXmax();
  Double_t *newbins = new Double_t[bins+1];
  newbins[0] = from;
  Double_t factor = TMath::Power(to/from, 1./bins);
  for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
  axisafter->Set(bins, newbins);
  delete [] newbins;
}

//_________________________________________________________________________________
void AliAnalysisTaskOmegaMCStudies::SetLogBinningXTH2(TH2* histoRebin){
  TAxis *axisafter = histoRebin->GetXaxis();
  Int_t bins = axisafter->GetNbins();
  Double_t from = axisafter->GetXmin();
  Double_t to = axisafter->GetXmax();
  Double_t *newbins = new Double_t[bins+1];
  newbins[0] = from;
  Double_t factor = TMath::Power(to/from, 1./bins);
  for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
  axisafter->Set(bins, newbins);
  delete [] newbins;
}

//_________________________________________________________________________________
void AliAnalysisTaskOmegaMCStudies::SetLogBinningYTH2(TH2* histoRebin){
  TAxis *axisafter = histoRebin->GetYaxis();
  Int_t bins = axisafter->GetNbins();
  Double_t from = axisafter->GetXmin();
  Double_t to = axisafter->GetXmax();
  Double_t *newbins = new Double_t[bins+1];
  newbins[0] = from;
  Double_t factor = TMath::Power(to/from, 1./bins);
  for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
  axisafter->Set(bins, newbins);
  delete [] newbins;
}
