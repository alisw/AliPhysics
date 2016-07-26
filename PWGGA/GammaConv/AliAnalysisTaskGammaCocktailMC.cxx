/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.  *
*                                                                         *
* Author: Friederike Bock, Mike Sas                                       *
* Version 1.0                                                             *
*                                                                         *
*                                                                         *
* Permission to use, copy, modify and distribute this software and its    *
* documentation strictly for non-commercial purposes is hereby granted    *
* without fee, provided that the above copyright notice appears in all    *
* copies and that both the copyright notice and this permission notice    *
* appear in the supporting documentation. The authors make no claims      *
* about the suitability of this software for any purpose. It is           *
* provided "as is" without express or implied warranty.                   *
**************************************************************************/

//////////////////////////////////////////////////////////////////
//----------------------------------------------------------------
// Class used to do analysis on electromagnetic cocktail output
//----------------------------------------------------------------
//////////////////////////////////////////////////////////////////
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
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
#include "AliStack.h"
#include "AliAnalysisTaskGammaCocktailMC.h"
#include "AliVParticle.h"
#include "AliEventplane.h"
#include "AliInputEventHandler.h"
#include <vector>
#include <map>

ClassImp(AliAnalysisTaskGammaCocktailMC)

//________________________________________________________________________
AliAnalysisTaskGammaCocktailMC::AliAnalysisTaskGammaCocktailMC(): AliAnalysisTaskSE(),
  fOutputContainer(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fMCStack(NULL),
  fHistNEvents(NULL),
  fHistPtYGamma(NULL),
  fHistPhiGamma(NULL),
  fHistPhiInput(NULL),
  fHistPtYInput(NULL),
  fHistPtYGammaSource(NULL),
  fHistPtYInputRest(NULL),
  fHistPtYGammaSourceRest(NULL),
  fParticleList(NULL),
  fParticleListNames(NULL),
  fIsMC(1),
  fMaxY(2)
{
  
}

//________________________________________________________________________
AliAnalysisTaskGammaCocktailMC::AliAnalysisTaskGammaCocktailMC(const char *name):
  AliAnalysisTaskSE(name),
  fOutputContainer(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fMCStack(NULL),
  fHistNEvents(NULL),
  fHistPtYGamma(NULL),
  fHistPhiGamma(NULL),
  fHistPhiInput(NULL),
  fHistPtYInput(NULL),
  fHistPtYGammaSource(NULL),
  fHistPtYInputRest(NULL),
  fHistPtYGammaSourceRest(NULL),
  fParticleList(NULL),
  fParticleListNames(NULL),
  fIsMC(1),
  fMaxY(2)
{
  // Define output slots here
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskGammaCocktailMC::~AliAnalysisTaskGammaCocktailMC()
{

}

//________________________________________________________________________
void AliAnalysisTaskGammaCocktailMC::UserCreateOutputObjects(){
  
  // Create histograms
  if(fOutputContainer != NULL){
    delete fOutputContainer;
    fOutputContainer          = NULL;
  }
  if(fOutputContainer == NULL){
    fOutputContainer          = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }
  
  fHistNEvents                = new TH1F("NEvents", "NEvents", 1, 0, 1);
  fHistNEvents->Sumw2();
  fOutputContainer->Add(fHistNEvents);
  
  fHistPtYGamma                 = new TH2F("Pt_Y_Gamma","Pt_Y_Gamma", 500,0, 50, 400, -2.0, 2.0);
  fHistPtYGamma->Sumw2();
  fOutputContainer->Add(fHistPtYGamma);
  
  fHistPhiGamma = new TH2F("Pt_Phi_Gamma","Pt_Phi_Gamma", 500,0, 50, 100,0,7);
  fHistPhiGamma->Sumw2();
  fOutputContainer->Add(fHistPhiGamma);
  
//   "22" Gamma
//   "111" Pi0
//   "221" Eta
//   "331" EtaPrim
//   "223" Omega
//   "211" Pi+
//   "-211" Pi-
//   "113" rho0
//   "213" rho+
//   "-213" rho-
//   "333" phi
//   "443" J/psi
//   "1114" Delta-
//   "2114" Delta0
//   "2214" Delta+
//   "2224" Delta++
//   "3212" Sigma0
  const Int_t nInputParticles = 14;
  Int_t   fParticleList_local[] = {111,221,331,223,113,213,-213,333,443,1114,2114,2214,2224,3212};
  TString fParticleListNames_local[] = {"Pi0","Eta","EtaPrim","omega","rho0","rho+","rho-","phi","J/psi","Delta-","Delta0","Delta+","Delta++","Sigma0"};  
  fParticleList       = fParticleList_local;
  fParticleListNames  = fParticleListNames_local;
  fHistPtYInput       = new TH2F*[nInputParticles];
  fHistPtYGammaSource = new TH2F*[nInputParticles];
  fHistPhiInput       = new TH2F*[nInputParticles];
  
  for(Int_t i=0; i<nInputParticles; i++){
    fHistPtYInput[i] = new TH2F(Form("Pt_Y_%s",fParticleListNames[i].Data()),Form("Pt_Y_%s",fParticleListNames[i].Data()), 500,0, 50, 400, -2.0, 2.0);
    fHistPtYInput[i]->Sumw2();
    fOutputContainer->Add(fHistPtYInput[i]);
    
    //Gammas from certain mother
    fHistPtYGammaSource[i] = new TH2F(Form("Pt_Y_Gamma_From_%s",fParticleListNames[i].Data()),Form("Pt_Y_Gamma_From_%s",fParticleListNames[i].Data()), 500,0, 50, 400, -2.0, 2.0);
    fHistPtYGammaSource[i]->Sumw2();
    fOutputContainer->Add(fHistPtYGammaSource[i]);
    
    //phi distributions
    fHistPhiInput[i] = new TH2F(Form("Pt_Phi_%s",fParticleListNames[i].Data()),Form("Pt_Phi_%s",fParticleListNames[i].Data()), 500,0, 50, 100,0,7);
    fHistPhiInput[i]->Sumw2();
    fOutputContainer->Add(fHistPhiInput[i]);
  }
  
  fHistPtYInputRest = new TH1I("Pdg_primary_rest","Pdg_primary_rest", 5000,0, 5000);
  fHistPtYInputRest->Sumw2();
  fOutputContainer->Add(fHistPtYInputRest);
  
  //Gammas from certain mother
  fHistPtYGammaSourceRest = new TH1I("Pdg_Gamma_From_rest","Pdg_Gamma_From_rest", 5000,0, 5000);
  fHistPtYGammaSourceRest->Sumw2();
  fOutputContainer->Add(fHistPtYGammaSourceRest);

  PostData(1, fOutputContainer);
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaCocktailMC::UserExec(Option_t *)
{

  fInputEvent = InputEvent();
//   cout << "I found an Event" << endl;
  
  fMCEvent = MCEvent();
  if(fMCEvent == NULL) fIsMC = 0;
  
  if (fIsMC==0) return;
//   cout << "I found an MC header" << endl;
    
  fMCStack = fMCEvent->Stack();
  if(fMCStack == NULL) fIsMC = 0;
  if (fIsMC==0) return;
  
  fHistNEvents->Fill(0.5);
//   cout << "the stack is intact" << endl;
  ProcessMCParticles();

  PostData(1, fOutputContainer);
}


//________________________________________________________________________
void AliAnalysisTaskGammaCocktailMC::ProcessMCParticles(){

  // Loop over all primary MC particle  
  for(UInt_t i = 0; i < fMCStack->GetNtrack(); i++) {
    // fill primary histograms
    TParticle* particle         = NULL;
    particle                    = (TParticle *)fMCStack->Particle(i);
    if (!particle) continue;
    Bool_t hasMother            = kFALSE;
    Bool_t particleIsPrimary    = kTRUE;
//     cout << i << "\t"<< particle->GetMother(0) << endl;
    if (particle->GetMother(0)>-1){
      hasMother = kTRUE;
      particleIsPrimary = kFALSE;
    }
    TParticle* motherParticle   = NULL;
    if( hasMother ) motherParticle = (TParticle *)fMCStack->Particle(particle->GetMother(0));
    if (motherParticle){
      hasMother                 = kTRUE;
    }else{
      hasMother                 = kFALSE;
    }
    
    Bool_t motherIsPrimary    = kFALSE;
    if(hasMother){
      if(motherParticle->GetMother(0)>-1)motherIsPrimary = kTRUE;
    }

//     if (!(abs(particle->GetPdgCode()) == 111 || abs(particle->GetPdgCode()) == 221 || abs(particle->GetPdgCode()) == 331 ||
//       abs(particle->GetPdgCode()) == 223 || abs(particle->GetPdgCode()) == 211 )  )
//       continue;
    
    if (!(fabs(particle->Energy()-particle->Pz())>0.)) continue;
    Double_t yPre = (particle->Energy()+particle->Pz())/(particle->Energy()-particle->Pz());
//     cout << i << "\t"<< particle->GetPdgCode() << "\t"<< particle->Pz() << "\t" << particle->Energy()<< "\t" << particle->Energy()-particle->Pz() << "\t"<< yPre << endl;
    if (yPre == 0.) continue;
    
    Double_t y = 0.5*TMath::Log(yPre); 
    if (fabs(y) > fMaxY) continue;
    
    if(particle->GetPdgCode()==22 && hasMother==kTRUE){
      if(motherIsPrimary || !IsMotherInList(motherParticle)){
        fHistPtYGamma->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
        fHistPhiGamma->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
        switch(motherParticle->GetPdgCode()){
        case 111:
          fHistPtYGammaSource[0]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          break;
        case 221:
          fHistPtYGammaSource[1]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          break;
        case 331:
          fHistPtYGammaSource[2]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          break;
        case 223:
          fHistPtYGammaSource[3]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          break;
        case 113:
          fHistPtYGammaSource[4]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          break;
        case 213:
          fHistPtYGammaSource[5]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          break;
        case -213:
          fHistPtYGammaSource[6]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          break;
        case 333:
          fHistPtYGammaSource[7]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          break;
        case 443:
          fHistPtYGammaSource[8]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          break;
        case 1114:
          fHistPtYGammaSource[9]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          break;
        case 2114:
          fHistPtYGammaSource[10]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          break;
        case 2214:
          fHistPtYGammaSource[11]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          break;
        case 2224:
          fHistPtYGammaSource[12]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          break;
        case 3212:
          fHistPtYGammaSource[13]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          break;
          default:
            fHistPtYGammaSourceRest->Fill(motherParticle->GetPdgCode());
            break;
        }
      }
    }
    
//     fParticleList_local[] = {111,221,331,223,113,213,-213,333,443,1114,2114,2214,2224,3212};
//     fParticleListNames_local[] = {"Pi0","Eta","EtaPrim","omega","rho0","rho+","rho-","phi","J/psi","Delta-","Delta0","Delta+","Delta++","Sigma0"};
    if(particle->GetPdgCode()!=22 && particleIsPrimary){
      switch(particle->GetPdgCode()){
        case 111:
          fHistPtYInput[0]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPhiInput[0]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          break;
        case 221:
          fHistPtYInput[1]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPhiInput[1]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          break;
        case 331:
          fHistPtYInput[2]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPhiInput[2]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          break;
        case 223:
          fHistPtYInput[3]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPhiInput[3]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          break;
        case 113:
          fHistPtYInput[4]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPhiInput[4]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          break;
        case 213:
          fHistPtYInput[5]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPhiInput[5]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          break;
        case -213:
          fHistPtYInput[6]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPhiInput[6]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          break;
        case 333:
          fHistPtYInput[7]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPhiInput[7]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          break;
        case 443:
          fHistPtYInput[8]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPhiInput[8]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          break;
        case 1114:
          fHistPtYInput[9]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPhiInput[9]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          break;
        case 2114:
          fHistPtYInput[10]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPhiInput[10]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          break;
        case 2214:
          fHistPtYInput[11]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPhiInput[11]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          break;
        case 2224:
          fHistPtYInput[12]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPhiInput[12]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          break;
        case 3212:
          fHistPtYInput[13]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPhiInput[13]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          break;
        default:
          fHistPtYInputRest->Fill(particle->GetPdgCode());
          break;
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskGammaCocktailMC::Terminate(const Option_t *)
{
  
  //fOutputContainer->Print(); // Will crash on GRID
}


//_________________________________________________________________________________
void AliAnalysisTaskGammaCocktailMC::SetLogBinningXTH1(TH1* histoRebin){
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
void AliAnalysisTaskGammaCocktailMC::SetLogBinningXTH2(TH2* histoRebin){
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
Bool_t AliAnalysisTaskGammaCocktailMC::IsMotherInList(TParticle* mother){
  
  Int_t PdgMother = mother->GetPdgCode();
  for(Int_t i=0;i<6;i++){
    if(PdgMother==fParticleList[i]) return kTRUE;
  }
  return kFALSE;
}
