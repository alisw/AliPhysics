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
#include <map>
#include <vector>
#include <algorithm>

ClassImp(AliAnalysisTaskGammaCocktailMC)

//________________________________________________________________________
AliAnalysisTaskGammaCocktailMC::AliAnalysisTaskGammaCocktailMC(): AliAnalysisTaskSE(),
  fOutputContainer(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fMCStack(NULL),
  fDoLightOutput(kFALSE),
  fHistNEvents(NULL),
  fHistPtYGamma(NULL),
  fHistPtPhiGamma(NULL),
  fHistPtPhiGammaSource(NULL),
  fHistPtPhiInput(NULL),
  fHistPtYInput(NULL),
  fHistPtYGammaSource(NULL),
  fHistPtAlphaInput(NULL),
  fHistPtDeltaPhiInput(NULL),
  fHistDecayChannelsInput(NULL),
  fHistPtGammaSourceInput(NULL),
  fHistPhiGammaSourceInput(NULL),
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
  fDoLightOutput(kFALSE),
  fHistNEvents(NULL),
  fHistPtYGamma(NULL),
  fHistPtPhiGamma(NULL),
  fHistPtPhiGammaSource(NULL),
  fHistPtPhiInput(NULL),
  fHistPtYInput(NULL),
  fHistPtYGammaSource(NULL),
  fHistPtAlphaInput(NULL),
  fHistPtDeltaPhiInput(NULL),
  fHistDecayChannelsInput(NULL),
  fHistPtGammaSourceInput(NULL),
  fHistPhiGammaSourceInput(NULL),
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

  fHistNEvents = (TH1F*)SetHist1D(fHistNEvents,"f","NEvents","","N_{evt}",1,0,1,kTRUE);
  fOutputContainer->Add(fHistNEvents);
  
  fHistPtYGamma = (TH2F*)SetHist2D(fHistPtYGamma,"f","Pt_Y_Gamma","#it{p}_{T}","Y",500,0, 50,400,-2.0,2.0,kTRUE);
  fOutputContainer->Add(fHistPtYGamma);
  
  fHistPtPhiGamma = (TH2F*)SetHist2D(fHistPtPhiGamma,"f","Pt_Phi_Gamma","#it{p}_{T}","#phi",500,0,50,100,0,7,kTRUE);
  fOutputContainer->Add(fHistPtPhiGamma);
  
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
  fParticleList            = fParticleList_local;
  fParticleListNames       = fParticleListNames_local;
  fHistPtYInput            = new TH2F*[nInputParticles];
  fHistPtYGammaSource      = new TH2F*[nInputParticles];
  fHistPtAlphaInput        = new TH2F*[nInputParticles];
  fHistPtDeltaPhiInput     = new TH2F*[nInputParticles];
  fHistDecayChannelsInput  = new TH1F*[nInputParticles];
  fHistPtPhiGammaSource    = new TH2F*[nInputParticles];
  fHistPtPhiInput          = new TH2F*[nInputParticles];
  fHistPtGammaSourceInput  = new TH2F*[nInputParticles];
  fHistPhiGammaSourceInput = new TH2F*[nInputParticles];
  
  // delta phi binning
  Double_t *binsDeltaPhi = new Double_t[83];
  binsDeltaPhi[0] = 0;
  Double_t binWidth = 0.02;
  for (Int_t i=1; i<83; i++) {
    if (i<21)
      binWidth = 0.02;
    else if (i>=21 && i<43)
      binWidth = 0.05;
    else if (i>=43 && i<68)
      binWidth = 0.1;
    else
      binWidth = 0.2;
    binsDeltaPhi[i] = binsDeltaPhi[i-1]+binWidth;
  }
    
  for(Int_t i=0; i<nInputParticles; i++){
    
    fHistPtYInput[i] = (TH2F*)SetHist2D(fHistPtYInput[i],"f",Form("Pt_Y_%s",fParticleListNames[i].Data()),"#it{p}_{T}","Y",500,0,50,400,-2.0,2.0,kTRUE);
    fOutputContainer->Add(fHistPtYInput[i]);
    
    //Gammas from certain mother
    fHistPtYGammaSource[i] = (TH2F*)SetHist2D(fHistPtYGammaSource[i],"f",Form("Pt_Y_Gamma_From_%s",fParticleListNames[i].Data()),"#it{p}_{T}","Y",500,0,50,400,-2.0,2.0,kTRUE);
    fOutputContainer->Add(fHistPtYGammaSource[i]);
    
    //phi distributions
    fHistPtPhiInput[i] = (TH2F*)SetHist2D(fHistPtPhiInput[i],"f",Form("Pt_Phi_%s",fParticleListNames[i].Data()),"#it{p}_{T}","#phi",500,0,50,100,0,7,kTRUE);
    fOutputContainer->Add(fHistPtPhiInput[i]);

    fHistPtPhiGammaSource[i] = (TH2F*)SetHist2D(fHistPtPhiGammaSource[i],"f",Form("Pt_Phi_Gamma_From_%s",fParticleListNames[i].Data()),"#it{p}_{T}","#phi",500,0,50,100,0,7,kTRUE);
    fOutputContainer->Add(fHistPtPhiGammaSource[i]);
      
    // correlation gamma from certain mother to mother
    fHistPtGammaSourceInput[i] = (TH2F*)SetHist2D(fHistPtGammaSourceInput[i],"f",Form("PtGamma_PtMother_%s",fParticleListNames[i].Data()),"#it{p}_{T,daughter}","#it{p}_{T,mother}",500,0,50,500,0,50,kTRUE);
    fOutputContainer->Add(fHistPtGammaSourceInput[i]);

    fHistPhiGammaSourceInput[i] = (TH2F*)SetHist2D(fHistPhiGammaSourceInput[i],"f",Form("PhiGamma_PhiMother_%s",fParticleListNames[i].Data()),"#phi_{daughter}","#phi_{mother}",100,0,7,100,0,7,kTRUE);
    fOutputContainer->Add(fHistPhiGammaSourceInput[i]);
    
    // lightweight output
    if (!fDoLightOutput || (fDoLightOutput && (i < 4 || i == 7))) {
      // decay channels mother
      fHistDecayChannelsInput[i] = (TH1F*)SetHist1D(fHistDecayChannelsInput[i],"f",Form("DecayChannels_%s",fParticleListNames[i].Data()),"","", 20,-0.5,19.5,kTRUE);
      fOutputContainer->Add(fHistDecayChannelsInput[i]);
      
      // gamma delta phi
      fHistPtDeltaPhiInput[i] = (TH2F*)SetHist2D(fHistPtDeltaPhiInput[i],"f",Form("Pt_DeltaPhi_%s",fParticleListNames[i].Data()),"#it{p}_{T}","#Delta#phi_{#gamma_{1}#gamma_{2}}",500,0,50,82,binsDeltaPhi,kTRUE);
      fOutputContainer->Add(fHistPtDeltaPhiInput[i]);
      
      // alpha mother
      fHistPtAlphaInput[i] = (TH2F*)SetHist2D(fHistPtAlphaInput[i],"f",Form("Pt_Alpha_%s",fParticleListNames[i].Data()),"#it{p}_{T}","#alpha",500,0,50,100,0,1,kTRUE);
      fOutputContainer->Add(fHistPtAlphaInput[i]);
    } else {
      fHistDecayChannelsInput[i] = NULL;
      fHistPtDeltaPhiInput[i] = NULL;
      fHistPtAlphaInput[i] = NULL;
    }
  }
  InitializeDecayChannelHist();
  delete[] binsDeltaPhi;
  
  fHistPtYInputRest = (TH1I*)SetHist1D(fHistPtYInputRest,"f","Pdg_primary_rest","PDG code","",5000,0,5000,kTRUE);
  fOutputContainer->Add(fHistPtYInputRest);
  
  //Gammas from certain mother
  fHistPtYGammaSourceRest = (TH1I*)SetHist1D(fHistPtYGammaSourceRest,"f","Pdg_Gamma_From_rest","PDG code mother","",5000,0,5000,kTRUE);
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
      
    TParticle* daughter0 = NULL;
    TParticle* daughter1 = NULL;

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
        fHistPtPhiGamma->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
        switch(motherParticle->GetPdgCode()){
        case 111:
          fHistPtYGammaSource[0]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[0]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourceInput[0]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourceInput[0]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 221:
          fHistPtYGammaSource[1]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[1]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourceInput[1]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourceInput[1]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 331:
          fHistPtYGammaSource[2]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[2]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourceInput[2]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourceInput[2]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 223:
          fHistPtYGammaSource[3]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[3]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourceInput[3]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourceInput[3]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 113:
          fHistPtYGammaSource[4]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[4]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourceInput[4]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourceInput[4]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 213:
          fHistPtYGammaSource[5]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[5]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourceInput[5]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourceInput[5]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case -213:
          fHistPtYGammaSource[6]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[6]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourceInput[6]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourceInput[6]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 333:
          fHistPtYGammaSource[7]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[7]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourceInput[7]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourceInput[7]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 443:
          fHistPtYGammaSource[8]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[8]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourceInput[8]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourceInput[8]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 1114:
          fHistPtYGammaSource[9]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[9]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourceInput[9]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourceInput[9]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 2114:
          fHistPtYGammaSource[10]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[10]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourceInput[10]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourceInput[10]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 2214:
          fHistPtYGammaSource[11]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[11]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourceInput[11]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourceInput[11]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 2224:
          fHistPtYGammaSource[12]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[12]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourceInput[12]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourceInput[12]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 3212:
          fHistPtYGammaSource[13]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[13]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourceInput[13]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourceInput[13]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
          default:
            fHistPtYGammaSourceRest->Fill(motherParticle->GetPdgCode());
            break;
        }
      }
    }
    
//     fParticleList_local[] = {111,221,331,223,113,213,-213,333,443,1114,2114,2214,2224,3212};
//     fParticleListNames_local[] = {"Pi0","Eta","EtaPrim","omega","rho0","rho+","rho-","phi","J/psi","Delta-","Delta0","Delta+","Delta++","Sigma0"};
    Double_t alpha    = -1;
    Double_t deltaPhi = -1;
    if(particle->GetPdgCode()!=22 && particleIsPrimary){
        
      if (particle->GetNDaughters() == 2) {
        daughter0 = (TParticle*)fMCStack->Particle(particle->GetFirstDaughter());
        daughter1 = (TParticle*)fMCStack->Particle(particle->GetLastDaughter());
                
        if (daughter0->GetPdgCode()==22 || daughter1->GetPdgCode()==22) {
          alpha    = TMath::Abs((daughter0->Energy() - daughter1->Energy())/(daughter0->Energy() + daughter1->Energy()));
          deltaPhi = TMath::Abs(daughter0->Phi() - daughter1->Phi());
        }
      }
        
      switch(particle->GetPdgCode()){
        case 111:
          fHistPtYInput[0]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[0]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtAlphaInput[0]->Fill(particle->Pt(), alpha, particle->GetWeight());
          fHistPtDeltaPhiInput[0]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          fHistDecayChannelsInput[0]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[0]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          break;
        case 221:
          fHistPtYInput[1]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[1]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtAlphaInput[1]->Fill(particle->Pt(), alpha, particle->GetWeight());
          fHistPtDeltaPhiInput[1]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          fHistDecayChannelsInput[1]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[1]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          break;
        case 331:
          fHistPtYInput[2]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[2]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtAlphaInput[2]->Fill(particle->Pt(), alpha, particle->GetWeight());
          fHistPtDeltaPhiInput[2]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          fHistDecayChannelsInput[2]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[2]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          break;
        case 223:
          fHistPtYInput[3]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[3]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtAlphaInput[3]->Fill(particle->Pt(), alpha, particle->GetWeight());
          fHistPtDeltaPhiInput[3]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          fHistDecayChannelsInput[3]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[3]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          break;
        case 113:
          fHistPtYInput[4]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[4]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          if (fHistPtAlphaInput[4]) fHistPtAlphaInput[4]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (fHistPtDeltaPhiInput[4]) fHistPtDeltaPhiInput[4]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          if (fHistDecayChannelsInput[4]) fHistDecayChannelsInput[4]->Fill(0., particle->GetWeight());
          if (fHistDecayChannelsInput[4]) fHistDecayChannelsInput[4]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          break;
        case 213:
          fHistPtYInput[5]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[5]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          if (fHistPtAlphaInput[5]) fHistPtAlphaInput[5]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (fHistPtDeltaPhiInput[5]) fHistPtDeltaPhiInput[5]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          if (fHistDecayChannelsInput[5]) fHistDecayChannelsInput[5]->Fill(0., particle->GetWeight());
          if (fHistDecayChannelsInput[5]) fHistDecayChannelsInput[5]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          break;
        case -213:
          fHistPtYInput[6]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[6]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          if (fHistPtAlphaInput[6]) fHistPtAlphaInput[6]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (fHistPtDeltaPhiInput[6]) fHistPtDeltaPhiInput[6]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          if (fHistDecayChannelsInput[6]) fHistDecayChannelsInput[6]->Fill(0., particle->GetWeight());
          if (fHistDecayChannelsInput[6]) fHistDecayChannelsInput[6]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          break;
        case 333:
          fHistPtYInput[7]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[7]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtAlphaInput[7]->Fill(particle->Pt(), alpha, particle->GetWeight());
          fHistPtDeltaPhiInput[7]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          fHistDecayChannelsInput[7]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[7]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          break;
        case 443:
          fHistPtYInput[8]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[8]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          if (fHistPtAlphaInput[8]) fHistPtAlphaInput[8]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (fHistPtDeltaPhiInput[8]) fHistPtDeltaPhiInput[8]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          if (fHistDecayChannelsInput[8]) fHistDecayChannelsInput[8]->Fill(0., particle->GetWeight());
          if (fHistDecayChannelsInput[8]) fHistDecayChannelsInput[8]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          break;
        case 1114:
          fHistPtYInput[9]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[9]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          if (fHistPtAlphaInput[9]) fHistPtAlphaInput[9]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (fHistPtDeltaPhiInput[9]) fHistPtDeltaPhiInput[9]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          if (fHistDecayChannelsInput[9]) fHistDecayChannelsInput[9]->Fill(0., particle->GetWeight());
          if (fHistDecayChannelsInput[9]) fHistDecayChannelsInput[9]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          break;
        case 2114:
          fHistPtYInput[10]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[10]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          if (fHistPtAlphaInput[10]) fHistPtAlphaInput[10]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (fHistPtDeltaPhiInput[10]) fHistPtDeltaPhiInput[10]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          if (fHistDecayChannelsInput[10]) fHistDecayChannelsInput[10]->Fill(0., particle->GetWeight());
          if (fHistDecayChannelsInput[10]) fHistDecayChannelsInput[10]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          break;
        case 2214:
          fHistPtYInput[11]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[11]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          if (fHistPtAlphaInput[11]) fHistPtAlphaInput[11]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (fHistPtDeltaPhiInput[11]) fHistPtDeltaPhiInput[11]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          if (fHistDecayChannelsInput[11]) fHistDecayChannelsInput[11]->Fill(0., particle->GetWeight());
          if (fHistDecayChannelsInput[11]) fHistDecayChannelsInput[11]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          break;
        case 2224:
          fHistPtYInput[12]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[12]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          if (fHistPtAlphaInput[12]) fHistPtAlphaInput[12]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (fHistPtDeltaPhiInput[12]) fHistPtDeltaPhiInput[12]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          if (fHistDecayChannelsInput[12]) fHistDecayChannelsInput[12]->Fill(0., particle->GetWeight());
          if (fHistDecayChannelsInput[12]) fHistDecayChannelsInput[12]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          break;
        case 3212:
          fHistPtYInput[13]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[13]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          if (fHistPtAlphaInput[13]) fHistPtAlphaInput[13]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (fHistPtDeltaPhiInput[13]) fHistPtDeltaPhiInput[13]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          if (fHistDecayChannelsInput[13]) fHistDecayChannelsInput[13]->Fill(0., particle->GetWeight());
          if (fHistDecayChannelsInput[13]) fHistDecayChannelsInput[13]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
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

//_________________________________________________________________________________
void AliAnalysisTaskGammaCocktailMC::InitializeDecayChannelHist() {

  fHistDecayChannelsInput[0]->GetXaxis()->SetBinLabel(1,"all");
  fHistDecayChannelsInput[0]->GetXaxis()->SetBinLabel(2,"2#gamma");
  fHistDecayChannelsInput[0]->GetXaxis()->SetBinLabel(3,"e^{+}e^{-}#gamma");
  fHistDecayChannelsInput[0]->GetXaxis()->SetBinLabel(4,"2e^{+}2e^{-}");
  fHistDecayChannelsInput[0]->GetXaxis()->SetBinLabel(5,"e^{+}e^{-}");
  fHistDecayChannelsInput[0]->GetXaxis()->SetBinLabel(20,"rest");

  fHistDecayChannelsInput[1]->GetXaxis()->SetBinLabel(1,"all");
  fHistDecayChannelsInput[1]->GetXaxis()->SetBinLabel(2,"2#gamma");
  fHistDecayChannelsInput[1]->GetXaxis()->SetBinLabel(3,"3#pi^{0}");
  fHistDecayChannelsInput[1]->GetXaxis()->SetBinLabel(4,"#pi^{0}2#gamma");
  fHistDecayChannelsInput[1]->GetXaxis()->SetBinLabel(5,"#pi^{+}#pi^{-}#gamma");
  fHistDecayChannelsInput[1]->GetXaxis()->SetBinLabel(6,"e^{+}e^{-}#gamma");
  fHistDecayChannelsInput[1]->GetXaxis()->SetBinLabel(7,"#mu^{+}#mu^{-}#gamma");
  fHistDecayChannelsInput[1]->GetXaxis()->SetBinLabel(8,"2e^{+}2e^{-}");
  fHistDecayChannelsInput[1]->GetXaxis()->SetBinLabel(9,"#pi^{+}#pi^{-}2#gamma");
  fHistDecayChannelsInput[1]->GetXaxis()->SetBinLabel(20,"rest");

  fHistDecayChannelsInput[2]->GetXaxis()->SetBinLabel(1,"all");
  fHistDecayChannelsInput[2]->GetXaxis()->SetBinLabel(2,"#pi^{+}#pi^{-}#eta");
  fHistDecayChannelsInput[2]->GetXaxis()->SetBinLabel(3,"#rho^{0}#gamma");
  fHistDecayChannelsInput[2]->GetXaxis()->SetBinLabel(4,"#pi^{+}#pi^{-}#gamma");
  fHistDecayChannelsInput[2]->GetXaxis()->SetBinLabel(5,"#omega#gamma");
  fHistDecayChannelsInput[2]->GetXaxis()->SetBinLabel(6,"2#gamma");
  fHistDecayChannelsInput[2]->GetXaxis()->SetBinLabel(7,"#mu^{+}#mu^{-}#gamma");
  fHistDecayChannelsInput[2]->GetXaxis()->SetBinLabel(20,"rest");

  fHistDecayChannelsInput[3]->GetXaxis()->SetBinLabel(1,"all");
  fHistDecayChannelsInput[3]->GetXaxis()->SetBinLabel(2,"#pi^{+}#pi^{-}#pi^{0}");
  fHistDecayChannelsInput[3]->GetXaxis()->SetBinLabel(3,"#pi^{0}#gamma");
  fHistDecayChannelsInput[3]->GetXaxis()->SetBinLabel(4,"#eta#gamma");
  fHistDecayChannelsInput[3]->GetXaxis()->SetBinLabel(5,"#pi^{0}e^{+}e^{-}");
  fHistDecayChannelsInput[3]->GetXaxis()->SetBinLabel(6,"e^{+}e^{-}");
  fHistDecayChannelsInput[3]->GetXaxis()->SetBinLabel(7,"2#pi^{0}#gamma");
  fHistDecayChannelsInput[3]->GetXaxis()->SetBinLabel(20,"rest");
  
  fHistDecayChannelsInput[7]->GetXaxis()->SetBinLabel(1,"all");
  fHistDecayChannelsInput[7]->GetXaxis()->SetBinLabel(2,"K^{+}K^{-}");
  fHistDecayChannelsInput[7]->GetXaxis()->SetBinLabel(3,"K^{0}_{L}K^{0}_{S}");
  fHistDecayChannelsInput[7]->GetXaxis()->SetBinLabel(4,"#eta#gamma");
  fHistDecayChannelsInput[7]->GetXaxis()->SetBinLabel(5,"#pi^{0}#gamma");
  fHistDecayChannelsInput[7]->GetXaxis()->SetBinLabel(6,"e^{+}e^{-}");
  fHistDecayChannelsInput[7]->GetXaxis()->SetBinLabel(7,"#eta e^{+}e^{-}");
  fHistDecayChannelsInput[7]->GetXaxis()->SetBinLabel(8,"#pi^{+}#pi^{-}#gamma");
  fHistDecayChannelsInput[7]->GetXaxis()->SetBinLabel(9,"f_{0}(980)#gamma");
  fHistDecayChannelsInput[7]->GetXaxis()->SetBinLabel(10,"2#pi^{0}#gamma");
  fHistDecayChannelsInput[7]->GetXaxis()->SetBinLabel(11,"#pi^{0}e^{+}e^{-}");
  fHistDecayChannelsInput[7]->GetXaxis()->SetBinLabel(12,"#pi^{0}#eta#gamma");
  fHistDecayChannelsInput[7]->GetXaxis()->SetBinLabel(13,"a_{0}(980)#gamma");
  fHistDecayChannelsInput[7]->GetXaxis()->SetBinLabel(14,"#eta'#gamma");
  fHistDecayChannelsInput[7]->GetXaxis()->SetBinLabel(15,"#mu^{+}#mu^{-}#gamma");
  fHistDecayChannelsInput[7]->GetXaxis()->SetBinLabel(20,"rest");

  if (!fDoLightOutput) {
    fHistDecayChannelsInput[4]->GetXaxis()->SetBinLabel(1,"all");
    fHistDecayChannelsInput[4]->GetXaxis()->SetBinLabel(2,"2#pi^{0}");
    fHistDecayChannelsInput[4]->GetXaxis()->SetBinLabel(3,"#pi^{+}#pi^{-}#gamma");
    fHistDecayChannelsInput[4]->GetXaxis()->SetBinLabel(4,"#pi^{0}#gamma");
    fHistDecayChannelsInput[4]->GetXaxis()->SetBinLabel(5,"#eta#gamma");
    fHistDecayChannelsInput[4]->GetXaxis()->SetBinLabel(6,"2#pi^{0}#gamma");
    fHistDecayChannelsInput[4]->GetXaxis()->SetBinLabel(7,"e^{+}e^{-}");
    fHistDecayChannelsInput[4]->GetXaxis()->SetBinLabel(20,"rest");

    fHistDecayChannelsInput[5]->GetXaxis()->SetBinLabel(1,"all");
    fHistDecayChannelsInput[5]->GetXaxis()->SetBinLabel(2,"#pi^{+}#pi^{0}");
    fHistDecayChannelsInput[5]->GetXaxis()->SetBinLabel(3,"#pi^{+}#gamma");
    fHistDecayChannelsInput[5]->GetXaxis()->SetBinLabel(20,"rest");

    fHistDecayChannelsInput[6]->GetXaxis()->SetBinLabel(1,"all");
    fHistDecayChannelsInput[6]->GetXaxis()->SetBinLabel(2,"#pi^{-}#pi^{0}");
    fHistDecayChannelsInput[6]->GetXaxis()->SetBinLabel(3,"#pi^{-}#gamma");
    fHistDecayChannelsInput[6]->GetXaxis()->SetBinLabel(20,"rest");

    fHistDecayChannelsInput[8]->GetXaxis()->SetBinLabel(1,"all");
    fHistDecayChannelsInput[8]->GetXaxis()->SetBinLabel(2,"3g");
    fHistDecayChannelsInput[8]->GetXaxis()->SetBinLabel(3,"2g#gamma");
    fHistDecayChannelsInput[8]->GetXaxis()->SetBinLabel(4,"e^{+}e^{-}");
    fHistDecayChannelsInput[8]->GetXaxis()->SetBinLabel(5,"e^{+}e^{-}#gamma");
    fHistDecayChannelsInput[8]->GetXaxis()->SetBinLabel(20,"rest");

    fHistDecayChannelsInput[9]->GetXaxis()->SetBinLabel(1,"all");
    fHistDecayChannelsInput[9]->GetXaxis()->SetBinLabel(2,"n#pi^{-}");
    fHistDecayChannelsInput[9]->GetXaxis()->SetBinLabel(20,"rest");

    fHistDecayChannelsInput[10]->GetXaxis()->SetBinLabel(1,"all");
    fHistDecayChannelsInput[10]->GetXaxis()->SetBinLabel(2,"n#pi^{0}");
    fHistDecayChannelsInput[10]->GetXaxis()->SetBinLabel(3,"p#pi^{-}");
    fHistDecayChannelsInput[10]->GetXaxis()->SetBinLabel(4,"n#gamma");
    fHistDecayChannelsInput[10]->GetXaxis()->SetBinLabel(20,"rest");

    fHistDecayChannelsInput[11]->GetXaxis()->SetBinLabel(1,"all");
    fHistDecayChannelsInput[11]->GetXaxis()->SetBinLabel(2,"n#pi^{+}");
    fHistDecayChannelsInput[11]->GetXaxis()->SetBinLabel(3,"p#pi^{0}");
    fHistDecayChannelsInput[11]->GetXaxis()->SetBinLabel(4,"p#gamma");
    fHistDecayChannelsInput[11]->GetXaxis()->SetBinLabel(20,"rest");

    fHistDecayChannelsInput[12]->GetXaxis()->SetBinLabel(1,"all");
    fHistDecayChannelsInput[12]->GetXaxis()->SetBinLabel(2,"p#pi^{+}");
    fHistDecayChannelsInput[12]->GetXaxis()->SetBinLabel(20,"rest");

    fHistDecayChannelsInput[13]->GetXaxis()->SetBinLabel(1,"all");
    fHistDecayChannelsInput[13]->GetXaxis()->SetBinLabel(2,"#Lambda#gamma");
    fHistDecayChannelsInput[13]->GetXaxis()->SetBinLabel(3,"#Lambda e^{+}e^{-}");
    fHistDecayChannelsInput[13]->GetXaxis()->SetBinLabel(20,"rest");
  }
}

//_________________________________________________________________________________
Float_t AliAnalysisTaskGammaCocktailMC::GetDecayChannel(AliStack* stack, TParticle* part) {
    
  Int_t nDaughters = part->GetNDaughters();
  if (nDaughters > 10) return 19.;
  
  std::vector<Long64_t> *PdgDaughter = new std::vector<Long64_t>(nDaughters);
  Long64_t tempPdgCode = 0;
  for (Int_t i=0; i<nDaughters; i++) {
    tempPdgCode = (Long64_t)((TParticle*)stack->Particle(part->GetFirstDaughter()+i))->GetPdgCode();
    if (TMath::Abs(tempPdgCode) == 111 || TMath::Abs(tempPdgCode) == 113 || TMath::Abs(tempPdgCode) == 130 || TMath::Abs(tempPdgCode) == 310 || TMath::Abs(tempPdgCode) == 223 || TMath::Abs(tempPdgCode) == 221 || TMath::Abs(tempPdgCode) == 331 || TMath::Abs(tempPdgCode) == 2112 || TMath::Abs(tempPdgCode) == 3122 || TMath::Abs(tempPdgCode) == 9000111 || TMath::Abs(tempPdgCode) == 9010221)
      tempPdgCode = TMath::Abs(tempPdgCode);
    PdgDaughter->at(i) = tempPdgCode;
  }
  std::sort(PdgDaughter->begin(), PdgDaughter->end());
  
  switch (part->GetPdgCode()) {
    case 111:
      if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 22)
        return 1.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11 && PdgDaughter->at(2) == 22)
        return 2.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == -11 && PdgDaughter->at(2) == 11 && PdgDaughter->at(3) == 11)
        return 3.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11)
        return 4.;
      else
        return 19.;
      break;
    case 221:
      if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 22)
        return 1.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 111)
        return 2.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 22 && PdgDaughter->at(2) == 111)
        return 3.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 22 && PdgDaughter->at(2) == 211)
        return 4.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11 && PdgDaughter->at(2) == 22)
        return 5.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -13 && PdgDaughter->at(1) == 13 && PdgDaughter->at(2) == 22)
        return 6.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == -11 && PdgDaughter->at(2) == 11 && PdgDaughter->at(3) == 11)
        return 7.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 22 && PdgDaughter->at(2) == 22 && PdgDaughter->at(3) == 211)
        return 8.;
      else
        return 19.;
      break;
    case 331:
      if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 211 && PdgDaughter->at(2) == 221)
        return 1.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 113)
        return 2.;
      if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 22 && PdgDaughter->at(2) == 211)
        return 3.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 223)
        return 4.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 22)
        return 5.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -13 && PdgDaughter->at(1) == 13 && PdgDaughter->at(2) == 22)
        return 6.;
      else
        return 19.;
      break;
    case 223:
      if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 211)
        return 1.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 111)
        return 2.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 221)
        return 3.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11 && PdgDaughter->at(2) == 111)
        return 4.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11)
        return 5.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 111)
        return 6.;
      else
        return 19.;
      break;
    case 113:
      if (nDaughters == 2 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 111)
        return 1.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 22 && PdgDaughter->at(2) == 211)
        return 2.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 111)
        return 3.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 221)
        return 4.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 111)
        return 5.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11)
        return 4.;
      else
        return 19.;
      break;
    case 213:
      if (nDaughters == 2 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 211)
        return 1.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 211)
        return 2.;
      else
        return 19.;
      break;
    case -213:
      if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 111)
        return 1.;
      else if (nDaughters == 2 && PdgDaughter->at(1) == -211 && PdgDaughter->at(0) == 22)
        return 2.;
      else
        return 19.;
      break;
    case 333:
      if (nDaughters == 2 && PdgDaughter->at(0) == -321 && PdgDaughter->at(1) == 321)
        return 1.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 130 && PdgDaughter->at(1) == 310)
        return 2.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 221)
        return 3.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 111)
        return 4.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11)
        return 5.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11 && PdgDaughter->at(2) == 221)
        return 6.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 22 && PdgDaughter->at(2) == 211)
        return 7.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 9010221)
        return 8.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 111)
        return 9.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11 && PdgDaughter->at(2) == 111)
        return 10.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 221)
        return 11.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 9000111)
        return 12.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 331)
        return 13.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -13 && PdgDaughter->at(1) == 13 && PdgDaughter->at(2) == 22)
        return 14.;
      else
        return 19.;
      break;
    case 443:
      if (nDaughters == 3 && (PdgDaughter->at(0) == 21 || PdgDaughter->at(0) == 9) && (PdgDaughter->at(1) == 21 || PdgDaughter->at(1) == 9) && (PdgDaughter->at(2) == 21 || PdgDaughter->at(2) == 9))
        return 1.;
      else if (nDaughters == 3 && (PdgDaughter->at(0) == 21 || PdgDaughter->at(0) == 9) && (PdgDaughter->at(1) == 21 || PdgDaughter->at(1) == 9) && PdgDaughter->at(2) == 22)
        return 2.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11)
        return 3.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11 && PdgDaughter->at(2) == 22)
        return 4.;
      else
        return 19.;
      break;
    case 1114:
      if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 2112)
        return 1.;
      else
        return 19.;
      break;
    case 2114:
      if (nDaughters == 2 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 2112)
        return 1.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 2212)
        return 2.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 2112)
        return 3.;
      else
        return 19.;
      break;
    case 2214:
      if (nDaughters == 2 && PdgDaughter->at(0) == 211 && PdgDaughter->at(1) == 2112)
        return 1.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 2212)
        return 2.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 2212)
        return 3.;
      else
        return 19.;
      break;
    case 2224:
      if (nDaughters == 2 && PdgDaughter->at(0) == 211 && PdgDaughter->at(1) == 2212)
        return 1.;
      else
        return 19.;
      break;
    case 3212:
      if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 3122)
        return 1.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11 && PdgDaughter->at(2) == 3122)
        return 2.;
      else
        return 19.;
      break;
    default:
      return -1.;
      break;
  }

  delete PdgDaughter;
}

//_________________________________________________________________________________
TH1* AliAnalysisTaskGammaCocktailMC::SetHist1D(TH1* hist, TString histType, TString histName, TString xTitle, TString yTitle, Int_t nBinsX, Double_t xMin, Double_t xMax, Bool_t optSumw2) {

  if (histType.CompareTo("f") == 0 || histType.CompareTo("F") == 0)
    hist = new TH1F(histName, histName, nBinsX, xMin, xMax);
  if (histType.CompareTo("i") == 0 || histType.CompareTo("I") == 0)
    hist = new TH1I(histName, histName, nBinsX, xMin, xMax);
  
  hist->GetXaxis()->SetTitle(xTitle);
  hist->GetYaxis()->SetTitle(yTitle);
  
  if (optSumw2)
   hist->Sumw2();
  
  return hist;
}

//_________________________________________________________________________________
TH2* AliAnalysisTaskGammaCocktailMC::SetHist2D(TH2* hist, TString histType, TString histName, TString xTitle, TString yTitle, Int_t nBinsX, Double_t xMin, Double_t xMax, Int_t nBinsY, Double_t yMin, Double_t yMax, Bool_t optSumw2) {

  if (histType.CompareTo("f") == 0 || histType.CompareTo("F") == 0)
    hist = new TH2F(histName, histName, nBinsX, xMin, xMax, nBinsY, yMin, yMax);
  if (histType.CompareTo("i") == 0 || histType.CompareTo("I") == 0)
    hist = new TH2I(histName, histName, nBinsX, xMin, xMax, nBinsY, yMin, yMax);
  
  hist->GetXaxis()->SetTitle(xTitle);
  hist->GetYaxis()->SetTitle(yTitle);
  
  if (optSumw2)
    hist->Sumw2();
  
  return hist;
}

//_________________________________________________________________________________
TH2* AliAnalysisTaskGammaCocktailMC::SetHist2D(TH2* hist, TString histType, TString histName, TString xTitle, TString yTitle, Int_t nBinsX, Double_t xMin, Double_t xMax, Int_t nBinsY, Double_t* binsY, Bool_t optSumw2) {
  
  if (histType.CompareTo("f") == 0 || histType.CompareTo("F") == 0)
    hist = new TH2F(histName, histName, nBinsX, xMin, xMax, nBinsY, binsY);
  if (histType.CompareTo("i") == 0 || histType.CompareTo("I") == 0)
    hist = new TH2I(histName, histName, nBinsX, xMin, xMax, nBinsY, binsY);
  
  hist->GetXaxis()->SetTitle(xTitle);
  hist->GetYaxis()->SetTitle(yTitle);
  
  if (optSumw2)
    hist->Sumw2();
  
  return hist;
}

