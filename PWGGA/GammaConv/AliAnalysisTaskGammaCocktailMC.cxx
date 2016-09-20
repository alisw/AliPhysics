/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.  *
*                                                                         *
* Author: Friederike Bock, Mike Sas, Lucas Altenkämper                    *
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
#include "TF1.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TObject.h"
#include "TObjArray.h"
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
#include "AliMCGenHandler.h"
#include "AliGenEMCocktailV2.h"
#include "AliGenerator.h"
#include "AliPythia6.h"
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
  fHasPtParametrization(kFALSE),
  fHasMother{kFALSE},
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
  fHistPythiaBR(NULL),
  fHistPtGammaSourceInput(NULL),
  fHistPhiGammaSourceInput(NULL),
  fHistPdgInputRest(NULL),
  fHistPdgGammaSourceRest(NULL),
  fParticleList(NULL),
  fParticleListNames(NULL),
  fPtParametrization{NULL},
  fCocktailSettings{NULL},
  fMtScalingFactors(NULL),
  fUserInfo(NULL),
  fOutputTree(NULL),
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
  fHasPtParametrization(kFALSE),
  fHasMother{kFALSE},
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
  fHistPythiaBR(NULL),
  fHistPtGammaSourceInput(NULL),
  fHistPhiGammaSourceInput(NULL),
  fHistPdgInputRest(NULL),
  fHistPdgGammaSourceRest(NULL),
  fParticleList(NULL),
  fParticleListNames(NULL),
  fPtParametrization{NULL},
  fCocktailSettings{NULL},
  fMtScalingFactors(NULL),
  fUserInfo(NULL),
  fOutputTree(NULL),
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
  
  // tree + user info list to protect contents from merging
  fOutputTree = new TTree("cocktailSettings", "cocktailSettings");
  fUserInfo   = (TList*)fOutputTree->GetUserInfo();
  
  AliMCGenHandler* mcGenHandler           = (AliMCGenHandler*)AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler();
  const AliGenerator* mcGenerator         = mcGenHandler->GetGenerator();
  TString mcGeneratorClassName;
  if (mcGenerator)  mcGeneratorClassName  = mcGenerator->ClassName();
  else              mcGeneratorClassName  = "";
  
  if (mcGenerator && mcGeneratorClassName.CompareTo("AliGenEMCocktailV2") == 0) {
    
    AliGenEMCocktailV2* mcCocktailGen = (AliGenEMCocktailV2*)mcGenerator;
    
    // has mother i
    SetHasMother((UInt_t)mcCocktailGen->GetSelectedMothers());
    
    // pt parametrizations
    GetAndSetPtParametrizations(mcCocktailGen, fHasPtParametrization);
    for (Int_t i=0; i<14; i++) {
      if (fHasMother[i]) fUserInfo->Add(fPtParametrization[i]);
    }
    
    // cocktail settings
    Double_t ptMin, ptMax;
    mcCocktailGen->GetPtRange(ptMin, ptMax);
    fCocktailSettings[0] = new TObjString(Form("collSys_%d",  mcCocktailGen->GetCollisionSystem()));
    fCocktailSettings[1] = new TObjString(Form("cent_%d",     mcCocktailGen->GetCentrality()));
    fCocktailSettings[2] = new TObjString(Form("decayMode_%.0f", mcCocktailGen->GetDecayMode()));
    fCocktailSettings[3] = new TObjString(Form("selectMothers_%d", mcCocktailGen->GetSelectedMothers()));
    fCocktailSettings[4] = new TObjString(Form("paramFile_%s", (mcCocktailGen->GetParametrizationFile()).Data()));
    fCocktailSettings[5] = new TObjString(Form("nParticles_%d", mcCocktailGen->GetNumberOfParticles()));
    fCocktailSettings[6] = new TObjString(Form("ptMin_%.2f", ptMin));
    fCocktailSettings[7] = new TObjString(Form("ptMax_%.2f", ptMax));
    fCocktailSettings[8] = new TObjString(Form("weightMode_%.0f", mcCocktailGen->GetWeightingMode()));
    for (Int_t i=0; i<9; i++) fUserInfo->Add(fCocktailSettings[i]);

    // mt scaling params
    TH1D* mtFactorHisto = (TH1D*)mcCocktailGen->GetMtScalingFactors();
    fMtScalingFactors   = new TH1D(*mtFactorHisto);
    fUserInfo->Add(fMtScalingFactors);

  } else {
    for (Int_t i=0; i<14; i++)
      fHasMother[i] = kTRUE;
  }
  
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
  const Int_t nInputParticles         = 14;
  Int_t   fParticleList_local[]       = {111,221,331,223,113,213,-213,333,443,1114,2114,2214,2224,3212};
  TString fParticleListNames_local[]  = {"Pi0","Eta","EtaPrim","omega","rho0","rho+","rho-","phi","J/psi","Delta-","Delta0","Delta+","Delta++","Sigma0"};
  fParticleList            = fParticleList_local;
  fParticleListNames       = fParticleListNames_local;
  fHistPtYInput            = new TH2F*[nInputParticles];
  fHistPtYGammaSource      = new TH2F*[nInputParticles];
  fHistPtAlphaInput        = new TH2F*[nInputParticles];
  fHistPtDeltaPhiInput     = new TH2F*[nInputParticles];
  fHistDecayChannelsInput  = new TH1F*[nInputParticles];
  fHistPythiaBR            = new TH1F*[nInputParticles];
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
    if (fHasMother[i]) fOutputContainer->Add(fHistPtYInput[i]);
    
    //Gammas from certain mother
    fHistPtYGammaSource[i] = (TH2F*)SetHist2D(fHistPtYGammaSource[i],"f",Form("Pt_Y_Gamma_From_%s",fParticleListNames[i].Data()),"#it{p}_{T}","Y",500,0,50,400,-2.0,2.0,kTRUE);
    if (fHasMother[i]) fOutputContainer->Add(fHistPtYGammaSource[i]);
    
    //phi distributions
    fHistPtPhiInput[i] = (TH2F*)SetHist2D(fHistPtPhiInput[i],"f",Form("Pt_Phi_%s",fParticleListNames[i].Data()),"#it{p}_{T}","#phi",500,0,50,100,0,7,kTRUE);
    if (fHasMother[i]) fOutputContainer->Add(fHistPtPhiInput[i]);

    fHistPtPhiGammaSource[i] = (TH2F*)SetHist2D(fHistPtPhiGammaSource[i],"f",Form("Pt_Phi_Gamma_From_%s",fParticleListNames[i].Data()),"#it{p}_{T}","#phi",500,0,50,100,0,7,kTRUE);
    if (fHasMother[i]) fOutputContainer->Add(fHistPtPhiGammaSource[i]);
      
    // correlation gamma from certain mother to mother
    fHistPtGammaSourceInput[i] = (TH2F*)SetHist2D(fHistPtGammaSourceInput[i],"f",Form("PtGamma_PtMother_%s",fParticleListNames[i].Data()),"#it{p}_{T,daughter}","#it{p}_{T,mother}",500,0,50,500,0,50,kTRUE);
    if (fHasMother[i]) fOutputContainer->Add(fHistPtGammaSourceInput[i]);

    fHistPhiGammaSourceInput[i] = (TH2F*)SetHist2D(fHistPhiGammaSourceInput[i],"f",Form("PhiGamma_PhiMother_%s",fParticleListNames[i].Data()),"#phi_{daughter}","#phi_{mother}",100,0,7,100,0,7,kTRUE);
    if (fHasMother[i]) fOutputContainer->Add(fHistPhiGammaSourceInput[i]);
    
    // decay channels mother
    fHistDecayChannelsInput[i] = (TH1F*)SetHist1D(fHistDecayChannelsInput[i],"f",Form("DecayChannels_%s",fParticleListNames[i].Data()),"","", 20,-0.5,19.5,kTRUE);
    InitializeDecayChannelHist(fHistDecayChannelsInput[i], i);
    if (fHasMother[i]) fOutputContainer->Add(fHistDecayChannelsInput[i]);
    
    // BR from pythia
    fHistPythiaBR[i] = (TH1F*)SetHist1D(fHistPythiaBR[i],"f",Form("PythiaBR_%s",fParticleListNames[i].Data()),"","", 20,-0.5,19.5,kTRUE);
    InitializeDecayChannelHist(fHistPythiaBR[i], i);
    FillPythiaBranchingRatio(fHistPythiaBR[i], i);
    if (fHasMother[i]) fUserInfo->Add(fHistPythiaBR[i]);
    
    // lightweight output
    if (!fDoLightOutput || (fDoLightOutput && (i < 4 || i == 7))) {
      // gamma delta phi
      fHistPtDeltaPhiInput[i] = (TH2F*)SetHist2D(fHistPtDeltaPhiInput[i],"f",Form("Pt_DeltaPhi_%s",fParticleListNames[i].Data()),"#it{p}_{T}","#Delta#phi_{#gamma_{1}#gamma_{2}}",500,0,50,82,binsDeltaPhi,kTRUE);
      if (fHasMother[i]) fOutputContainer->Add(fHistPtDeltaPhiInput[i]);
      
      // alpha mother
      fHistPtAlphaInput[i] = (TH2F*)SetHist2D(fHistPtAlphaInput[i],"f",Form("Pt_Alpha_%s",fParticleListNames[i].Data()),"#it{p}_{T}","#alpha",500,0,50,100,-1,1,kTRUE);
      if (fHasMother[i]) fOutputContainer->Add(fHistPtAlphaInput[i]);
    } else {
      fHistPtDeltaPhiInput[i] = NULL;
      fHistPtAlphaInput[i] = NULL;
    }
  }
  delete[] binsDeltaPhi;
  
  fHistPdgInputRest = (TH1I*)SetHist1D(fHistPdgInputRest,"f","Pdg_primary_rest","PDG code","",5000,0,5000,kTRUE);
  fOutputContainer->Add(fHistPdgInputRest);
  
  //Gammas from certain mother
  fHistPdgGammaSourceRest = (TH1I*)SetHist1D(fHistPdgGammaSourceRest,"f","Pdg_Gamma_From_rest","PDG code mother","",5000,0,5000,kTRUE);
  fOutputContainer->Add(fHistPdgGammaSourceRest);
  
  fOutputContainer->Add(fOutputTree);

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

//_____________________________________________________________________________
void AliAnalysisTaskGammaCocktailMC::GetAndSetPtParametrizations(AliGenEMCocktailV2* mcCocktailGen, Bool_t setParams)
{
  if (setParams || !mcCocktailGen) return;

  for (Int_t i=0; i<14; i++) fPtParametrization[i] = NULL;
  
  TF1* fct        = NULL;
  TString fctName = "";
  for (Int_t i=0; i<16; i++) {
    fct = (TF1*)mcCocktailGen->GetPtParametrization(i);
    if (fct) {
      fctName = fct->GetName();
      if (fctName.BeginsWith("111_pt"))  fPtParametrization[0]   = new TF1(*fct);
      if (fctName.BeginsWith("221_pt"))  fPtParametrization[1]   = new TF1(*fct);
      if (fctName.BeginsWith("331_pt"))  fPtParametrization[2]   = new TF1(*fct);
      if (fctName.BeginsWith("223_pt"))  fPtParametrization[3]   = new TF1(*fct);
      if (fctName.BeginsWith("113_pt"))  fPtParametrization[4]   = new TF1(*fct);
      if (fctName.BeginsWith("213_pt"))  fPtParametrization[5]   = new TF1(*fct);
      if (fctName.BeginsWith("-213_pt")) fPtParametrization[6]   = new TF1(*fct);
      if (fctName.BeginsWith("333_pt"))  fPtParametrization[7]   = new TF1(*fct);
      if (fctName.BeginsWith("443_pt"))  fPtParametrization[8]   = new TF1(*fct);
      if (fctName.BeginsWith("1114_pt")) fPtParametrization[9]   = new TF1(*fct);
      if (fctName.BeginsWith("2114_pt")) fPtParametrization[10]  = new TF1(*fct);
      if (fctName.BeginsWith("2214_pt")) fPtParametrization[11]  = new TF1(*fct);
      if (fctName.BeginsWith("2224_pt")) fPtParametrization[12]  = new TF1(*fct);
      if (fctName.BeginsWith("3212_pt")) fPtParametrization[13]  = new TF1(*fct);
    }
  }

  fHasPtParametrization = kTRUE;
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaCocktailMC::SetHasMother(UInt_t selectedMothers) {
  
  for (Int_t i=0; i<14; i++) fHasMother[i] = kFALSE;
  
  if (selectedMothers&AliGenEMCocktailV2::kGenPizero)    fHasMother[0] = kTRUE;
  if (selectedMothers&AliGenEMCocktailV2::kGenEta)       fHasMother[1] = kTRUE;
  if (selectedMothers&AliGenEMCocktailV2::kGenEtaprime)  fHasMother[2] = kTRUE;
  if (selectedMothers&AliGenEMCocktailV2::kGenOmega)     fHasMother[3] = kTRUE;
  if (selectedMothers&AliGenEMCocktailV2::kGenRho0)      fHasMother[4] = kTRUE;
  if (selectedMothers&AliGenEMCocktailV2::kGenRhoPl)     fHasMother[5] = kTRUE;
  if (selectedMothers&AliGenEMCocktailV2::kGenRhoMi)     fHasMother[6] = kTRUE;
  if (selectedMothers&AliGenEMCocktailV2::kGenPhi)       fHasMother[7] = kTRUE;
  if (selectedMothers&AliGenEMCocktailV2::kGenJpsi)      fHasMother[8] = kTRUE;
  if (selectedMothers&AliGenEMCocktailV2::kGenDeltaMi)   fHasMother[9] = kTRUE;
  if (selectedMothers&AliGenEMCocktailV2::kGenDeltaZero) fHasMother[10] = kTRUE;
  if (selectedMothers&AliGenEMCocktailV2::kGenDeltaPl)   fHasMother[11] = kTRUE;
  if (selectedMothers&AliGenEMCocktailV2::kGenDeltaPlPl) fHasMother[12] = kTRUE;
  if (selectedMothers&AliGenEMCocktailV2::kGenSigma0)    fHasMother[13] = kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskGammaCocktailMC::ProcessMCParticles(){

  // Loop over all primary MC particle
  for(Long_t i = 0; i < fMCStack->GetNtrack(); i++) {
    // fill primary histograms
    TParticle* particle         = NULL;
    particle                    = (TParticle *)fMCStack->Particle(i);
    if (!particle) continue;
    Bool_t hasMother            = kFALSE;
    Bool_t particleIsPrimary    = kTRUE;

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
    
    Bool_t motherIsPrimary      = kFALSE;
    if(hasMother){
      if(motherParticle->GetMother(0)>-1)motherIsPrimary = kFALSE;
      else motherIsPrimary                               = kTRUE;
    }
      
    TParticle* daughter0 = NULL;
    TParticle* daughter1 = NULL;

    if (!(fabs(particle->Energy()-particle->Pz())>0.)) continue;
    Double_t yPre = (particle->Energy()+particle->Pz())/(particle->Energy()-particle->Pz());
//     cout << i << "\t"<< particle->GetPdgCode() << "\t"<< particle->Pz() << "\t" << particle->Energy()<< "\t" << particle->Energy()-particle->Pz() << "\t"<< yPre << endl;
    if (yPre == 0.) continue;
    
    Double_t y = 0.5*TMath::Log(yPre); 
    if (fabs(y) > fMaxY) continue;
    
    if(particle->GetPdgCode()==22 && hasMother==kTRUE){
      if(motherIsPrimary){
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
          fHistPdgGammaSourceRest->Fill(motherParticle->GetPdgCode());
          break;
        }
      }
    }
    
//     fParticleList_local[] = {111,221,331,223,113,213,-213,333,443,1114,2114,2214,2224,3212};
//     fParticleListNames_local[] = {"Pi0","Eta","EtaPrim","omega","rho0","rho+","rho-","phi","J/psi","Delta-","Delta0","Delta+","Delta++","Sigma0"};
    if(particle->GetPdgCode()!=22 && particleIsPrimary){

      Double_t alpha    = -999;
      Double_t deltaPhi = -999;
      if (particle->GetNDaughters() == 2) {
        daughter0 = (TParticle*)fMCStack->Particle(particle->GetFirstDaughter());
        daughter1 = (TParticle*)fMCStack->Particle(particle->GetLastDaughter());
                
        if (daughter0->GetPdgCode()==22 || daughter1->GetPdgCode()==22) {
          Double_t firstEnergy, secondEnergy;
          if (daughter0->GetPdgCode()==22) {
            firstEnergy = daughter0->Energy();
            secondEnergy = daughter1->Energy();
          } else {
            firstEnergy = daughter1->Energy();
            secondEnergy = daughter0->Energy();
          }
          alpha = (firstEnergy - secondEnergy)/(firstEnergy + secondEnergy);
        }
        
        if (daughter0->GetPdgCode()==22 && daughter1->GetPdgCode()==22) {
          deltaPhi = TMath::Abs(daughter0->Phi() - daughter1->Phi());
        }
      }
        
      switch(particle->GetPdgCode()){
        case 111:
          fHistPtYInput[0]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[0]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[0]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[0]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          if (alpha>=-1) fHistPtAlphaInput[0]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0) fHistPtDeltaPhiInput[0]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 221:
          fHistPtYInput[1]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[1]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[1]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[1]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          if (alpha>=-1) fHistPtAlphaInput[1]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0) fHistPtDeltaPhiInput[1]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 331:
          fHistPtYInput[2]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[2]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[2]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[2]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          if (alpha>=-1) fHistPtAlphaInput[2]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0) fHistPtDeltaPhiInput[2]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 223:
          fHistPtYInput[3]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[3]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[3]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[3]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          if (alpha>=-1) fHistPtAlphaInput[3]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0) fHistPtDeltaPhiInput[3]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 113:
          fHistPtYInput[4]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[4]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[4]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[4]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          if (alpha>=-1 && fHistPtAlphaInput[4]) fHistPtAlphaInput[4]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0 && fHistPtDeltaPhiInput[4]) fHistPtDeltaPhiInput[4]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 213:
          fHistPtYInput[5]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[5]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[5]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[5]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          if (alpha>=-1 && fHistPtAlphaInput[5]) fHistPtAlphaInput[5]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0 && fHistPtDeltaPhiInput[5]) fHistPtDeltaPhiInput[5]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case -213:
          fHistPtYInput[6]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[6]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[6]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[6]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          if (alpha>=-1 && fHistPtAlphaInput[6]) fHistPtAlphaInput[6]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0 && fHistPtDeltaPhiInput[6]) fHistPtDeltaPhiInput[6]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 333:
          fHistPtYInput[7]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[7]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[7]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[7]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          if (alpha >= -1) fHistPtAlphaInput[7]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0) fHistPtDeltaPhiInput[7]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 443:
          fHistPtYInput[8]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[8]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[8]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[8]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          if (alpha>=-1 && fHistPtAlphaInput[8]) fHistPtAlphaInput[8]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0 && fHistPtDeltaPhiInput[8]) fHistPtDeltaPhiInput[8]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 1114:
          fHistPtYInput[9]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[9]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[9]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[9]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          if (alpha>=-1 && fHistPtAlphaInput[9]) fHistPtAlphaInput[9]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0 && fHistPtDeltaPhiInput[9]) fHistPtDeltaPhiInput[9]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 2114:
          fHistPtYInput[10]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[10]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[10]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[10]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          if (alpha>=-1 && fHistPtAlphaInput[10]) fHistPtAlphaInput[10]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0 && fHistPtDeltaPhiInput[10]) fHistPtDeltaPhiInput[10]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 2214:
          fHistPtYInput[11]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[11]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[11]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[11]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          if (alpha>=-1 && fHistPtAlphaInput[11]) fHistPtAlphaInput[11]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0 && fHistPtDeltaPhiInput[11]) fHistPtDeltaPhiInput[11]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 2224:
          fHistPtYInput[12]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[12]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[12]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[12]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          if (alpha>=-1 && fHistPtAlphaInput[12]) fHistPtAlphaInput[12]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0 && fHistPtDeltaPhiInput[12]) fHistPtDeltaPhiInput[12]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 3212:
          fHistPtYInput[13]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[13]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[13]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[13]->Fill(GetDecayChannel(fMCStack, particle), particle->GetWeight());
          if (alpha>=-1 && fHistPtAlphaInput[13]) fHistPtAlphaInput[13]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0 && fHistPtDeltaPhiInput[13]) fHistPtDeltaPhiInput[13]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        default:
          fHistPdgInputRest->Fill(particle->GetPdgCode());
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
void AliAnalysisTaskGammaCocktailMC::InitializeDecayChannelHist(TH1F* hist, Int_t np) {
  
  switch (np) {
    case 0:
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#gamma #gamma");
      hist->GetXaxis()->SetBinLabel(3,"e^{+} e^{-} #gamma");
      hist->GetXaxis()->SetBinLabel(4,"e^{+} e^{-} e^{+} e^{-}");
      hist->GetXaxis()->SetBinLabel(5,"e^{+} e^{-}");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;

    case 1:
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#gamma #gamma");
      hist->GetXaxis()->SetBinLabel(3,"#pi^{0} #pi^{0} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(4,"#pi^{0} #gamma #gamma");
      hist->GetXaxis()->SetBinLabel(5,"#pi^{+} #pi^{-} #gamma");
      hist->GetXaxis()->SetBinLabel(6,"e^{+} e^{-} #gamma");
      hist->GetXaxis()->SetBinLabel(7,"#mu^{+} #mu^{-} #gamma");
      hist->GetXaxis()->SetBinLabel(8,"e^{+} e^{-} e^{+} e^{-}");
      hist->GetXaxis()->SetBinLabel(9,"#pi^{+} #pi^{-} #gamma #gamma");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;
      
    case 2:
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#pi^{+} #pi^{-} #eta");
      hist->GetXaxis()->SetBinLabel(3,"#rho^{0} #gamma");
      hist->GetXaxis()->SetBinLabel(4,"#pi^{+} #pi^{-} #gamma");
      hist->GetXaxis()->SetBinLabel(5,"#omega #gamma");
      hist->GetXaxis()->SetBinLabel(6,"#gamma #gamma");
      hist->GetXaxis()->SetBinLabel(7,"#mu^{+} #mu^{-} #gamma");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;
      
    case 3:
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#pi^{+} #pi^{-} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(3,"#pi^{0} #gamma");
      hist->GetXaxis()->SetBinLabel(4,"#eta #gamma");
      hist->GetXaxis()->SetBinLabel(5,"#pi^{0} e^{+} e^{-}");
      hist->GetXaxis()->SetBinLabel(6,"e^{+} e^{-}");
      hist->GetXaxis()->SetBinLabel(7,"#pi^{0} #pi^{0} #gamma");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;
      
    case 4:
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#pi^{+} #pi^{-}");
      hist->GetXaxis()->SetBinLabel(3,"#pi^{+} #pi^{-} #gamma");
      hist->GetXaxis()->SetBinLabel(4,"#pi^{0} #gamma");
      hist->GetXaxis()->SetBinLabel(5,"#eta #gamma");
      hist->GetXaxis()->SetBinLabel(6,"#pi^{0} #pi^{0} #gamma");
      hist->GetXaxis()->SetBinLabel(7,"e^{+} e^{-}");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;
      
    case 5:
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#pi^{+} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(3,"#pi^{+} #gamma");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;
      
    case 6:
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#pi^{-} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(3,"#pi^{-} #gamma");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;
      
    case 7:
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"K^{+} K^{-}");
      hist->GetXaxis()->SetBinLabel(3,"K^{0}_{L} K^{0}_{S}");
      hist->GetXaxis()->SetBinLabel(4,"#eta #gamma");
      hist->GetXaxis()->SetBinLabel(5,"#pi^{0} #gamma");
      hist->GetXaxis()->SetBinLabel(6,"e^{+} e^{-}");
      hist->GetXaxis()->SetBinLabel(7,"#eta e^{+} e^{-}");
      hist->GetXaxis()->SetBinLabel(8,"#pi^{+} #pi^{-} #gamma");
      hist->GetXaxis()->SetBinLabel(9,"f_{0}(980) #gamma");
      hist->GetXaxis()->SetBinLabel(10,"#pi^{0} #pi^{0} #gamma");
      hist->GetXaxis()->SetBinLabel(11,"#pi^{0} e^{+} e^{-}");
      hist->GetXaxis()->SetBinLabel(12,"#pi^{0} #eta #gamma");
      hist->GetXaxis()->SetBinLabel(13,"a_{0}(980) #gamma");
      hist->GetXaxis()->SetBinLabel(14,"#eta' #gamma");
      hist->GetXaxis()->SetBinLabel(15,"#mu^{+} #mu^{-} #gamma");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;
      
    case 8:
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"ggg");
      hist->GetXaxis()->SetBinLabel(3,"gg #gamma");
      hist->GetXaxis()->SetBinLabel(4,"e^{+} e^{-}");
      hist->GetXaxis()->SetBinLabel(5,"e^{+} e^{-} #gamma");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;
      
    case 9:
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"n #pi^{-}");
      hist->GetXaxis()->SetBinLabel(3,"X #gamma");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;
      
    case 10:
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"n #pi^{0}");
      hist->GetXaxis()->SetBinLabel(3,"p #pi^{-}");
      hist->GetXaxis()->SetBinLabel(4,"n #gamma");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;

    case 11:
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"n #pi^{+}");
      hist->GetXaxis()->SetBinLabel(3,"p #pi^{0}");
      hist->GetXaxis()->SetBinLabel(4,"p #gamma");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;
      
    case 12:
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"p #pi^{+}");
      hist->GetXaxis()->SetBinLabel(3,"X #gamma");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;
      
    case 13:
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#Lambda #gamma");
      hist->GetXaxis()->SetBinLabel(3,"#Lambda e^{+} e^{-}");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;
      
    default:
      break;
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
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 22 && PdgDaughter->at(2) == 211)
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
      if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 211)
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
        return 6.;
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
      else if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 22)
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
      else if (std::find(PdgDaughter->begin(), PdgDaughter->end(), 22) != PdgDaughter->end())
        return 2.;
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
      else if (std::find(PdgDaughter->begin(), PdgDaughter->end(), 22) != PdgDaughter->end())
        return 2.;
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
void AliAnalysisTaskGammaCocktailMC::FillPythiaBranchingRatio(TH1F* histo, Int_t np) {
  
  Int_t kc, kfdp, nPart, firstChannel, lastChannel;
  Double_t BR, BRtot;
  std::vector<Int_t> pdgCodes;

  switch (np) {
    case 0:
      kc            = (AliPythia6::Instance())->Pycomp(111);
      firstChannel  = (AliPythia6::Instance())->GetMDCY(kc,2);
      lastChannel   = firstChannel + (AliPythia6::Instance())->GetMDCY(kc,3) - 1;
      BRtot         = 0.;
      for (Int_t channel=firstChannel; channel<=lastChannel; channel++) {
        BR          = (AliPythia6::Instance())->GetBRAT(channel);
        BRtot       = BRtot + BR;
        nPart       = 0;
        for (Int_t i=1; i<=5; i++) {
          if ((AliPythia6::Instance())->GetKFDP(channel,i)) {
            pdgCodes.push_back((AliPythia6::Instance())->GetKFDP(channel,i));
            nPart++;
          }
        }
        std::sort(pdgCodes.begin(), pdgCodes.end());
        if (nPart == 2 && pdgCodes[0] == 22 && pdgCodes[1] == 22)
          histo->SetBinContent(2, BR);
        else if (nPart == 3 && pdgCodes[0] == -11 && pdgCodes[1] == 11 && pdgCodes[2] == 22)
          histo->SetBinContent(3, BR);
        else if (nPart == 4 && pdgCodes[0] == -11 && pdgCodes[1] == -11 && pdgCodes[2] == 11 && pdgCodes[3] == 11)
          histo->SetBinContent(4, BR);
        else if (nPart == 2 && pdgCodes[0] == -11 && pdgCodes[1] == 11)
          histo->SetBinContent(5, BR);
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;

    case 1:
      kc            = (AliPythia6::Instance())->Pycomp(221);
      firstChannel  = (AliPythia6::Instance())->GetMDCY(kc,2);
      lastChannel   = firstChannel + (AliPythia6::Instance())->GetMDCY(kc,3) - 1;
      BRtot         = 0.;
      for (Int_t channel=firstChannel; channel<=lastChannel; channel++) {
        BR          = (AliPythia6::Instance())->GetBRAT(channel);
        BRtot       = BRtot + BR;
        nPart       = 0;
        for (Int_t i=1; i<=5; i++) {
          if ((AliPythia6::Instance())->GetKFDP(channel,i)) {
            pdgCodes.push_back((AliPythia6::Instance())->GetKFDP(channel,i));
            nPart++;
          }
        }
        std::sort(pdgCodes.begin(), pdgCodes.end());
        if (nPart == 2 && pdgCodes[0] == 22 && pdgCodes[1] == 22)
          histo->SetBinContent(2, BR);
        else if (nPart == 3 && pdgCodes[0] == 111 && pdgCodes[1] == 111 && pdgCodes[2] == 111)
          histo->SetBinContent(3, BR);
        else if (nPart == 3 && pdgCodes[0] == 22 && pdgCodes[1] == 22 && pdgCodes[2] == 111)
          histo->SetBinContent(4, BR);
        else if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == 22 && pdgCodes[2] == 211)
          histo->SetBinContent(5, BR);
        else if (nPart == 3 && pdgCodes[0] == -11 && pdgCodes[1] == 11 && pdgCodes[2] == 22)
          histo->SetBinContent(6, BR);
        else if (nPart == 3 && pdgCodes[0] == -13 && pdgCodes[1] == 13 && pdgCodes[2] == 22)
          histo->SetBinContent(7, BR);
        else if (nPart == 4 && pdgCodes[0] == -11 && pdgCodes[1] == -11 && pdgCodes[2] == 11 && pdgCodes[3] == 11)
          histo->SetBinContent(8, BR);
        else if (nPart == 4 && pdgCodes[0] == -211 && pdgCodes[1] == 22 && pdgCodes[2] == 22 && pdgCodes[3] == 211)
          histo->SetBinContent(9, BR);
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;

    case 2:
      kc            = (AliPythia6::Instance())->Pycomp(331);
      firstChannel  = (AliPythia6::Instance())->GetMDCY(kc,2);
      lastChannel   = firstChannel + (AliPythia6::Instance())->GetMDCY(kc,3) - 1;
      BRtot         = 0.;
      for (Int_t channel=firstChannel; channel<=lastChannel; channel++) {
        BR          = (AliPythia6::Instance())->GetBRAT(channel);
        BRtot       = BRtot + BR;
        nPart       = 0;
        for (Int_t i=1; i<=5; i++) {
          if ((AliPythia6::Instance())->GetKFDP(channel,i)) {
            pdgCodes.push_back((AliPythia6::Instance())->GetKFDP(channel,i));
            nPart++;
          }
        }
        std::sort(pdgCodes.begin(), pdgCodes.end());
        if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == 211 && pdgCodes[2] == 221)
          histo->SetBinContent(2, BR);
        else if (nPart == 2 && pdgCodes[0] == 22 && pdgCodes[1] == 113)
          histo->SetBinContent(3, BR);
        else if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == 22 && pdgCodes[2] == 211)
          histo->SetBinContent(4, BR);
        else if (nPart == 2 && pdgCodes[0] == 22 && pdgCodes[1] == 223)
          histo->SetBinContent(5, BR);
        else if (nPart == 2 && pdgCodes[0] == 22 && pdgCodes[1] == 22)
          histo->SetBinContent(6, BR);
        else if (nPart == 3 && pdgCodes[0] == -13 && pdgCodes[1] == 13 && pdgCodes[2] == 22)
          histo->SetBinContent(7, BR);
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;
      
    case 3:
      kc            = (AliPythia6::Instance())->Pycomp(223);
      firstChannel  = (AliPythia6::Instance())->GetMDCY(kc,2);
      lastChannel   = firstChannel + (AliPythia6::Instance())->GetMDCY(kc,3) - 1;
      BRtot         = 0.;
      for (Int_t channel=firstChannel; channel<=lastChannel; channel++) {
        BR          = (AliPythia6::Instance())->GetBRAT(channel);
        BRtot       = BRtot + BR;
        nPart       = 0;
        for (Int_t i=1; i<=5; i++) {
          if ((AliPythia6::Instance())->GetKFDP(channel,i)) {
            pdgCodes.push_back((AliPythia6::Instance())->GetKFDP(channel,i));
            nPart++;
          }
        }
        std::sort(pdgCodes.begin(), pdgCodes.end());
        if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == 111 && pdgCodes[2] == 211)
          histo->SetBinContent(2, BR);
        else if (nPart == 2 && pdgCodes[0] == 22 && pdgCodes[1] == 111)
          histo->SetBinContent(3, BR);
        else if (nPart == 2 && pdgCodes[0] == 22 && pdgCodes[1] == 221)
          histo->SetBinContent(4, BR);
        else if (nPart == 3 && pdgCodes[0] == -11 && pdgCodes[1] == 11 && pdgCodes[2] == 111)
          histo->SetBinContent(5, BR);
        else if (nPart == 2 && pdgCodes[0] == -11 && pdgCodes[1] == 11)
          histo->SetBinContent(6, BR);
        else if (nPart == 3 && pdgCodes[0] == 22 && pdgCodes[1] == 111 && pdgCodes[2] == 111)
          histo->SetBinContent(7, BR);
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;

    case 4:
      kc            = (AliPythia6::Instance())->Pycomp(113);
      firstChannel  = (AliPythia6::Instance())->GetMDCY(kc,2);
      lastChannel   = firstChannel + (AliPythia6::Instance())->GetMDCY(kc,3) - 1;
      BRtot         = 0.;
      for (Int_t channel=firstChannel; channel<=lastChannel; channel++) {
        BR          = (AliPythia6::Instance())->GetBRAT(channel);
        BRtot       = BRtot + BR;
        nPart       = 0;
        for (Int_t i=1; i<=5; i++) {
          if ((AliPythia6::Instance())->GetKFDP(channel,i)) {
            pdgCodes.push_back((AliPythia6::Instance())->GetKFDP(channel,i));
            nPart++;
          }
        }
        std::sort(pdgCodes.begin(), pdgCodes.end());
        if (nPart == 2 && pdgCodes[0] == -211 && pdgCodes[1] == 211)
          histo->SetBinContent(2, BR);
        else if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == 22 && pdgCodes[2] == 211)
          histo->SetBinContent(3, BR);
        else if (nPart == 2 && pdgCodes[0] == 22 && pdgCodes[1] == 111)
          histo->SetBinContent(4, BR);
        else if (nPart == 2 && pdgCodes[0] == 22 && pdgCodes[1] == 221)
          histo->SetBinContent(5, BR);
        else if (nPart == 3 && pdgCodes[0] == 22 && pdgCodes[1] == 111 && pdgCodes[2] == 111)
          histo->SetBinContent(6, BR);
        else if (nPart == 2 && pdgCodes[0] == -11 && pdgCodes[1] == 11)
          histo->SetBinContent(7, BR);
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;
      
    case 5:
      kc            = (AliPythia6::Instance())->Pycomp(213);
      firstChannel  = (AliPythia6::Instance())->GetMDCY(kc,2);
      lastChannel   = firstChannel + (AliPythia6::Instance())->GetMDCY(kc,3) - 1;
      BRtot         = 0.;
      for (Int_t channel=firstChannel; channel<=lastChannel; channel++) {
        BR          = (AliPythia6::Instance())->GetBRAT(channel);
        BRtot       = BRtot + BR;
        nPart       = 0;
        for (Int_t i=1; i<=5; i++) {
          if ((AliPythia6::Instance())->GetKFDP(channel,i)) {
            pdgCodes.push_back((AliPythia6::Instance())->GetKFDP(channel,i));
            nPart++;
          }
        }
        std::sort(pdgCodes.begin(), pdgCodes.end());
        if (nPart == 2 && pdgCodes[0] == 111 && pdgCodes[1] == 211)
          histo->SetBinContent(2, BR);
        else if (nPart == 2 && pdgCodes[0] == 22 && pdgCodes[1] == 211)
          histo->SetBinContent(3, BR);
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;
      
    case 6:
      kc            = (AliPythia6::Instance())->Pycomp(213);      // is rho- (-213), but Pycomp handels like rho+ (213)
      firstChannel  = (AliPythia6::Instance())->GetMDCY(kc,2);
      lastChannel   = firstChannel + (AliPythia6::Instance())->GetMDCY(kc,3) - 1;
      BRtot         = 0.;
      for (Int_t channel=firstChannel; channel<=lastChannel; channel++) {
        BR          = (AliPythia6::Instance())->GetBRAT(channel);
        BRtot       = BRtot + BR;
        nPart       = 0;
        for (Int_t i=1; i<=5; i++) {
          if ((AliPythia6::Instance())->GetKFDP(channel,i)) {
            pdgCodes.push_back((AliPythia6::Instance())->GetKFDP(channel,i));
            nPart++;
          }
        }
        std::sort(pdgCodes.begin(), pdgCodes.end());
        if (nPart == 2 && pdgCodes[0] == 111 && pdgCodes[1] == 211)
          histo->SetBinContent(2, BR);
        else if (nPart == 2 && pdgCodes[0] == 22 && pdgCodes[1] == 211)
          histo->SetBinContent(3, BR);
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;
      
    case 7:
      kc            = (AliPythia6::Instance())->Pycomp(333);
      firstChannel  = (AliPythia6::Instance())->GetMDCY(kc,2);
      lastChannel   = firstChannel + (AliPythia6::Instance())->GetMDCY(kc,3) - 1;
      BRtot         = 0.;
      for (Int_t channel=firstChannel; channel<=lastChannel; channel++) {
        BR          = (AliPythia6::Instance())->GetBRAT(channel);
        BRtot       = BRtot + BR;
        nPart       = 0;
        for (Int_t i=1; i<=5; i++) {
          if ((AliPythia6::Instance())->GetKFDP(channel,i)) {
            pdgCodes.push_back((AliPythia6::Instance())->GetKFDP(channel,i));
            nPart++;
          }
        }
        std::sort(pdgCodes.begin(), pdgCodes.end());
        if (nPart == 2 && pdgCodes[0] == -321 && pdgCodes[1] == 321)
          histo->SetBinContent(2, BR);
        else if (nPart == 2 && pdgCodes[0] == 130 && pdgCodes[1] == 310)
          histo->SetBinContent(3, BR);
        else if (nPart == 2 && pdgCodes[0] == 22 && pdgCodes[1] == 221)
          histo->SetBinContent(4, BR);
        else if (nPart == 2 && pdgCodes[0] == 22 && pdgCodes[1] == 111)
          histo->SetBinContent(5, BR);
        else if (nPart == 2 && pdgCodes[0] == -11 && pdgCodes[1] == 11)
          histo->SetBinContent(6, BR);
        else if (nPart == 3 && pdgCodes[0] == -11 && pdgCodes[1] == 11 && pdgCodes[2] == 221)
          histo->SetBinContent(7, BR);
        else if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == 22 && pdgCodes[2] == 211)
          histo->SetBinContent(8, BR);
        else if (nPart == 2 && pdgCodes[0] == 22 && pdgCodes[1] == 9010221)
          histo->SetBinContent(9, BR);
        else if (nPart == 3 && pdgCodes[0] == 22 && pdgCodes[1] == 111 && pdgCodes[2] == 111)
          histo->SetBinContent(10, BR);
        else if (nPart == 3 && pdgCodes[0] == -11 && pdgCodes[1] == 11 && pdgCodes[2] == 111)
          histo->SetBinContent(11, BR);
        else if (nPart == 3 && pdgCodes[0] == 22 && pdgCodes[1] == 111 && pdgCodes[2] == 221)
          histo->SetBinContent(12, BR);
        else if (nPart == 2 && pdgCodes[0] == 22 && pdgCodes[1] == 9000111)
          histo->SetBinContent(13, BR);
        else if (nPart == 2 && pdgCodes[0] == 22 && pdgCodes[1] == 331)
          histo->SetBinContent(14, BR);
        else if (nPart == 3 && pdgCodes[0] == -13 && pdgCodes[1] == 13 && pdgCodes[2] == 22)
          histo->SetBinContent(15, BR);
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;
      
    case 8:
      kc            = (AliPythia6::Instance())->Pycomp(443);
      firstChannel  = (AliPythia6::Instance())->GetMDCY(kc,2);
      lastChannel   = firstChannel + (AliPythia6::Instance())->GetMDCY(kc,3) - 1;
      BRtot         = 0.;
      for (Int_t channel=firstChannel; channel<=lastChannel; channel++) {
        BR          = (AliPythia6::Instance())->GetBRAT(channel);
        BRtot       = BRtot + BR;
        nPart       = 0;
        for (Int_t i=1; i<=5; i++) {
          if ((AliPythia6::Instance())->GetKFDP(channel,i)) {
            pdgCodes.push_back((AliPythia6::Instance())->GetKFDP(channel,i));
            nPart++;
          }
        }
        std::sort(pdgCodes.begin(), pdgCodes.end());
        if (nPart == 3 && (pdgCodes[0] == 21 || pdgCodes[0] == 9) && (pdgCodes[1] == 21 || pdgCodes[1] == 9) && (pdgCodes[2] == 21 || pdgCodes[2] == 9))
          histo->SetBinContent(2, BR);
        else if (nPart == 3 && (pdgCodes[0] == 21 || pdgCodes[0] == 9) && (pdgCodes[1] == 21 || pdgCodes[1] == 9) && pdgCodes[2] == 22)
          histo->SetBinContent(3, BR);
        else if (nPart == 2 && pdgCodes[0] == -11 && pdgCodes[1] == 11)
          histo->SetBinContent(4, BR);
        else if (nPart == 3 && pdgCodes[0] == -11 && pdgCodes[1] == 11 && pdgCodes[2] == 22)
          histo->SetBinContent(5, BR);
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;
      
    case 9:
      kc            = (AliPythia6::Instance())->Pycomp(1114);
      firstChannel  = (AliPythia6::Instance())->GetMDCY(kc,2);
      lastChannel   = firstChannel + (AliPythia6::Instance())->GetMDCY(kc,3) - 1;
      BRtot         = 0.;
      for (Int_t channel=firstChannel; channel<=lastChannel; channel++) {
        BR          = (AliPythia6::Instance())->GetBRAT(channel);
        BRtot       = BRtot + BR;
        nPart       = 0;
        for (Int_t i=1; i<=5; i++) {
          if ((AliPythia6::Instance())->GetKFDP(channel,i)) {
            pdgCodes.push_back((AliPythia6::Instance())->GetKFDP(channel,i));
            nPart++;
          }
        }
        std::sort(pdgCodes.begin(), pdgCodes.end());
        if (nPart == 2 && pdgCodes[0] == -211 && pdgCodes[0] == 2112)
          histo->SetBinContent(2, BR);
        else if (std::find(pdgCodes.begin(), pdgCodes.end(), 22) != pdgCodes.end())
          histo->SetBinContent(3, BR+histo->GetBinContent(3));
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;
      
    case 10:
      kc            = (AliPythia6::Instance())->Pycomp(2114);
      firstChannel  = (AliPythia6::Instance())->GetMDCY(kc,2);
      lastChannel   = firstChannel + (AliPythia6::Instance())->GetMDCY(kc,3) - 1;
      BRtot         = 0.;
      for (Int_t channel=firstChannel; channel<=lastChannel; channel++) {
        BR          = (AliPythia6::Instance())->GetBRAT(channel);
        BRtot       = BRtot + BR;
        nPart       = 0;
        for (Int_t i=1; i<=5; i++) {
          if ((AliPythia6::Instance())->GetKFDP(channel,i)) {
            pdgCodes.push_back((AliPythia6::Instance())->GetKFDP(channel,i));
            nPart++;
          }
        }
        std::sort(pdgCodes.begin(), pdgCodes.end());
        if (nPart == 2 && pdgCodes[0] == 111 && pdgCodes[1] == 2112)
          histo->SetBinContent(2, BR);
        else if (nPart == 2 && pdgCodes[0] == -211 && pdgCodes[1] == 2212)
          histo->SetBinContent(3, BR);
        else if (nPart == 2 && pdgCodes[0] == 22 && pdgCodes[1] == 2112)
          histo->SetBinContent(4, BR);
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;
      
    case 11:
      kc            = (AliPythia6::Instance())->Pycomp(2214);
      firstChannel  = (AliPythia6::Instance())->GetMDCY(kc,2);
      lastChannel   = firstChannel + (AliPythia6::Instance())->GetMDCY(kc,3) - 1;
      BRtot         = 0.;
      for (Int_t channel=firstChannel; channel<=lastChannel; channel++) {
        BR          = (AliPythia6::Instance())->GetBRAT(channel);
        BRtot       = BRtot + BR;
        nPart       = 0;
        for (Int_t i=1; i<=5; i++) {
          if ((AliPythia6::Instance())->GetKFDP(channel,i)) {
            pdgCodes.push_back((AliPythia6::Instance())->GetKFDP(channel,i));
            nPart++;
          }
        }
        std::sort(pdgCodes.begin(), pdgCodes.end());
        if (nPart == 2 && pdgCodes[0] == 211 && pdgCodes[1] == 2112)
          histo->SetBinContent(2, BR);
        else if (nPart == 2 && pdgCodes[0] == 111 && pdgCodes[1] == 2212)
          histo->SetBinContent(3, BR);
        else if (nPart == 2 && pdgCodes[0] == 22 && pdgCodes[1] == 2212)
          histo->SetBinContent(4, BR);
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;
      
    case 12:
      kc            = (AliPythia6::Instance())->Pycomp(2224);
      firstChannel  = (AliPythia6::Instance())->GetMDCY(kc,2);
      lastChannel   = firstChannel + (AliPythia6::Instance())->GetMDCY(kc,3) - 1;
      BRtot         = 0.;
      for (Int_t channel=firstChannel; channel<=lastChannel; channel++) {
        BR          = (AliPythia6::Instance())->GetBRAT(channel);
        BRtot       = BRtot + BR;
        nPart       = 0;
        for (Int_t i=1; i<=5; i++) {
          if ((AliPythia6::Instance())->GetKFDP(channel,i)) {
            pdgCodes.push_back((AliPythia6::Instance())->GetKFDP(channel,i));
            nPart++;
          }
        }
        std::sort(pdgCodes.begin(), pdgCodes.end());
        if (nPart == 2 && pdgCodes[0] == 211 && pdgCodes[0] == 2212)
          histo->SetBinContent(2, BR);
        else if (std::find(pdgCodes.begin(), pdgCodes.end(), 22) != pdgCodes.end())
          histo->SetBinContent(3, BR+histo->GetBinContent(3));
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;
      
    case 13:
      kc            = (AliPythia6::Instance())->Pycomp(3212);
      firstChannel  = (AliPythia6::Instance())->GetMDCY(kc,2);
      lastChannel   = firstChannel + (AliPythia6::Instance())->GetMDCY(kc,3) - 1;
      BRtot         = 0.;
      for (Int_t channel=firstChannel; channel<=lastChannel; channel++) {
        BR          = (AliPythia6::Instance())->GetBRAT(channel);
        BRtot       = BRtot + BR;
        nPart       = 0;
        for (Int_t i=1; i<=5; i++) {
          if ((AliPythia6::Instance())->GetKFDP(channel,i)) {
            pdgCodes.push_back((AliPythia6::Instance())->GetKFDP(channel,i));
            nPart++;
          }
        }
        std::sort(pdgCodes.begin(), pdgCodes.end());
        if (nPart == 2 && pdgCodes[0] == 22 && pdgCodes[1] == 3122)
          histo->SetBinContent(2, BR);
        else if (nPart == 3 && pdgCodes[0] == -11 && pdgCodes[1] == 11 && pdgCodes[2] == 3122)
          histo->SetBinContent(3, BR);
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;
      
    default:
      break;
  }
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

