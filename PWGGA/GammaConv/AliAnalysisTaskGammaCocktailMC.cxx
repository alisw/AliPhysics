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
  fMCGenHandler(NULL),
  fMCGenerator(NULL),
  fMCCocktailGen(NULL),
  fDoLightOutput(kFALSE),
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
  fHistPtGammaSourcePtInput(NULL),
  fHistPhiGammaSourcePhiInput(NULL),
  fHistPdgInputRest(NULL),
  fHistPdgGammaSourceRest(NULL),
  fParticleList(NULL),
  fParticleListNames(NULL),
  fPtParametrization{NULL},
  fPtParametrizationProton(NULL),
  fCocktailSettings{NULL},
  fMtScalingFactors(NULL),
  fPtYDistributions{NULL},
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
  fMCGenHandler(NULL),
  fMCGenerator(NULL),
  fMCCocktailGen(NULL),
  fDoLightOutput(kFALSE),
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
  fHistPtGammaSourcePtInput(NULL),
  fHistPhiGammaSourcePhiInput(NULL),
  fHistPdgInputRest(NULL),
  fHistPdgGammaSourceRest(NULL),
  fParticleList(NULL),
  fParticleListNames(NULL),
  fPtParametrization{NULL},
  fPtParametrizationProton(NULL),
  fCocktailSettings{NULL},
  fMtScalingFactors(NULL),
  fPtYDistributions{NULL},
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
  for (Int_t i=0; i<10; i++) {
    if (fCocktailSettings[i]) delete fCocktailSettings[i];
  }
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
  
  fMCGenHandler                           = (AliMCGenHandler*)AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler();
  fMCGenerator                            = fMCGenHandler->GetGenerator();
  TString mcGeneratorClassName            = "";
  if (fMCGenerator) mcGeneratorClassName  = fMCGenerator->ClassName();
  
  if (mcGeneratorClassName.CompareTo("AliGenEMCocktailV2") == 0) {
    
    fMCCocktailGen = (AliGenEMCocktailV2*)fMCGenerator;
    
    // has mother i
    SetHasMother((UInt_t)fMCCocktailGen->GetSelectedMothers());
    
    // pt parametrizations
    GetAndSetPtParametrizations(fMCCocktailGen);
    for (Int_t i=0; i<17; i++) {
      if (fHasMother[i]) fUserInfo->Add(fPtParametrization[i]);
    }
    if (fPtParametrizationProton) fUserInfo->Add(fPtParametrizationProton);
    
    // cocktail settings
    Double_t ptMin, ptMax;
    fMCCocktailGen->GetPtRange(ptMin, ptMax);
    fCocktailSettings[0] = new TObjString(Form("collSys_%d",fMCCocktailGen->GetCollisionSystem()));
    fCocktailSettings[1] = new TObjString(Form("cent_%d",fMCCocktailGen->GetCentrality()));
    fCocktailSettings[2] = new TObjString(Form("decayMode_%.0f",fMCCocktailGen->GetDecayMode()));
    fCocktailSettings[3] = new TObjString(Form("selectMothers_%d",fMCCocktailGen->GetSelectedMothers()));
    fCocktailSettings[4] = new TObjString(Form("paramFile_%s",(fMCCocktailGen->GetParametrizationFile()).Data()));
    fCocktailSettings[5] = new TObjString(Form("paramDir_%s",(fMCCocktailGen->GetParametrizationFileDirectory()).Data()));
    fCocktailSettings[6] = new TObjString(Form("nParticles_%d",fMCCocktailGen->GetNumberOfParticles()));
    fCocktailSettings[7] = new TObjString(Form("ptMin_%.2f",ptMin));
    fCocktailSettings[8] = new TObjString(Form("ptMax_%.2f",ptMax));
    fCocktailSettings[9] = new TObjString(Form("weightMode_%.0f",fMCCocktailGen->GetWeightingMode()));
    fCocktailSettings[10] = new TObjString(Form("dynamicalPtRang_%d",fMCCocktailGen->GetDynamicalPtRangeOption()));
    fCocktailSettings[11] = new TObjString(Form("yWeights_%d",fMCCocktailGen->GetYWeightOption()));
    for (Int_t i=0; i<12; i++) fUserInfo->Add(fCocktailSettings[i]);
    
    // mt scaling params
    fMtScalingFactors = (TH1D*)fMCCocktailGen->GetMtScalingFactors();
    fUserInfo->Add(fMtScalingFactors);

    // pt-y distributions
    GetAndSetPtYDistributions(fMCCocktailGen);
    for (Int_t i=0; i<17; i++) {
      if (fHasMother[i]) fUserInfo->Add(fPtYDistributions[i]);
    }
  } else {
    for (Int_t i=0; i<17; i++) fHasMother[i] = kTRUE;
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
  const Int_t nInputParticles         = 17;
  Int_t   fParticleList_local[]       = {111,221,331,223,113,213,-213,333,443,1114,2114,2214,2224,3212,310,130,3122};
  TString fParticleListNames_local[]  = {"Pi0","Eta","EtaPrim","omega","rho0","rho+","rho-","phi","J/psi","Delta-","Delta0","Delta+","Delta++","Sigma0","K0s","K0l","Lambda"};
  fParticleList                   = fParticleList_local;
  fParticleListNames              = fParticleListNames_local;
  fHistPtYInput                   = new TH2F*[nInputParticles];
  fHistPtPhiInput                 = new TH2F*[nInputParticles];
  fHistDecayChannelsInput         = new TH1F*[nInputParticles];
  fHistPythiaBR                   = new TH1F*[nInputParticles];
  fHistPtYGammaSource             = new TH2F*[nInputParticles];
  fHistPtPhiGammaSource           = new TH2F*[nInputParticles];
  fHistPtAlphaInput               = new TH2F*[nInputParticles];
  fHistPtDeltaPhiInput            = new TH2F*[nInputParticles];
  fHistPtGammaSourcePtInput       = new TH2F*[nInputParticles];
  fHistPhiGammaSourcePhiInput     = new TH2F*[nInputParticles];
  
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
    
    if (fHasMother[i]) {
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
      fHistPtGammaSourcePtInput[i] = (TH2F*)SetHist2D(fHistPtGammaSourcePtInput[i],"f",Form("PtGamma_PtMother_%s",fParticleListNames[i].Data()),"#it{p}_{T,daughter}","#it{p}_{T,mother}",500,0,50,500,0,50,kTRUE);
      fOutputContainer->Add(fHistPtGammaSourcePtInput[i]);
      
      fHistPhiGammaSourcePhiInput[i] = (TH2F*)SetHist2D(fHistPhiGammaSourcePhiInput[i],"f",Form("PhiGamma_PhiMother_%s",fParticleListNames[i].Data()),"#phi_{daughter}","#phi_{mother}",100,0,7,100,0,7,kTRUE);
      fOutputContainer->Add(fHistPhiGammaSourcePhiInput[i]);
      
      // decay channels mother
      fHistDecayChannelsInput[i] = (TH1F*)SetHist1D(fHistDecayChannelsInput[i],"f",Form("DecayChannels_%s",fParticleListNames[i].Data()),"","", 20,-0.5,19.5,kTRUE);
      InitializeDecayChannelHist(fHistDecayChannelsInput[i], i);
      fOutputContainer->Add(fHistDecayChannelsInput[i]);
      
      // BR from pythia
      fHistPythiaBR[i] = (TH1F*)SetHist1D(fHistPythiaBR[i],"f",Form("PythiaBR_%s",fParticleListNames[i].Data()),"","", 20,-0.5,19.5,kTRUE);
      InitializeDecayChannelHist(fHistPythiaBR[i], i);
      FillPythiaBranchingRatio(fHistPythiaBR[i], i);
      fUserInfo->Add(fHistPythiaBR[i]);
    } else {
      fHistPtYInput[i]                = NULL;
      fHistPtYGammaSource[i]          = NULL;
      fHistPtPhiInput[i]              = NULL;
      fHistPtPhiGammaSource[i]        = NULL;
      fHistPtGammaSourcePtInput[i]    = NULL;
      fHistPhiGammaSourcePhiInput[i]  = NULL;
      fHistDecayChannelsInput[i]      = NULL;
      fHistPythiaBR[i]                = NULL;
    }
    
    // lightweight output
    if (!fDoLightOutput || (fDoLightOutput && (i < 4 || i == 7))) {
      
      if (fHasMother[i]) {
        // gamma delta phi
        fHistPtDeltaPhiInput[i] = (TH2F*)SetHist2D(fHistPtDeltaPhiInput[i],"f",Form("Pt_DeltaPhi_%s",fParticleListNames[i].Data()),"#it{p}_{T}","#Delta#phi_{#gamma_{1}#gamma_{2}}",500,0,50,82,binsDeltaPhi,kTRUE);
        fOutputContainer->Add(fHistPtDeltaPhiInput[i]);
        
        // alpha mother
        fHistPtAlphaInput[i] = (TH2F*)SetHist2D(fHistPtAlphaInput[i],"f",Form("Pt_Alpha_%s",fParticleListNames[i].Data()),"#it{p}_{T}","#alpha",500,0,50,100,-1,1,kTRUE);
        fOutputContainer->Add(fHistPtAlphaInput[i]);
      } else {
        fHistPtDeltaPhiInput[i] = NULL;
        fHistPtAlphaInput[i]    = NULL;
      }
    } else {
      fHistPtDeltaPhiInput[i] = NULL;
      fHistPtAlphaInput[i]    = NULL;
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
  
  fMCEvent = MCEvent();
  if(fMCEvent == NULL) fIsMC = 0;
  
  if (fIsMC==0) return;
  
  fHistNEvents->Fill(0.5);
  ProcessMCParticles();

  PostData(1, fOutputContainer);
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaCocktailMC::GetAndSetPtParametrizations(AliGenEMCocktailV2* fMCCocktailGen)
{
  if (!fMCCocktailGen) return;

  for (Int_t i=0; i<17; i++) fPtParametrization[i] = NULL;
  fPtParametrizationProton = NULL;
  
  TF1* fct        = NULL;
  TString fctName = "";
  for (Int_t i=0; i<19; i++) {
    fct = (TF1*)fMCCocktailGen->GetPtParametrization(i);
    if (fct) {
      fctName = fct->GetName();
      if (fctName.BeginsWith("111_pt")  && fHasMother[0])  fPtParametrization[0]   = fct;
      if (fctName.BeginsWith("221_pt")  && fHasMother[1])  fPtParametrization[1]   = fct;
      if (fctName.BeginsWith("331_pt")  && fHasMother[2])  fPtParametrization[2]   = fct;
      if (fctName.BeginsWith("223_pt")  && fHasMother[3])  fPtParametrization[3]   = fct;
      if (fctName.BeginsWith("113_pt")  && fHasMother[4])  fPtParametrization[4]   = fct;
      if (fctName.BeginsWith("213_pt")  && fHasMother[5])  fPtParametrization[5]   = fct;
      if (fctName.BeginsWith("-213_pt") && fHasMother[6])  fPtParametrization[6]   = fct;
      if (fctName.BeginsWith("333_pt")  && fHasMother[7])  fPtParametrization[7]   = fct;
      if (fctName.BeginsWith("443_pt")  && fHasMother[8])  fPtParametrization[8]   = fct;
      if (fctName.BeginsWith("1114_pt") && fHasMother[9])  fPtParametrization[9]   = fct;
      if (fctName.BeginsWith("2114_pt") && fHasMother[10]) fPtParametrization[10]  = fct;
      if (fctName.BeginsWith("2214_pt") && fHasMother[11]) fPtParametrization[11]  = fct;
      if (fctName.BeginsWith("2224_pt") && fHasMother[12]) fPtParametrization[12]  = fct;
      if (fctName.BeginsWith("3212_pt") && fHasMother[13]) fPtParametrization[13]  = fct;
      if (fctName.BeginsWith("310_pt")  && fHasMother[14]) fPtParametrization[14]  = fct;
      if (fctName.BeginsWith("130_pt")  && fHasMother[15]) fPtParametrization[15]  = fct;
      if (fctName.BeginsWith("3122_pt") && fHasMother[16]) fPtParametrization[16]  = fct;
      if (fctName.BeginsWith("2212_pt")) fPtParametrizationProton = fct;
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaCocktailMC::GetAndSetPtYDistributions(AliGenEMCocktailV2* fMCCocktailGen)
{
  if (!fMCCocktailGen) return;

  for (Int_t i=0; i<17; i++) fPtYDistributions[i] = NULL;

  TH2F* tempPtY = NULL;
  TString tempPtYName = "";
  for (Int_t i=0; i<18; i++) {
    tempPtY = (TH2F*)fMCCocktailGen->GetPtYDistribution(i);
    if (tempPtY) {
      tempPtYName = tempPtY->GetName();
      if (tempPtYName.BeginsWith("111_pt_y")  && fHasMother[0])  fPtYDistributions[0]   = tempPtY;
      if (tempPtYName.BeginsWith("221_pt_y")  && fHasMother[1])  fPtYDistributions[1]   = tempPtY;
      if (tempPtYName.BeginsWith("331_pt_y")  && fHasMother[2])  fPtYDistributions[2]   = tempPtY;
      if (tempPtYName.BeginsWith("223_pt_y")  && fHasMother[3])  fPtYDistributions[3]   = tempPtY;
      if (tempPtYName.BeginsWith("113_pt_y")  && fHasMother[4])  fPtYDistributions[4]   = tempPtY;
      if (tempPtYName.BeginsWith("213_pt_y")  && fHasMother[5])  fPtYDistributions[5]   = tempPtY;
      if (tempPtYName.BeginsWith("-213_pt_y") && fHasMother[6])  fPtYDistributions[6]   = tempPtY;
      if (tempPtYName.BeginsWith("333_pt_y")  && fHasMother[7])  fPtYDistributions[7]   = tempPtY;
      if (tempPtYName.BeginsWith("443_pt_y")  && fHasMother[8])  fPtYDistributions[8]   = tempPtY;
      if (tempPtYName.BeginsWith("1114_pt_y") && fHasMother[9])  fPtYDistributions[9]   = tempPtY;
      if (tempPtYName.BeginsWith("2114_pt_y") && fHasMother[10]) fPtYDistributions[10]  = tempPtY;
      if (tempPtYName.BeginsWith("2214_pt_y") && fHasMother[11]) fPtYDistributions[11]  = tempPtY;
      if (tempPtYName.BeginsWith("2224_pt_y") && fHasMother[12]) fPtYDistributions[12]  = tempPtY;
      if (tempPtYName.BeginsWith("3212_pt_y") && fHasMother[13]) fPtYDistributions[13]  = tempPtY;
      if (tempPtYName.BeginsWith("310_pt_y")  && fHasMother[14]) fPtYDistributions[14]  = tempPtY;
      if (tempPtYName.BeginsWith("130_pt_y")  && fHasMother[15]) fPtYDistributions[15]  = tempPtY;
      if (tempPtYName.BeginsWith("3122_pt_y") && fHasMother[16]) fPtYDistributions[16]  = tempPtY;
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaCocktailMC::SetHasMother(UInt_t selectedMothers) {
  
  for (Int_t i=0; i<17; i++) fHasMother[i] = kFALSE;
  
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
  if (selectedMothers&AliGenEMCocktailV2::kGenK0s)       fHasMother[14] = kTRUE;
  if (selectedMothers&AliGenEMCocktailV2::kGenK0l)       fHasMother[15] = kTRUE;
  if (selectedMothers&AliGenEMCocktailV2::kGenLambda)    fHasMother[16] = kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskGammaCocktailMC::ProcessMCParticles(){

  // Loop over all primary MC particle
  for(Long_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {
    // fill primary histograms
    TParticle* particle         = NULL;
    particle                    = (TParticle *)fMCEvent->Particle(i);
    if (!particle) continue;
    Bool_t hasMother            = kFALSE;
    Bool_t particleIsPrimary    = kTRUE;

    if (particle->GetMother(0)>-1){
      hasMother = kTRUE;
      particleIsPrimary = kFALSE;
    }
    TParticle* motherParticle   = NULL;
    if( hasMother ) motherParticle = (TParticle *)fMCEvent->Particle(particle->GetMother(0));
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

    if (!(TMath::Abs(particle->Energy()-particle->Pz())>0.)) continue;
    Double_t yPre = (particle->Energy()+particle->Pz())/(particle->Energy()-particle->Pz());
//     cout << i << "\t"<< particle->GetPdgCode() << "\t"<< particle->Pz() << "\t" << particle->Energy()<< "\t" << particle->Energy()-particle->Pz() << "\t"<< yPre << endl;
    if (yPre == 0.) continue;
    
    Double_t y = 0.5*TMath::Log(yPre); 
    if (TMath::Abs(y) > fMaxY) continue;
    
    if(particle->GetPdgCode()==22 && hasMother==kTRUE){
      if(motherIsPrimary && fHasMother[GetParticlePosLocal(motherParticle->GetPdgCode())]){
        fHistPtYGamma->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
        fHistPtPhiGamma->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
        switch(motherParticle->GetPdgCode()){
        case 111:
          fHistPtYGammaSource[0]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[0]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourcePtInput[0]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourcePhiInput[0]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 221:
          fHistPtYGammaSource[1]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[1]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourcePtInput[1]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourcePhiInput[1]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 331:
          fHistPtYGammaSource[2]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[2]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourcePtInput[2]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourcePhiInput[2]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 223:
          fHistPtYGammaSource[3]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[3]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourcePtInput[3]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourcePhiInput[3]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 113:
          fHistPtYGammaSource[4]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[4]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourcePtInput[4]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourcePhiInput[4]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 213:
          fHistPtYGammaSource[5]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[5]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourcePtInput[5]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourcePhiInput[5]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case -213:
          fHistPtYGammaSource[6]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[6]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourcePtInput[6]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourcePhiInput[6]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 333:
          fHistPtYGammaSource[7]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[7]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourcePtInput[7]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourcePhiInput[7]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 443:
          fHistPtYGammaSource[8]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[8]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourcePtInput[8]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourcePhiInput[8]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 1114:
          fHistPtYGammaSource[9]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[9]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourcePtInput[9]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourcePhiInput[9]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 2114:
          fHistPtYGammaSource[10]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[10]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourcePtInput[10]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourcePhiInput[10]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 2214:
          fHistPtYGammaSource[11]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[11]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourcePtInput[11]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourcePhiInput[11]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 2224:
          fHistPtYGammaSource[12]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[12]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourcePtInput[12]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourcePhiInput[12]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 3212:
          fHistPtYGammaSource[13]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[13]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourcePtInput[13]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourcePhiInput[13]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 310:
          fHistPtYGammaSource[14]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[14]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourcePtInput[14]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourcePhiInput[14]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 130:
          fHistPtYGammaSource[15]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[15]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourcePtInput[15]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourcePhiInput[15]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        case 3122:
          fHistPtYGammaSource[16]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiGammaSource[16]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistPtGammaSourcePtInput[16]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
          fHistPhiGammaSourcePhiInput[16]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
          break;
        default:
          fHistPdgGammaSourceRest->Fill(motherParticle->GetPdgCode());
          break;
        }
      }
    }
    
    if(particle->GetPdgCode()!=22 && particleIsPrimary && fHasMother[GetParticlePosLocal(particle->GetPdgCode())]){

      Double_t alpha    = -999;
      Double_t deltaPhi = -999;
      if (particle->GetNDaughters() == 2) {
        daughter0 = (TParticle*)fMCEvent->Particle(particle->GetFirstDaughter());
        daughter1 = (TParticle*)fMCEvent->Particle(particle->GetLastDaughter());
                
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
          fHistDecayChannelsInput[0]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          if (alpha>=-1) fHistPtAlphaInput[0]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0) fHistPtDeltaPhiInput[0]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 221:
          fHistPtYInput[1]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[1]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[1]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[1]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          if (alpha>=-1) fHistPtAlphaInput[1]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0) fHistPtDeltaPhiInput[1]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 331:
          fHistPtYInput[2]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[2]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[2]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[2]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          if (alpha>=-1) fHistPtAlphaInput[2]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0) fHistPtDeltaPhiInput[2]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 223:
          fHistPtYInput[3]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[3]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[3]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[3]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          if (alpha>=-1) fHistPtAlphaInput[3]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0) fHistPtDeltaPhiInput[3]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 113:
          fHistPtYInput[4]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[4]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[4]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[4]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          if (alpha>=-1 && fHistPtAlphaInput[4]) fHistPtAlphaInput[4]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0 && fHistPtDeltaPhiInput[4]) fHistPtDeltaPhiInput[4]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 213:
          fHistPtYInput[5]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[5]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[5]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[5]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          if (alpha>=-1 && fHistPtAlphaInput[5]) fHistPtAlphaInput[5]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0 && fHistPtDeltaPhiInput[5]) fHistPtDeltaPhiInput[5]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case -213:
          fHistPtYInput[6]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[6]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[6]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[6]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          if (alpha>=-1 && fHistPtAlphaInput[6]) fHistPtAlphaInput[6]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0 && fHistPtDeltaPhiInput[6]) fHistPtDeltaPhiInput[6]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 333:
          fHistPtYInput[7]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[7]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[7]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[7]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          if (alpha >= -1) fHistPtAlphaInput[7]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0) fHistPtDeltaPhiInput[7]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 443:
          fHistPtYInput[8]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[8]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[8]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[8]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          if (alpha>=-1 && fHistPtAlphaInput[8]) fHistPtAlphaInput[8]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0 && fHistPtDeltaPhiInput[8]) fHistPtDeltaPhiInput[8]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 1114:
          fHistPtYInput[9]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[9]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[9]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[9]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          if (alpha>=-1 && fHistPtAlphaInput[9]) fHistPtAlphaInput[9]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0 && fHistPtDeltaPhiInput[9]) fHistPtDeltaPhiInput[9]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 2114:
          fHistPtYInput[10]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[10]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[10]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[10]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          if (alpha>=-1 && fHistPtAlphaInput[10]) fHistPtAlphaInput[10]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0 && fHistPtDeltaPhiInput[10]) fHistPtDeltaPhiInput[10]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 2214:
          fHistPtYInput[11]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[11]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[11]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[11]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          if (alpha>=-1 && fHistPtAlphaInput[11]) fHistPtAlphaInput[11]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0 && fHistPtDeltaPhiInput[11]) fHistPtDeltaPhiInput[11]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 2224:
          fHistPtYInput[12]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[12]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[12]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[12]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          if (alpha>=-1 && fHistPtAlphaInput[12]) fHistPtAlphaInput[12]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0 && fHistPtDeltaPhiInput[12]) fHistPtDeltaPhiInput[12]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 3212:
          fHistPtYInput[13]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[13]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[13]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[13]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          if (alpha>=-1 && fHistPtAlphaInput[13]) fHistPtAlphaInput[13]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0 && fHistPtDeltaPhiInput[13]) fHistPtDeltaPhiInput[13]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 310:
          fHistPtYInput[14]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[14]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[14]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[14]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          if (alpha>=-1 && fHistPtAlphaInput[14]) fHistPtAlphaInput[14]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0 && fHistPtDeltaPhiInput[14]) fHistPtDeltaPhiInput[14]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 130:
          fHistPtYInput[15]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[15]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[15]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[15]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          if (alpha>=-1 && fHistPtAlphaInput[15]) fHistPtAlphaInput[15]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0 && fHistPtDeltaPhiInput[15]) fHistPtDeltaPhiInput[15]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
          break;
        case 3122:
          fHistPtYInput[16]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[16]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[16]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[16]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          if (alpha>=-1 && fHistPtAlphaInput[16]) fHistPtAlphaInput[16]->Fill(particle->Pt(), alpha, particle->GetWeight());
          if (deltaPhi>=0 && fHistPtDeltaPhiInput[16]) fHistPtDeltaPhiInput[16]->Fill(particle->Pt(), deltaPhi, particle->GetWeight());
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

    case 14:
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#pi^{0} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(3,"#pi^{+} #pi^{-}");
      hist->GetXaxis()->SetBinLabel(4,"#pi^{+} #pi^{-} #gamma");
      hist->GetXaxis()->SetBinLabel(5,"#pi^{0} #gamma #gamma");
      hist->GetXaxis()->SetBinLabel(6,"#gamma #gamma");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;

    case 15:
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#pi^{0} #pi^{0} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(3,"#pi^{+} #pi^{-} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(4,"#pi^{#pm} e^{#mp} #nu #gamma");
      hist->GetXaxis()->SetBinLabel(5,"#pi^{#pm} #mu^{#mp} #nu #gamma");
      hist->GetXaxis()->SetBinLabel(6,"#pi^{+} #pi^{-} #gamma");
      hist->GetXaxis()->SetBinLabel(7,"#pi^{0} #gamma #gamma");
      hist->GetXaxis()->SetBinLabel(8,"#pi^{0} e^{+} e^{-} #gamma");
      hist->GetXaxis()->SetBinLabel(9,"#gamma #gamma");
      hist->GetXaxis()->SetBinLabel(10,"e^{+} e^{-} #gamma");
      hist->GetXaxis()->SetBinLabel(11,"e^{+} e^{-} #gamma #gamma");
      hist->GetXaxis()->SetBinLabel(12,"#mu^{+} #mu^{-} #gamma");
      hist->GetXaxis()->SetBinLabel(13,"#mu^{+} #mu^{-} #gamma #gamma");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;
      
    case 16:
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"p #pi^{-}");
      hist->GetXaxis()->SetBinLabel(3,"n #pi^{0}");
      hist->GetXaxis()->SetBinLabel(4,"n #gamma");
      hist->GetXaxis()->SetBinLabel(5,"p #pi^{-} #gamma");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;

    default:
      break;
  }
}

//_________________________________________________________________________________
Float_t AliAnalysisTaskGammaCocktailMC::GetDecayChannel(AliMCEvent *mcEvent, TParticle* part) {
    
  Int_t nDaughters = part->GetNDaughters();
  if (nDaughters > 10) return 19.;
  
  std::vector<Long64_t> *PdgDaughter = new std::vector<Long64_t>(nDaughters);
  Long64_t tempPdgCode = 0;
  for (Int_t i=0; i<nDaughters; i++) {
    tempPdgCode = (Long64_t)((TParticle*)mcEvent->Particle(part->GetFirstDaughter()+i))->GetPdgCode();
    if (TMath::Abs(tempPdgCode) == 111 || TMath::Abs(tempPdgCode) == 113 || TMath::Abs(tempPdgCode) == 130 || TMath::Abs(tempPdgCode) == 310 || TMath::Abs(tempPdgCode) == 223 || TMath::Abs(tempPdgCode) == 221 || TMath::Abs(tempPdgCode) == 331 || TMath::Abs(tempPdgCode) == 2112 || TMath::Abs(tempPdgCode) == 3122 || TMath::Abs(tempPdgCode) == 9000111 || TMath::Abs(tempPdgCode) == 9010221)
      tempPdgCode = TMath::Abs(tempPdgCode);
    PdgDaughter->at(i) = tempPdgCode;
  }
  std::sort(PdgDaughter->begin(), PdgDaughter->end());
  
  Double_t returnVal  = -1.;
  
  switch (part->GetPdgCode()) {
    case 111:
      if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 22)
        returnVal = 1.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11 && PdgDaughter->at(2) == 22)
        returnVal = 2.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == -11 && PdgDaughter->at(2) == 11 && PdgDaughter->at(3) == 11)
        returnVal = 3.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11)
        returnVal = 4.;
      else
        returnVal = 19.;
      break;
      
    case 221:
      if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 22)
        returnVal = 1.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 111)
        returnVal = 2.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 22 && PdgDaughter->at(2) == 111)
        returnVal = 3.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 22 && PdgDaughter->at(2) == 211)
        returnVal = 4.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11 && PdgDaughter->at(2) == 22)
        returnVal = 5.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -13 && PdgDaughter->at(1) == 13 && PdgDaughter->at(2) == 22)
        returnVal = 6.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == -11 && PdgDaughter->at(2) == 11 && PdgDaughter->at(3) == 11)
        returnVal = 7.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 22 && PdgDaughter->at(2) == 22 && PdgDaughter->at(3) == 211)
        returnVal = 8.;
      else
        returnVal = 19.;
      break;
      
    case 331:
      if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 211 && PdgDaughter->at(2) == 221)
        returnVal = 1.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 113)
        returnVal = 2.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 22 && PdgDaughter->at(2) == 211)
        returnVal = 3.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 223)
        returnVal = 4.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 22)
        returnVal = 5.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -13 && PdgDaughter->at(1) == 13 && PdgDaughter->at(2) == 22)
        returnVal = 6.;
      else
        returnVal = 19.;
      break;
      
    case 223:
      if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 211)
        returnVal = 1.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 111)
        returnVal = 2.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 221)
        returnVal = 3.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11 && PdgDaughter->at(2) == 111)
        returnVal = 4.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11)
        returnVal = 5.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 111)
        returnVal = 6.;
      else
        returnVal = 19.;
      break;
      
    case 113:
      if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 211)
        returnVal = 1.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 22 && PdgDaughter->at(2) == 211)
        returnVal = 2.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 111)
        returnVal = 3.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 221)
        returnVal = 4.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 111)
        returnVal = 5.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11)
        returnVal = 6.;
      else
        returnVal = 19.;
      break;
      
    case 213:
      if (nDaughters == 2 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 211)
        returnVal = 1.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 211)
        returnVal = 2.;
      else
        returnVal = 19.;
      break;
      
    case -213:
      if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 111)
        returnVal = 1.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 22)
        returnVal = 2.;
      else
        returnVal = 19.;
      break;
      
    case 333:
      if (nDaughters == 2 && PdgDaughter->at(0) == -321 && PdgDaughter->at(1) == 321)
        returnVal = 1.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 130 && PdgDaughter->at(1) == 310)
        returnVal = 2.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 221)
        returnVal = 3.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 111)
        returnVal = 4.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11)
        returnVal = 5.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11 && PdgDaughter->at(2) == 221)
        returnVal = 6.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 22 && PdgDaughter->at(2) == 211)
        returnVal = 7.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 9010221)
        returnVal = 8.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 111)
        returnVal = 9.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11 && PdgDaughter->at(2) == 111)
        returnVal = 10.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 221)
        returnVal = 11.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 9000111)
        returnVal = 12.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 331)
        returnVal = 13.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -13 && PdgDaughter->at(1) == 13 && PdgDaughter->at(2) == 22)
        returnVal = 14.;
      else
        returnVal = 19.;
      break;
      
    case 443:
      if (nDaughters == 3 && (PdgDaughter->at(0) == 21 || PdgDaughter->at(0) == 9) && (PdgDaughter->at(1) == 21 || PdgDaughter->at(1) == 9) && (PdgDaughter->at(2) == 21 || PdgDaughter->at(2) == 9))
        returnVal = 1.;
      else if (nDaughters == 3 && (PdgDaughter->at(0) == 21 || PdgDaughter->at(0) == 9) && (PdgDaughter->at(1) == 21 || PdgDaughter->at(1) == 9) && PdgDaughter->at(2) == 22)
        returnVal = 2.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11)
        returnVal = 3.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11 && PdgDaughter->at(2) == 22)
        returnVal = 4.;
      else
        returnVal = 19.;
      break;
      
    case 1114:
      if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 2112)
        returnVal = 1.;
      else if (std::find(PdgDaughter->begin(), PdgDaughter->end(), 22) != PdgDaughter->end())
        returnVal = 2.;
      else
        returnVal = 19.;
      break;
      
    case 2114:
      if (nDaughters == 2 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 2112)
        returnVal = 1.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 2212)
        returnVal = 2.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 2112)
        returnVal = 3.;
      else
        returnVal = 19.;
      break;
      
    case 2214:
      if (nDaughters == 2 && PdgDaughter->at(0) == 211 && PdgDaughter->at(1) == 2112)
        returnVal = 1.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 2212)
        returnVal = 2.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 2212)
        returnVal = 3.;
      else
        returnVal = 19.;
      break;
      
    case 2224:
      if (nDaughters == 2 && PdgDaughter->at(0) == 211 && PdgDaughter->at(1) == 2212)
        returnVal = 1.;
      else if (std::find(PdgDaughter->begin(), PdgDaughter->end(), 22) != PdgDaughter->end())
        returnVal = 2.;
      else
        returnVal = 19.;
      break;
      
    case 3212:
      if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 3122)
        returnVal = 1.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11 && PdgDaughter->at(2) == 3122)
        returnVal = 2.;
      else
        returnVal = 19.;
      break;
      
    case 310:
      if (nDaughters == 2 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 111)
        returnVal = 1.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 211)
        returnVal = 2.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 22 && PdgDaughter->at(2) == 211)
        returnVal = 3.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 111)
        returnVal = 4.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 22)
        returnVal = 5.;
      else
        returnVal = 19.;
      break;
      
    case 130:
      if (nDaughters == 3 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 111)
        returnVal = 1.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 211)
        returnVal = 2.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -12 && PdgDaughter->at(1) == 11 && PdgDaughter->at(2) == 22 && PdgDaughter->at(3) == 211)
        returnVal = 3.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == -11 && PdgDaughter->at(2) == 12 && PdgDaughter->at(3) == 22)
        returnVal = 3.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -14 && PdgDaughter->at(1) == 13 && PdgDaughter->at(2) == 22 && PdgDaughter->at(3) == 211)
        returnVal = 4.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == -13 && PdgDaughter->at(2) == 14 && PdgDaughter->at(3) == 22)
        returnVal = 4.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 22 && PdgDaughter->at(2) == 211)
        returnVal = 5.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 22 && PdgDaughter->at(2) == 111)
        returnVal = 6.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11 && PdgDaughter->at(2) == 22 && PdgDaughter->at(3) == 111)
        returnVal = 7.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 22)
        returnVal = 8.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11 && PdgDaughter->at(2) == 22)
        returnVal = 9.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11 && PdgDaughter->at(2) == 22 && PdgDaughter->at(3) == 22)
        returnVal = 10.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -13 && PdgDaughter->at(1) == 13 && PdgDaughter->at(2) == 22)
        returnVal = 11.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -13 && PdgDaughter->at(1) == 13 && PdgDaughter->at(2) == 22 && PdgDaughter->at(3) == 22)
        returnVal = 12.;
      else
        returnVal = 19.;
      break;
      
    case 3122:
      if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 2212)
        returnVal = 1.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 2112)
        returnVal = 2.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 2112)
        returnVal = 3.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 22 && PdgDaughter->at(2) == 2212)
        returnVal = 4.;
      else
        returnVal = 19.;
      break;
      
    default:
      returnVal = -1.;
      break;
  }

  delete PdgDaughter;
  return returnVal;
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
      
    case 14:
      kc            = (AliPythia6::Instance())->Pycomp(310);
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
        if (nPart == 2 && pdgCodes[0] == 111 && pdgCodes[1] == 111)
          histo->SetBinContent(2, BR);
        else if (nPart == 2 && pdgCodes[0] == -211 && pdgCodes[1] == 211)
          histo->SetBinContent(3, BR);
        else if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == 22 && pdgCodes[2] == 211)
          histo->SetBinContent(4, BR);
        else if (nPart == 3 && pdgCodes[0] == 22 && pdgCodes[1] == 22 && pdgCodes[2] == 111)
          histo->SetBinContent(5, BR);
        else if (nPart == 2 && pdgCodes[0] == 22 && pdgCodes[1] == 22)
          histo->SetBinContent(6, BR);
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;

    case 15:
      kc            = (AliPythia6::Instance())->Pycomp(130);
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
        if (nPart == 3 && pdgCodes[0] == 111 && pdgCodes[1] == 111 && pdgCodes[2] == 111)
          histo->SetBinContent(2, BR);
        else if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == 111 && pdgCodes[2] == 211)
          histo->SetBinContent(3, BR);
        else if (nPart == 4 && pdgCodes[0] == -12 && pdgCodes[1] == 11 && pdgCodes[2] == 22 && pdgCodes[3] == 211)
          histo->SetBinContent(4, BR+histo->GetBinContent(4));
        else if (nPart == 4 && pdgCodes[0] == -211 && pdgCodes[1] == -11 && pdgCodes[2] == 12 && pdgCodes[3] == 22)
          histo->SetBinContent(4, BR+histo->GetBinContent(4));
        else if (nPart == 4 && pdgCodes[0] == -14 && pdgCodes[1] == 13 && pdgCodes[2] == 22 && pdgCodes[3] == 211)
          histo->SetBinContent(5, BR+histo->GetBinContent(5));
        else if (nPart == 4 && pdgCodes[0] == -211 && pdgCodes[1] == -13 && pdgCodes[2] == 14 && pdgCodes[3] == 22)
          histo->SetBinContent(5, BR+histo->GetBinContent(5));
        else if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == 22 && pdgCodes[2] == 211)
          histo->SetBinContent(6, BR);
        else if (nPart == 3 && pdgCodes[0] == 22 && pdgCodes[1] == 22 && pdgCodes[2] == 111)
          histo->SetBinContent(7, BR);
        else if (nPart == 4 && pdgCodes[0] == -11 && pdgCodes[1] == 11 && pdgCodes[2] == 22 && pdgCodes[3] == 111)
          histo->SetBinContent(8, BR);
        else if (nPart == 2 && pdgCodes[0] == 22 && pdgCodes[1] == 22)
          histo->SetBinContent(9, BR);
        else if (nPart == 3 && pdgCodes[0] == -11 && pdgCodes[1] == 11 && pdgCodes[2] == 22)
          histo->SetBinContent(10, BR);
        else if (nPart == 4 && pdgCodes[0] == -11 && pdgCodes[1] == 11 && pdgCodes[2] == 22 && pdgCodes[3] == 22)
          histo->SetBinContent(11, BR);
        else if (nPart == 3 && pdgCodes[0] == -13 && pdgCodes[1] == 13 && pdgCodes[2] == 22)
          histo->SetBinContent(12, BR);
        else if (nPart == 4 && pdgCodes[0] == -13 && pdgCodes[1] == 13 && pdgCodes[2] == 22 && pdgCodes[3] == 22)
          histo->SetBinContent(13, BR);
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;

    case 16:
      kc            = (AliPythia6::Instance())->Pycomp(3122);
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
        if (nPart == 2 && pdgCodes[0] == -211 && pdgCodes[1] == 2212)
          histo->SetBinContent(2, BR);
        else if (nPart == 2 && pdgCodes[0] == 111 && pdgCodes[1] == 2112)
          histo->SetBinContent(3, BR);
        else if (nPart == 2 && pdgCodes[0] == 22 && pdgCodes[1] == 2112)
          histo->SetBinContent(4, BR);
        else if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == 22 && pdgCodes[2] == 2212)
          histo->SetBinContent(4, BR);
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
Int_t AliAnalysisTaskGammaCocktailMC::GetParticlePosLocal(Int_t pdg) {

  Int_t returnVal = -9999;
  
  switch (pdg) {
    case 111:
      returnVal = 0;
      break;
    case 221:
      returnVal = 1;
      break;
    case 331:
      returnVal = 2;
      break;
    case 223:
      returnVal = 3;
      break;
    case 113:
      returnVal = 4;
      break;
    case 213:
      returnVal = 5;
      break;
    case -213:
      returnVal = 6;
      break;
    case 333:
      returnVal = 7;
      break;
    case 443:
      returnVal = 8;
      break;
    case 1114:
      returnVal = 9;
      break;
    case 2114:
      returnVal = 10;
      break;
    case 2214:
      returnVal = 11;
      break;
    case 2224:
      returnVal = 12;
      break;
    case 3212:
      returnVal = 13;
      break;
    case 310:
      returnVal = 14;
      break;
    case 130:
      returnVal = 15;
      break;
    case 3122:
      returnVal = 16;
      break;
    default:
      break;
  }
  
  return returnVal;
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
