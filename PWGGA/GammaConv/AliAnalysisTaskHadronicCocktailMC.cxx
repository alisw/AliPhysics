/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.  *
 *                                                                         *
 * Author: Friederike Bock, Lucas Altenkämper                              *
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
// Class used to do analysis on hadronic cocktail output
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
#include "AliAnalysisTaskHadronicCocktailMC.h"
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

ClassImp(AliAnalysisTaskHadronicCocktailMC)

//________________________________________________________________________
AliAnalysisTaskHadronicCocktailMC::AliAnalysisTaskHadronicCocktailMC(): AliAnalysisTaskSE(),
  fOutputContainer(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fMCGenHandler(NULL),
  fMCGenerator(NULL),
  fMCCocktailGen(NULL),
  fAnalyzedMeson(0),
  fAnalyzeNeutralPi(kFALSE),
  fAnalyzeChargedPi(kFALSE),
  fDoLightOutput(kFALSE),
  fHasMother{kFALSE},
  fHistNEvents(NULL),
  fHistPtPhiDaughterSource(NULL),
  fHistPtPhiInput(NULL),
  fHistPtYInput(NULL),
  fHistPtYDaughterSource(NULL),
  fHistDecayChannelsInput(NULL),
  fHistPythiaBR(NULL),
  fHistPtDaughterPtSourceInput(NULL),
  fHistPhiDaughterPhiSourceInput(NULL),
  fHistPdgInputRest(NULL),
  fHistPdgDaughterSourceRest(NULL),
  fHistPtYGammaFromXFromInput(NULL),
  fHistPtPhiGammaFromXFromInput(NULL),
  fHistPtYGammaFromPi0FromInput(NULL),
  fHistPtPhiGammaFromPi0FromInput(NULL),
  fParticleList(NULL),
  fParticleListNames(NULL),
  fPtParametrization{NULL},
  fPtParametrizationProton(NULL),
  fPtParametrizationPi0(NULL),
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
AliAnalysisTaskHadronicCocktailMC::AliAnalysisTaskHadronicCocktailMC(const char *name):
  AliAnalysisTaskSE(name),
  fOutputContainer(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fMCGenHandler(NULL),
  fMCGenerator(NULL),
  fMCCocktailGen(NULL),
  fAnalyzedMeson(0),
  fAnalyzeNeutralPi(kFALSE),
  fAnalyzeChargedPi(kFALSE),
  fDoLightOutput(kFALSE),
  fHasMother{kFALSE},
  fHistNEvents(NULL),
  fHistPtPhiDaughterSource(NULL),
  fHistPtPhiInput(NULL),
  fHistPtYInput(NULL),
  fHistPtYDaughterSource(NULL),
  fHistDecayChannelsInput(NULL),
  fHistPythiaBR(NULL),
  fHistPtDaughterPtSourceInput(NULL),
  fHistPhiDaughterPhiSourceInput(NULL),
  fHistPdgInputRest(NULL),
  fHistPdgDaughterSourceRest(NULL),
  fHistPtYGammaFromXFromInput(NULL),
  fHistPtPhiGammaFromXFromInput(NULL),
  fHistPtYGammaFromPi0FromInput(NULL),
  fHistPtPhiGammaFromPi0FromInput(NULL),
  fParticleList(NULL),
  fParticleListNames(NULL),
  fPtParametrization{NULL},
  fPtParametrizationProton(NULL),
  fPtParametrizationPi0(NULL),
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

AliAnalysisTaskHadronicCocktailMC::~AliAnalysisTaskHadronicCocktailMC()
{
  for (Int_t i=0; i<12; i++) {
    if (fCocktailSettings[i]) delete fCocktailSettings[i];
  }
}

//________________________________________________________________________
void AliAnalysisTaskHadronicCocktailMC::UserCreateOutputObjects(){
  
  // Create histograms
  if(fOutputContainer != NULL){
    delete fOutputContainer;
    fOutputContainer          = NULL;
  }
  if(fOutputContainer == NULL){
    fOutputContainer          = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }
  
  TString                     fAnalyzedParticle = "";
  if (fAnalyzedMeson==0)      fAnalyzedParticle = "Pi0";
  else if (fAnalyzedMeson==1) fAnalyzedParticle = "Eta";
  else if (fAnalyzedMeson==2) fAnalyzedParticle = "Pi+-";
  else AliFatalClass("Analyzed particle not recognized!");

  // tree + user info list to protect contents from merging
  fOutputTree = new TTree("cocktailSettings", "cocktailSettings");
  fUserInfo   = (TList*)fOutputTree->GetUserInfo();
  
  fMCGenHandler                           = (AliMCGenHandler*)AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler();
  fMCGenerator                            = fMCGenHandler->GetGenerator();
  TString mcGeneratorClassName            = "";
  if (fMCGenerator)  mcGeneratorClassName = fMCGenerator->ClassName();
  
  if (mcGeneratorClassName.CompareTo("AliGenEMCocktailV2") == 0) {
    
    fMCCocktailGen = (AliGenEMCocktailV2*)fMCGenerator;
    
    // has mother i
    SetHasMother((UInt_t)fMCCocktailGen->GetSelectedMothers());
    
    // pt parametrizations
    GetAndSetPtParametrizations(fMCCocktailGen);
    for (Int_t i=0; i<24; i++) {
      if (fHasMother[i]) fUserInfo->Add(fPtParametrization[i]);
    }
    if (fPtParametrizationProton) fUserInfo->Add(fPtParametrizationProton);
    if (fPtParametrizationPi0)    fUserInfo->Add(fPtParametrizationPi0);
    
    // cocktail settings
    Double_t ptMin, ptMax;
    fMCCocktailGen->GetPtRange(ptMin, ptMax);
    fCocktailSettings[0] = new TObjString(Form("collSys_%d",          fMCCocktailGen->GetCollisionSystem()));
    fCocktailSettings[1] = new TObjString(Form("cent_%d",             fMCCocktailGen->GetCentrality()));
    fCocktailSettings[2] = new TObjString(Form("decayMode_%.0f",      fMCCocktailGen->GetDecayMode()));
    fCocktailSettings[3] = new TObjString(Form("selectMothers_%d",    fMCCocktailGen->GetSelectedMothers()));
    fCocktailSettings[4] = new TObjString(Form("paramFile_%s",        (fMCCocktailGen->GetParametrizationFile()).Data()));
    fCocktailSettings[5] = new TObjString(Form("paramDir_%s",         (fMCCocktailGen->GetParametrizationFileDirectory()).Data()));
    fCocktailSettings[6] = new TObjString(Form("nParticles_%d",       fMCCocktailGen->GetNumberOfParticles()));
    fCocktailSettings[7] = new TObjString(Form("ptMin_%.2f",          ptMin));
    fCocktailSettings[8] = new TObjString(Form("ptMax_%.2f",          ptMax));
    fCocktailSettings[9] = new TObjString(Form("weightMode_%.0f",     fMCCocktailGen->GetWeightingMode()));
    fCocktailSettings[10] = new TObjString(Form("dynamicalPtRange_%d",fMCCocktailGen->GetDynamicalPtRangeOption()));
    fCocktailSettings[11] = new TObjString(Form("yWeights_%d",        fMCCocktailGen->GetYWeightOption()));
    for (Int_t i=0; i<12; i++) fUserInfo->Add(fCocktailSettings[i]);
    
    // mt scaling params
    fMtScalingFactors = (TH1D*)fMCCocktailGen->GetMtScalingFactors();
    fUserInfo->Add(fMtScalingFactors);

    // pt-y distributions
    GetAndSetPtYDistributions(fMCCocktailGen);
    for (Int_t i=0; i<24; i++) {
      if (fHasMother[i]) fUserInfo->Add(fPtYDistributions[i]);
    }
  } else {
    for (Int_t i=0; i<24; i++) fHasMother[i] = kTRUE;
  }
  
  fHistNEvents = (TH1F*)SetHist1D(fHistNEvents,"f","NEvents","","N_{evt}",1,0,1,kTRUE);
  fOutputContainer->Add(fHistNEvents);
  
  const Int_t nInputParticles         = 24;
  Int_t   fParticleList_local[]       = {221,310,130,3122,113,331,223,213,-213,333,443,2114,2214,1114,2224,321,-321,-3334,3334,-3312,3312,3224,3114,313};
  TString fParticleListNames_local[]  = {"Eta","K0s","K0l","Lambda","rho0","EtaPrim","omega","rho+","rho-","phi","J/psi","Delta0","Delta+","Delta-","Delta++",
                                          "K+","K-","Omega+","Omega-","Xi+","Xi-","Sigma(1385)+", "Sigma(1385)-","K*(892)0"};
  
  // pi0/eta/pi+- from X
  fParticleList                   = fParticleList_local;
  fParticleListNames              = fParticleListNames_local;
  fHistPtYInput                   = new TH2F*[nInputParticles];
  fHistPtYDaughterSource          = new TH2F*[nInputParticles];
  fHistDecayChannelsInput         = new TH1F*[nInputParticles];
  fHistPythiaBR                   = new TH1F*[nInputParticles];
  fHistPtPhiDaughterSource        = new TH2F*[nInputParticles];
  fHistPtPhiInput                 = new TH2F*[nInputParticles];
  fHistPtDaughterPtSourceInput    = new TH2F*[nInputParticles];
  fHistPhiDaughterPhiSourceInput  = new TH2F*[nInputParticles];
  for(Int_t i=0; i<nInputParticles; i++){
    if (fHasMother[i]) {
      fHistPtYInput[i] = (TH2F*)SetHist2D(fHistPtYInput[i],"f",Form("Pt_Y_%s",fParticleListNames[i].Data()),"#it{p}_{T}","Y",500,0,50,400,-2.0,2.0,kTRUE);
      fOutputContainer->Add(fHistPtYInput[i]);
      
      // pi0/eta/pi+- from certain mother
      fHistPtYDaughterSource[i] = (TH2F*)SetHist2D(fHistPtYDaughterSource[i],"f",Form("Pt_Y_%s_From_%s",fAnalyzedParticle.Data(),fParticleListNames[i].Data()),"#it{p}_{T}","Y",500,0,50,400,-2.0,2.0,kTRUE);
      fOutputContainer->Add(fHistPtYDaughterSource[i]);
      
      // phi distributions
      fHistPtPhiInput[i] = (TH2F*)SetHist2D(fHistPtPhiInput[i],"f",Form("Pt_Phi_%s",fParticleListNames[i].Data()),"#it{p}_{T}","#phi",500,0,50,100,0,7,kTRUE);
      fOutputContainer->Add(fHistPtPhiInput[i]);
      
      fHistPtPhiDaughterSource[i] = (TH2F*)SetHist2D(fHistPtPhiDaughterSource[i],"f",Form("Pt_Phi_%s_From_%s",fAnalyzedParticle.Data(),fParticleListNames[i].Data()),"#it{p}_{T}","#phi",500,0,50,100,0,7,kTRUE);
      fOutputContainer->Add(fHistPtPhiDaughterSource[i]);
      
      // correlation gamma from certain mother to mother
      fHistPtDaughterPtSourceInput[i] = (TH2F*)SetHist2D(fHistPtDaughterPtSourceInput[i],"f",Form("Pt%s_PtMother_%s",fAnalyzedParticle.Data(),fParticleListNames[i].Data()),"#it{p}_{T,daughter}","#it{p}_{T,mother}",500,0,50,500,0,50,kTRUE);
      fOutputContainer->Add(fHistPtDaughterPtSourceInput[i]);
      
      fHistPhiDaughterPhiSourceInput[i] = (TH2F*)SetHist2D(fHistPhiDaughterPhiSourceInput[i],"f",Form("Phi%s_PhiMother_%s",fAnalyzedParticle.Data(),fParticleListNames[i].Data()),"#phi_{daughter}","#phi_{mother}",100,0,7,100,0,7,kTRUE);
      fOutputContainer->Add(fHistPhiDaughterPhiSourceInput[i]);
      
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
      fHistPtYInput[i]                  = NULL;
      fHistPtYDaughterSource[i]         = NULL;
      fHistPtPhiInput[i]                = NULL;
      fHistPtPhiDaughterSource[i]       = NULL;
      fHistPtDaughterPtSourceInput[i]   = NULL;
      fHistPhiDaughterPhiSourceInput[i] = NULL;
      fHistDecayChannelsInput[i]        = NULL;
      fHistPythiaBR[i]                  = NULL;
    }
  }
  
  // gamma from X (or pi0) from X
  fHistPtYGammaFromXFromInput     = new TH2F*[3];
  fHistPtPhiGammaFromXFromInput   = new TH2F*[3];
  fHistPtYGammaFromPi0FromInput   = new TH2F*[3];
  fHistPtPhiGammaFromPi0FromInput = new TH2F*[3];
  for (Int_t i = 0; i<3; i++) {
    if (fHasMother[i+1]) {
      
      fHistPtYGammaFromXFromInput[i] = (TH2F*)SetHist2D(fHistPtYGammaFromXFromInput[i],"f",Form("Pt_Y_Gamma_From_X_From_%s",fParticleListNames[i+1].Data()),"#it{p}_{T}","Y",500,0,50,400,-2.0,2.0,kTRUE);
      fOutputContainer->Add(fHistPtYGammaFromXFromInput[i]);

      fHistPtYGammaFromPi0FromInput[i] = (TH2F*)SetHist2D(fHistPtYGammaFromPi0FromInput[i],"f",Form("Pt_Y_Gamma_From_Pi0_From_%s",fParticleListNames[i+1].Data()),"#it{p}_{T}","Y",500,0,50,400,-2.0,2.0,kTRUE);
      fOutputContainer->Add(fHistPtYGammaFromPi0FromInput[i]);

      fHistPtPhiGammaFromXFromInput[i] = (TH2F*)SetHist2D(fHistPtPhiGammaFromXFromInput[i],"f",Form("Pt_Phi_Gamma_From_X_From_%s",fParticleListNames[i+1].Data()),"#it{p}_{T}","#phi",500,0,50,100,0,7,kTRUE);
      fOutputContainer->Add(fHistPtPhiGammaFromXFromInput[i]);

      fHistPtPhiGammaFromPi0FromInput[i] = (TH2F*)SetHist2D(fHistPtPhiGammaFromPi0FromInput[i],"f",Form("Pt_Phi_Gamma_From_Pi0_From_%s",fParticleListNames[i+1].Data()),"#it{p}_{T}","#phi",500,0,50,100,0,7,kTRUE);
      fOutputContainer->Add(fHistPtPhiGammaFromPi0FromInput[i]);
      
    } else {
      fHistPtYGammaFromXFromInput[i]      = NULL;
      fHistPtPhiGammaFromXFromInput[i]    = NULL;
      fHistPtYGammaFromPi0FromInput[i]    = NULL;
      fHistPtPhiGammaFromPi0FromInput[i]  = NULL;
    }
  }
  
  fHistPdgInputRest = (TH1I*)SetHist1D(fHistPdgInputRest,"f","Pdg_primary_rest","PDG code","",5000,0,5000,kTRUE);
  fOutputContainer->Add(fHistPdgInputRest);
  
  fHistPdgDaughterSourceRest = (TH1I*)SetHist1D(fHistPdgDaughterSourceRest,"f",Form("Pdg_%s_From_rest",fAnalyzedParticle.Data()),"PDG code mother","",5000,0,5000,kTRUE);
  fOutputContainer->Add(fHistPdgDaughterSourceRest);
  
  fOutputContainer->Add(fOutputTree);
  
  PostData(1, fOutputContainer);
}

//_____________________________________________________________________________
void AliAnalysisTaskHadronicCocktailMC::UserExec(Option_t *)
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
void AliAnalysisTaskHadronicCocktailMC::GetAndSetPtParametrizations(AliGenEMCocktailV2* fMCCocktailGen)
{
  if (!fMCCocktailGen) return;
  
  for (Int_t i=0; i<24; i++) fPtParametrization[i] = NULL;
  fPtParametrizationProton = NULL;
  fPtParametrizationPi0    = NULL;
  
  TF1* fct        = NULL;
  TString fctName = "";
  for (Int_t i=0; i<27; i++) {
    fct = (TF1*)fMCCocktailGen->GetPtParametrization(i);
    if (fct) {
      fctName = fct->GetName();
      if (fctName.BeginsWith("221_pt")    && fHasMother[0])  fPtParametrization[0]   = fct;
      if (fctName.BeginsWith("310_pt")    && fHasMother[1])  fPtParametrization[1]   = fct;
      if (fctName.BeginsWith("130_pt")    && fHasMother[2])  fPtParametrization[2]   = fct;
      if (fctName.BeginsWith("3122_pt")   && fHasMother[3])  fPtParametrization[3]   = fct;
      if (fctName.BeginsWith("113_pt")    && fHasMother[4])  fPtParametrization[4]   = fct;
      if (fctName.BeginsWith("331_pt")    && fHasMother[5])  fPtParametrization[5]   = fct;
      if (fctName.BeginsWith("223_pt")    && fHasMother[6])  fPtParametrization[6]   = fct;
      if (fctName.BeginsWith("213_pt")    && fHasMother[7])  fPtParametrization[7]   = fct;
      if (fctName.BeginsWith("-213_pt")   && fHasMother[8])  fPtParametrization[8]   = fct;
      if (fctName.BeginsWith("333_pt")    && fHasMother[9])  fPtParametrization[9]   = fct;
      if (fctName.BeginsWith("443_pt")    && fHasMother[10]) fPtParametrization[10]  = fct;
      if (fctName.BeginsWith("2114_pt")   && fHasMother[11]) fPtParametrization[11]  = fct;
      if (fctName.BeginsWith("2214_pt")   && fHasMother[12]) fPtParametrization[12]  = fct;
      if (fctName.BeginsWith("1114_pt")   && fHasMother[13]) fPtParametrization[13]  = fct;
      if (fctName.BeginsWith("2224_pt")   && fHasMother[14]) fPtParametrization[14]  = fct;
      if (fctName.BeginsWith("321_pt")    && fHasMother[15]) fPtParametrization[15]  = fct;
      if (fctName.BeginsWith("-321_pt")   && fHasMother[16]) fPtParametrization[16]  = fct;
      if (fctName.BeginsWith("-3334_pt")  && fHasMother[17]) fPtParametrization[17]  = fct;
      if (fctName.BeginsWith("3334_pt")   && fHasMother[18]) fPtParametrization[18]  = fct;
      if (fctName.BeginsWith("-3312_pt")  && fHasMother[19]) fPtParametrization[19]  = fct;
      if (fctName.BeginsWith("3312_pt")   && fHasMother[20]) fPtParametrization[20]  = fct;
      if (fctName.BeginsWith("3224_pt")   && fHasMother[21]) fPtParametrization[21]  = fct;
      if (fctName.BeginsWith("3114_pt")   && fHasMother[22]) fPtParametrization[22]  = fct;
      if (fctName.BeginsWith("313_pt")    && fHasMother[23]) fPtParametrization[23]  = fct;
      if (fctName.BeginsWith("2212_pt")) fPtParametrizationProton = fct;
      if (fctName.BeginsWith("111_pt"))  fPtParametrizationPi0    = fct;
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskHadronicCocktailMC::GetAndSetPtYDistributions(AliGenEMCocktailV2* fMCCocktailGen)
{
  if (!fMCCocktailGen) return;

  for (Int_t i=0; i<24; i++) fPtYDistributions[i] = NULL;

  TH2F* tempPtY = NULL;
  TString tempPtYName = "";
  for (Int_t i=0; i<26; i++) {
    tempPtY = (TH2F*)fMCCocktailGen->GetPtYDistribution(i);
    if (tempPtY) {
      tempPtYName = tempPtY->GetName();
      if (tempPtYName.BeginsWith("221_pt_y")    && fHasMother[0])  fPtYDistributions[0]   = tempPtY;
      if (tempPtYName.BeginsWith("310_pt_y")    && fHasMother[1])  fPtYDistributions[1]   = tempPtY;
      if (tempPtYName.BeginsWith("130_pt_y")    && fHasMother[2])  fPtYDistributions[2]   = tempPtY;
      if (tempPtYName.BeginsWith("3122_pt_y")   && fHasMother[3])  fPtYDistributions[3]   = tempPtY;
      if (tempPtYName.BeginsWith("113_pt_y")    && fHasMother[4])  fPtYDistributions[4]   = tempPtY;
      if (tempPtYName.BeginsWith("331_pt_y")    && fHasMother[5])  fPtYDistributions[5]   = tempPtY;
      if (tempPtYName.BeginsWith("223_pt_y")    && fHasMother[6])  fPtYDistributions[6]   = tempPtY;
      if (tempPtYName.BeginsWith("213_pt_y")    && fHasMother[7])  fPtYDistributions[7]   = tempPtY;
      if (tempPtYName.BeginsWith("-213_pt_y")   && fHasMother[8])  fPtYDistributions[8]   = tempPtY;
      if (tempPtYName.BeginsWith("333_pt_y")    && fHasMother[9])  fPtYDistributions[9]   = tempPtY;
      if (tempPtYName.BeginsWith("443_pt_y")    && fHasMother[10]) fPtYDistributions[10]  = tempPtY;
      if (tempPtYName.BeginsWith("2114_pt_y")   && fHasMother[11]) fPtYDistributions[11]  = tempPtY;
      if (tempPtYName.BeginsWith("2214_pt_y")   && fHasMother[12]) fPtYDistributions[12]  = tempPtY;
      if (tempPtYName.BeginsWith("1114_pt_y")   && fHasMother[13]) fPtYDistributions[13]  = tempPtY;
      if (tempPtYName.BeginsWith("2224_pt_y")   && fHasMother[14]) fPtYDistributions[14]  = tempPtY;
      if (tempPtYName.BeginsWith("321_pt_y")    && fHasMother[15]) fPtYDistributions[15]  = tempPtY;
      if (tempPtYName.BeginsWith("-321_pt_y")   && fHasMother[16]) fPtYDistributions[16]  = tempPtY;
      if (tempPtYName.BeginsWith("-3334_pt_y")  && fHasMother[17]) fPtYDistributions[17]  = tempPtY;
      if (tempPtYName.BeginsWith("3334_pt_y")   && fHasMother[18]) fPtYDistributions[18]  = tempPtY;
      if (tempPtYName.BeginsWith("-3312_pt_y")  && fHasMother[19]) fPtYDistributions[19]  = tempPtY;
      if (tempPtYName.BeginsWith("3312_pt_y")   && fHasMother[20]) fPtYDistributions[20]  = tempPtY;
      if (tempPtYName.BeginsWith("3224_pt_y")   && fHasMother[21]) fPtYDistributions[21]  = tempPtY;
      if (tempPtYName.BeginsWith("3114_pt_y")   && fHasMother[22]) fPtYDistributions[22]  = tempPtY;
      if (tempPtYName.BeginsWith("313_pt_y")    && fHasMother[23]) fPtYDistributions[23]  = tempPtY;
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskHadronicCocktailMC::SetHasMother(UInt_t selectedMothers) {
  
  for (Int_t i=0; i<24; i++) fHasMother[i] = kFALSE;
  
  // selects mother particles according to choice and possible decays (i.e. into pi/pi+-/eta)
  if (                   (selectedMothers&AliGenEMCocktailV2::kGenEta)      && (fAnalyzeNeutralPi || fAnalyzeChargedPi))    fHasMother[0] = kTRUE;
  if (                   (selectedMothers&AliGenEMCocktailV2::kGenK0s)      && (fAnalyzeNeutralPi || fAnalyzeChargedPi))    fHasMother[1] = kTRUE;
  if (                   (selectedMothers&AliGenEMCocktailV2::kGenK0l)      && (fAnalyzeNeutralPi || fAnalyzeChargedPi))    fHasMother[2] = kTRUE;
  if (                   (selectedMothers&AliGenEMCocktailV2::kGenLambda)   && (fAnalyzeNeutralPi || fAnalyzeChargedPi))    fHasMother[3] = kTRUE;
  if (                    selectedMothers&AliGenEMCocktailV2::kGenRho0)                                                     fHasMother[4] = kTRUE;
  if (!fDoLightOutput && (selectedMothers&AliGenEMCocktailV2::kGenEtaprime))                                                fHasMother[5] = kTRUE;
  if (!fDoLightOutput && (selectedMothers&AliGenEMCocktailV2::kGenOmega))                                                   fHasMother[6] = kTRUE;
  if (!fDoLightOutput && (selectedMothers&AliGenEMCocktailV2::kGenRhoPl))                                                   fHasMother[7] = kTRUE;
  if (!fDoLightOutput && (selectedMothers&AliGenEMCocktailV2::kGenRhoMi))                                                   fHasMother[8] = kTRUE;
  if (!fDoLightOutput && (selectedMothers&AliGenEMCocktailV2::kGenPhi))                                                     fHasMother[9] = kTRUE;
  if (!fDoLightOutput && (selectedMothers&AliGenEMCocktailV2::kGenJpsi))                                                    fHasMother[10] = kTRUE;
  if (!fDoLightOutput && (selectedMothers&AliGenEMCocktailV2::kGenDeltaZero)  && (fAnalyzeNeutralPi || fAnalyzeChargedPi))  fHasMother[11] = kTRUE;
  if (!fDoLightOutput && (selectedMothers&AliGenEMCocktailV2::kGenDeltaPl)    && (fAnalyzeNeutralPi || fAnalyzeChargedPi))  fHasMother[12] = kTRUE;
  if (!fDoLightOutput && (selectedMothers&AliGenEMCocktailV2::kGenDeltaMi)    && fAnalyzeChargedPi)                         fHasMother[13] = kTRUE;
  if (!fDoLightOutput && (selectedMothers&AliGenEMCocktailV2::kGenDeltaPlPl)  && fAnalyzeChargedPi)                         fHasMother[14] = kTRUE;
  if (!fDoLightOutput && (selectedMothers&AliGenEMCocktailV2::kGenKPl)        && (fAnalyzeNeutralPi || fAnalyzeChargedPi))  fHasMother[15] = kTRUE;
  if (!fDoLightOutput && (selectedMothers&AliGenEMCocktailV2::kGenKMi)        && (fAnalyzeNeutralPi || fAnalyzeChargedPi))  fHasMother[16] = kTRUE;
  if (!fDoLightOutput && (selectedMothers&AliGenEMCocktailV2::kGenOmegaPl)    && (fAnalyzeNeutralPi || fAnalyzeChargedPi))  fHasMother[17] = kTRUE;
  if (!fDoLightOutput && (selectedMothers&AliGenEMCocktailV2::kGenOmegaMi)    && (fAnalyzeNeutralPi || fAnalyzeChargedPi))  fHasMother[18] = kTRUE;
  if (!fDoLightOutput && (selectedMothers&AliGenEMCocktailV2::kGenXiPl)       && fAnalyzeChargedPi)                         fHasMother[19] = kTRUE;
  if (!fDoLightOutput && (selectedMothers&AliGenEMCocktailV2::kGenXiMi)       && fAnalyzeChargedPi)                         fHasMother[20] = kTRUE;
  if (!fDoLightOutput && (selectedMothers&AliGenEMCocktailV2::kGenSigmaPl)    && fAnalyzeChargedPi)                         fHasMother[21] = kTRUE;
  if (!fDoLightOutput && (selectedMothers&AliGenEMCocktailV2::kGenSigmaMi)    && fAnalyzeChargedPi)                         fHasMother[22] = kTRUE;
  if (!fDoLightOutput && (selectedMothers&AliGenEMCocktailV2::kGenK0star)     && (fAnalyzeNeutralPi || fAnalyzeChargedPi))  fHasMother[23] = kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskHadronicCocktailMC::ProcessMCParticles(){
  
  // Loop over all primary MC particle
  for(Long_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {
    // fill primary histograms
    TParticle* particle         = NULL;
    particle                    = (TParticle *)fMCEvent->Particle(i);
    if (!particle) continue;
    Bool_t hasMother            = kFALSE;
    Bool_t particleIsPrimary    = kTRUE;
    
    if (particle->GetMother(0)>-1){
      hasMother         = kTRUE;
      particleIsPrimary = kFALSE;
    }
    TParticle*      motherParticle  = NULL;
    if (hasMother)  motherParticle  = (TParticle*)fMCEvent->Particle(particle->GetMother(0));
    if (motherParticle) hasMother   = kTRUE;
    else                hasMother   = kFALSE;
    
    Bool_t motherIsPrimary                                = kFALSE;
    if(hasMother){
      if(motherParticle->GetMother(0)>-1) motherIsPrimary = kFALSE;
      else                                motherIsPrimary = kTRUE;
    }
    
    TParticle* grandMotherParticle  = NULL;
    Bool_t motherHasMother          = kFALSE;
    if (hasMother && !motherIsPrimary) {
      grandMotherParticle           = (TParticle*)fMCEvent->Particle(motherParticle->GetMother(0));
      motherHasMother               = kTRUE;
    }

    Bool_t grandMotherIsPrimary                                      = kFALSE;
    if (motherHasMother) {
      if(grandMotherParticle->GetMother(0)>-1)  grandMotherIsPrimary = kFALSE;
      else                                      grandMotherIsPrimary = kTRUE;
    }

    if (!(TMath::Abs(particle->Energy()-particle->Pz())>0.)) continue;
    Double_t yPre = (particle->Energy()+particle->Pz())/(particle->Energy()-particle->Pz());
    if (yPre == 0.) continue;
    
    Double_t y = 0.5*TMath::Log(yPre);
    if (TMath::Abs(y) > fMaxY) continue;
    
    Int_t                       PdgAnalyzedParticle = 0;
    if (fAnalyzedMeson==0)      PdgAnalyzedParticle = 111;
    else if (fAnalyzedMeson==1) PdgAnalyzedParticle = 221;
    else if (fAnalyzedMeson==2) PdgAnalyzedParticle = 211; // set to pi+ pdg code, but both are accepted, no distinction between pi+ and pi-
    
    // pi0/eta/pi+- from source
    if(TMath::Abs(particle->GetPdgCode())==PdgAnalyzedParticle && hasMother==kTRUE){
      if(motherIsPrimary && fHasMother[GetParticlePosLocal(motherParticle->GetPdgCode())]){

        switch(motherParticle->GetPdgCode()){
          case 221:
            fHistPtYDaughterSource[0]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiDaughterSource[0]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            fHistPtDaughterPtSourceInput[0]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
            fHistPhiDaughterPhiSourceInput[0]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
            break;
          case 310:
            fHistPtYDaughterSource[1]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiDaughterSource[1]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            fHistPtDaughterPtSourceInput[1]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
            fHistPhiDaughterPhiSourceInput[1]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
            break;
          case 130:
            fHistPtYDaughterSource[2]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiDaughterSource[2]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            fHistPtDaughterPtSourceInput[2]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
            fHistPhiDaughterPhiSourceInput[2]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
            break;
          case 3122:
            fHistPtYDaughterSource[3]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiDaughterSource[3]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            fHistPtDaughterPtSourceInput[3]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
            fHistPhiDaughterPhiSourceInput[3]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
            break;
          case 113:
            fHistPtYDaughterSource[4]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiDaughterSource[4]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            fHistPtDaughterPtSourceInput[4]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
            fHistPhiDaughterPhiSourceInput[4]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
            break;
          case 331:
            fHistPtYDaughterSource[5]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiDaughterSource[5]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            fHistPtDaughterPtSourceInput[5]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
            fHistPhiDaughterPhiSourceInput[5]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
            break;
          case 223:
            fHistPtYDaughterSource[6]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiDaughterSource[6]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            fHistPtDaughterPtSourceInput[6]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
            fHistPhiDaughterPhiSourceInput[6]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
            break;
          case 213:
            fHistPtYDaughterSource[7]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiDaughterSource[7]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            fHistPtDaughterPtSourceInput[7]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
            fHistPhiDaughterPhiSourceInput[7]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
            break;
          case -213:
            fHistPtYDaughterSource[8]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiDaughterSource[8]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            fHistPtDaughterPtSourceInput[8]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
            fHistPhiDaughterPhiSourceInput[8]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
            break;
          case 333:
            fHistPtYDaughterSource[9]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiDaughterSource[9]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            fHistPtDaughterPtSourceInput[9]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
            fHistPhiDaughterPhiSourceInput[9]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
            break;
          case 443:
            fHistPtYDaughterSource[10]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiDaughterSource[10]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            fHistPtDaughterPtSourceInput[10]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
            fHistPhiDaughterPhiSourceInput[10]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
            break;
          case 2114:
            fHistPtYDaughterSource[11]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiDaughterSource[11]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            fHistPtDaughterPtSourceInput[11]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
            fHistPhiDaughterPhiSourceInput[11]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
            break;
          case 2214:
            fHistPtYDaughterSource[12]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiDaughterSource[12]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            fHistPtDaughterPtSourceInput[12]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
            fHistPhiDaughterPhiSourceInput[12]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
            break;
          case 1114:
            fHistPtYDaughterSource[13]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiDaughterSource[13]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            fHistPtDaughterPtSourceInput[13]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
            fHistPhiDaughterPhiSourceInput[13]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
            break;
          case 2224:
            fHistPtYDaughterSource[14]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiDaughterSource[14]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            fHistPtDaughterPtSourceInput[14]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
            fHistPhiDaughterPhiSourceInput[14]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
            break;
          case 321:
            fHistPtYDaughterSource[15]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiDaughterSource[15]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            fHistPtDaughterPtSourceInput[15]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
            fHistPhiDaughterPhiSourceInput[15]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
            break;
          case -321:
            fHistPtYDaughterSource[16]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiDaughterSource[16]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            fHistPtDaughterPtSourceInput[16]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
            fHistPhiDaughterPhiSourceInput[16]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
            break;
          case -3334:
            fHistPtYDaughterSource[17]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiDaughterSource[17]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            fHistPtDaughterPtSourceInput[17]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
            fHistPhiDaughterPhiSourceInput[17]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
            break;
          case 3334:
            fHistPtYDaughterSource[18]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiDaughterSource[18]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            fHistPtDaughterPtSourceInput[18]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
            fHistPhiDaughterPhiSourceInput[18]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
            break;
          case -3312:
            fHistPtYDaughterSource[19]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiDaughterSource[19]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            fHistPtDaughterPtSourceInput[19]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
            fHistPhiDaughterPhiSourceInput[19]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
            break;
          case 3312:
            fHistPtYDaughterSource[20]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiDaughterSource[20]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            fHistPtDaughterPtSourceInput[20]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
            fHistPhiDaughterPhiSourceInput[20]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
            break;
          case 3224:
            fHistPtYDaughterSource[21]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiDaughterSource[21]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            fHistPtDaughterPtSourceInput[21]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
            fHistPhiDaughterPhiSourceInput[21]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
            break;
          case 3114:
            fHistPtYDaughterSource[22]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiDaughterSource[22]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            fHistPtDaughterPtSourceInput[22]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
            fHistPhiDaughterPhiSourceInput[22]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
            break;
          case 313:
            fHistPtYDaughterSource[23]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiDaughterSource[23]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            fHistPtDaughterPtSourceInput[23]->Fill(particle->Pt(), motherParticle->Pt(), particle->GetWeight());
            fHistPhiDaughterPhiSourceInput[23]->Fill(particle->Phi(), motherParticle->Phi(), particle->GetWeight());
            break;
          default:
            fHistPdgDaughterSourceRest->Fill(motherParticle->GetPdgCode());
            break;
        }
      }
    }

    // source
    if(particle->GetPdgCode()!=PdgAnalyzedParticle && particleIsPrimary && fHasMother[GetParticlePosLocal(particle->GetPdgCode())]){
      
      switch(particle->GetPdgCode()){
        case 221:
          fHistPtYInput[0]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[0]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[0]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[0]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          break;
        case 310:
          fHistPtYInput[1]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[1]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[1]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[1]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          break;
        case 130:
          fHistPtYInput[2]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[2]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[2]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[2]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          break;
        case 3122:
          fHistPtYInput[3]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[3]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[3]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[3]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          break;
        case 113:
          fHistPtYInput[4]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[4]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[4]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[4]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          break;
        case 331:
          fHistPtYInput[5]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[5]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[5]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[5]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          break;
        case 223:
          fHistPtYInput[6]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[6]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[6]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[6]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          break;
        case 213:
          fHistPtYInput[7]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[7]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[7]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[7]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          break;
        case -213:
          fHistPtYInput[8]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[8]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[8]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[8]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          break;
        case 333:
          fHistPtYInput[9]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[9]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[9]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[9]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          break;
        case 443:
          fHistPtYInput[10]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[10]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[10]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[10]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          break;
        case 2114:
          fHistPtYInput[11]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[11]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[11]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[11]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          break;
        case 2214:
          fHistPtYInput[12]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[12]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[12]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[12]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          break;
        case 1114:
          fHistPtYInput[13]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[13]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[13]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[13]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          break;
        case 2224:
          fHistPtYInput[14]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[14]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[14]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[14]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          break;
        case 321:
          fHistPtYInput[15]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[15]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[15]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[15]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          break;
        case -321:
          fHistPtYInput[16]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[16]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[16]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[16]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          break;
        case -3334:
          fHistPtYInput[17]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[17]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[17]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[17]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          break;
        case 3334:
          fHistPtYInput[18]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[18]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[18]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[18]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          break;
        case -3312:
          fHistPtYInput[19]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[19]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[19]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[19]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          break;
        case 3312:
          fHistPtYInput[20]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[20]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[20]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[20]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          break;
        case 3224:
          fHistPtYInput[21]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[21]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[21]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[21]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          break;
        case 3114:
          fHistPtYInput[22]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[22]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[22]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[22]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          break;
        case 313:
          fHistPtYInput[23]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
          fHistPtPhiInput[23]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
          fHistDecayChannelsInput[23]->Fill(0., particle->GetWeight());
          fHistDecayChannelsInput[23]->Fill(GetDecayChannel(fMCEvent, particle), particle->GetWeight());
          break;
        default:
          fHistPdgInputRest->Fill(particle->GetPdgCode());
          break;
      }
    }
    
    // gamma from X/pi0 from source
    if (particle->GetPdgCode()==22 && motherHasMother) {
      if (grandMotherIsPrimary && fHasMother[GetParticlePosLocal(grandMotherParticle->GetPdgCode())]) {
        
        switch(grandMotherParticle->GetPdgCode()){
          case 310:
            fHistPtYGammaFromXFromInput[0]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiGammaFromXFromInput[0]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            if(motherParticle->GetPdgCode() == 111) {
              fHistPtYGammaFromPi0FromInput[0]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
              fHistPtPhiGammaFromPi0FromInput[0]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            }
            break;
          case 130:
            fHistPtYGammaFromXFromInput[1]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiGammaFromXFromInput[1]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            if(motherParticle->GetPdgCode() == 111) {
              fHistPtYGammaFromPi0FromInput[1]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
              fHistPtPhiGammaFromPi0FromInput[1]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            }
            break;
          case 3122:
            fHistPtYGammaFromXFromInput[2]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
            fHistPtPhiGammaFromXFromInput[2]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            if(motherParticle->GetPdgCode() == 111) {
              fHistPtYGammaFromPi0FromInput[2]->Fill(particle->Pt(), particle->Y(), particle->GetWeight());
              fHistPtPhiGammaFromPi0FromInput[2]->Fill(particle->Pt(), particle->Phi(), particle->GetWeight());
            }
            break;
        }
      }
    }
  }
}

//_________________________________________________________________________________
void AliAnalysisTaskHadronicCocktailMC::Terminate(const Option_t *)
{
  //fOutputContainer->Print(); // Will crash on GRID
}

//_________________________________________________________________________________
void AliAnalysisTaskHadronicCocktailMC::SetAnalyzedParticle(Int_t flag){
  if (flag>=0 && flag<3) {
    fAnalyzedMeson = flag;
    // neutral pion flag
    if (fAnalyzedMeson==0)  fAnalyzeNeutralPi = kTRUE;
    else                    fAnalyzeNeutralPi = kFALSE;
    // charged pion flag
    if (fAnalyzedMeson==2)  fAnalyzeChargedPi = kTRUE;
    else                    fAnalyzeChargedPi = kFALSE;
  } else {
    AliFatalClass(Form("Requested particle type flag (%d) for analysis not recognized!",flag));
  }
}

//_________________________________________________________________________________
void AliAnalysisTaskHadronicCocktailMC::SetLogBinningXTH1(TH1* histoRebin){
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
void AliAnalysisTaskHadronicCocktailMC::SetLogBinningXTH2(TH2* histoRebin){
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
void AliAnalysisTaskHadronicCocktailMC::InitializeDecayChannelHist(TH1F* hist, Int_t np) {
  
  switch (np) {

    case 0: // eta
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#gamma #gamma");
      hist->GetXaxis()->SetBinLabel(3,"#pi^{0} #pi^{0} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(4,"#pi^{0} #gamma #gamma");
      hist->GetXaxis()->SetBinLabel(5,"#pi^{+} #pi^{-} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(6,"#pi^{+} #pi^{-} #gamma");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;
      
    case 1: // K0s
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#pi^{0} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(3,"#pi^{+} #pi^{-}");
      hist->GetXaxis()->SetBinLabel(4,"#pi^{+} #pi^{-} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(5,"#pi^{0} #gamma #gamma");
      hist->GetXaxis()->SetBinLabel(6,"#pi^{0} e^{+} e^{-}");
      hist->GetXaxis()->SetBinLabel(7,"#pi^{0} #mu^{+} #mu^{-}");
      hist->GetXaxis()->SetBinLabel(8,"#pi^{+} #pi^{-} #gamma");
      hist->GetXaxis()->SetBinLabel(9,"#pi^{+} #pi^{-} e^{+} e^{-}");
      hist->GetXaxis()->SetBinLabel(10,"#pi^{#pm} e^{#mp} #nu");
      hist->GetXaxis()->SetBinLabel(11,"#pi^{#pm} #mu^{#mp} #nu");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;

    case 2: // K0l
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#pi^{0} #pi^{0} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(3,"#pi^{+} #pi^{-} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(4,"#pi^{0} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(5,"#pi^{0} #gamma #gamma");
      hist->GetXaxis()->SetBinLabel(6,"#pi^{0} e^{+} e^{-} #gamma");
      hist->GetXaxis()->SetBinLabel(7,"#pi^{0} #pi^{#pm} e^{#mp} #nu");
      hist->GetXaxis()->SetBinLabel(8,"#pi^{+} #pi^{-}");
      hist->GetXaxis()->SetBinLabel(9,"#pi^{#pm} e^{#mp} #nu");
      hist->GetXaxis()->SetBinLabel(10,"#pi^{#pm} #mu^{#mp} #nu");
      hist->GetXaxis()->SetBinLabel(11,"#pi^{#pm} e^{#mp} #nu #e{+} e^{-}");
      hist->GetXaxis()->SetBinLabel(12,"#pi^{#pm} e^{#mp} #nu #gamma");
      hist->GetXaxis()->SetBinLabel(13,"#pi^{#pm} #mu^{#mp} #nu #gamma");
      hist->GetXaxis()->SetBinLabel(14,"#pi^{+} #pi^{-} #gamma");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;
      
    case 3: // Lambda
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"p #pi^{-}");
      hist->GetXaxis()->SetBinLabel(3,"n #pi^{0}");
      hist->GetXaxis()->SetBinLabel(4,"p #pi^{-} #gamma");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;

    case 4: // rho0
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#pi^{+} #pi^{-}");
      hist->GetXaxis()->SetBinLabel(3,"#pi^{0} #gamma");
      hist->GetXaxis()->SetBinLabel(4,"#eta #gamma");
      hist->GetXaxis()->SetBinLabel(5,"#pi^{0} #pi^{0} #gamma");
      hist->GetXaxis()->SetBinLabel(6,"#pi^{+} #pi^{-} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(7,"#pi^{+} #pi^{-} #pi^{0} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(8,"#pi^{+} #pi^{-} #gamma");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;

    case 5: // eta'
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#pi^{+} #pi^{-} #eta");
      hist->GetXaxis()->SetBinLabel(3,"#pi^{0} #pi^{0} #eta");
      hist->GetXaxis()->SetBinLabel(4,"#pi^{0} #pi^{0} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(5,"#pi^{+} #pi^{-} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(6,"#pi^{#pm} #rho^{#mp}");
      hist->GetXaxis()->SetBinLabel(7,"#pi^{+} #pi^{-} #pi^{0} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(8,"#pi^{+} #pi^{-} e^{+} e^{-}");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;

    case 6: // omega
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#pi^{+} #pi^{-} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(3,"#pi^{0} #gamma");
      hist->GetXaxis()->SetBinLabel(4,"#eta #gamma");
      hist->GetXaxis()->SetBinLabel(5,"#pi^{0} e^{+} e^{-}");
      hist->GetXaxis()->SetBinLabel(6,"#pi^{0} #mu^{+} #mu^{-}");
      hist->GetXaxis()->SetBinLabel(7,"#pi^{0} #pi^{0} #gamma");
      hist->GetXaxis()->SetBinLabel(8,"#pi^{+} #pi^{-}");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;

    case 7: // rho+
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#pi^{+} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(3,"#pi^{+} #gamma");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;

    case 8: // rho-
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#pi^{-} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(3,"#pi^{-} #gamma");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;

    case 9: // phi
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"K^{+} K^{-}");
      hist->GetXaxis()->SetBinLabel(3,"K^{0}_{L} K^{0}_{S}");
      hist->GetXaxis()->SetBinLabel(4,"#eta #gamma");
      hist->GetXaxis()->SetBinLabel(5,"#pi^{0} #gamma");
      hist->GetXaxis()->SetBinLabel(6,"#eta e^{+} e^{-}");
      hist->GetXaxis()->SetBinLabel(7,"#omega #pi^{0}");
      hist->GetXaxis()->SetBinLabel(8,"#pi^{0} #pi^{0} #gamma");
      hist->GetXaxis()->SetBinLabel(9,"#pi^{0} e^{+} e^{-}");
      hist->GetXaxis()->SetBinLabel(10,"#pi^{0} #eta #gamma");
      hist->GetXaxis()->SetBinLabel(11,"#pi^{+} #pi^{-} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(12,"#rho^{0} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(13,"#rho^{#pm} #pi^{#mp}");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;

    case 10: // J/psi
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#pi^{0} X");
      hist->GetXaxis()->SetBinLabel(3,"#eta X");
      hist->GetXaxis()->SetBinLabel(4,"#pi^{+} X");
      hist->GetXaxis()->SetBinLabel(5,"#pi^{-} X");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;
      
    case 11: //Delta0
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"n #pi^{0}");
      hist->GetXaxis()->SetBinLabel(3,"p #pi^{-}");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;
      
    case 12: // Delta+
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"n #pi^{+}");
      hist->GetXaxis()->SetBinLabel(3,"p #pi^{0}");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;

    case 13: // Delta-
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"n #pi^{-}");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;
      
    case 14: //Delta++
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"p #pi^{+}");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;

    case 15: // K+
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#mu^{+} #nu");
      hist->GetXaxis()->SetBinLabel(3,"#mu^{+} #nu #gamma");
      hist->GetXaxis()->SetBinLabel(4,"#pi^{0} e^{+} #nu");
      hist->GetXaxis()->SetBinLabel(5,"#pi^{0} #mu^{+} #nu");
      hist->GetXaxis()->SetBinLabel(6,"#pi^{0} #pi^{0} e^{+} #nu");
      hist->GetXaxis()->SetBinLabel(7,"#pi^{+} #pi^{-} e^{+} #nu");
      hist->GetXaxis()->SetBinLabel(8,"#pi^{+} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(9,"#pi^{+} #pi^{0} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(10,"#pi^{+} #pi^{+} #pi^{-}");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;

    case 16: // K-
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#mu^{-} #nu");
      hist->GetXaxis()->SetBinLabel(3,"#mu^{-} #nu #gamma");
      hist->GetXaxis()->SetBinLabel(4,"#pi^{0} e^{-} #nu");
      hist->GetXaxis()->SetBinLabel(5,"#pi^{0} #mu^{-} #nu");
      hist->GetXaxis()->SetBinLabel(6,"#pi^{0} #pi^{0} e^{-} #nu");
      hist->GetXaxis()->SetBinLabel(7,"#pi^{+} #pi^{-} e^{-} #nu");
      hist->GetXaxis()->SetBinLabel(8,"#pi^{-} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(9,"#pi^{-} #pi^{0} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(10,"#pi^{+} #pi^{-} #pi^{-}");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;

    case 17: // Omega+
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#Lambda K^{+}");
      hist->GetXaxis()->SetBinLabel(3,"#Xi^{0} #pi^{+}");
      hist->GetXaxis()->SetBinLabel(4,"#Xi^{+} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(5,"#Xi^{+} #pi^{+} #pi^{-}");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;

    case 18: // Omega-
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#Lambda K^{-}");
      hist->GetXaxis()->SetBinLabel(3,"#Xi^{0} #pi^{-}");
      hist->GetXaxis()->SetBinLabel(4,"#Xi^{-} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(5,"#Xi^{-} #pi^{+} #pi^{-}");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;

    case 19: // Xi+
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#Lambda #pi^{+}");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;

    case 20: // Xi-
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#Lambda #pi^{-}");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;

    case 21: // Sigma(1385)+
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#Lambda #pi^{+}");
      hist->GetXaxis()->SetBinLabel(3,"#Sigma^{0} #pi^{+}");
      hist->GetXaxis()->SetBinLabel(4,"#Sigma^{+} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;

    case 22: // Sigma(1385)-
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"#Lambda #pi^{-}");
      hist->GetXaxis()->SetBinLabel(3,"#Sigma^{0} #pi^{-}");
      hist->GetXaxis()->SetBinLabel(4,"#Sigma^{-} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;

    case 23: // K*(892)0
      hist->GetXaxis()->SetBinLabel(1,"all");
      hist->GetXaxis()->SetBinLabel(2,"K^{+} #pi^{-}");
      hist->GetXaxis()->SetBinLabel(3,"K^{-} #pi^{+}");
      hist->GetXaxis()->SetBinLabel(4,"K^{0} #pi^{0}");
      hist->GetXaxis()->SetBinLabel(20,"rest");
      break;

    default:
      break;
  }
}

//_________________________________________________________________________________
Float_t AliAnalysisTaskHadronicCocktailMC::GetDecayChannel(AliMCEvent* mcEvent, TParticle* part) {

  Int_t nDaughters = part->GetNDaughters();
  if (nDaughters > 10) return 19.;
  
  std::vector<Long64_t> *PdgDaughter = new std::vector<Long64_t>(nDaughters);
  Long64_t tempPdgCode = 0;
  for (Int_t i=0; i<nDaughters; i++) {
    tempPdgCode = (Long64_t)((TParticle*)mcEvent->Particle(part->GetFirstDaughter()+i))->GetPdgCode();
    if (TMath::Abs(tempPdgCode) == 111 || TMath::Abs(tempPdgCode) == 113 || TMath::Abs(tempPdgCode) == 130 || TMath::Abs(tempPdgCode) == 310 || TMath::Abs(tempPdgCode) == 223 || TMath::Abs(tempPdgCode) == 221 || TMath::Abs(tempPdgCode) == 331 || TMath::Abs(tempPdgCode) == 2112 || TMath::Abs(tempPdgCode) == 3122 || TMath::Abs(tempPdgCode) == 9000111 || TMath::Abs(tempPdgCode) == 9010221 || TMath::Abs(tempPdgCode) == 3322)
      tempPdgCode = TMath::Abs(tempPdgCode);
    PdgDaughter->at(i) = tempPdgCode;
  }
  std::sort(PdgDaughter->begin(), PdgDaughter->end());
  
  Double_t returnVal = -1.;
  
  switch (part->GetPdgCode()) {
    case 221:
      if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 22)
        returnVal = 1.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 111)
        returnVal = 2.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 22 && PdgDaughter->at(2) == 111)
        returnVal = 3.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 211)
        returnVal = 4.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 22 && PdgDaughter->at(2) == 211)
        returnVal = 5.;
      else
        returnVal = 19.;
      break;
      
    case 310:
      if (nDaughters == 2 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 111)
        returnVal = 1.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 211)
        returnVal = 2.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 211)
        returnVal = 3.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 22 && PdgDaughter->at(2) == 111)
        returnVal = 4.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11 && PdgDaughter->at(2) == 111)
        returnVal = 5.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -13 && PdgDaughter->at(1) == 13 && PdgDaughter->at(2) == 111)
        returnVal = 6.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 22 && PdgDaughter->at(2) == 211)
        returnVal = 7.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == -11  && PdgDaughter->at(2) == 11 && PdgDaughter->at(3) == 211)
        returnVal = 8.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == -11 && PdgDaughter->at(2) == 12)
        returnVal = 9.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -12 && PdgDaughter->at(1) == 11 && PdgDaughter->at(2) == 211)
        returnVal = 9.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == -13 && PdgDaughter->at(2) == 14)
        returnVal = 10.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -14 && PdgDaughter->at(1) == 13 && PdgDaughter->at(2) == 211)
        returnVal = 10.;
      else
        returnVal = 19.;
      break;
      
    case 130:
      if (nDaughters == 3 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 111)
        returnVal = 1.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 211)
        returnVal = 2.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 111)
        returnVal = 3.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 22 && PdgDaughter->at(2) == 111)
        returnVal = 4.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11  && PdgDaughter->at(2) == 22 && PdgDaughter->at(3) == 111)
        returnVal = 5.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == -11  && PdgDaughter->at(2) == 12 && PdgDaughter->at(3) == 111)
        returnVal = 6.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -12 && PdgDaughter->at(1) == 11  && PdgDaughter->at(2) == 111 && PdgDaughter->at(3) == 211)
        returnVal = 6.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 211)
        returnVal = 7.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == -11  && PdgDaughter->at(2) == 12)
        returnVal = 8.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -12 && PdgDaughter->at(1) == 11 && PdgDaughter->at(2) == 211)
        returnVal = 8.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == -13  && PdgDaughter->at(2) == 14)
        returnVal = 9.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -14 && PdgDaughter->at(1) == 13 && PdgDaughter->at(2) == 211)
        returnVal = 9.;
      else if (nDaughters == 5 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == -11 && PdgDaughter->at(2) == -11 && PdgDaughter->at(3) == 11 && PdgDaughter->at(4) == 12)
        returnVal = 10.;
      else if (nDaughters == 5 && PdgDaughter->at(0) == -12 && PdgDaughter->at(1) == -11 && PdgDaughter->at(2) == 11 && PdgDaughter->at(3) == 11 && PdgDaughter->at(4) == 211)
        returnVal = 10.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == -11  && PdgDaughter->at(2) == 12 && PdgDaughter->at(3) == 22)
        returnVal = 11.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -12 && PdgDaughter->at(1) == 11 && PdgDaughter->at(2) == 22 && PdgDaughter->at(3) == 211)
        returnVal = 11.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == -13  && PdgDaughter->at(2) == 14 && PdgDaughter->at(3) == 22)
        returnVal = 12.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -14 && PdgDaughter->at(1) == 13 && PdgDaughter->at(2) == 22 && PdgDaughter->at(3) == 211)
        returnVal = 12.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 22 && PdgDaughter->at(2) == 211)
        returnVal = 13.;
      else
        returnVal = 19.;
      break;
      
    case 3122:
      if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 2212)
        returnVal = 1.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 2112)
        returnVal = 2.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 22 && PdgDaughter->at(2) == 2212)
        returnVal = 3.;
      else
        returnVal = 19.;
      break;
      
    case 113:
      if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 211)
        returnVal = 1.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 111)
        returnVal = 2.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 221)
        returnVal = 3.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 111)
        returnVal = 4.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 211)
        returnVal = 5.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 111 && PdgDaughter->at(3) == 211)
        returnVal = 6.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 22 && PdgDaughter->at(2) == 211)
        returnVal = 7.;
      else
        returnVal = 19.;
      break;

    case 331:
      if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 211 && PdgDaughter->at(2) == 221)
        returnVal = 1.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 221)
        returnVal = 2.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 111)
        returnVal = 3.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 211)
        returnVal = 4.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 213)
        returnVal = 5.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -213 && PdgDaughter->at(1) == 211)
        returnVal = 5.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 111 && PdgDaughter->at(3) == 211)
        returnVal = 6.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == -11 && PdgDaughter->at(2) == 11 && PdgDaughter->at(3) == 211)
        returnVal = 7.;
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
      else if (nDaughters == 3 && PdgDaughter->at(0) == -13 && PdgDaughter->at(1) == 13 && PdgDaughter->at(2) == 111)
        returnVal = 5.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 111)
        returnVal = 6.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 211)
        returnVal = 7.;
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
        returnVal = 1.;
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
      else if (nDaughters == 3 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11 && PdgDaughter->at(2) == 221)
        returnVal = 5.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 223)
        returnVal = 6.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 111)
        returnVal = 7.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 11 && PdgDaughter->at(2) == 111)
        returnVal = 8.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == 22 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 221)
        returnVal = 9.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 211)
        returnVal = 10.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 113)
        returnVal = 11.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 213)
        returnVal = 12.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -213 && PdgDaughter->at(1) == 211)
        returnVal = 12.;
      else
        returnVal = 19.;
      break;
      
    case 443:
      if (std::find(PdgDaughter->begin(), PdgDaughter->end(), 111) != PdgDaughter->end())
        returnVal = 1.;
      else if (std::find(PdgDaughter->begin(), PdgDaughter->end(), 221) != PdgDaughter->end())
        returnVal = 2.;
      else if (std::find(PdgDaughter->begin(), PdgDaughter->end(), 211) != PdgDaughter->end())
        returnVal = 3.;
      else if (std::find(PdgDaughter->begin(), PdgDaughter->end(), -211) != PdgDaughter->end())
        returnVal = 4.;
      else
        returnVal = 19.;
      break;
      
    case 2114:
      if (nDaughters == 2 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 2112)
        returnVal = 1.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 2212)
        returnVal = 2.;
      else
        returnVal = 19.;
      break;
      
    case 2214:
      if (nDaughters == 2 && PdgDaughter->at(0) == 211 && PdgDaughter->at(1) == 2112)
        returnVal = 1.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 2212)
        returnVal = 2.;
      else
        returnVal = 19.;
      break;

    case 1114:
      if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 2112)
        returnVal = 1.;
      else
        returnVal = 19.;
      break;

    case 2224:
      if (nDaughters == 2 && PdgDaughter->at(0) == 211 && PdgDaughter->at(1) == 2212)
        returnVal = 1.;
      else
        returnVal = 19.;
      break;
      
    case 321:
      if (nDaughters == 2 && PdgDaughter->at(0) == -13 && PdgDaughter->at(1) == 14)
        returnVal = 1.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -13 && PdgDaughter->at(1) == 14  && PdgDaughter->at(2) == 22)
        returnVal = 2.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 12  && PdgDaughter->at(2) == 111)
        returnVal = 3.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -13 && PdgDaughter->at(1) == 14  && PdgDaughter->at(2) == 111)
        returnVal = 4.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -11 && PdgDaughter->at(1) == 12  && PdgDaughter->at(2) == 111  && PdgDaughter->at(3) == 111)
        returnVal = 5.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == -11 && PdgDaughter->at(2) == 12  && PdgDaughter->at(3) == 211)
        returnVal = 6.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 211)
        returnVal = 7.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 211)
        returnVal = 8.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 211 && PdgDaughter->at(2) == 211)
        returnVal = 9.;
      else
        returnVal = 19.;
      break;

    case -321:
      if (nDaughters == 2 && PdgDaughter->at(0) == -14 && PdgDaughter->at(1) == 13)
        returnVal = 1.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -14 && PdgDaughter->at(1) == 13  && PdgDaughter->at(2) == 22)
        returnVal = 2.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -12 && PdgDaughter->at(1) == 11  && PdgDaughter->at(2) == 111)
        returnVal = 3.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -14 && PdgDaughter->at(1) == 13  && PdgDaughter->at(2) == 111)
        returnVal = 4.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -12 && PdgDaughter->at(1) == 11  && PdgDaughter->at(2) == 111  && PdgDaughter->at(3) == 111)
        returnVal = 5.;
      else if (nDaughters == 4 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == -12 && PdgDaughter->at(2) == 11  && PdgDaughter->at(3) == 211)
        returnVal = 6.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 111)
        returnVal = 7.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 111 && PdgDaughter->at(2) == 111)
        returnVal = 8.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == -211 && PdgDaughter->at(2) == 211)
        returnVal = 9.;
      else
        returnVal = 19.;
      break;

    case -3334:
      if (nDaughters == 2 && PdgDaughter->at(0) == 321 && PdgDaughter->at(1) == 3122)
        returnVal = 1.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 211 && PdgDaughter->at(1) == 3322)
        returnVal = 2.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -3312 && PdgDaughter->at(1) == 111)
        returnVal = 3.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -3312 && PdgDaughter->at(1) == -211 && PdgDaughter->at(2) == 211)
        returnVal = 4.;
      else
        returnVal = 19.;
      break;

    case 3334:
      if (nDaughters == 2 && PdgDaughter->at(0) == -321 && PdgDaughter->at(1) == 3122)
        returnVal = 1.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 3322)
        returnVal = 2.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 3312)
        returnVal = 3.;
      else if (nDaughters == 3 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 211 && PdgDaughter->at(2) == 3312)
        returnVal = 4.;
      else
        returnVal = 19.;
      break;

    case -3312:
      if (nDaughters == 2 && PdgDaughter->at(0) == 211 && PdgDaughter->at(1) == 3122)
        returnVal = 1.;
      else
        returnVal = 19.;
      break;

    case 3312:
      if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 3122)
        returnVal = 1.;
      else
        returnVal = 19.;
      break;

    case 3224:
      if (nDaughters == 2 && PdgDaughter->at(0) == 211 && PdgDaughter->at(1) == 3122)
        returnVal = 1.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 211 && PdgDaughter->at(1) == 3212)
        returnVal = 2.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 3222)
        returnVal = 3.;
      else
        returnVal = 19.;
      break;

    case 3114:
      if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 3122)
        returnVal = 1.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 3212)
        returnVal = 2.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 3112)
        returnVal = 3.;
      else
        returnVal = 19.;
      break;

    case 313:
      if (nDaughters == 2 && PdgDaughter->at(0) == -211 && PdgDaughter->at(1) == 321)
        returnVal = 1.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == -321 && PdgDaughter->at(1) == 211)
        returnVal = 2.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 130)
        returnVal = 3.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 310)
        returnVal = 3.;
      else if (nDaughters == 2 && PdgDaughter->at(0) == 111 && PdgDaughter->at(1) == 311)
        returnVal = 3.;
      else
        returnVal = 19.;
      break;

    default:
      return -1.;
      break;
  }

  delete PdgDaughter;

  return returnVal;
}

//_________________________________________________________________________________
void AliAnalysisTaskHadronicCocktailMC::FillPythiaBranchingRatio(TH1F* histo, Int_t np) {

  Int_t kc, kfdp, nPart, firstChannel, lastChannel;
  Double_t BR, BRtot;
  std::vector<Int_t> pdgCodes;
  
  switch (np) {
      
    case 0:
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
        else if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == 111 && pdgCodes[2] == 211)
          histo->SetBinContent(5, BR);
        else if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == 22 && pdgCodes[2] == 211)
          histo->SetBinContent(6, BR);
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;
      
    case 1:
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
        else if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == 111 && pdgCodes[2] == 211)
          histo->SetBinContent(4, BR);
        else if (nPart == 3 && pdgCodes[0] == 22 && pdgCodes[1] == 22 && pdgCodes[2] == 111)
          histo->SetBinContent(5, BR);
        else if (nPart == 3 && pdgCodes[0] == -11 && pdgCodes[1] == 11 && pdgCodes[2] == 111)
          histo->SetBinContent(6, BR);
        else if (nPart == 3 && pdgCodes[0] == -13 && pdgCodes[1] == 13 && pdgCodes[2] == 111)
          histo->SetBinContent(7, BR);
        else if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == 22 && pdgCodes[2] == 211)
          histo->SetBinContent(8, BR);
        else if (nPart == 4 && pdgCodes[0] == -211 && pdgCodes[1] == -11 && pdgCodes[2] == 11 && pdgCodes[3] == 211)
          histo->SetBinContent(9, BR);
        else if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == -11 && pdgCodes[2] == 12)
          histo->SetBinContent(10, BR+histo->GetBinContent(10));
        else if (nPart == 3 && pdgCodes[0] == -12 && pdgCodes[1] == 11 && pdgCodes[2] == 211)
          histo->SetBinContent(10, BR+histo->GetBinContent(10));
        else if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == -13 && pdgCodes[2] == 14)
          histo->SetBinContent(11, BR+histo->GetBinContent(11));
        else if (nPart == 3 && pdgCodes[0] == -14 && pdgCodes[1] == 13 && pdgCodes[2] == 211)
          histo->SetBinContent(11, BR+histo->GetBinContent(11));
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;
      
    case 2:
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
        else if (nPart == 2 && pdgCodes[0] == 111 && pdgCodes[1] == 111)
          histo->SetBinContent(4, BR);
        else if (nPart == 3 && pdgCodes[0] == 22 && pdgCodes[1] == 22 && pdgCodes[2] == 111)
          histo->SetBinContent(5, BR);
        else if (nPart == 4 && pdgCodes[0] == -11 && pdgCodes[1] == 11 && pdgCodes[2] == 22 && pdgCodes[3] == 111)
          histo->SetBinContent(6, BR);
        else if (nPart == 4 && pdgCodes[0] == -12 && pdgCodes[1] == 11 && pdgCodes[2] == 111 && pdgCodes[3] == 211)
          histo->SetBinContent(7, BR+histo->GetBinContent(7));
        else if (nPart == 4 && pdgCodes[0] == -211 && pdgCodes[1] == -11 && pdgCodes[2] == 12 && pdgCodes[3] == 111)
          histo->SetBinContent(7, BR+histo->GetBinContent(7));
        else if (nPart == 2 && pdgCodes[0] == -211 && pdgCodes[1] == 211)
          histo->SetBinContent(8, BR);
        else if (nPart == 3 && pdgCodes[0] == -12 && pdgCodes[1] == 11 && pdgCodes[2] == 211)
          histo->SetBinContent(9, BR+histo->GetBinContent(9));
        else if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == -11 && pdgCodes[2] == 12)
          histo->SetBinContent(9, BR+histo->GetBinContent(9));
        else if (nPart == 3 && pdgCodes[0] == -14 && pdgCodes[1] == 13 && pdgCodes[2] == 211)
          histo->SetBinContent(10, BR+histo->GetBinContent(10));
        else if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == -13 && pdgCodes[2] == 14)
          histo->SetBinContent(10, BR+histo->GetBinContent(10));
        else if (nPart == 5 && pdgCodes[0] == -12 && pdgCodes[1] == -11  && pdgCodes[2] == 11  && pdgCodes[3] == 11 && pdgCodes[4] == 211)
          histo->SetBinContent(11, BR+histo->GetBinContent(11));
        else if (nPart == 5 && pdgCodes[0] == -211 && pdgCodes[1] == -11 && pdgCodes[2] == -11  && pdgCodes[3] == 11 && pdgCodes[4] == 12)
          histo->SetBinContent(11, BR+histo->GetBinContent(11));
        else if (nPart == 4 && pdgCodes[0] == -12 && pdgCodes[1] == 11 && pdgCodes[2] == 22 && pdgCodes[3] == 211)
          histo->SetBinContent(12, BR+histo->GetBinContent(12));
        else if (nPart == 4 && pdgCodes[0] == -211 && pdgCodes[1] == -11 && pdgCodes[2] == 12 && pdgCodes[3] == 22)
          histo->SetBinContent(12, BR+histo->GetBinContent(12));
        else if (nPart == 4 && pdgCodes[0] == -14 && pdgCodes[1] == 13 && pdgCodes[2] == 22 && pdgCodes[3] == 211)
          histo->SetBinContent(13, BR+histo->GetBinContent(13));
        else if (nPart == 4 && pdgCodes[0] == -211 && pdgCodes[1] == -13 && pdgCodes[2] == 14 && pdgCodes[3] == 22)
          histo->SetBinContent(13, BR+histo->GetBinContent(13));
        else if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == 22 && pdgCodes[2] == 211)
          histo->SetBinContent(14, BR);
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;
      
    case 3:
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
        else if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == 22 && pdgCodes[2] == 2212)
          histo->SetBinContent(4, BR);
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
        else if (nPart == 2 && pdgCodes[0] == 22 && pdgCodes[1] == 111)
          histo->SetBinContent(3, BR);
        else if (nPart == 2 && pdgCodes[0] == 22 && pdgCodes[1] == 221)
          histo->SetBinContent(4, BR);
        else if (nPart == 3 && pdgCodes[0] == 22 && pdgCodes[1] == 111 && pdgCodes[2] == 111)
          histo->SetBinContent(5, BR);
        else if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == 111 && pdgCodes[2] == 211)
          histo->SetBinContent(6, BR);
        else if (nPart == 4 && pdgCodes[0] == -211 && pdgCodes[1] == 111 && pdgCodes[2] == 111 && pdgCodes[3] == 211)
          histo->SetBinContent(7, BR);
        else if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == 22 && pdgCodes[2] == 211)
          histo->SetBinContent(8, BR);
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;
      
    case 5:
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
        else if (nPart == 3 && pdgCodes[0] == 111 && pdgCodes[1] == 111 && pdgCodes[2] == 221)
          histo->SetBinContent(3, BR);
        else if (nPart == 3 && pdgCodes[0] == 111 && pdgCodes[1] == 111 && pdgCodes[2] == 111)
          histo->SetBinContent(4, BR);
        else if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == 111 && pdgCodes[2] == 211)
          histo->SetBinContent(5, BR);
        else if (nPart == 2 && pdgCodes[0] == -211 && pdgCodes[1] == 213)
          histo->SetBinContent(6, BR+histo->GetBinContent(6));
        else if (nPart == 2 && pdgCodes[0] == -213 && pdgCodes[1] == 211)
          histo->SetBinContent(6, BR+histo->GetBinContent(6));
        else if (nPart == 4 && pdgCodes[0] == -211 && pdgCodes[1] == 111 && pdgCodes[2] == 111 && pdgCodes[3] == 211)
          histo->SetBinContent(7, BR+histo->GetBinContent(7));
        else if (nPart == 4 && pdgCodes[0] == -211 && pdgCodes[1] == -11 && pdgCodes[2] == 11 && pdgCodes[3] == 211)
          histo->SetBinContent(8, BR+histo->GetBinContent(8));
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;
      
    case 6:
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
        else if (nPart == 3 && pdgCodes[0] == -13 && pdgCodes[1] == 13 && pdgCodes[2] == 111)
          histo->SetBinContent(6, BR);
        else if (nPart == 3 && pdgCodes[0] == 22 && pdgCodes[1] == 111 && pdgCodes[2] == 111)
          histo->SetBinContent(7, BR);
        else if (nPart == 2 && pdgCodes[0] == -211 && pdgCodes[1] == 211)
          histo->SetBinContent(8, BR);
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;
      
    case 7: case 8: // PYTHIA doesn't distinguish between rho+/rho-
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

    case 9:
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
        else if (nPart == 3 && pdgCodes[0] == -11 && pdgCodes[1] == 11 && pdgCodes[2] == 221)
          histo->SetBinContent(6, BR);
        else if (nPart == 2 && pdgCodes[0] == 111 && pdgCodes[1] == 223)
          histo->SetBinContent(7, BR);
        else if (nPart == 3 && pdgCodes[0] == 22 && pdgCodes[1] == 111 && pdgCodes[2] == 111)
          histo->SetBinContent(8, BR);
        else if (nPart == 3 && pdgCodes[0] == -11 && pdgCodes[1] == 11 && pdgCodes[2] == 111)
          histo->SetBinContent(9, BR);
        else if (nPart == 3 && pdgCodes[0] == 22 && pdgCodes[1] == 111 && pdgCodes[2] == 221)
          histo->SetBinContent(10, BR);
        else if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == 111 && pdgCodes[2] == 211)
          histo->SetBinContent(11, BR);
        else if (nPart == 2 && pdgCodes[0] == 111 && pdgCodes[1] == 113)
          histo->SetBinContent(12, BR);
        else if (nPart == 2 && pdgCodes[0] == -211 && pdgCodes[1] == 213)
          histo->SetBinContent(13, BR+histo->GetBinContent(13));
        else if (nPart == 2 && pdgCodes[0] == -213 && pdgCodes[1] == 211)
          histo->SetBinContent(13, BR+histo->GetBinContent(13));
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;
      
    case 10:
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
        if (std::find(pdgCodes.begin(), pdgCodes.end(), 111) != pdgCodes.end())
          histo->SetBinContent(2, BR+histo->GetBinContent(2));
        else if (std::find(pdgCodes.begin(), pdgCodes.end(), 221) != pdgCodes.end())
          histo->SetBinContent(3, BR+histo->GetBinContent(3));
        else if (std::find(pdgCodes.begin(), pdgCodes.end(), 211) != pdgCodes.end())
          histo->SetBinContent(4, BR+histo->GetBinContent(4));
        else if (std::find(pdgCodes.begin(), pdgCodes.end(), -211) != pdgCodes.end())
          histo->SetBinContent(5, BR+histo->GetBinContent(5));
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;
      
    case 11:
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
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;
      
    case 12:
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
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;

    case 13:
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
        if (nPart == 2 && pdgCodes[0] == -211 && pdgCodes[1] == 2112)
          histo->SetBinContent(2, BR);
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;

    case 14:
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
        if (nPart == 2 && pdgCodes[0] == 211 && pdgCodes[1] == 2212)
          histo->SetBinContent(2, BR);
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;

    case 15: case 16: // PYTHIA doesn't distinguish between K+/K-
      kc            = (AliPythia6::Instance())->Pycomp(321);
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
        if (nPart == 2 && pdgCodes[0] == -13 && pdgCodes[1] == 14)
          histo->SetBinContent(2, BR);
        else if (nPart == 3 && pdgCodes[0] == -13 && pdgCodes[1] == 14 && pdgCodes[2] == 22)
          histo->SetBinContent(3, BR);
        else if (nPart == 3 && pdgCodes[0] == -11 && pdgCodes[1] == 12 && pdgCodes[2] == 111)
          histo->SetBinContent(4, BR);
        else if (nPart == 3 && pdgCodes[0] == -13 && pdgCodes[1] == 14 && pdgCodes[2] == 111)
          histo->SetBinContent(5, BR);
        else if (nPart == 4 && pdgCodes[0] == -11 && pdgCodes[1] == 12 && pdgCodes[2] == 111 && pdgCodes[3] == 111)
          histo->SetBinContent(6, BR);
        else if (nPart == 4 && pdgCodes[0] == -211 && pdgCodes[1] == -11 && pdgCodes[2] == 12 && pdgCodes[3] == 211)
          histo->SetBinContent(7, BR);
        else if (nPart == 2 && pdgCodes[0] == 111 && pdgCodes[1] == 211)
          histo->SetBinContent(8, BR);
        else if (nPart == 3 && pdgCodes[0] == 111 && pdgCodes[1] == 111 && pdgCodes[2] == 211)
          histo->SetBinContent(9, BR);
        else if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == 211 && pdgCodes[2] == 211)
          histo->SetBinContent(10, BR);
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;

    case 17: case 18: // PYTHIA doesn't distinguish between Omega-/Omega+
      kc            = (AliPythia6::Instance())->Pycomp(3334);
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
        if (nPart == 2 && pdgCodes[0] == -321 && pdgCodes[1] == 3122)
          histo->SetBinContent(2, BR);
        else if (nPart == 2 && pdgCodes[0] == -211 && pdgCodes[1] == 3322)
          histo->SetBinContent(3, BR);
        else if (nPart == 2 && pdgCodes[0] == 111 && pdgCodes[1] == 3312)
          histo->SetBinContent(4, BR);
        else if (nPart == 3 && pdgCodes[0] == -211 && pdgCodes[1] == 211 && pdgCodes[2] == 3312)
          histo->SetBinContent(5, BR);
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;

    case 19: case 20: // PYTHIA doesn't distinguish between Xi-/Xi+
      kc            = (AliPythia6::Instance())->Pycomp(3312);
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
        if (nPart == 2 && pdgCodes[0] == -211 && pdgCodes[1] == 3122)
          histo->SetBinContent(2, BR);
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;
      
    case 21:
      kc            = (AliPythia6::Instance())->Pycomp(3224);
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
        if (nPart == 2 && pdgCodes[0] == 211 && pdgCodes[1] == 3122)
          histo->SetBinContent(2, BR);
        else if (nPart == 2 && pdgCodes[0] == 211 && pdgCodes[1] == 3212)
          histo->SetBinContent(3, BR);
        else if (nPart == 2 && pdgCodes[0] == 111 && pdgCodes[1] == 3222)
          histo->SetBinContent(4, BR);
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;

    case 22:
      kc            = (AliPythia6::Instance())->Pycomp(3114);
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
        if (nPart == 2 && pdgCodes[0] == -211 && pdgCodes[1] == 3122)
          histo->SetBinContent(2, BR);
        else if (nPart == 2 && pdgCodes[0] == -211 && pdgCodes[1] == 3212)
          histo->SetBinContent(3, BR);
        else if (nPart == 2 && pdgCodes[0] == 111 && pdgCodes[1] == 3112)
          histo->SetBinContent(4, BR);
        else
          histo->SetBinContent(20, BR+histo->GetBinContent(20));
        pdgCodes.clear();
      }
      histo->SetBinContent(1, BRtot);
      pdgCodes.clear();
      break;

    case 23:
      kc            = (AliPythia6::Instance())->Pycomp(313);
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
        if (nPart == 2 && pdgCodes[0] == -211 && pdgCodes[1] == 321)
          histo->SetBinContent(2, BR);
        else if (nPart == 2 && pdgCodes[0] == -321 && pdgCodes[1] == 211)
          histo->SetBinContent(3, BR);
        else if (nPart == 2 && pdgCodes[0] == 111 && pdgCodes[1] == 130)
          histo->SetBinContent(4, BR+histo->GetBinContent(4));
        else if (nPart == 2 && pdgCodes[0] == 111 && pdgCodes[1] == 310)
          histo->SetBinContent(4, BR+histo->GetBinContent(4));
        else if (nPart == 2 && pdgCodes[0] == 111 && pdgCodes[1] == 311)
          histo->SetBinContent(4, BR+histo->GetBinContent(4));
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
Int_t AliAnalysisTaskHadronicCocktailMC::GetParticlePosLocal(Int_t pdg) {
  
  Int_t returnVal = -9999;
  
  switch (pdg) {
    case 221:
      returnVal = 0;
      break;
    case 310:
      returnVal = 1;
      break;
    case 130:
      returnVal = 2;
      break;
    case 3122:
      returnVal = 3;
      break;
    case 113:
      returnVal = 4;
      break;
    case 331:
      returnVal = 5;
      break;
    case 223:
      returnVal = 6;
      break;
    case 213:
      returnVal = 7;
      break;
    case -213:
      returnVal = 8;
      break;
    case 333:
      returnVal = 9;
      break;
    case 443:
      returnVal = 10;
      break;
    case 2114:
      returnVal = 11;
      break;
    case 2214:
      returnVal = 12;
      break;
    case 1114:
      returnVal = 13;
      break;
    case 2224:
      returnVal = 14;
      break;
    case 321:
      returnVal = 15;
      break;
    case -321:
      returnVal = 16;
      break;
    case -3334:
      returnVal = 17;
      break;
    case 3334:
      returnVal = 18;
      break;
    case -3312:
      returnVal = 19;
      break;
    case 3312:
      returnVal = 20;
      break;
    case 3224:
      returnVal = 21;
      break;
    case 3114:
      returnVal = 22;
      break;
    case 313:
      returnVal = 23;
      break;
    default:
      break;
  }
  
  return returnVal;
}

//_________________________________________________________________________________
TH1* AliAnalysisTaskHadronicCocktailMC::SetHist1D(TH1* hist, TString histType, TString histName, TString xTitle, TString yTitle, Int_t nBinsX, Double_t xMin, Double_t xMax, Bool_t optSumw2) {
  
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
TH2* AliAnalysisTaskHadronicCocktailMC::SetHist2D(TH2* hist, TString histType, TString histName, TString xTitle, TString yTitle, Int_t nBinsX, Double_t xMin, Double_t xMax, Int_t nBinsY, Double_t yMin, Double_t yMax, Bool_t optSumw2) {
  
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
TH2* AliAnalysisTaskHadronicCocktailMC::SetHist2D(TH2* hist, TString histType, TString histName, TString xTitle, TString yTitle, Int_t nBinsX, Double_t xMin, Double_t xMax, Int_t nBinsY, Double_t* binsY, Bool_t optSumw2) {
  
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

