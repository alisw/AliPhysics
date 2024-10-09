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
// Class used to do analysis on conversion photons + calo photons
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
#include "AliAnalysisTaskGammaPureMC.h"
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

#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequenceArea.hh>

ClassImp(AliAnalysisTaskGammaPureMC)

//________________________________________________________________________
AliAnalysisTaskGammaPureMC::AliAnalysisTaskGammaPureMC(): AliAnalysisTaskSE(),
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
  fHistPtYGamma(nullptr),
  fHistPtYDirGamma(nullptr),
  fHistPtYPi0FromEta(nullptr),
  fHistPtYPi0FromLambda(nullptr),
  fHistPtYPi0FromK(nullptr),
  fHistPtYPiPlFromK(nullptr),
  fHistPtYPiMiFromK(nullptr),
  fHistPtYPi0GG(nullptr),
  fHistPtYPi0GGPCMAcc(nullptr),
  fHistPtYPi0GGEMCAcc(nullptr),
  fHistPtYPi0GGPHOAcc(nullptr),
  fHistPtYPi0GGPCMEMCAcc(nullptr),
  fHistPtYPi0GGPCMPHOAcc(nullptr),
  fHistPtAlphaPi0GGPCMAcc(nullptr),
  fHistPtAlphaPi0GGEMCAcc(nullptr),
  fHistPtAlphaPi0GGPHOAcc(nullptr),
  fHistPtAlphaPi0GGPCMEMCAcc(nullptr),
  fHistPtAlphaPi0GGPCMPHOAcc(nullptr),
  fHistPtYEtaGG(nullptr),
  fHistPtYEtaGGPCMAcc(nullptr),
  fHistPtYEtaGGEMCAcc(nullptr),
  fHistPtYEtaGGPHOAcc(nullptr),
  fHistPtYEtaGGPCMEMCAcc(nullptr),
  fHistPtYEtaGGPCMPHOAcc(nullptr),
  fHistPtAlphaEtaGGPCMAcc(nullptr),
  fHistPtAlphaEtaGGEMCAcc(nullptr),
  fHistPtAlphaEtaGGPHOAcc(nullptr),
  fHistPtAlphaEtaGGPCMEMCAcc(nullptr),
  fHistPtAlphaEtaGGPCMPHOAcc(nullptr),
  fHistPtYEtaPrimeGG(nullptr),
  fHistPtYEtaPrimeGGPCMAcc(nullptr),
  fHistPtYEtaPrimeGGEMCAcc(nullptr),
  fHistPtYEtaPrimeGGPHOAcc(nullptr),
  fHistPtYEtaPrimeGGPCMEMCAcc(nullptr),
  fHistPtYEtaPrimeGGPCMPHOAcc(nullptr),
  fHistPtAlphaEtaPrimeGGPCMAcc(nullptr),
  fHistPtAlphaEtaPrimeGGEMCAcc(nullptr),
  fHistPtAlphaEtaPrimeGGPHOAcc(nullptr),
  fHistPtAlphaEtaPrimeGGPCMEMCAcc(nullptr),
  fHistPtAlphaEtaPrimeGGPCMPHOAcc(nullptr),
  fHistPtYPi0FromKGG(nullptr),
  fHistPtYPi0FromKGGPCMAcc(nullptr),
  fHistPtYPi0FromKGGEMCAcc(nullptr),
  fHistPtYPi0FromKGGPCMEMCAcc(nullptr),
  fHistPtYPi0FromKGGEMCPCMAcc(nullptr),
  fHistPtYPi0FromKGGEMCAccSamePi0(nullptr),
  fHistPtYPi0FromKGGEMCAccDiffPi0(nullptr),
  fHistPtAlphaPi0FromKGG(nullptr),
  fHistPtAlphaPi0FromKGGPCMAcc(nullptr),
  fHistPtAlphaPi0FromKGGEMCAcc(nullptr),
  fHistPtAlphaPi0FromKGGPCMEMCAcc(nullptr),
  fHistPtAlphaPi0FromKGGEMCPCMAcc(nullptr),
  fHistPtAlphaPi0FromKGGEMCAccSamePi0(nullptr),
  fHistPtAlphaPi0FromKGGEMCAccDiffPi0(nullptr),
  fHistV0Mult(nullptr),
  fHistPtV0MultPi0GG(nullptr),
  fHistPtV0MultEtaGG(nullptr),
  fHistPtV0MultEtaPrimeGG(nullptr),
  fHistPtV0MultPi0GGPrompt(nullptr),
  fHistPtV0MultPi0GGFromEta(nullptr),
  fHistPtV0MultPi0GGFromOmega(nullptr),
  fHistPtV0MultPi0GGFromRest(nullptr),
  fHistPtV0MultPi0GGFromRho(nullptr),
  fHistPtV0MultEtaGGPrompt(nullptr),
  fHistPtV0MultEtaGGFromEtaPrim(nullptr),
  fHistPtV0MultEtaGGFromRest(nullptr),
  fHistPtV0MultGamma(nullptr),
  fHistPtV0MultDirGamma(nullptr),
  fHistPtV0MultChargedPi(nullptr),
  fHistPi0PtJetPt(nullptr),
  fHistEtaPtJetPt(nullptr),
  fHistPi0ZJetPt(nullptr),
  fHistEtaZJetPt(nullptr),
  fHistJetPtY(nullptr),
  fHistJetEta(nullptr),
  fHistJetPhi(nullptr),
  fIsK0(1),
  fIsMC(1),
  fMaxpT(100),
  fDoMultStudies(0),
  fNTracksInV0Acc(0),
  fIsEvtINELgtZERO(0),
  fDoFeedDownStudies(false),
  fHistPtPi0FromDecay(nullptr),
  fHistPtEtaFromDecay(nullptr),
  fHistPtOmegaFromDecay(nullptr),
  fHistPtYPi0Primordial(nullptr),
  fHistPtYEtaPrimordial(nullptr),
  fHistPtYOmegaPrimordial(nullptr),
  fDoJetStudies(false),
  fJetRadius(0.4),
  fJetMinE(1.),
  fJetAccEta(0.8),
  fJetParticleAcc(0.4),
  fJetParticleAccFF(1.2),
  fJetAlgorithm(fastjet::antikt_algorithm),
  fJetStrategy(fastjet::Best),
  fJetAreaType(fastjet::active_area),
  fJetRecombScheme(fastjet::BIpt_scheme),
  fJetGhostArea(0.01),
  fGhostEtaMax(1.5),
  fActiveAreaRepeats(1),
  fAreaType(fastjet::active_area),
  fVecJets({})
{

}

//________________________________________________________________________
AliAnalysisTaskGammaPureMC::AliAnalysisTaskGammaPureMC(const char *name):
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
  fHistPtYGamma(nullptr),
  fHistPtYDirGamma(nullptr),
  fHistPtYPi0FromEta(nullptr),
  fHistPtYPi0FromLambda(nullptr),
  fHistPtYPi0FromK(nullptr),
  fHistPtYPiPlFromK(nullptr),
  fHistPtYPiMiFromK(nullptr),
  fHistPtYPi0GG(nullptr),
  fHistPtYPi0GGPCMAcc(nullptr),
  fHistPtYPi0GGEMCAcc(nullptr),
  fHistPtYPi0GGPHOAcc(nullptr),
  fHistPtYPi0GGPCMEMCAcc(nullptr),
  fHistPtYPi0GGPCMPHOAcc(nullptr),
  fHistPtAlphaPi0GGPCMAcc(nullptr),
  fHistPtAlphaPi0GGEMCAcc(nullptr),
  fHistPtAlphaPi0GGPHOAcc(nullptr),
  fHistPtAlphaPi0GGPCMEMCAcc(nullptr),
  fHistPtAlphaPi0GGPCMPHOAcc(nullptr),
  fHistPtYEtaGG(nullptr),
  fHistPtYEtaGGPCMAcc(nullptr),
  fHistPtYEtaGGEMCAcc(nullptr),
  fHistPtYEtaGGPHOAcc(nullptr),
  fHistPtYEtaGGPCMEMCAcc(nullptr),
  fHistPtYEtaGGPCMPHOAcc(nullptr),
  fHistPtAlphaEtaGGPCMAcc(nullptr),
  fHistPtAlphaEtaGGEMCAcc(nullptr),
  fHistPtAlphaEtaGGPHOAcc(nullptr),
  fHistPtAlphaEtaGGPCMEMCAcc(nullptr),
  fHistPtAlphaEtaGGPCMPHOAcc(nullptr),
  fHistPtYEtaPrimeGG(nullptr),
  fHistPtYEtaPrimeGGPCMAcc(nullptr),
  fHistPtYEtaPrimeGGEMCAcc(nullptr),
  fHistPtYEtaPrimeGGPHOAcc(nullptr),
  fHistPtYEtaPrimeGGPCMEMCAcc(nullptr),
  fHistPtYEtaPrimeGGPCMPHOAcc(nullptr),
  fHistPtAlphaEtaPrimeGGPCMAcc(nullptr),
  fHistPtAlphaEtaPrimeGGEMCAcc(nullptr),
  fHistPtAlphaEtaPrimeGGPHOAcc(nullptr),
  fHistPtAlphaEtaPrimeGGPCMEMCAcc(nullptr),
  fHistPtAlphaEtaPrimeGGPCMPHOAcc(nullptr),
  fHistPtYPi0FromKGG(nullptr),
  fHistPtYPi0FromKGGPCMAcc(nullptr),
  fHistPtYPi0FromKGGEMCAcc(nullptr),
  fHistPtYPi0FromKGGPCMEMCAcc(nullptr),
  fHistPtYPi0FromKGGEMCPCMAcc(nullptr),
  fHistPtYPi0FromKGGEMCAccSamePi0(nullptr),
  fHistPtYPi0FromKGGEMCAccDiffPi0(nullptr),
  fHistPtAlphaPi0FromKGG(nullptr),
  fHistPtAlphaPi0FromKGGPCMAcc(nullptr),
  fHistPtAlphaPi0FromKGGEMCAcc(nullptr),
  fHistPtAlphaPi0FromKGGPCMEMCAcc(nullptr),
  fHistPtAlphaPi0FromKGGEMCPCMAcc(nullptr),
  fHistPtAlphaPi0FromKGGEMCAccSamePi0(nullptr),
  fHistPtAlphaPi0FromKGGEMCAccDiffPi0(nullptr),
  fHistV0Mult(nullptr),
  fHistPtV0MultPi0GG(nullptr),
  fHistPtV0MultEtaGG(nullptr),
  fHistPtV0MultEtaPrimeGG(nullptr),
  fHistPtV0MultPi0GGPrompt(nullptr),
  fHistPtV0MultPi0GGFromEta(nullptr),
  fHistPtV0MultPi0GGFromOmega(nullptr),
  fHistPtV0MultPi0GGFromRest(nullptr),
  fHistPtV0MultPi0GGFromRho(nullptr),
  fHistPtV0MultEtaGGPrompt(nullptr),
  fHistPtV0MultEtaGGFromEtaPrim(nullptr),
  fHistPtV0MultEtaGGFromRest(nullptr),
  fHistPtV0MultGamma(nullptr),
  fHistPtV0MultDirGamma(nullptr),
  fHistPtV0MultChargedPi(nullptr),
  fHistPi0PtJetPt(nullptr),
  fHistEtaPtJetPt(nullptr),
  fHistPi0ZJetPt(nullptr),
  fHistEtaZJetPt(nullptr),
  fHistJetPtY(nullptr),
  fHistJetEta(nullptr),
  fHistJetPhi(nullptr),
  fIsK0(1),
  fIsMC(1),
  fMaxpT(100),
  fDoMultStudies(0),
  fNTracksInV0Acc(0),
  fIsEvtINELgtZERO(0),
  fDoFeedDownStudies(false),
  fHistPtPi0FromDecay(nullptr),
  fHistPtEtaFromDecay(nullptr),
  fHistPtOmegaFromDecay(nullptr),
  fHistPtYPi0Primordial(nullptr),
  fHistPtYEtaPrimordial(nullptr),
  fHistPtYOmegaPrimordial(nullptr),
  fDoJetStudies(false),
  fJetRadius(0.4),
  fJetMinE(1.),
  fJetAccEta(0.8),
  fJetParticleAcc(0.4),
  fJetParticleAccFF(1.2),
  fJetAlgorithm(fastjet::antikt_algorithm),
  fJetStrategy(fastjet::Best),
  fJetAreaType(fastjet::active_area),
  fJetRecombScheme(fastjet::BIpt_scheme),
  fJetGhostArea(0.01),
  fGhostEtaMax(1.5),
  fActiveAreaRepeats(1),
  fAreaType(fastjet::active_area),
  fVecJets({})
{
  // Define output slots here
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskGammaPureMC::~AliAnalysisTaskGammaPureMC()
{

}

//________________________________________________________________________
void AliAnalysisTaskGammaPureMC::UserCreateOutputObjects(){

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

  fHistXSection               		= new TH1D("XSection", "XSection", 1000000, 0, 1e4);

  //   SetLogBinningXTH1(fHistXSection);
  fHistXSection->Sumw2();
  fOutputContainer->Add(fHistXSection);

  fHistPtHard                 		= new TH1F("PtHard", "PtHard", 400, 0, 200);
  fHistPtHard->Sumw2();
  fOutputContainer->Add(fHistPtHard);

  fHistPtYPi0                		= new TH2F("Pt_Y_Pi0","Pt_Y_Pi0", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYPi0->Sumw2();
  fOutputContainer->Add(fHistPtYPi0);

  fHistPtYPiPl                		= new TH2F("Pt_Y_PiPl","Pt_Y_PiPl", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYPiPl->Sumw2();
  fOutputContainer->Add(fHistPtYPiPl);

  fHistPtYPiMi                		= new TH2F("Pt_Y_PiMi","Pt_Y_PiMi", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYPiMi->Sumw2();
  fOutputContainer->Add(fHistPtYPiMi);

  fHistPtYEta                 		= new TH2F("Pt_Y_Eta","Pt_Y_Eta", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYEta->Sumw2();
  fOutputContainer->Add(fHistPtYEta);

  fHistPtYEtaPrime             		= new TH2F("Pt_Y_EtaPrime","Pt_Y_EtaPrime", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYEtaPrime->Sumw2();
  fOutputContainer->Add(fHistPtYEtaPrime);

  fHistPtYOmega               		= new TH2F("Pt_Y_Omega","Pt_Y_Omega", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYOmega->Sumw2();
  fOutputContainer->Add(fHistPtYOmega);

  fHistPtYRho0                		= new TH2F("Pt_Y_Rho0","Pt_Y_Rho0", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYRho0->Sumw2();
  fOutputContainer->Add(fHistPtYRho0);

  fHistPtYRhoPl               		= new TH2F("Pt_Y_RhoPl","Pt_Y_RhoPl", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYRhoPl->Sumw2();
  fOutputContainer->Add(fHistPtYRhoPl);

  fHistPtYRhoMi               		= new TH2F("Pt_Y_RhoMi","Pt_Y_RhoMi", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYRhoMi->Sumw2();
  fOutputContainer->Add(fHistPtYRhoMi);

  fHistPtYPhi                 		= new TH2F("Pt_Y_Phi","Pt_Y_Phi", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYPhi->Sumw2();
  fOutputContainer->Add(fHistPtYPhi);

  fHistPtYJPsi                		= new TH2F("Pt_Y_JPsi","Pt_Y_JPsi", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYJPsi->Sumw2();
  fOutputContainer->Add(fHistPtYJPsi);

  fHistPtYSigma0              		= new TH2F("Pt_Y_Sigma0","Pt_Y_Sigma0", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYSigma0->Sumw2();
  fOutputContainer->Add(fHistPtYSigma0);

  fHistPtYK0s                 		= new TH2F("Pt_Y_K0s","Pt_Y_K0s", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYK0s->Sumw2();
  fOutputContainer->Add(fHistPtYK0s);

  fHistPtYK0l                 		= new TH2F("Pt_Y_K0l","Pt_Y_K0l", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYK0l->Sumw2();
  fOutputContainer->Add(fHistPtYK0l);

  fHistPtYK0star              		= new TH2F("Pt_Y_K0star","Pt_Y_K0star", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYK0star->Sumw2();
  fOutputContainer->Add(fHistPtYK0star);

  fHistPtYDeltaPlPl           		= new TH2F("Pt_Y_DeltaPlPl","Pt_Y_DeltaPlPl", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYDeltaPlPl->Sumw2();
  fOutputContainer->Add(fHistPtYDeltaPlPl);

  fHistPtYDeltaPl             		= new TH2F("Pt_Y_DeltaPl","Pt_Y_DeltaPl", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYDeltaPl->Sumw2();
  fOutputContainer->Add(fHistPtYDeltaPl);

  fHistPtYDeltaMi             		= new TH2F("Pt_Y_DeltaMi","Pt_Y_DeltaMi", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYDeltaMi->Sumw2();
  fOutputContainer->Add(fHistPtYDeltaMi);

  fHistPtYDelta0              		= new TH2F("Pt_Y_Delta0","Pt_Y_Delta0", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYDelta0->Sumw2();
  fOutputContainer->Add(fHistPtYDelta0);

  fHistPtYLambda              		= new TH2F("Pt_Y_Lambda","Pt_Y_Lambda", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYLambda->Sumw2();
  fOutputContainer->Add(fHistPtYLambda);

  fHistPtYKPl              		= new TH2F("Pt_Y_KPl","Pt_Y_KPl", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYKPl->Sumw2();
  fOutputContainer->Add(fHistPtYKPl);

  fHistPtYKMi              		= new TH2F("Pt_Y_KMi","Pt_Y_KMi", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYKMi->Sumw2();
  fOutputContainer->Add(fHistPtYKMi);

  fHistPtYGamma              		= new TH2F("Pt_Y_Gamma","Pt_Y_Gamma", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYGamma->Sumw2();
  fOutputContainer->Add(fHistPtYGamma);

  fHistPtYDirGamma              		= new TH2F("Pt_Y_GammaDir","Pt_Y_GammaDir", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYDirGamma->Sumw2();
  fOutputContainer->Add(fHistPtYDirGamma);

  fHistPtYPi0FromEta       		    = new TH2F("Pt_Y_Pi0FromEta","Pt_Y_Pi0FromEta", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYPi0FromEta->Sumw2();
  fOutputContainer->Add(fHistPtYPi0FromEta);

  fHistPtYPi0FromLambda       		= new TH2F("Pt_Y_Pi0FromLambda","Pt_Y_Pi0FromLambda", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYPi0FromLambda->Sumw2();
  fOutputContainer->Add(fHistPtYPi0FromLambda);

  fHistPtYPi0FromK            		= new TH2F("Pt_Y_Pi0FromK","Pt_Y_Pi0FromK", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYPi0FromK->Sumw2();
  fOutputContainer->Add(fHistPtYPi0FromK);

  fHistPtYPiPlFromK           		= new TH2F("Pt_Y_PiPlFromK","Pt_Y_PiPlFromK", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYPiPlFromK->Sumw2();
  fOutputContainer->Add(fHistPtYPiPlFromK);

  fHistPtYPiMiFromK           		= new TH2F("Pt_Y_PiMiFromK","Pt_Y_PiMiFromK", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYPiMiFromK->Sumw2();
  fOutputContainer->Add(fHistPtYPiMiFromK);


  fHistPtYPi0GG               		= new TH2F("Pt_Y_Pi0GG","Pt_Y_Pi0GG", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYPi0GG->Sumw2();
  fOutputContainer->Add(fHistPtYPi0GG);
  fHistPtYPi0GGPCMAcc         		= new TH2F("Pt_Y_Pi0GGPCMAcc","Pt_Y_Pi0GGPCMAcc", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYPi0GGPCMAcc->Sumw2();
  fOutputContainer->Add(fHistPtYPi0GGPCMAcc);
  fHistPtYPi0GGEMCAcc         		= new TH2F("Pt_Y_Pi0GGEMCAcc","Pt_Y_Pi0GGEMCAcc", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYPi0GGEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtYPi0GGEMCAcc);
  fHistPtYPi0GGPHOAcc         		= new TH2F("Pt_Y_Pi0GGPHOAcc","Pt_Y_Pi0GGPHOAcc", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYPi0GGPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtYPi0GGPHOAcc);
  fHistPtYPi0GGPCMEMCAcc      		= new TH2F("Pt_Y_Pi0GGPCMEMCAcc","Pt_Y_Pi0GGPCMEMCAcc", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYPi0GGPCMEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtYPi0GGPCMEMCAcc);
  fHistPtYPi0GGPCMPHOAcc      		= new TH2F("Pt_Y_Pi0GGPCMPHOAcc","Pt_Y_Pi0GGPCMPHOAcc", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYPi0GGPCMPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtYPi0GGPCMPHOAcc);

  fHistPtYEtaGG               		= new TH2F("Pt_Y_EtaGG","Pt_Y_EtaGG", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYEtaGG->Sumw2();
  fOutputContainer->Add(fHistPtYEtaGG);
  fHistPtYEtaGGPCMAcc         		= new TH2F("Pt_Y_EtaGGPCMAcc","Pt_Y_EtaGGPCMAcc", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYEtaGGPCMAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaGGPCMAcc);
  fHistPtYEtaGGEMCAcc         		= new TH2F("Pt_Y_EtaGGEMCAcc","Pt_Y_EtaGGEMCAcc", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYEtaGGEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaGGEMCAcc);
  fHistPtYEtaGGPHOAcc         		= new TH2F("Pt_Y_EtaGGPHOAcc","Pt_Y_EtaGGPHOAcc", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYEtaGGPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaGGPHOAcc);
  fHistPtYEtaGGPCMEMCAcc      		= new TH2F("Pt_Y_EtaGGPCMEMCAcc","Pt_Y_EtaGGPCMEMCAcc", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYEtaGGPCMEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaGGPCMEMCAcc);
  fHistPtYEtaGGPCMPHOAcc      		= new TH2F("Pt_Y_EtaGGPCMPHOAcc","Pt_Y_EtaGGPCMPHOAcc", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYEtaGGPCMPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaGGPCMPHOAcc);

  fHistPtYEtaPrimeGG           		= new TH2F("Pt_Y_EtaPrimeGG","Pt_Y_EtaPrimeGG", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYEtaPrimeGG->Sumw2();
  fOutputContainer->Add(fHistPtYEtaPrimeGG);
  fHistPtYEtaPrimeGGPCMAcc     		= new TH2F("Pt_Y_EtaPrimeGGPCMAcc","Pt_Y_EtaPrimeGGPCMAcc", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYEtaPrimeGGPCMAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaPrimeGGPCMAcc);
  fHistPtYEtaPrimeGGEMCAcc     		= new TH2F("Pt_Y_EtaPrimeGGEMCAcc","Pt_Y_EtaPrimeGGEMCAcc", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYEtaPrimeGGEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaPrimeGGEMCAcc);
  fHistPtYEtaPrimeGGPHOAcc     		= new TH2F("Pt_Y_EtaPrimeGGPHOAcc","Pt_Y_EtaPrimeGGPHOAcc", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYEtaPrimeGGPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaPrimeGGPHOAcc);
  fHistPtYEtaPrimeGGPCMEMCAcc  		= new TH2F("Pt_Y_EtaPrimeGGPCMEMCAcc","Pt_Y_EtaPrimeGGPCMEMCAcc", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYEtaPrimeGGPCMEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaPrimeGGPCMEMCAcc);
  fHistPtYEtaPrimeGGPCMPHOAcc  		= new TH2F("Pt_Y_EtaPrimeGGPCMPHOAcc","Pt_Y_EtaPrimeGGPCMPHOAcc", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
  fHistPtYEtaPrimeGGPCMPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaPrimeGGPCMPHOAcc);

  fHistPtAlphaPi0GGPCMAcc     		= new TH2F("Pt_Alpha_Pi0GGPCMAcc","Pt_Alpha_Pi0GGPCMAcc", fMaxpT*5, 0.1, fMaxpT/2, 100, 0., 1.);
  SetLogBinningXTH2(fHistPtAlphaPi0GGPCMAcc);
  fHistPtAlphaPi0GGPCMAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaPi0GGPCMAcc);
  fHistPtAlphaPi0GGEMCAcc     		= new TH2F("Pt_Alpha_Pi0GGEMCAcc","Pt_Alpha_Pi0GGEMCAcc", fMaxpT*5, 0.1, fMaxpT/2, 100, 0., 1.);
  SetLogBinningXTH2(fHistPtAlphaPi0GGEMCAcc);
  fHistPtAlphaPi0GGEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaPi0GGEMCAcc);
  fHistPtAlphaPi0GGPHOAcc     		= new TH2F("Pt_Alpha_Pi0GGPHOAcc","Pt_Alpha_Pi0GGPHOAcc", fMaxpT*5, 0.1, fMaxpT/2, 100, 0., 1.);
  SetLogBinningXTH2(fHistPtAlphaPi0GGPHOAcc);
  fHistPtAlphaPi0GGPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaPi0GGPHOAcc);
  fHistPtAlphaPi0GGPCMEMCAcc  		= new TH2F("Pt_Alpha_Pi0GGPCMEMCAcc","Pt_Alpha_Pi0GGPCMEMCAcc", fMaxpT*5, 0.1, fMaxpT/2, 200, -1., 1.);
  SetLogBinningXTH2(fHistPtAlphaPi0GGPCMEMCAcc);
  fHistPtAlphaPi0GGPCMEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaPi0GGPCMEMCAcc);
  fHistPtAlphaPi0GGPCMPHOAcc  		= new TH2F("Pt_Alpha_Pi0GGPCMPHOAcc","Pt_Alpha_Pi0GGPCMPHOAcc", fMaxpT*5, 0.1, fMaxpT/2, 200, -1., 1.);
  SetLogBinningXTH2(fHistPtAlphaPi0GGPCMPHOAcc);
  fHistPtAlphaPi0GGPCMPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaPi0GGPCMPHOAcc);

  fHistPtAlphaEtaGGPCMAcc     		= new TH2F("Pt_Alpha_EtaGGPCMAcc","Pt_Alpha_EtaGGPCMAcc", fMaxpT*5, 0.1, fMaxpT/2, 100, 0., 1.);
  SetLogBinningXTH2(fHistPtAlphaEtaGGPCMAcc);
  fHistPtAlphaEtaGGPCMAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaEtaGGPCMAcc);
  fHistPtAlphaEtaGGEMCAcc     		= new TH2F("Pt_Alpha_EtaGGEMCAcc","Pt_Alpha_EtaGGEMCAcc", fMaxpT*5, 0.1, fMaxpT/2, 100, 0., 1.);
  SetLogBinningXTH2(fHistPtAlphaEtaGGEMCAcc);
  fHistPtAlphaEtaGGEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaEtaGGEMCAcc);
  fHistPtAlphaEtaGGPHOAcc     		= new TH2F("Pt_Alpha_EtaGGPHOAcc","Pt_Alpha_EtaGGPHOAcc", fMaxpT*5, 0.1, fMaxpT/2, 100, 0., 1.);
  SetLogBinningXTH2(fHistPtAlphaEtaGGPHOAcc);
  fHistPtAlphaEtaGGPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaEtaGGPHOAcc);
  fHistPtAlphaEtaGGPCMEMCAcc  		= new TH2F("Pt_Alpha_EtaGGPCMEMCAcc","Pt_Alpha_EtaGGPCMEMCAcc", fMaxpT*5, 0.1, fMaxpT/2, 200, -1., 1.);
  SetLogBinningXTH2(fHistPtAlphaEtaGGPCMEMCAcc);
  fHistPtAlphaEtaGGPCMEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaEtaGGPCMEMCAcc);
  fHistPtAlphaEtaGGPCMPHOAcc  		= new TH2F("Pt_Alpha_EtaGGPCMPHOAcc","Pt_Alpha_EtaGGPCMPHOAcc", fMaxpT*5, 0.1, fMaxpT/2, 200, -1., 1.);
  SetLogBinningXTH2(fHistPtAlphaEtaGGPCMPHOAcc);
  fHistPtAlphaEtaGGPCMPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaEtaGGPCMPHOAcc);

  fHistPtAlphaEtaPrimeGGPCMAcc     		= new TH2F("Pt_Alpha_EtaPrimeGGPCMAcc","Pt_Alpha_EtaPrimeGGPCMAcc", fMaxpT*5, 0.1, fMaxpT/2, 100, 0., 1.);
  SetLogBinningXTH2(fHistPtAlphaEtaPrimeGGPCMAcc);
  fHistPtAlphaEtaPrimeGGPCMAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaEtaPrimeGGPCMAcc);
  fHistPtAlphaEtaPrimeGGEMCAcc     		= new TH2F("Pt_Alpha_EtaPrimeGGEMCAcc","Pt_Alpha_EtaPrimeGGEMCAcc", fMaxpT*5, 0.1, fMaxpT/2, 100, 0., 1.);
  SetLogBinningXTH2(fHistPtAlphaEtaPrimeGGEMCAcc);
  fHistPtAlphaEtaPrimeGGEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaEtaPrimeGGEMCAcc);
  fHistPtAlphaEtaPrimeGGPHOAcc     		= new TH2F("Pt_Alpha_EtaPrimeGGPHOAcc","Pt_Alpha_EtaPrimeGGPHOAcc", fMaxpT*5, 0.1, fMaxpT/2, 100, 0., 1.);
  SetLogBinningXTH2(fHistPtAlphaEtaPrimeGGPHOAcc);
  fHistPtAlphaEtaPrimeGGPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaEtaPrimeGGPHOAcc);
  fHistPtAlphaEtaPrimeGGPCMEMCAcc  		= new TH2F("Pt_Alpha_EtaPrimeGGPCMEMCAcc","Pt_Alpha_EtaPrimeGGPCMEMCAcc", fMaxpT*5, 0.1, fMaxpT/2, 200, -1., 1.);
  SetLogBinningXTH2(fHistPtAlphaEtaPrimeGGPCMEMCAcc);
  fHistPtAlphaEtaPrimeGGPCMEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaEtaPrimeGGPCMEMCAcc);
  fHistPtAlphaEtaPrimeGGPCMPHOAcc  		= new TH2F("Pt_Alpha_EtaPrimeGGPCMPHOAcc","Pt_Alpha_EtaPrimeGGPCMPHOAcc", fMaxpT*5, 0.1, fMaxpT/2, 200, -1., 1.);
  SetLogBinningXTH2(fHistPtAlphaEtaPrimeGGPCMPHOAcc);
  fHistPtAlphaEtaPrimeGGPCMPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaEtaPrimeGGPCMPHOAcc);

  if (fIsK0 == 1){
        fHistPtYPi0FromKGG          		= new TH2F("Pt_Y_Pi0FromKGG","Pt_Y_Pi0FromKGG", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
        fHistPtYPi0FromKGG->Sumw2();
        fOutputContainer->Add(fHistPtYPi0FromKGG);
        fHistPtYPi0FromKGGPCMAcc    		= new TH2F("Pt_Y_Pi0FromKGGPCMAcc","Pt_Y_Pi0FromKGGPCMAcc", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
        fHistPtYPi0FromKGGPCMAcc->Sumw2();
        fOutputContainer->Add(fHistPtYPi0FromKGGPCMAcc);
        fHistPtYPi0FromKGGEMCAcc    		= new TH2F("Pt_Y_Pi0FromKGGEMCAcc","Pt_Y_Pi0FromKGGEMCAcc", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
        fHistPtYPi0FromKGGEMCAcc->Sumw2();
        fOutputContainer->Add(fHistPtYPi0FromKGGEMCAcc);
        fHistPtYPi0FromKGGPCMEMCAcc 		= new TH2F("Pt_Y_Pi0FromKGGPCMEMCAcc","Pt_Y_Pi0FromKGGPCMEMCAcc", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
        fHistPtYPi0FromKGGPCMEMCAcc->Sumw2();
        fOutputContainer->Add(fHistPtYPi0FromKGGPCMEMCAcc);
        fHistPtYPi0FromKGGEMCPCMAcc 		= new TH2F("Pt_Y_Pi0FromKGGEMCPCMAcc","Pt_Y_Pi0FromKGGEMCPCMAcc", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
        fHistPtYPi0FromKGGEMCPCMAcc->Sumw2();
        fOutputContainer->Add(fHistPtYPi0FromKGGEMCPCMAcc);
        fHistPtYPi0FromKGGEMCAccSamePi0    = new TH2F("Pt_Y_Pi0FromKGGEMCAccSamePi0","Pt_Y_Pi0FromKGGEMCAccSamePi0", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
        fHistPtYPi0FromKGGEMCAccSamePi0->Sumw2();
        fOutputContainer->Add(fHistPtYPi0FromKGGEMCAccSamePi0);
        fHistPtYPi0FromKGGEMCAccDiffPi0  	= new TH2F("Pt_Y_Pi0FromKGGEMCAccDiffPi0","Pt_Y_Pi0FromKGGEMCAccDiffPi0", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
        fHistPtYPi0FromKGGEMCAccDiffPi0->Sumw2();
        fOutputContainer->Add(fHistPtYPi0FromKGGEMCAccDiffPi0);
        fHistPtAlphaPi0FromKGG              = new TH2F("Pt_Alpha_Pi0FromKGG","Pt_Alpha_Pi0FromKGG", fMaxpT*5, 0.1, fMaxpT/2, 200, -1., 1.);
        SetLogBinningXTH2(fHistPtAlphaPi0FromKGG);
        fHistPtAlphaPi0FromKGG->Sumw2();
        fOutputContainer->Add(fHistPtAlphaPi0FromKGG);
        fHistPtAlphaPi0FromKGGPCMAcc        = new TH2F("Pt_Alpha_Pi0FromKGGPCMAcc","Pt_Alpha_Pi0FromKGGPCMAcc", fMaxpT*5, 0.1, fMaxpT/2, 200, -1., 1.);
        SetLogBinningXTH2(fHistPtAlphaPi0FromKGGPCMAcc);
        fHistPtAlphaPi0FromKGGPCMAcc->Sumw2();
        fOutputContainer->Add(fHistPtAlphaPi0FromKGGPCMAcc);
        fHistPtAlphaPi0FromKGGEMCAcc      = new TH2F("Pt_Alpha_Pi0FromKGGEMCAcc","Pt_Alpha_Pi0FromKGGEMCAcc", fMaxpT*5, 0.1, fMaxpT/2, 200, -1., 1.);
        SetLogBinningXTH2(fHistPtAlphaPi0FromKGGEMCAcc);
        fHistPtAlphaPi0FromKGGEMCAcc->Sumw2();
        fOutputContainer->Add(fHistPtAlphaPi0FromKGGEMCAcc);
        fHistPtAlphaPi0FromKGGPCMEMCAcc      = new TH2F("Pt_Alpha_Pi0FromKGGPCMEMCAcc","Pt_Alpha_Pi0FromKGGPCMEMCAcc", fMaxpT*5, 0.1, fMaxpT/2, 200, -1., 1.);
        SetLogBinningXTH2(fHistPtAlphaPi0FromKGGPCMEMCAcc);
        fHistPtAlphaPi0FromKGGPCMEMCAcc->Sumw2();
        fOutputContainer->Add(fHistPtAlphaPi0FromKGGPCMEMCAcc);
        fHistPtAlphaPi0FromKGGEMCPCMAcc      = new TH2F("Pt_Alpha_Pi0FromKGGEMCPCMAcc","Pt_Alpha_Pi0FromKGGEMCPCMAcc", fMaxpT*5, 0.1, fMaxpT/2, 200, -1., 1.);
    	  SetLogBinningXTH2(fHistPtAlphaPi0FromKGGEMCPCMAcc);
        fHistPtAlphaPi0FromKGGEMCPCMAcc->Sumw2();
        fOutputContainer->Add(fHistPtAlphaPi0FromKGGEMCPCMAcc);
        fHistPtAlphaPi0FromKGGEMCAccSamePi0      = new TH2F("Pt_Alpha_Pi0FromKGGEMCAccSamePi0","Pt_Alpha_Pi0FromKGGEMCAccSamePi0", fMaxpT*5, 0.1, fMaxpT/2, 200, -1., 1.);
    	  SetLogBinningXTH2(fHistPtAlphaPi0FromKGGEMCAccSamePi0);
        fHistPtAlphaPi0FromKGGEMCAccSamePi0->Sumw2();
        fOutputContainer->Add(fHistPtAlphaPi0FromKGGEMCAccSamePi0);
        fHistPtAlphaPi0FromKGGEMCAccDiffPi0      = new TH2F("Pt_Alpha_Pi0FromKGGEMCAccDiffPi0","Pt_Alpha_Pi0FromKGGEMCAccDiffPi0", fMaxpT*5, 0.1, fMaxpT/2, 200, -1., 1.);
    	  SetLogBinningXTH2(fHistPtAlphaPi0FromKGGEMCAccDiffPi0);
        fHistPtAlphaPi0FromKGGEMCAccDiffPi0->Sumw2();
        fOutputContainer->Add(fHistPtAlphaPi0FromKGGEMCAccDiffPi0);
  }

  if(fDoMultStudies){
    fHistV0Mult = new TH1D("V0Multiplicity", "V0Multiplicity", 1000, -0.5, 1000 - 0.5);
    fHistV0Mult->Sumw2();
    fOutputContainer->Add(fHistV0Mult);

    fHistPtV0MultPi0GG      = new TH2F("Pt_V0Mult_Pi0FromGG","Pt_V0Mult_Pi0FromGG", fMaxpT*10, 0., fMaxpT, 500, -0.5, 500 - 0.5);
    fHistPtV0MultPi0GG->Sumw2();
    fOutputContainer->Add(fHistPtV0MultPi0GG);

    fHistPtV0MultEtaGG      = new TH2F("Pt_V0Mult_EtaFromGG","Pt_V0Mult_EtaFromGG", fMaxpT*10, 0., fMaxpT, 500, -0.5, 500 - 0.5);
    fHistPtV0MultEtaGG->Sumw2();
    fOutputContainer->Add(fHistPtV0MultEtaGG);

    fHistPtV0MultEtaPrimeGG      = new TH2F("Pt_V0Mult_EtaPrimeFromGG","Pt_V0Mult_EtaPrimeFromGG", fMaxpT*10, 0., fMaxpT/2, 500, -0.5, 500 - 0.5);
    fHistPtV0MultEtaPrimeGG->Sumw2();
    fOutputContainer->Add(fHistPtV0MultEtaPrimeGG);

    // investigate prompt and feed down pi0
    fHistPtV0MultPi0GGPrompt      = new TH2F("Pt_V0Mult_Pi0FromGG_Prompt","Pt_V0Mult_Pi0FromGG_Prompt", fMaxpT*10, 0., fMaxpT, 500, -0.5, 500 - 0.5);
    fHistPtV0MultPi0GGPrompt->Sumw2();
    fOutputContainer->Add(fHistPtV0MultPi0GGPrompt);

    fHistPtV0MultPi0GGFromEta      = new TH2F("Pt_V0Mult_Pi0FromGG_FromEta","Pt_V0Mult_Pi0FromGG_FromEta", fMaxpT*10, 0., fMaxpT, 500, -0.5, 500 - 0.5);
    fHistPtV0MultPi0GGFromEta->Sumw2();
    fOutputContainer->Add(fHistPtV0MultPi0GGFromEta);

    fHistPtV0MultPi0GGFromOmega      = new TH2F("Pt_V0Mult_Pi0FromGG_FromOmega","Pt_V0Mult_Pi0FromGG_FromOmega", fMaxpT*10, 0., fMaxpT, 500, -0.5, 500 - 0.5);
    fHistPtV0MultPi0GGFromOmega->Sumw2();
    fOutputContainer->Add(fHistPtV0MultPi0GGFromOmega);

    fHistPtV0MultPi0GGFromRest      = new TH2F("Pt_V0Mult_Pi0FromGG_FromRest","Pt_V0Mult_Pi0FromGG_FromRest", fMaxpT*10, 0., fMaxpT, 500, -0.5, 500 - 0.5);
    fHistPtV0MultPi0GGFromRest->Sumw2();
    fOutputContainer->Add(fHistPtV0MultPi0GGFromRest);

    fHistPtV0MultPi0GGFromRho      = new TH2F("Pt_V0Mult_Pi0FromGG_FromRho","Pt_V0Mult_Pi0FromGG_FromRho", fMaxpT*10, 0., fMaxpT, 500, -0.5, 500 - 0.5);
    fHistPtV0MultPi0GGFromRho->Sumw2();
    fOutputContainer->Add(fHistPtV0MultPi0GGFromRho);

    // investigate prompt and feed down eta
    fHistPtV0MultEtaGGPrompt      = new TH2F("Pt_V0Mult_EtaFromGG_Prompt","Pt_V0Mult_EtaFromGG_Prompt", fMaxpT*10, 0., fMaxpT, 500, -0.5, 500 - 0.5);
    fHistPtV0MultEtaGGPrompt->Sumw2();
    fOutputContainer->Add(fHistPtV0MultEtaGGPrompt);

    fHistPtV0MultEtaGGFromEtaPrim      = new TH2F("Pt_V0Mult_EtaFromGG_FromEtaPrim","Pt_V0Mult_EtaFromGG_FromEtaPrim", fMaxpT*10, 0., fMaxpT, 500, -0.5, 500 - 0.5);
    fHistPtV0MultEtaGGFromEtaPrim->Sumw2();
    fOutputContainer->Add(fHistPtV0MultEtaGGFromEtaPrim);

    fHistPtV0MultEtaGGFromRest      = new TH2F("Pt_V0Mult_EtaFromGG_FromRest","Pt_V0Mult_EtaFromGG_FromRest", fMaxpT*10, 0., fMaxpT, 500, -0.5, 500 - 0.5);
    fHistPtV0MultEtaGGFromRest->Sumw2();
    fOutputContainer->Add(fHistPtV0MultEtaGGFromRest);

    fHistPtV0MultGamma      = new TH2F("Pt_V0Mult_Gamma","Pt_V0Mult_Gamma", fMaxpT*10, 0., fMaxpT, 500, -0.5, 500 - 0.5);
    fHistPtV0MultGamma->Sumw2();
    fOutputContainer->Add(fHistPtV0MultGamma);

    fHistPtV0MultDirGamma      = new TH2F("Pt_V0Mult_GammaDir","Pt_V0Mult_GammaDir", fMaxpT*10, 0., fMaxpT, 500, -0.5, 500 - 0.5);
    fHistPtV0MultDirGamma->Sumw2();
    fOutputContainer->Add(fHistPtV0MultDirGamma);

    fHistPtV0MultChargedPi   = new TH2F("Pt_V0Mult_ChargedPi","Pt_V0Mult_ChargedPi", fMaxpT*10, 0., fMaxpT, 500, -0.5, 500 - 0.5);
    fHistPtV0MultChargedPi->Sumw2();
    fOutputContainer->Add(fHistPtV0MultChargedPi);

  }

  if(fDoFeedDownStudies){
    fHistPtYPi0Primordial = new TH2F("Pt_Y_Pi0Primordial","Pt_Y_Pi0Primordial", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
    fHistPtYPi0Primordial->Sumw2();
    fOutputContainer->Add(fHistPtYPi0Primordial);
    fHistPtYEtaPrimordial = new TH2F("Pt_Y_EtaPrimordial","Pt_Y_EtaPrimordial", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
    fHistPtYEtaPrimordial->Sumw2();
    fOutputContainer->Add(fHistPtYEtaPrimordial);
    fHistPtYOmegaPrimordial = new TH2F("Pt_Y_OmegaPrimordial","Pt_Y_OmegaPrimordial", fMaxpT*10, 0, fMaxpT, 200, -1.0, 1.0);
    fHistPtYOmegaPrimordial->Sumw2();
    fOutputContainer->Add(fHistPtYOmegaPrimordial);

    const Int_t NFeedDownMothers = 20;
    TString FeedDownMotherNames[NFeedDownMothers] = {"u","d","s","c","b","t","qq","g","#rho","#omega","#eta","#eta'","K^{*}","#Delta","#Lambda","#Sigma","#phi","D","#Xi^{0}","Other"};

    fHistPtPi0FromDecay = new TH2F("Pt_Pi0FromDecay","Pt_Pi0FromDecay",20, 0.5, 20.5, fMaxpT*10, 0, fMaxpT);
    for (Int_t iMother = 0; iMother < NFeedDownMothers; iMother++)
      fHistPtPi0FromDecay->GetXaxis()->SetBinLabel(iMother+1,FeedDownMotherNames[iMother]);
    fHistPtPi0FromDecay->Sumw2();
    fOutputContainer->Add(fHistPtPi0FromDecay);

    fHistPtEtaFromDecay = new TH2F("Pt_EtaFromDecay","Pt_EtaFromDecay",20, 0.5, 20.5, fMaxpT*10, 0, fMaxpT);
    for (Int_t iMother = 0; iMother < NFeedDownMothers; iMother++)
      fHistPtEtaFromDecay->GetXaxis()->SetBinLabel(iMother+1,FeedDownMotherNames[iMother]);
    fHistPtEtaFromDecay->Sumw2();
    fOutputContainer->Add(fHistPtEtaFromDecay);

    fHistPtOmegaFromDecay = new TH2F("Pt_OmegaFromDecay","Pt_OmegaFromDecay",20, 0.5, 20.5, fMaxpT*10, 0, fMaxpT);
    for (Int_t iMother = 0; iMother < NFeedDownMothers; iMother++)
      fHistPtOmegaFromDecay->GetXaxis()->SetBinLabel(iMother+1,FeedDownMotherNames[iMother]);
    fHistPtOmegaFromDecay->Sumw2();
    fOutputContainer->Add(fHistPtOmegaFromDecay);
  }

  if(fDoJetStudies){
    std::vector<double> vecJetPt;
    double jetPt = 0.;
    double maxJetPt = 500.;
    double epsilon = 1e-20;
    while(jetPt < maxJetPt){
      vecJetPt.push_back(jetPt);
      if(jetPt < 5. - epsilon) jetPt +=1.;
      else if(jetPt < 50. - epsilon) jetPt +=5.;
      else if(jetPt < 100. - epsilon) jetPt +=10.;
      else if(jetPt < 200. - epsilon) jetPt +=50.;
      else if(jetPt < 500. - epsilon) jetPt +=100.;
      else vecJetPt.push_back(maxJetPt);
    }
    
    fHistJetPtY = new TH2F("JetPtY", "JetPtY", vecJetPt.size()-1, vecJetPt.data(), 100, -2, 2);
    fHistJetPtY->Sumw2();
    fOutputContainer->Add(fHistJetPtY);
    fHistJetEta = new TH1D("JetEta", "JetEta", 100, -1., 1.);
    fHistJetEta->Sumw2();
    fOutputContainer->Add(fHistJetEta);
    fHistJetPhi = new TH1D("JetPhi", "JetPhi", 100, 0., 6.5);
    fHistJetPhi->Sumw2();
    fOutputContainer->Add(fHistJetPhi);

    fHistPi0PtJetPt = new TH2F("Pi0PtVsJetPt_Eta08", "Pi0PtVsJetPt_Eta08", fMaxpT*10, 0., fMaxpT, vecJetPt.size()-1, vecJetPt.data());
    fHistPi0PtJetPt->Sumw2();
    fOutputContainer->Add(fHistPi0PtJetPt);

    fHistEtaPtJetPt = new TH2F("EtaPtVsJetPt_Eta08", "EtaPtVsJetPt_Eta08", fMaxpT*10, 0., fMaxpT, vecJetPt.size()-1, vecJetPt.data());
    fHistEtaPtJetPt->Sumw2();
    fOutputContainer->Add(fHistEtaPtJetPt);

    fHistPi0ZJetPt = new TH2F("Pi0ZVsJetPt_Eta08", "Pi0ZVsJetPt_Eta08", 100 ,0, 1, vecJetPt.size()-1, vecJetPt.data());
    fHistPi0ZJetPt->Sumw2();
    fOutputContainer->Add(fHistPi0ZJetPt);

    fHistEtaZJetPt = new TH2F("EtaZVsJetPt_Eta08", "EtaZVsJetPt_Eta08", 100 ,0, 1, vecJetPt.size()-1, vecJetPt.data());
    fHistEtaZJetPt->Sumw2();
    fOutputContainer->Add(fHistEtaZJetPt);
  }

  PostData(1, fOutputContainer);
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaPureMC::UserExec(Option_t *)
{

  fInputEvent = InputEvent();
    // std::cout << "I found an Event" << std::endl;

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

  if(fDoMultStudies)
  {
    ProcessMultiplicity();
  }

  if(fDoJetStudies){
    ProcessJets();
  }

  ProcessMCParticles();


  PostData(1, fOutputContainer);
}

void AliAnalysisTaskGammaPureMC::ProcessMultiplicity()
{
  // set number of tracks in V0 acceptance to 0
  fNTracksInV0Acc = 0;
  // set INEL>0 to false
  fIsEvtINELgtZERO = false;

  // Loop over all primary MC particle
  for(Long_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {
    AliVParticle* particle     = nullptr;
    particle                    = (AliVParticle *)fMCEvent->GetTrack(i);
    if (!particle) continue;

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
  if(fDoMultStudies == 2){
    if(fIsEvtINELgtZERO == true){
      fHistV0Mult->Fill(fNTracksInV0Acc);
    }
  } else {
    fHistV0Mult->Fill(fNTracksInV0Acc);
  }
}


int AliAnalysisTaskGammaPureMC::ReturnFeedDownBinFromPDG(int pdgcode) {
  switch (TMath::Abs(pdgcode))
  {
  case 2:
    return 1; // u
  case 1:
    return 2; // d
  case 3:
    return 3; // s
  case 4:
    return 4; // c
  case 5:
    return 5; // b
  case 6:
    return 6; // t
  case 1103:
  case 2101:
  case 2103:
  case 2203:
  case 3101:
  case 3103:
  case 3201:
  case 3203:
  case 3303:
    return 7; // diquark
  case 21:
  case 9:
    return 8; // g
  case 113:
  case 213:
    return 9; // charged rho
  case 223:
    return 10; // omega
  case 221:
    return 11; // eta
  case 331:
    return 12; // eta'
  case 313:
  case 323:
  case 10311:
    return 13; // K* 0 and +
  case 2224:
  case 2214:
  case 2114:
  case 1114:
    return 14; // Delta
  case 3122:
  case 4122: // Lambda+C
    return 15; // Lambda
  case 3222:
  case 3212:
  case 3112:
  case 3214:
  case 3114:
  case 3224:
    return 16; // Sigma
  case 333:
    return 17; // phi
  case 421:
  case 411:
  case 413:
  case 423:
  case 431:
    return 18; // D
  case 3322:
  case 3312:
  case 3324:
  case 3314:
    return 19; // Xi
  default:
    // std::cout << "Unknown PDG code of mother: " << pdgcode << std::endl;
    return 20; // Other
  }
}

//________________________________________________________________________
bool AliAnalysisTaskGammaPureMC::isPhotonFromDecay(int pdgCodeMother){
  switch (std::abs(pdgCodeMother))
  {
    case 111:     // pi0
    case 221:     // eta
    case 223:     // omega
    case 331:     // eta prime
    case 333:     // phi
    case 113:     // rho 0
    case 213:     // rho+-
      return true;
    default:
      return false;
  }

}

//________________________________________________________________________
void AliAnalysisTaskGammaPureMC::ProcessMCParticles()
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

    int absmotherpdg = 0;
    if (hasMother)
      absmotherpdg = TMath::Abs(motherParticle->PdgCode());

    const std::array<int, 20> kAcceptPdgCodes = {kPdgPi0, kPdgEta, kPdgEtaPrime, kPdgOmega, kPdgPiPlus, kPdgRho0, kPdgPhi, kPdgJPsi, kPdgSigma0, kPdgK0Short, kPdgDeltaPlus, kPdgDeltaPlusPlus, kPdgDeltaMinus, kPdgDelta0, kPdgRhoPlus, kPdgKStar, kPdgK0Long, kPdgLambda, kPdgKPlus, kGamma};
    if(std::find(kAcceptPdgCodes.begin(), kAcceptPdgCodes.end(), TMath::Abs(particle->PdgCode())) ==  kAcceptPdgCodes.end()) continue;  // species not supported

    if (!(TMath::Abs(particle->E()-particle->Pz())>0.)) continue;
    Double_t yPre = (particle->E()+particle->Pz())/(particle->E()-particle->Pz());
    //     cout << i << "\t"<< particle->GetPdgCode() << "\t"<< particle->Pz() << "\t" << particle->Energy()<< "\t" << particle->Energy()-particle->Pz() << "\t"<< yPre << endl;
    if( yPre <= 0 ) continue;

    Double_t y = 0.5*TMath::Log(yPre);


    if (y > 1.000) continue;
    switch(particle->PdgCode()){
    case kPdgPi0:
      fHistPtYPi0->Fill(particle->Pt(), particle->Y());
      if (hasMother){
        if (TMath::Abs(motherParticle->PdgCode()) == kPdgK0Short ||
            TMath::Abs(motherParticle->PdgCode()) == kPdgK0Long ||
            TMath::Abs(motherParticle->PdgCode()) == kPdgKPlus
        )
          fHistPtYPi0FromK->Fill(particle->Pt(), particle->Y());
        if (TMath::Abs(motherParticle->PdgCode()) == kPdgLambda)
          fHistPtYPi0FromLambda->Fill(particle->Pt(), particle->Y());
        if (motherParticle->PdgCode() == kPdgEta)
          fHistPtYPi0FromEta->Fill(particle->Pt(), particle->Y());
        if(fDoFeedDownStudies){
          fHistPtPi0FromDecay->Fill(ReturnFeedDownBinFromPDG(absmotherpdg),particle->Pt());
          if(ReturnFeedDownBinFromPDG(absmotherpdg) < 9)
            fHistPtYPi0Primordial->Fill(particle->Pt(), particle->Y());
        }
      }

      // fill primary pi0s in eta > 0.8
      if(fDoMultStudies){
        if(std::abs(particle->Y()) <= 0.8){
          if(!IsSecondary(motherParticle)){
            if((fDoMultStudies == 2 && fIsEvtINELgtZERO == true) || fDoMultStudies == 1){ // select only INEL>0 events for multiplicity
              int sourcePi0 = ReturnFeedDownBinFromPDG(absmotherpdg);
              fHistPtV0MultPi0GG->Fill(particle->Pt(), fNTracksInV0Acc);
              if(sourcePi0 < 9 || (absmotherpdg == 2212 && particle->GetMother() < 2)){ // (absmotherpdg == 2212 && particle->GetMother() < 2) is for EPOS as EPOS does not have quarks
                fHistPtV0MultPi0GGPrompt->Fill(particle->Pt(), fNTracksInV0Acc);
              } else if(sourcePi0 == 9){
                fHistPtV0MultPi0GGFromRho->Fill(particle->Pt(), fNTracksInV0Acc);
              } else if(sourcePi0 == 10){
                fHistPtV0MultPi0GGFromEta->Fill(particle->Pt(), fNTracksInV0Acc);
              } else if(sourcePi0 == 11){
                fHistPtV0MultPi0GGFromOmega->Fill(particle->Pt(), fNTracksInV0Acc);
              } else {
                fHistPtV0MultPi0GGFromRest->Fill(particle->Pt(), fNTracksInV0Acc);
              }                 
            }
          }
        }
      }
      if(fDoJetStudies){
        if(!IsSecondary(motherParticle)){
          int index = -1;
          double R = 0;
          if(IsParticleInJet(fVecJets, particle->Eta(), particle->Phi(), index, R)){
            if(std::abs(particle->Y()) <= fJetParticleAcc){
              fHistPi0PtJetPt->Fill(particle->Pt(), fVecJets[index].pt());
            }
            if(std::abs(particle->Y()) <= fJetParticleAccFF){
              fHistPi0ZJetPt->Fill(GetFrag(fVecJets[index], particle), fVecJets[index].pt());
            }
          }
        }
      }
      break;
    case kPdgEta:
      fHistPtYEta->Fill(particle->Pt(), particle->Y());
      if(fDoFeedDownStudies){
        if (hasMother)
          fHistPtEtaFromDecay->Fill(ReturnFeedDownBinFromPDG(absmotherpdg),particle->Pt());
        if(ReturnFeedDownBinFromPDG(absmotherpdg) < 9)
          fHistPtYEtaPrimordial->Fill(particle->Pt(), particle->Y());
      }
      // fill primary etas in eta > 0.8
      if(fDoMultStudies){
        if(std::abs(particle->Y()) <= 0.8){
          if(!IsSecondary(motherParticle)){
            if((fDoMultStudies == 2 && fIsEvtINELgtZERO == true) || fDoMultStudies == 1){ // select only INEL>0 events for multiplicity
              fHistPtV0MultEtaGG->Fill(particle->Pt(), fNTracksInV0Acc);
              int sourceEta = ReturnFeedDownBinFromPDG(absmotherpdg);
              if(sourceEta < 9 || (absmotherpdg == 2212 && particle->GetMother() < 2)){ // (absmotherpdg == 2212 && particle->GetMother() < 2) is for EPOS as EPOS does not have quarks
                fHistPtV0MultEtaGGPrompt->Fill(particle->Pt(), fNTracksInV0Acc);
              } else if(sourceEta == 12){
                fHistPtV0MultEtaGGFromEtaPrim->Fill(particle->Pt(), fNTracksInV0Acc);
              } else {
                fHistPtV0MultEtaGGFromRest->Fill(particle->Pt(), fNTracksInV0Acc);
              }
            }
          }
        }
      }
      if(fDoJetStudies){
        if(!IsSecondary(motherParticle)){
          int index = -1;
          double R = 0;
          if(IsParticleInJet(fVecJets, particle->Eta(), particle->Phi(), index, R)){
            if(std::abs(particle->Y()) <= fJetParticleAcc){
              fHistEtaPtJetPt->Fill(particle->Pt(), fVecJets[index].pt());
            }
            if(std::abs(particle->Y()) <= fJetParticleAccFF){
              fHistEtaZJetPt->Fill(GetFrag(fVecJets[index], particle), fVecJets[index].pt());
            }
          }
        }
      }
      break;
    case kPdgEtaPrime:
      fHistPtYEtaPrime->Fill(particle->Pt(), particle->Y());
      if(fDoMultStudies){
        if(std::abs(particle->Y()) <= 0.8){
          if(!IsSecondary(motherParticle)){
            if(fDoMultStudies == 2){ // select only INEL>0 events for multiplicity
              if(fIsEvtINELgtZERO == true){
                fHistPtV0MultEtaPrimeGG->Fill(particle->Pt(), fNTracksInV0Acc);
              }
            } else {
              fHistPtV0MultEtaPrimeGG->Fill(particle->Pt(), fNTracksInV0Acc);
            }
          }
        }
      }
      break;
    case kPdgOmega:
      fHistPtYOmega->Fill(particle->Pt(), particle->Y());
      if(fDoFeedDownStudies){
        if (hasMother)
          fHistPtOmegaFromDecay->Fill(ReturnFeedDownBinFromPDG(absmotherpdg),particle->Pt());
        if(ReturnFeedDownBinFromPDG(absmotherpdg) < 9)
          fHistPtYOmegaPrimordial->Fill(particle->Pt(), particle->Y());
      }
      break;
    case kPdgPiPlus:
      fHistPtYPiPl->Fill(particle->Pt(), particle->Y());
      if (hasMother){
        if (TMath::Abs(motherParticle->PdgCode()) == kPdgK0Short ||
            TMath::Abs(motherParticle->PdgCode()) == kPdgK0Long ||
            TMath::Abs(motherParticle->PdgCode()) == kPdgKPlus
        )
          fHistPtYPiPlFromK->Fill(particle->Pt(), particle->Y());
      }
      if(fDoMultStudies){
        if(std::abs(particle->Y()) <= 0.8){
          if(!IsSecondary(motherParticle)){
            if((fDoMultStudies == 2 && fIsEvtINELgtZERO == true) || fDoMultStudies == 1){ // select only INEL>0 events for multiplicity
              fHistPtV0MultChargedPi->Fill(particle->Pt(), fNTracksInV0Acc);
            }
          }
        }
      }
      break;
    case kPdgPiMinus:
      fHistPtYPiMi->Fill(particle->Pt(), particle->Y());
      if (hasMother){
        if (TMath::Abs(motherParticle->PdgCode()) == kPdgK0Short ||
            TMath::Abs(motherParticle->PdgCode()) == kPdgK0Long ||
            TMath::Abs(motherParticle->PdgCode()) == kPdgKPlus
        )
          fHistPtYPiMiFromK->Fill(particle->Pt(), particle->Y());
      }
      if(fDoMultStudies){
        if(std::abs(particle->Y()) <= 0.8){
          if(!IsSecondary(motherParticle)){
            if((fDoMultStudies == 2 && fIsEvtINELgtZERO == true) || fDoMultStudies == 1){ // select only INEL>0 events for multiplicity
              fHistPtV0MultChargedPi->Fill(particle->Pt(), fNTracksInV0Acc);
            }
          }
        }
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
    case kGamma:
      bool isDecayPhoton = isPhotonFromDecay(absmotherpdg);
      fHistPtYGamma->Fill(particle->Pt(), particle->Y());
      if(!isDecayPhoton){ // only direct photons
        fHistPtYDirGamma->Fill(particle->Pt(), particle->Y());
      }
      if(fDoMultStudies){
        if(std::abs(particle->Y()) <= 0.8){
          // For gammas, also check grandmother for secondary
          bool isSecondaryGrandmother = false;
          if (motherParticle->GetMother()>-1){
            hasMother                 = kTRUE;
            AliVParticle* grandMotherParticle = (AliVParticle *)fMCEvent->GetTrack(motherParticle->GetMother());
            isSecondaryGrandmother = IsSecondary(grandMotherParticle);
          }
          if(!IsSecondary(motherParticle) && !isSecondaryGrandmother){
            if((fDoMultStudies == 2 && fIsEvtINELgtZERO == true) || fDoMultStudies == 1){ // select only INEL>0 events for multiplicity
              fHistPtV0MultGamma->Fill(particle->Pt(), fNTracksInV0Acc);
              bool isDecayPhoton = isPhotonFromDecay(absmotherpdg);
              if(!isDecayPhoton){ // only direct photons
                fHistPtV0MultDirGamma->Fill(particle->Pt(), fNTracksInV0Acc);
              }                  
            }
          }
        }
      }
    }

    // from here on, we are only intested in particles considered primaries in ALICE
    if ((particle->PdgCode()== kPdgPi0 || particle->PdgCode()== kPdgEta) && hasMother){
      if (TMath::Abs(motherParticle->PdgCode()) == kPdgK0Short ||
          TMath::Abs(motherParticle->PdgCode()) == kPdgK0Long  ||
          TMath::Abs(motherParticle->PdgCode()) == kPdgKPlus  ||
          TMath::Abs(motherParticle->PdgCode()) == kPdgLambda ||
          TMath::Abs(motherParticle->PdgCode()) == kPdgSigma0 ||
          TMath::Abs(motherParticle->PdgCode()) == kPdgSigmaPlus ||
          TMath::Abs(motherParticle->PdgCode()) == kPdgSigmaMinus ||
          TMath::Abs(motherParticle->PdgCode()) == kPdgXi0 ||
          TMath::Abs(motherParticle->PdgCode()) == kPdgXiMinus
      )
        continue;
    }

    // just looking at pi0, etas, etaprims
    if (particle->PdgCode()==kPdgPi0 || particle->PdgCode()==kPdgEta || particle->PdgCode() == kPdgEtaPrime){
      if (particle->GetNDaughters() != 2) continue;   // only the two particle decays
      UChar_t acceptanceGamma[2] = {0,0};
      Double_t energyGamma[2] = {0,0};
      Bool_t allOK[2] = {kFALSE,kFALSE};


      for(Int_t j=0;j<2;j++){
        AliVParticle *daughter = (AliVParticle*) fMCEvent->GetTrack(particle->GetDaughterLabel(j));
        if (!daughter) continue;

        // Is Daughter a Photon?
        if(daughter->PdgCode() == 22) allOK[j] =kTRUE;
        if(IsInPCMAcceptance(daughter))  SETBIT(acceptanceGamma[j], kPCMAcceptance);
        if(IsInPHOSAcceptance(daughter)) SETBIT(acceptanceGamma[j], kPHOSAcceptance);
        if(IsInEMCalAcceptance(daughter)) SETBIT(acceptanceGamma[j], kEMCALAcceptance);
        energyGamma[j] = daughter->E();


      }

      if (!(allOK[0] && allOK[1])) continue;

      Double_t alpha = (energyGamma[0]-energyGamma[1])/(energyGamma[0]+energyGamma[1]);

      if (particle->PdgCode()==kPdgPi0){
        fHistPtYPi0GG->Fill(particle->Pt(), particle->Y());
        if (TESTBIT(acceptanceGamma[0], kPCMAcceptance) && TESTBIT(acceptanceGamma[1], kPCMAcceptance)){
          fHistPtYPi0GGPCMAcc->Fill(particle->Pt(), particle->Y());
          fHistPtAlphaPi0GGPCMAcc->Fill(particle->Pt(), TMath::Abs(alpha));
        }
        if (TESTBIT(acceptanceGamma[0], kEMCALAcceptance) && TESTBIT(acceptanceGamma[1], kEMCALAcceptance)){
          fHistPtYPi0GGEMCAcc->Fill(particle->Pt(), particle->Y());
          fHistPtAlphaPi0GGEMCAcc->Fill(particle->Pt(), TMath::Abs(alpha));
        }
        if (TESTBIT(acceptanceGamma[0], kPHOSAcceptance) && TESTBIT(acceptanceGamma[1], kPHOSAcceptance)){
          fHistPtYPi0GGPHOAcc->Fill(particle->Pt(), particle->Y());
          fHistPtAlphaPi0GGPHOAcc->Fill(particle->Pt(), TMath::Abs(alpha));
        }
        if ( (TESTBIT(acceptanceGamma[0], kPCMAcceptance)  && TESTBIT(acceptanceGamma[1], kEMCALAcceptance)) ||
            (TESTBIT(acceptanceGamma[1], kPCMAcceptance) && TESTBIT(acceptanceGamma[0], kPCMAcceptance))
        ){
          fHistPtYPi0GGPCMEMCAcc->Fill(particle->Pt(), particle->Y());
          if (!TESTBIT(acceptanceGamma[1], kEMCALAcceptance)) alpha = (energyGamma[1]-energyGamma[0])/(energyGamma[0]+energyGamma[1]);
          fHistPtAlphaPi0GGPCMEMCAcc->Fill(particle->Pt(), alpha);
        }
        if ( (TESTBIT(acceptanceGamma[0], kPCMAcceptance) && TESTBIT(acceptanceGamma[1], kPHOSAcceptance)) ||
            (TESTBIT(acceptanceGamma[1], kPCMAcceptance) && TESTBIT(acceptanceGamma[0], kPHOSAcceptance))
        ){
          fHistPtYPi0GGPCMPHOAcc->Fill(particle->Pt(), particle->Y());
          if (!TESTBIT(acceptanceGamma[1], kPHOSAcceptance)) alpha = (energyGamma[1]-energyGamma[0])/(energyGamma[0]+energyGamma[1]);
          fHistPtAlphaPi0GGPCMPHOAcc->Fill(particle->Pt(), alpha);
        }
      }
      if (particle->PdgCode()==kPdgEta){
        fHistPtYEtaGG->Fill(particle->Pt(), particle->Y());
        if (TESTBIT(acceptanceGamma[0], kPCMAcceptance) && TESTBIT(acceptanceGamma[1], kPCMAcceptance)){
          fHistPtYEtaGGPCMAcc->Fill(particle->Pt(), particle->Y());
          fHistPtAlphaEtaGGPCMAcc->Fill(particle->Pt(), TMath::Abs(alpha));
        }
        if (TESTBIT(acceptanceGamma[0], kEMCALAcceptance) && TESTBIT(acceptanceGamma[1], kEMCALAcceptance)){
          fHistPtYEtaGGEMCAcc->Fill(particle->Pt(), particle->Y());
          fHistPtAlphaEtaGGEMCAcc->Fill(particle->Pt(), TMath::Abs(alpha));
        }
        if (TESTBIT(acceptanceGamma[0], kPHOSAcceptance) && TESTBIT(acceptanceGamma[1], kPHOSAcceptance)){
          fHistPtYEtaGGPHOAcc->Fill(particle->Pt(), particle->Y());
          fHistPtAlphaEtaGGPHOAcc->Fill(particle->Pt(), TMath::Abs(alpha));
        }
        if ( (TESTBIT(acceptanceGamma[0], kPCMAcceptance) && TESTBIT(acceptanceGamma[1], kEMCALAcceptance)) ||
            (TESTBIT(acceptanceGamma[1], kPCMAcceptance) && TESTBIT(acceptanceGamma[0], kEMCALAcceptance))
        ){
          fHistPtYEtaGGPCMEMCAcc->Fill(particle->Pt(), particle->Y());
          if (!TESTBIT(acceptanceGamma[1], kEMCALAcceptance)) alpha = (energyGamma[1]-energyGamma[0])/(energyGamma[0]+energyGamma[1]);
          fHistPtAlphaEtaGGPCMEMCAcc->Fill(particle->Pt(), alpha);
        }
        if ( (TESTBIT(acceptanceGamma[0], kPCMAcceptance) && TESTBIT(acceptanceGamma[1], kPHOSAcceptance)) ||
            (TESTBIT(acceptanceGamma[1], kPCMAcceptance) && TESTBIT(acceptanceGamma[0], kPHOSAcceptance))
        ){
          fHistPtYEtaGGPCMPHOAcc->Fill(particle->Pt(), particle->Y());
          if (TESTBIT(!acceptanceGamma[1],kPHOSAcceptance)) alpha = (energyGamma[1]-energyGamma[0])/(energyGamma[0]+energyGamma[1]);
          fHistPtAlphaEtaGGPCMPHOAcc->Fill(particle->Pt(), alpha);
        }
      }
      if (particle->PdgCode()==kPdgEtaPrime){
        fHistPtYEtaPrimeGG->Fill(particle->Pt(), particle->Y());
        if (TESTBIT(acceptanceGamma[0], kPCMAcceptance) && TESTBIT(acceptanceGamma[1], kPCMAcceptance)){
          fHistPtYEtaPrimeGGPCMAcc->Fill(particle->Pt(), particle->Y());
          fHistPtAlphaEtaPrimeGGPCMAcc->Fill(particle->Pt(), TMath::Abs(alpha));
        }
        if (TESTBIT(acceptanceGamma[0], kEMCALAcceptance) && TESTBIT(acceptanceGamma[1], kEMCALAcceptance)){
          fHistPtYEtaPrimeGGEMCAcc->Fill(particle->Pt(), particle->Y());
          fHistPtAlphaEtaPrimeGGEMCAcc->Fill(particle->Pt(), TMath::Abs(alpha));
        }
        if (TESTBIT(acceptanceGamma[0], kPHOSAcceptance) && TESTBIT(acceptanceGamma[1], kPHOSAcceptance)){
          fHistPtYEtaPrimeGGPHOAcc->Fill(particle->Pt(), particle->Y());
          fHistPtAlphaEtaPrimeGGPHOAcc->Fill(particle->Pt(), TMath::Abs(alpha));
        }
        if ( (TESTBIT(acceptanceGamma[0], kPCMAcceptance)  && TESTBIT(acceptanceGamma[1], kEMCALAcceptance)) ||
            (TESTBIT(acceptanceGamma[1], kPCMAcceptance) && TESTBIT(acceptanceGamma[0], kEMCALAcceptance))
        ){
          fHistPtYEtaPrimeGGPCMEMCAcc->Fill(particle->Pt(), particle->Y());
          if (!TESTBIT(acceptanceGamma[1], kEMCALAcceptance)) alpha = (energyGamma[1]-energyGamma[0])/(energyGamma[0]+energyGamma[1]);
          fHistPtAlphaEtaPrimeGGPCMEMCAcc->Fill(particle->Pt(), alpha);
        }
        if ( (TESTBIT(acceptanceGamma[0], kPCMAcceptance) && TESTBIT(acceptanceGamma[1], kPHOSAcceptance)) ||
            (TESTBIT(acceptanceGamma[1], kPCMAcceptance) && TESTBIT(acceptanceGamma[0], kPHOSAcceptance))
        ){
          fHistPtYEtaPrimeGGPCMPHOAcc->Fill(particle->Pt(), particle->Y());
          if (TESTBIT(!acceptanceGamma[1],kPHOSAcceptance)) alpha = (energyGamma[1]-energyGamma[0])/(energyGamma[0]+energyGamma[1]);
          fHistPtAlphaEtaPrimeGGPCMPHOAcc->Fill(particle->Pt(), alpha);
        }
      }
    }


    if(fIsK0 == 0) continue;
    if( particle->PdgCode() == kPdgK0Short){
      if (particle->GetNDaughters() != 2) continue;
      //UChar_t acceptanceGamma[2] = {0,0};
      Double_t energyPi0[2] = {0,0};
      Bool_t allOK[2] = {kFALSE,kFALSE};
      UChar_t gdAcceptanceGamma[4] = {0,0,0,0};
      //Double_t gdEnergyGamma[4] = {0,0,0,0};
      Bool_t allGDOK[4] = {kFALSE, kFALSE, kFALSE,kFALSE};
      for(Int_t k=0;k<2;k++){
        AliVParticle *daughter = (AliVParticle*) fMCEvent->GetTrack(particle->GetDaughterLabel(k));
        if (!daughter) continue;

        // Is Daughter a pi0?
        if (daughter->PdgCode() == kPdgPi0){
          allOK[k] = kTRUE;
          if(daughter->GetNDaughters() != 2) continue;
          energyPi0[k] = daughter->E();
          for(Int_t h=0;h<2;h++){
            AliVParticle *granddaughter = (AliVParticle*) fMCEvent->GetTrack(daughter->GetDaughterLabel(k));
            if(granddaughter->PdgCode() == 22) allGDOK[2*k + h] = kTRUE;
            if(IsInPCMAcceptance(granddaughter))  SETBIT(gdAcceptanceGamma[2*k+h], kPCMAcceptance);
            if(IsInEMCalAcceptance(granddaughter)) SETBIT(gdAcceptanceGamma[2*k+h], kEMCALAcceptance);
            //gdEnergyGamma[2*k+h] = granddaughter->Energy();
          }
        }
      }

	  Double_t alpha_k0 = (energyPi0[0]-energyPi0[1])/(energyPi0[0]+energyPi0[1]);

      if(allOK[0] && allOK[1]){
            fHistPtYPi0FromKGG->Fill(particle->Pt(), particle->Y());
            fHistPtAlphaPi0FromKGG->Fill(particle->Pt(), alpha_k0);
      }

      if (!(allGDOK[0] && allGDOK[1] && allGDOK[2] && allGDOK[3])) continue;

      if (TESTBIT(gdAcceptanceGamma[0], kPCMAcceptance) && TESTBIT(gdAcceptanceGamma[1], kPCMAcceptance)
        && TESTBIT(gdAcceptanceGamma[2], kPCMAcceptance) && TESTBIT(gdAcceptanceGamma[3], kPCMAcceptance))
        {
          fHistPtYPi0FromKGGPCMAcc->Fill(particle->Pt(),particle->Y());
          fHistPtAlphaPi0FromKGGPCMAcc->Fill(particle->Pt(), alpha_k0);
        }

      if (TESTBIT(gdAcceptanceGamma[0], kEMCALAcceptance) && TESTBIT(gdAcceptanceGamma[1], kEMCALAcceptance)
        && TESTBIT(gdAcceptanceGamma[2], kEMCALAcceptance) && TESTBIT(gdAcceptanceGamma[3], kEMCALAcceptance) )
        {
          fHistPtYPi0FromKGGEMCAcc->Fill(particle->Pt(),particle->Y());
          fHistPtAlphaPi0FromKGGEMCAcc->Fill(particle->Pt(), alpha_k0);
      	}


      if ((TESTBIT(gdAcceptanceGamma[0], kEMCALAcceptance) && TESTBIT(gdAcceptanceGamma[1], kPCMAcceptance)
        && TESTBIT(gdAcceptanceGamma[2], kPCMAcceptance) && TESTBIT(gdAcceptanceGamma[3], kPCMAcceptance)) ||
        (TESTBIT(gdAcceptanceGamma[0], kPCMAcceptance) && TESTBIT(gdAcceptanceGamma[1], kEMCALAcceptance)
        && TESTBIT(gdAcceptanceGamma[2], kPCMAcceptance) && TESTBIT(gdAcceptanceGamma[3], kPCMAcceptance)) ||
        (TESTBIT(gdAcceptanceGamma[0], kPCMAcceptance) && TESTBIT(gdAcceptanceGamma[1], kPCMAcceptance)
        && TESTBIT(gdAcceptanceGamma[2], kEMCALAcceptance) && TESTBIT(gdAcceptanceGamma[3], kPCMAcceptance)) ||
        (TESTBIT(gdAcceptanceGamma[0], kPCMAcceptance) && TESTBIT(gdAcceptanceGamma[1], kPCMAcceptance)
        && TESTBIT(gdAcceptanceGamma[2], kPCMAcceptance) && TESTBIT(gdAcceptanceGamma[3], kEMCALAcceptance)))
        {
          fHistPtYPi0FromKGGPCMEMCAcc->Fill(particle->Pt(),particle->Y());
          fHistPtAlphaPi0FromKGGPCMEMCAcc->Fill(particle->Pt(), alpha_k0);
        }

      if ((TESTBIT(gdAcceptanceGamma[0], kPCMAcceptance) && TESTBIT(gdAcceptanceGamma[1], kEMCALAcceptance)
        && TESTBIT(gdAcceptanceGamma[2], kEMCALAcceptance) && TESTBIT(gdAcceptanceGamma[3], kEMCALAcceptance)) ||
        (TESTBIT(gdAcceptanceGamma[0], kEMCALAcceptance) && TESTBIT(gdAcceptanceGamma[1], kPCMAcceptance)
        && TESTBIT(gdAcceptanceGamma[2], kEMCALAcceptance) && TESTBIT(gdAcceptanceGamma[3], kEMCALAcceptance)) ||
        (TESTBIT(gdAcceptanceGamma[0], kEMCALAcceptance) && TESTBIT(gdAcceptanceGamma[1], kEMCALAcceptance)
        && TESTBIT(gdAcceptanceGamma[2], kPCMAcceptance) && TESTBIT(gdAcceptanceGamma[3], kEMCALAcceptance)) ||
        (TESTBIT(gdAcceptanceGamma[0], kEMCALAcceptance) && TESTBIT(gdAcceptanceGamma[1], kEMCALAcceptance)
        && TESTBIT(gdAcceptanceGamma[2], kEMCALAcceptance) && TESTBIT(gdAcceptanceGamma[3], kPCMAcceptance)))
	    {
          fHistPtYPi0FromKGGEMCPCMAcc->Fill(particle->Pt(),particle->Y());
          fHistPtAlphaPi0FromKGGEMCPCMAcc->Fill(particle->Pt(), alpha_k0);
		}

      if ((TESTBIT(gdAcceptanceGamma[0], kEMCALAcceptance) && TESTBIT(gdAcceptanceGamma[1], kEMCALAcceptance) &&
        TESTBIT(gdAcceptanceGamma[2], kPCMAcceptance) && TESTBIT(gdAcceptanceGamma[3], kPCMAcceptance))||
        (TESTBIT(gdAcceptanceGamma[0], kPCMAcceptance) && TESTBIT(gdAcceptanceGamma[1], kPCMAcceptance) &&
        TESTBIT(gdAcceptanceGamma[2], kEMCALAcceptance) && TESTBIT(gdAcceptanceGamma[3], kEMCALAcceptance)))
        {
          fHistPtYPi0FromKGGEMCAccSamePi0->Fill(particle->Pt(),particle->Y());
          fHistPtAlphaPi0FromKGGEMCAccSamePi0->Fill(particle->Pt(), alpha_k0);
        }

      if ((TESTBIT(gdAcceptanceGamma[0], kEMCALAcceptance) && TESTBIT(gdAcceptanceGamma[2], kEMCALAcceptance)&&
        TESTBIT(gdAcceptanceGamma[1], kPCMAcceptance) && TESTBIT(gdAcceptanceGamma[3], kPCMAcceptance))||
        (TESTBIT(gdAcceptanceGamma[0], kEMCALAcceptance) && TESTBIT(gdAcceptanceGamma[3], kEMCALAcceptance)&&
        TESTBIT(gdAcceptanceGamma[1], kPCMAcceptance) && TESTBIT(gdAcceptanceGamma[2], kPCMAcceptance))||
        (TESTBIT(gdAcceptanceGamma[1], kEMCALAcceptance) && TESTBIT(gdAcceptanceGamma[2], kEMCALAcceptance)&&
        TESTBIT(gdAcceptanceGamma[0], kPCMAcceptance) && TESTBIT(gdAcceptanceGamma[3], kPCMAcceptance))||
        (TESTBIT(gdAcceptanceGamma[1], kEMCALAcceptance) && TESTBIT(gdAcceptanceGamma[3], kEMCALAcceptance)&&
        TESTBIT(gdAcceptanceGamma[0], kPCMAcceptance) && TESTBIT(gdAcceptanceGamma[2], kPCMAcceptance)))
        {
          fHistPtYPi0FromKGGEMCAccDiffPi0->Fill(particle->Pt(),particle->Y());
          fHistPtAlphaPi0FromKGGEMCAccDiffPi0->Fill(particle->Pt(), alpha_k0);
        }
    }
  }
}

//________________________________________________________________________
bool AliAnalysisTaskGammaPureMC::IsInPCMAcceptance(AliVParticle* part) const {
  const Double_t kBoundaryEta = 0.900001;
  if (part->Pt() > 0.050 && TMath::Abs(part->Eta()) < kBoundaryEta) return true;

  return false;
}

//________________________________________________________________________
bool AliAnalysisTaskGammaPureMC::IsInPHOSAcceptance(AliVParticle* part) const {
  const Double_t kBoundaryEtaMin = -0.13;
  const Double_t kBoundaryEtaMax = 0.13;
  const Double_t kBoundaryPhiMin = 4.54;
  const Double_t kBoundaryPhiMax = 5.59;
  if (part->Pt() < 0.300) return false;
  if (part->Eta() > kBoundaryEtaMax || part->Eta() < kBoundaryEtaMin) return false;
  if (part->Phi() > kBoundaryPhiMax || part->Phi() < kBoundaryPhiMin) return false;
  return true;
}

//________________________________________________________________________
bool AliAnalysisTaskGammaPureMC::IsInEMCalAcceptance(AliVParticle* part) const {
  const Double_t kBoundaryEtaMin = -0.6687;
  const Double_t kBoundaryEtaMax = 0.66465;
  const Double_t kBoundaryPhiMin = 1.39626;
  const Double_t kBoundaryPhiMax = 3.15;
  if (part->Pt() < 0.400) return false;
  if (part->Eta() > kBoundaryEtaMax || part->Eta() < kBoundaryEtaMin) return false;
  if (part->Phi() > kBoundaryPhiMax || part->Phi() < kBoundaryPhiMin) return false;
  return true;
}

//________________________________________________________________________
bool AliAnalysisTaskGammaPureMC::IsInV0Acceptance(AliVParticle* part) const {
  const Double_t kBoundaryEtaMinV0A = 2.8;
  const Double_t kBoundaryEtaMaxV0A = 5.1;
  const Double_t kBoundaryEtaMinV0C = -3.7;
  const Double_t kBoundaryEtaMaxV0C = -1.7;
  if (part->Eta() < kBoundaryEtaMaxV0A && part->Eta() > kBoundaryEtaMinV0A) return true;
  if (part->Eta() < kBoundaryEtaMaxV0C && part->Eta() > kBoundaryEtaMinV0C) return true;
  return false;
}

//________________________________________________________________________
bool AliAnalysisTaskGammaPureMC::IsSecondary(AliVParticle* motherParticle) const {
  if (TMath::Abs(motherParticle->PdgCode()) == kPdgK0Short ||
      TMath::Abs(motherParticle->PdgCode()) == kPdgK0Long  ||
      TMath::Abs(motherParticle->PdgCode()) == kPdgKPlus  ||
      TMath::Abs(motherParticle->PdgCode()) == kPdgLambda ||
      TMath::Abs(motherParticle->PdgCode()) == kPdgSigma0 ||
      TMath::Abs(motherParticle->PdgCode()) == kPdgSigmaPlus ||
      TMath::Abs(motherParticle->PdgCode()) == kPdgSigmaMinus ||
      TMath::Abs(motherParticle->PdgCode()) == kPdgXi0 ||
      TMath::Abs(motherParticle->PdgCode()) == kPdgXiMinus
  ){
    return true;
  }
  return false;
}

//________________________________________________________________________
void AliAnalysisTaskGammaPureMC::Terminate(const Option_t *)
{

  //fOutputContainer->Print(); // Will crash on GRID
}


//_________________________________________________________________________________
void AliAnalysisTaskGammaPureMC::SetLogBinningXTH1(TH1* histoRebin){
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
void AliAnalysisTaskGammaPureMC::SetLogBinningXTH2(TH2* histoRebin){
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
void AliAnalysisTaskGammaPureMC::ProcessJets(){
  fVecJets.clear();

  std::vector<fastjet::PseudoJet> vecParticles;
  for(Long_t i = 0; i < fMCEvent->GetNumberOfTracks(); i++) {
    // fill primary histograms
    AliVParticle* particle     = nullptr;
    particle                    = (AliVParticle *)fMCEvent->GetTrack(i);
    if (!particle) continue;

    // only select stable particles for jet finder
    if(fDoJetStudies == 1){ // reject all unstable particles
      if(particle->GetDaughterFirst() > 0){
        continue;
      }
    } else if(fDoJetStudies == 2){ // set pi0 and eta to stable and reject photons from those decays
      if(particle->GetDaughterFirst() > 0 && (particle->PdgCode() != 111 && particle->PdgCode() != 221)){
        continue;
      } else if(particle->GetMother() > 0){
        if(static_cast<AliVParticle*>(fMCEvent->GetTrack(particle->GetMother()))->PdgCode() == 111 || static_cast<AliVParticle*>(fMCEvent->GetTrack(particle->GetMother()))->PdgCode() == 221){
          continue;
        }
      }
    }
  
    fastjet::PseudoJet jetPart(particle->Px(), particle->Py(), particle->Pz(), particle->P());
    jetPart.set_user_index(i);
    vecParticles.push_back(jetPart);
  }

  fastjet::Selector sel_jets = fastjet::SelectorEMin(fJetMinE) * fastjet::SelectorEtaRange(-fJetAccEta, fJetAccEta);
  fastjet::JetDefinition jet_def(fJetAlgorithm, fJetRadius, fJetRecombScheme, fJetStrategy);
  fastjet::GhostedAreaSpec ghostSpec(fGhostEtaMax, fActiveAreaRepeats, fJetGhostArea, 0.1, 1e-100);
  fastjet::AreaDefinition area_def = fastjet::AreaDefinition(fAreaType,ghostSpec);
  fastjet::ClusterSequenceArea clust_seq_full(vecParticles, jet_def, area_def);
  fVecJets = sel_jets(clust_seq_full.inclusive_jets());

  for(const auto & jet : fVecJets){
    fHistJetPtY->Fill(jet.pt(), jet.eta());
    fHistJetEta->Fill(jet.eta());
    fHistJetPhi->Fill(jet.phi());
  }

}

//_________________________________________________________________________________
bool AliAnalysisTaskGammaPureMC::IsParticleInJet(const std::vector<fastjet::PseudoJet>& vecJet, double eta, double phi, int& index, double& R){
  R = 100; // some high value that will get replaced
  for(unsigned int i = 0; i < vecJet.size(); ++i){
    double dEta = std::abs(vecJet[i].eta() - eta);
    double dPhi = std::abs(vecJet[i].phi() - phi);
    double Rtmp = sqrt(dEta*dEta + dPhi * dPhi);
    if(Rtmp < fJetRadius && Rtmp < R){
      R = Rtmp;
      index = i;
    }
  }
  if(R >= 100) return false;
  return true;
}

//_________________________________________________________________________________
double AliAnalysisTaskGammaPureMC::GetFrag(const fastjet::PseudoJet& jet, const AliVParticle* part){
  double ScalarProd = jet.px()*part->Px() + jet.py()*part->Py() + jet.pz()*part->Pz();
  double jetP2 = jet.px()*jet.px() + jet.py()*jet.py() + jet.pz()*jet.pz();
  double z = ScalarProd/jetP2;
  return z;
}
