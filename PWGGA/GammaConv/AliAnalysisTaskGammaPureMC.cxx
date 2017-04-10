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
#include "AliAnalysisTaskGammaPureMC.h"
#include "AliVParticle.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenDPMjetEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliEventplane.h"
#include "AliInputEventHandler.h"
#include <vector>
#include <map>

ClassImp(AliAnalysisTaskGammaPureMC)

//________________________________________________________________________
AliAnalysisTaskGammaPureMC::AliAnalysisTaskGammaPureMC(): AliAnalysisTaskSE(),
  fOutputContainer(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fMCStack(NULL),
  fHistNEvents(NULL),
  fHistXSection(NULL),
  fHistPtHard(NULL),
  fHistPtYPi0(NULL),
  fHistPtYPiPl(NULL),
  fHistPtYPiMi(NULL),
  fHistPtYEta(NULL),
  fHistPtYEtaPrim(NULL),
  fHistPtYOmega(NULL),
  fHistPtYRho0(NULL),
  fHistPtYRhoPl(NULL),
  fHistPtYRhoMi(NULL),
  fHistPtYPhi(NULL),
  fHistPtYJPsi(NULL),
  fHistPtYSigma0(NULL),
  fHistPtYK0s(NULL),
  fHistPtYK0l(NULL),
  fHistPtYK0star(NULL),
  fHistPtYDeltaPlPl(NULL),
  fHistPtYDeltaPl(NULL),
  fHistPtYDeltaMi(NULL),
  fHistPtYDelta0(NULL),
  fHistPtYLambda(NULL),
  fHistPtYPi0FromEta(NULL),
  fHistPtYPi0FromLambda(NULL),
  fHistPtYPi0FromK(NULL),
  fHistPtYPiPlFromK(NULL),
  fHistPtYPiMiFromK(NULL),
  fHistPtYPi0GG(NULL),
  fHistPtYPi0GGPCMAcc(NULL),
  fHistPtYPi0GGEMCAcc(NULL),
  fHistPtYPi0GGPHOAcc(NULL),
  fHistPtYPi0GGPCMEMCAcc(NULL),
  fHistPtYPi0GGPCMPHOAcc(NULL),
  fHistPtAlphaPi0GGPCMAcc(NULL),
  fHistPtAlphaPi0GGEMCAcc(NULL),
  fHistPtAlphaPi0GGPHOAcc(NULL),
  fHistPtAlphaPi0GGPCMEMCAcc(NULL),
  fHistPtAlphaPi0GGPCMPHOAcc(NULL),
  fHistPtYEtaGG(NULL),
  fHistPtYEtaGGPCMAcc(NULL),
  fHistPtYEtaGGEMCAcc(NULL),
  fHistPtYEtaGGPHOAcc(NULL),
  fHistPtYEtaGGPCMEMCAcc(NULL),
  fHistPtYEtaGGPCMPHOAcc(NULL),
  fHistPtAlphaEtaGGPCMAcc(NULL),
  fHistPtAlphaEtaGGEMCAcc(NULL),
  fHistPtAlphaEtaGGPHOAcc(NULL),
  fHistPtAlphaEtaGGPCMEMCAcc(NULL),
  fHistPtAlphaEtaGGPCMPHOAcc(NULL),
  fHistPtYEtaPrimGG(NULL),
  fHistPtYEtaPrimGGPCMAcc(NULL),
  fHistPtYEtaPrimGGEMCAcc(NULL),
  fHistPtYEtaPrimGGPHOAcc(NULL),
  fHistPtYEtaPrimGGPCMEMCAcc(NULL),
  fHistPtYEtaPrimGGPCMPHOAcc(NULL),
  fIsMC(1)
{
  
}

//________________________________________________________________________
AliAnalysisTaskGammaPureMC::AliAnalysisTaskGammaPureMC(const char *name):
  AliAnalysisTaskSE(name),
  fOutputContainer(NULL),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fMCStack(NULL),
  fHistNEvents(NULL),
  fHistXSection(NULL),
  fHistPtHard(NULL),
  fHistPtYPi0(NULL),
  fHistPtYPiPl(NULL),
  fHistPtYPiMi(NULL),
  fHistPtYEta(NULL),
  fHistPtYEtaPrim(NULL),
  fHistPtYOmega(NULL),
  fHistPtYRho0(NULL),
  fHistPtYRhoPl(NULL),
  fHistPtYRhoMi(NULL),
  fHistPtYPhi(NULL),
  fHistPtYJPsi(NULL),
  fHistPtYSigma0(NULL),
  fHistPtYK0s(NULL),
  fHistPtYK0l(NULL),
  fHistPtYK0star(NULL),
  fHistPtYDeltaPlPl(NULL),
  fHistPtYDeltaPl(NULL),
  fHistPtYDeltaMi(NULL),
  fHistPtYDelta0(NULL),
  fHistPtYLambda(NULL),
  fHistPtYPi0FromEta(NULL),
  fHistPtYPi0FromLambda(NULL),
  fHistPtYPi0FromK(NULL),
  fHistPtYPiPlFromK(NULL),
  fHistPtYPiMiFromK(NULL),
  fHistPtYPi0GG(NULL),
  fHistPtYPi0GGPCMAcc(NULL),
  fHistPtYPi0GGEMCAcc(NULL),
  fHistPtYPi0GGPHOAcc(NULL),
  fHistPtYPi0GGPCMEMCAcc(NULL),
  fHistPtYPi0GGPCMPHOAcc(NULL),
  fHistPtAlphaPi0GGPCMAcc(NULL),
  fHistPtAlphaPi0GGEMCAcc(NULL),
  fHistPtAlphaPi0GGPHOAcc(NULL),
  fHistPtAlphaPi0GGPCMEMCAcc(NULL),
  fHistPtAlphaPi0GGPCMPHOAcc(NULL),
  fHistPtYEtaGG(NULL),
  fHistPtYEtaGGPCMAcc(NULL),
  fHistPtYEtaGGEMCAcc(NULL),
  fHistPtYEtaGGPHOAcc(NULL),
  fHistPtYEtaGGPCMEMCAcc(NULL),
  fHistPtYEtaGGPCMPHOAcc(NULL),
  fHistPtAlphaEtaGGPCMAcc(NULL),
  fHistPtAlphaEtaGGEMCAcc(NULL),
  fHistPtAlphaEtaGGPHOAcc(NULL),
  fHistPtAlphaEtaGGPCMEMCAcc(NULL),
  fHistPtAlphaEtaGGPCMPHOAcc(NULL),
  fHistPtYEtaPrimGG(NULL),
  fHistPtYEtaPrimGGPCMAcc(NULL),
  fHistPtYEtaPrimGGEMCAcc(NULL),
  fHistPtYEtaPrimGGPHOAcc(NULL),
  fHistPtYEtaPrimGGPCMEMCAcc(NULL),
  fHistPtYEtaPrimGGPCMPHOAcc(NULL),
  fIsMC(1)
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
  if(fOutputContainer != NULL){
    delete fOutputContainer;
    fOutputContainer          = NULL;
  }
  if(fOutputContainer == NULL){
    fOutputContainer          = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }
  
  fHistNEvents                = new TH1F("NEvents", "NEvents", 3, -0.5, 2.5);
  fHistNEvents->Sumw2();
    fOutputContainer->Add(fHistNEvents);

  fHistXSection               = new TH1D("XSection", "XSection", 1000000, 0, 1e4);
//   SetLogBinningXTH1(fHistXSection);
  fHistXSection->Sumw2();
    fOutputContainer->Add(fHistXSection);
  
  fHistPtHard                 = new TH1F("PtHard", "PtHard", 400, 0, 200);
  fHistPtHard->Sumw2();
    fOutputContainer->Add(fHistPtHard);
  
  fHistPtYPi0                 = new TH2F("Pt_Y_Pi0","Pt_Y_Pi0", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPi0->Sumw2();
  fOutputContainer->Add(fHistPtYPi0);

  fHistPtYPiPl                = new TH2F("Pt_Y_PiPl","Pt_Y_PiPl", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPiPl->Sumw2();
  fOutputContainer->Add(fHistPtYPiPl);

  fHistPtYPiMi                = new TH2F("Pt_Y_PiMi","Pt_Y_PiMi", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPiMi->Sumw2();
  fOutputContainer->Add(fHistPtYPiMi);
  
  fHistPtYEta                 = new TH2F("Pt_Y_Eta","Pt_Y_Eta", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEta->Sumw2();
  fOutputContainer->Add(fHistPtYEta);

  fHistPtYEtaPrim             = new TH2F("Pt_Y_EtaPrim","Pt_Y_EtaPrim", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaPrim->Sumw2();
  fOutputContainer->Add(fHistPtYEtaPrim);

  fHistPtYOmega               = new TH2F("Pt_Y_Omega","Pt_Y_Omega", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYOmega->Sumw2();
  fOutputContainer->Add(fHistPtYOmega);
  
  fHistPtYRho0                = new TH2F("Pt_Y_Rho0","Pt_Y_Rho0", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYRho0->Sumw2();
  fOutputContainer->Add(fHistPtYRho0);

  fHistPtYRhoPl               = new TH2F("Pt_Y_RhoPl","Pt_Y_RhoPl", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYRhoPl->Sumw2();
  fOutputContainer->Add(fHistPtYRhoPl);

  fHistPtYRhoMi               = new TH2F("Pt_Y_RhoMi","Pt_Y_RhoMi", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYRhoMi->Sumw2();
  fOutputContainer->Add(fHistPtYRhoMi);

  fHistPtYPhi                 = new TH2F("Pt_Y_Phi","Pt_Y_Phi", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPhi->Sumw2();
  fOutputContainer->Add(fHistPtYPhi);

  fHistPtYJPsi                = new TH2F("Pt_Y_JPsi","Pt_Y_JPsi", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYJPsi->Sumw2();
  fOutputContainer->Add(fHistPtYJPsi);

  fHistPtYSigma0              = new TH2F("Pt_Y_Sigma0","Pt_Y_Sigma0", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYSigma0->Sumw2();
  fOutputContainer->Add(fHistPtYSigma0);

  fHistPtYK0s                 = new TH2F("Pt_Y_K0s","Pt_Y_K0s", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYK0s->Sumw2();
  fOutputContainer->Add(fHistPtYK0s);

  fHistPtYK0l                 = new TH2F("Pt_Y_K0l","Pt_Y_K0l", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYK0l->Sumw2();
  fOutputContainer->Add(fHistPtYK0l);

  fHistPtYK0star              = new TH2F("Pt_Y_K0star","Pt_Y_K0star", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYK0star->Sumw2();
  fOutputContainer->Add(fHistPtYK0star);

  fHistPtYDeltaPlPl           = new TH2F("Pt_Y_DeltaPlPl","Pt_Y_DeltaPlPl", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYDeltaPlPl->Sumw2();
  fOutputContainer->Add(fHistPtYDeltaPlPl);

  fHistPtYDeltaPl             = new TH2F("Pt_Y_DeltaPl","Pt_Y_DeltaPl", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYDeltaPl->Sumw2();
  fOutputContainer->Add(fHistPtYDeltaPl);

  fHistPtYDeltaMi             = new TH2F("Pt_Y_DeltaMi","Pt_Y_DeltaMi", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYDeltaMi->Sumw2();
  fOutputContainer->Add(fHistPtYDeltaMi);

  fHistPtYDelta0              = new TH2F("Pt_Y_Delta0","Pt_Y_Delta0", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYDelta0->Sumw2();
  fOutputContainer->Add(fHistPtYDelta0);

  fHistPtYLambda              = new TH2F("Pt_Y_Lambda","Pt_Y_Lambda", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYLambda->Sumw2();
  fOutputContainer->Add(fHistPtYLambda);

  fHistPtYPi0FromEta          = new TH2F("Pt_Y_Pi0FromEta","Pt_Y_Pi0FromEta", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPi0FromEta->Sumw2();
  fOutputContainer->Add(fHistPtYPi0FromEta);

  fHistPtYPi0FromLambda       = new TH2F("Pt_Y_Pi0FromLambda","Pt_Y_Pi0FromLambda", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPi0FromLambda->Sumw2();
  fOutputContainer->Add(fHistPtYPi0FromLambda);

  fHistPtYPi0FromK            = new TH2F("Pt_Y_Pi0FromK","Pt_Y_Pi0FromK", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPi0FromK->Sumw2();
  fOutputContainer->Add(fHistPtYPi0FromK);

  fHistPtYPiPlFromK           = new TH2F("Pt_Y_PiPlFromK","Pt_Y_PiPlFromK", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPiPlFromK->Sumw2();
  fOutputContainer->Add(fHistPtYPiPlFromK);

  fHistPtYPiMiFromK           = new TH2F("Pt_Y_PiMiFromK","Pt_Y_PiMiFromK", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPiMiFromK->Sumw2();
  fOutputContainer->Add(fHistPtYPiMiFromK);

  
  fHistPtYPi0GG               = new TH2F("Pt_Y_Pi0GG","Pt_Y_Pi0GG", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPi0GG->Sumw2();
  fOutputContainer->Add(fHistPtYPi0GG);
  fHistPtYPi0GGPCMAcc         = new TH2F("Pt_Y_Pi0GGPCMAcc","Pt_Y_Pi0GGPCMAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPi0GGPCMAcc->Sumw2();
  fOutputContainer->Add(fHistPtYPi0GGPCMAcc);
  fHistPtYPi0GGEMCAcc         = new TH2F("Pt_Y_Pi0GGEMCAcc","Pt_Y_Pi0GGEMCAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPi0GGEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtYPi0GGEMCAcc);
  fHistPtYPi0GGPHOAcc         = new TH2F("Pt_Y_Pi0GGPHOAcc","Pt_Y_Pi0GGPHOAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPi0GGPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtYPi0GGPHOAcc);
  fHistPtYPi0GGPCMEMCAcc      = new TH2F("Pt_Y_Pi0GGPCMEMCAcc","Pt_Y_Pi0GGPCMEMCAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPi0GGPCMEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtYPi0GGPCMEMCAcc);
  fHistPtYPi0GGPCMPHOAcc      = new TH2F("Pt_Y_Pi0GGPCMPHOAcc","Pt_Y_Pi0GGPCMPHOAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPi0GGPCMPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtYPi0GGPCMPHOAcc);

  fHistPtYEtaGG               = new TH2F("Pt_Y_EtaGG","Pt_Y_EtaGG", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaGG->Sumw2();
  fOutputContainer->Add(fHistPtYEtaGG);
  fHistPtYEtaGGPCMAcc         = new TH2F("Pt_Y_EtaGGPCMAcc","Pt_Y_EtaGGPCMAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaGGPCMAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaGGPCMAcc);
  fHistPtYEtaGGEMCAcc         = new TH2F("Pt_Y_EtaGGEMCAcc","Pt_Y_EtaGGEMCAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaGGEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaGGEMCAcc);
  fHistPtYEtaGGPHOAcc         = new TH2F("Pt_Y_EtaGGPHOAcc","Pt_Y_EtaGGPHOAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaGGPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaGGPHOAcc);
  fHistPtYEtaGGPCMEMCAcc      = new TH2F("Pt_Y_EtaGGPCMEMCAcc","Pt_Y_EtaGGPCMEMCAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaGGPCMEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaGGPCMEMCAcc);
  fHistPtYEtaGGPCMPHOAcc      = new TH2F("Pt_Y_EtaGGPCMPHOAcc","Pt_Y_EtaGGPCMPHOAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaGGPCMPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaGGPCMPHOAcc);

  fHistPtYEtaPrimGG           = new TH2F("Pt_Y_EtaPrimGG","Pt_Y_EtaPrimGG", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaPrimGG->Sumw2();
  fOutputContainer->Add(fHistPtYEtaPrimGG);
  fHistPtYEtaPrimGGPCMAcc     = new TH2F("Pt_Y_EtaPrimGGPCMAcc","Pt_Y_EtaPrimGGPCMAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaPrimGGPCMAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaPrimGGPCMAcc);
  fHistPtYEtaPrimGGEMCAcc     = new TH2F("Pt_Y_EtaPrimGGEMCAcc","Pt_Y_EtaPrimGGEMCAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaPrimGGEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaPrimGGEMCAcc);
  fHistPtYEtaPrimGGPHOAcc     = new TH2F("Pt_Y_EtaPrimGGPHOAcc","Pt_Y_EtaPrimGGPHOAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaPrimGGPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaPrimGGPHOAcc);
  fHistPtYEtaPrimGGPCMEMCAcc  = new TH2F("Pt_Y_EtaPrimGGPCMEMCAcc","Pt_Y_EtaPrimGGPCMEMCAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaPrimGGPCMEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaPrimGGPCMEMCAcc);
  fHistPtYEtaPrimGGPCMPHOAcc  = new TH2F("Pt_Y_EtaPrimGGPCMPHOAcc","Pt_Y_EtaPrimGGPCMPHOAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaPrimGGPCMPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaPrimGGPCMPHOAcc);

  fHistPtAlphaPi0GGPCMAcc     = new TH2F("Pt_Alpha_Pi0GGPCMAcc","Pt_Alpha_Pi0GGPCMAcc", 500,0.1, 50, 100, 0., 1.);
  SetLogBinningXTH2(fHistPtAlphaPi0GGPCMAcc);
  fHistPtAlphaPi0GGPCMAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaPi0GGPCMAcc);
  fHistPtAlphaPi0GGEMCAcc     = new TH2F("Pt_Alpha_Pi0GGEMCAcc","Pt_Alpha_Pi0GGEMCAcc", 500,0.1, 50, 100, 0., 1.);
  SetLogBinningXTH2(fHistPtAlphaPi0GGEMCAcc);
  fHistPtAlphaPi0GGEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaPi0GGEMCAcc);
  fHistPtAlphaPi0GGPHOAcc     = new TH2F("Pt_Alpha_Pi0GGPHOAcc","Pt_Alpha_Pi0GGPHOAcc", 500,0.1, 50, 100, 0., 1.);
  SetLogBinningXTH2(fHistPtAlphaPi0GGPHOAcc);
  fHistPtAlphaPi0GGPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaPi0GGPHOAcc);
  fHistPtAlphaPi0GGPCMEMCAcc  = new TH2F("Pt_Alpha_Pi0GGPCMEMCAcc","Pt_Alpha_Pi0GGPCMEMCAcc", 500,0.1, 50, 200, -1., 1.);
  SetLogBinningXTH2(fHistPtAlphaPi0GGPCMEMCAcc);
  fHistPtAlphaPi0GGPCMEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaPi0GGPCMEMCAcc);
  fHistPtAlphaPi0GGPCMPHOAcc  = new TH2F("Pt_Alpha_Pi0GGPCMPHOAcc","Pt_Alpha_Pi0GGPCMPHOAcc", 500,0.1, 50, 200, -1., 1.);
  SetLogBinningXTH2(fHistPtAlphaPi0GGPCMPHOAcc);
  fHistPtAlphaPi0GGPCMPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaPi0GGPCMPHOAcc);

  fHistPtAlphaEtaGGPCMAcc     = new TH2F("Pt_Alpha_EtaGGPCMAcc","Pt_Alpha_EtaGGPCMAcc", 500,0.1, 50, 100, 0., 1.);
  SetLogBinningXTH2(fHistPtAlphaEtaGGPCMAcc);
  fHistPtAlphaEtaGGPCMAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaEtaGGPCMAcc);
  fHistPtAlphaEtaGGEMCAcc     = new TH2F("Pt_Alpha_EtaGGEMCAcc","Pt_Alpha_EtaGGEMCAcc", 500,0.1, 50, 100, 0., 1.);
  SetLogBinningXTH2(fHistPtAlphaEtaGGEMCAcc);
  fHistPtAlphaEtaGGEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaEtaGGEMCAcc);
  fHistPtAlphaEtaGGPHOAcc     = new TH2F("Pt_Alpha_EtaGGPHOAcc","Pt_Alpha_EtaGGPHOAcc", 500,0.1, 50, 100, 0., 1.);
  SetLogBinningXTH2(fHistPtAlphaEtaGGPHOAcc);
  fHistPtAlphaEtaGGPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaEtaGGPHOAcc);
  fHistPtAlphaEtaGGPCMEMCAcc  = new TH2F("Pt_Alpha_EtaGGPCMEMCAcc","Pt_Alpha_EtaGGPCMEMCAcc", 500,0.1, 50, 200, -1., 1.);
  SetLogBinningXTH2(fHistPtAlphaEtaGGPCMEMCAcc);
  fHistPtAlphaEtaGGPCMEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaEtaGGPCMEMCAcc);
  fHistPtAlphaEtaGGPCMPHOAcc  = new TH2F("Pt_Alpha_EtaGGPCMPHOAcc","Pt_Alpha_EtaGGPCMPHOAcc", 500,0.1, 50, 200, -1., 1.);
  SetLogBinningXTH2(fHistPtAlphaEtaGGPCMPHOAcc);
  fHistPtAlphaEtaGGPCMPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaEtaGGPCMPHOAcc);
  
  PostData(1, fOutputContainer);
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaPureMC::UserExec(Option_t *)
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
  
//   cout << "the stack is intact" << endl;
  
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
void AliAnalysisTaskGammaPureMC::ProcessMCParticles()
{

  // Loop over all primary MC particle  
  for(Long_t i = 0; i < fMCStack->GetNtrack(); i++) {
    // fill primary histograms
    TParticle* particle         = NULL;
    particle                    = (TParticle *)fMCStack->Particle(i);
    if (!particle) continue;
    Bool_t hasMother            = kFALSE;
//     cout << i << "\t"<< particle->GetMother(0) << endl;
    if (particle->GetMother(0)>-1) 
      hasMother                 = kTRUE;
    TParticle* motherParticle   = NULL;
    if( hasMother ) 
      motherParticle            = (TParticle *)fMCStack->Particle(particle->GetMother(0));
    if (motherParticle) 
      hasMother                 = kTRUE;
    else 
      hasMother                 = kFALSE;

    if (!(TMath::Abs(particle->GetPdgCode()) == 111 || TMath::Abs(particle->GetPdgCode()) == 221 || TMath::Abs(particle->GetPdgCode()) == 331 ||
      TMath::Abs(particle->GetPdgCode()) == 223 || TMath::Abs(particle->GetPdgCode()) == 211 || TMath::Abs(particle->GetPdgCode()) == 113 || TMath::Abs(particle->GetPdgCode()) == 333 || TMath::Abs(particle->GetPdgCode()) == 443 || TMath::Abs(particle->GetPdgCode()) == 3212 || TMath::Abs(particle->GetPdgCode()) == 310 || TMath::Abs(particle->GetPdgCode()) == 2224 || TMath::Abs(particle->GetPdgCode()) == 2214 || TMath::Abs(particle->GetPdgCode()) == 1114 || TMath::Abs(particle->GetPdgCode()) == 2114 || TMath::Abs(particle->GetPdgCode()) == 213 || TMath::Abs(particle->GetPdgCode()) == 313 || TMath::Abs(particle->GetPdgCode()) == 130 || TMath::Abs(particle->GetPdgCode()) == 3122 )  )
      continue;
    
    if (!(TMath::Abs(particle->Energy()-particle->Pz())>0.)) continue;
    Double_t yPre = (particle->Energy()+particle->Pz())/(particle->Energy()-particle->Pz());
//     cout << i << "\t"<< particle->GetPdgCode() << "\t"<< particle->Pz() << "\t" << particle->Energy()<< "\t" << particle->Energy()-particle->Pz() << "\t"<< yPre << endl;
    if (yPre == 0.) continue;
    
    Double_t y = 0.5*TMath::Log(yPre);
    
    
    if (y > 1.000) continue;
    switch(particle->GetPdgCode()){
      case 111:
        fHistPtYPi0->Fill(particle->Pt(), particle->Y());
        if (hasMother){
          if (TMath::Abs(motherParticle->GetPdgCode()) == 310 ||  // K0s
            TMath::Abs(motherParticle->GetPdgCode()) == 130 ||  // K0l
            TMath::Abs(motherParticle->GetPdgCode()) == 321  // K+/-
            )
            fHistPtYPi0FromK->Fill(particle->Pt(), particle->Y());
          if (TMath::Abs(motherParticle->GetPdgCode()) == 3122   // Lambda
            ) 
            fHistPtYPi0FromLambda->Fill(particle->Pt(), particle->Y());
          if (motherParticle->GetPdgCode() == 221)  // eta
            fHistPtYPi0FromEta->Fill(particle->Pt(), particle->Y());
        }
        break;
      case 221:
        fHistPtYEta->Fill(particle->Pt(), particle->Y());
        break;
      case 331:
        fHistPtYEtaPrim->Fill(particle->Pt(), particle->Y());
        break;
      case 223:
        fHistPtYOmega->Fill(particle->Pt(), particle->Y());
        break;
      case 211:
        fHistPtYPiPl->Fill(particle->Pt(), particle->Y());
        if (hasMother){
          if (TMath::Abs(motherParticle->GetPdgCode()) == 310 ||  // K0s
            TMath::Abs(motherParticle->GetPdgCode()) == 130 ||  // K0l
            TMath::Abs(motherParticle->GetPdgCode()) == 321  // K+/-
            )
            fHistPtYPiPlFromK->Fill(particle->Pt(), particle->Y());
        }
        break;
      case -211:
        fHistPtYPiMi->Fill(particle->Pt(), particle->Y());
        if (hasMother){
          if (TMath::Abs(motherParticle->GetPdgCode()) == 310 ||  // K0s
            TMath::Abs(motherParticle->GetPdgCode()) == 130 ||  // K0l
            TMath::Abs(motherParticle->GetPdgCode()) == 321  // K+/-
            )
            fHistPtYPiMiFromK->Fill(particle->Pt(), particle->Y());
        }
        break;
      case 113:
        fHistPtYRho0->Fill(particle->Pt(), particle->Y());
        break;
      case 213:
        fHistPtYRhoPl->Fill(particle->Pt(), particle->Y());
        break;
      case -213:
        fHistPtYRhoMi->Fill(particle->Pt(), particle->Y());
        break;
      case 333:
        fHistPtYPhi->Fill(particle->Pt(), particle->Y());
        break;
      case 443:
        fHistPtYJPsi->Fill(particle->Pt(), particle->Y());
        break;
      case 3212:
        fHistPtYSigma0->Fill(particle->Pt(), particle->Y());
        break;
      case 310:
        fHistPtYK0s->Fill(particle->Pt(), particle->Y());
        break;
      case 130:
        fHistPtYK0l->Fill(particle->Pt(), particle->Y());
        break;
      case 313:
        fHistPtYK0star->Fill(particle->Pt(), particle->Y());
        break;
      case 2224:
        fHistPtYDeltaPlPl->Fill(particle->Pt(), particle->Y());
        break;
      case 2214:
        fHistPtYDeltaPl->Fill(particle->Pt(), particle->Y());
        break;
      case 1114:
        fHistPtYDeltaMi->Fill(particle->Pt(), particle->Y());
        break;
      case 2114:
        fHistPtYDelta0->Fill(particle->Pt(), particle->Y());
        break;
      case 3122:
        fHistPtYLambda->Fill(particle->Pt(), particle->Y());
        break;
    }
    
    // from here on, we are only intested in particles considered primaries in ALICE
    if ((particle->GetPdgCode()==111 || particle->GetPdgCode()==221) && hasMother){
      if (TMath::Abs(motherParticle->GetPdgCode()) == 310 ||  // K0s
        TMath::Abs(motherParticle->GetPdgCode()) == 130  ||  // K0l
        TMath::Abs(motherParticle->GetPdgCode()) == 321  ||  // K+/-
        TMath::Abs(motherParticle->GetPdgCode()) == 3122 ||   // Lambdas
        TMath::Abs(motherParticle->GetPdgCode()) == 3212 ||   // Sigma0
        TMath::Abs(motherParticle->GetPdgCode()) == 3222 ||   // Sigmas
        TMath::Abs(motherParticle->GetPdgCode()) == 3112 ||   // Sigmas
        TMath::Abs(motherParticle->GetPdgCode()) == 3322 ||   // Cascades
        TMath::Abs(motherParticle->GetPdgCode()) == 3312    // Cascades
        )
        continue;
    }
    
    // just looking at pi0, etas, etaprims
    if (particle->GetPdgCode()==111 || particle->GetPdgCode()==221 || particle->GetPdgCode() == 331){
      if (particle->GetNDaughters() != 2) continue;   // only the two particle decays
      Int_t acceptanceGamma[2] = {0,0};
      Double_t energyGamma[2] = {0,0};
      Bool_t allOK[2] = {kFALSE,kFALSE};
      for(Int_t i=0;i<2;i++){
        TParticle *daughter=fMCStack->Particle(particle->GetDaughter(i));
        if (!daughter) continue;
        // Is Daughter a Photon?
        if(daughter->GetPdgCode() == 22) allOK[i] =kTRUE;
        IsInPCMAcceptance(daughter, acceptanceGamma[i]);
        IsInPHOSAcceptance(daughter, acceptanceGamma[i]);
        IsInEMCalAcceptance(daughter, acceptanceGamma[i]);
        energyGamma[i] = daughter->Energy();
      }
      if (!(allOK[0] && allOK[1])) continue;
      
      Double_t alpha = (energyGamma[0]-energyGamma[1])/(energyGamma[0]+energyGamma[1]);
      
      if (particle->GetPdgCode()==111){
        fHistPtYPi0GG->Fill(particle->Pt(), particle->Y());
        if (acceptanceGamma[0] > 0 && acceptanceGamma[1] > 0 ){
          fHistPtYPi0GGPCMAcc->Fill(particle->Pt(), particle->Y());
          fHistPtAlphaPi0GGPCMAcc->Fill(particle->Pt(), TMath::Abs(alpha));
        }
        if (acceptanceGamma[0] == 2 && acceptanceGamma[1] == 2 ){
          fHistPtYPi0GGEMCAcc->Fill(particle->Pt(), particle->Y());
          fHistPtAlphaPi0GGEMCAcc->Fill(particle->Pt(), TMath::Abs(alpha));
        }
        if (acceptanceGamma[0] == 3 && acceptanceGamma[1] == 3 ){
          fHistPtYPi0GGPHOAcc->Fill(particle->Pt(), particle->Y());
          fHistPtAlphaPi0GGPHOAcc->Fill(particle->Pt(), TMath::Abs(alpha));
        }
        if ( (acceptanceGamma[0] > 0 && acceptanceGamma[1] == 2) || 
          (acceptanceGamma[1] > 0 && acceptanceGamma[0] == 2)
        ){
          fHistPtYPi0GGPCMEMCAcc->Fill(particle->Pt(), particle->Y());
          if (acceptanceGamma[1]!=2) alpha = (energyGamma[1]-energyGamma[0])/(energyGamma[0]+energyGamma[1]);
          fHistPtAlphaPi0GGPCMEMCAcc->Fill(particle->Pt(), alpha);
        }
        if ( (acceptanceGamma[0] > 0 && acceptanceGamma[1] == 3) || 
          (acceptanceGamma[1] > 0 && acceptanceGamma[0] == 3)
        ){
          fHistPtYPi0GGPCMPHOAcc->Fill(particle->Pt(), particle->Y());
          if (acceptanceGamma[1]!=3) alpha = (energyGamma[1]-energyGamma[0])/(energyGamma[0]+energyGamma[1]);
          fHistPtAlphaPi0GGPCMPHOAcc->Fill(particle->Pt(), alpha);
        }
      }
      if (particle->GetPdgCode()==221){
        fHistPtYEtaGG->Fill(particle->Pt(), particle->Y());
        if (acceptanceGamma[0] > 0 && acceptanceGamma[1] > 0 ){
          fHistPtYEtaGGPCMAcc->Fill(particle->Pt(), particle->Y());
          fHistPtAlphaEtaGGPCMAcc->Fill(particle->Pt(), TMath::Abs(alpha));
        }
        if (acceptanceGamma[0] == 2 && acceptanceGamma[1] == 2 ){
          fHistPtYEtaGGEMCAcc->Fill(particle->Pt(), particle->Y());
          fHistPtAlphaEtaGGEMCAcc->Fill(particle->Pt(), TMath::Abs(alpha));
        }
        if (acceptanceGamma[0] == 3 && acceptanceGamma[1] == 3 ){
          fHistPtYEtaGGPHOAcc->Fill(particle->Pt(), particle->Y());
          fHistPtAlphaEtaGGPHOAcc->Fill(particle->Pt(), TMath::Abs(alpha));
        }
        if ( (acceptanceGamma[0] > 0 && acceptanceGamma[1] == 2) || 
          (acceptanceGamma[1] > 0 && acceptanceGamma[0] == 2)
        ){
          fHistPtYEtaGGPCMEMCAcc->Fill(particle->Pt(), particle->Y());
          if (acceptanceGamma[1]!=2) alpha = (energyGamma[1]-energyGamma[0])/(energyGamma[0]+energyGamma[1]);
          fHistPtAlphaEtaGGPCMEMCAcc->Fill(particle->Pt(), alpha);
        }
        if ( (acceptanceGamma[0] > 0 && acceptanceGamma[1] == 3) || 
          (acceptanceGamma[1] > 0 && acceptanceGamma[0] == 3)
        ){
          fHistPtYEtaGGPCMPHOAcc->Fill(particle->Pt(), particle->Y());
          if (acceptanceGamma[1]!=3) alpha = (energyGamma[1]-energyGamma[0])/(energyGamma[0]+energyGamma[1]);
          fHistPtAlphaEtaGGPCMPHOAcc->Fill(particle->Pt(), alpha);
        }
      }
      if (particle->GetPdgCode()==331){
        fHistPtYEtaPrimGG->Fill(particle->Pt(), particle->Y());
        if (acceptanceGamma[0] > 0 && acceptanceGamma[1] > 0 )
          fHistPtYEtaPrimGGPCMAcc->Fill(particle->Pt(), particle->Y());
        if (acceptanceGamma[0] == 2 && acceptanceGamma[1] == 2 )
          fHistPtYEtaPrimGGEMCAcc->Fill(particle->Pt(), particle->Y());
        if (acceptanceGamma[0] == 3 && acceptanceGamma[1] == 3 )
          fHistPtYEtaPrimGGPHOAcc->Fill(particle->Pt(), particle->Y());
        if ( (acceptanceGamma[0] > 0 && acceptanceGamma[1] == 2) || 
          (acceptanceGamma[1] > 0 && acceptanceGamma[0] == 2)  
        )
          fHistPtYEtaPrimGGPCMEMCAcc->Fill(particle->Pt(), particle->Y());
        if ( (acceptanceGamma[0] > 0 && acceptanceGamma[1] == 3) || 
          (acceptanceGamma[1] > 0 && acceptanceGamma[0] == 3)
        )
          fHistPtYEtaPrimGGPCMPHOAcc->Fill(particle->Pt(), particle->Y());
      }
      
      
    }
    
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskGammaPureMC::IsInPCMAcceptance(TParticle* part, Int_t& acceptance){
  Double_t boundaryEta = 0.900001;
  if (part->Pt() > 0.050 && TMath::Abs(part->Eta()) < boundaryEta)
    acceptance = 1;

  return;  
}

//________________________________________________________________________
void AliAnalysisTaskGammaPureMC::IsInPHOSAcceptance(TParticle* part, Int_t& acceptance){
  Double_t boundaryEtaMin = -0.13;
  Double_t boundaryEtaMax = 0.13;
  Double_t boundaryPhiMin = 4.54;
  Double_t boundaryPhiMax = 5.59;
  if (part->Pt() < 0.300) return;
  if (part->Eta() > boundaryEtaMax || part->Eta() < boundaryEtaMin) return;
  if (part->Phi() > boundaryPhiMax || part->Phi() < boundaryPhiMin) return;  
    acceptance = 3;
  return;  
}

//________________________________________________________________________
void AliAnalysisTaskGammaPureMC::IsInEMCalAcceptance(TParticle* part, Int_t& acceptance){
  Double_t boundaryEtaMin = -0.6687;
  Double_t boundaryEtaMax = 0.66465;
  Double_t boundaryPhiMin = 1.39626;
  Double_t boundaryPhiMax = 3.15;
  if (part->Pt() < 0.400) return;
  if (part->Eta() > boundaryEtaMax || part->Eta() < boundaryEtaMin) return;
  if (part->Phi() > boundaryPhiMax || part->Phi() < boundaryPhiMin) return;  
    acceptance = 2;
  return;  
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
