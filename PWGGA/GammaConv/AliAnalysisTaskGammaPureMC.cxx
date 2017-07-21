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
  fHistPtYEtaPrim(nullptr),
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
  fHistPtYEtaPrimGG(nullptr),
  fHistPtYEtaPrimGGPCMAcc(nullptr),
  fHistPtYEtaPrimGGEMCAcc(nullptr),
  fHistPtYEtaPrimGGPHOAcc(nullptr),
  fHistPtYEtaPrimGGPCMEMCAcc(nullptr),
  fHistPtYEtaPrimGGPCMPHOAcc(nullptr),
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
  fIsK0(1),
  fIsMC(1)
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
  fHistPtYEtaPrim(nullptr),
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
  fHistPtYEtaPrimGG(nullptr),
  fHistPtYEtaPrimGGPCMAcc(nullptr),
  fHistPtYEtaPrimGGEMCAcc(nullptr),
  fHistPtYEtaPrimGGPHOAcc(nullptr),
  fHistPtYEtaPrimGGPCMEMCAcc(nullptr),
  fHistPtYEtaPrimGGPCMPHOAcc(nullptr),
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
  fIsK0(1),
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

  fHistPtYPi0                		= new TH2F("Pt_Y_Pi0","Pt_Y_Pi0", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPi0->Sumw2();
  fOutputContainer->Add(fHistPtYPi0);

  fHistPtYPiPl                		= new TH2F("Pt_Y_PiPl","Pt_Y_PiPl", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPiPl->Sumw2();
  fOutputContainer->Add(fHistPtYPiPl);

  fHistPtYPiMi                		= new TH2F("Pt_Y_PiMi","Pt_Y_PiMi", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPiMi->Sumw2();
  fOutputContainer->Add(fHistPtYPiMi);

  fHistPtYEta                 		= new TH2F("Pt_Y_Eta","Pt_Y_Eta", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEta->Sumw2();
  fOutputContainer->Add(fHistPtYEta);

  fHistPtYEtaPrim             		= new TH2F("Pt_Y_EtaPrim","Pt_Y_EtaPrim", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaPrim->Sumw2();
  fOutputContainer->Add(fHistPtYEtaPrim);

  fHistPtYOmega               		= new TH2F("Pt_Y_Omega","Pt_Y_Omega", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYOmega->Sumw2();
  fOutputContainer->Add(fHistPtYOmega);

  fHistPtYRho0                		= new TH2F("Pt_Y_Rho0","Pt_Y_Rho0", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYRho0->Sumw2();
  fOutputContainer->Add(fHistPtYRho0);

  fHistPtYRhoPl               		= new TH2F("Pt_Y_RhoPl","Pt_Y_RhoPl", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYRhoPl->Sumw2();
  fOutputContainer->Add(fHistPtYRhoPl);

  fHistPtYRhoMi               		= new TH2F("Pt_Y_RhoMi","Pt_Y_RhoMi", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYRhoMi->Sumw2();
  fOutputContainer->Add(fHistPtYRhoMi);

  fHistPtYPhi                 		= new TH2F("Pt_Y_Phi","Pt_Y_Phi", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPhi->Sumw2();
  fOutputContainer->Add(fHistPtYPhi);

  fHistPtYJPsi                		= new TH2F("Pt_Y_JPsi","Pt_Y_JPsi", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYJPsi->Sumw2();
  fOutputContainer->Add(fHistPtYJPsi);

  fHistPtYSigma0              		= new TH2F("Pt_Y_Sigma0","Pt_Y_Sigma0", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYSigma0->Sumw2();
  fOutputContainer->Add(fHistPtYSigma0);

  fHistPtYK0s                 		= new TH2F("Pt_Y_K0s","Pt_Y_K0s", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYK0s->Sumw2();
  fOutputContainer->Add(fHistPtYK0s);

  fHistPtYK0l                 		= new TH2F("Pt_Y_K0l","Pt_Y_K0l", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYK0l->Sumw2();
  fOutputContainer->Add(fHistPtYK0l);

  fHistPtYK0star              		= new TH2F("Pt_Y_K0star","Pt_Y_K0star", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYK0star->Sumw2();
  fOutputContainer->Add(fHistPtYK0star);

  fHistPtYDeltaPlPl           		= new TH2F("Pt_Y_DeltaPlPl","Pt_Y_DeltaPlPl", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYDeltaPlPl->Sumw2();
  fOutputContainer->Add(fHistPtYDeltaPlPl);

  fHistPtYDeltaPl             		= new TH2F("Pt_Y_DeltaPl","Pt_Y_DeltaPl", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYDeltaPl->Sumw2();
  fOutputContainer->Add(fHistPtYDeltaPl);

  fHistPtYDeltaMi             		= new TH2F("Pt_Y_DeltaMi","Pt_Y_DeltaMi", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYDeltaMi->Sumw2();
  fOutputContainer->Add(fHistPtYDeltaMi);

  fHistPtYDelta0              		= new TH2F("Pt_Y_Delta0","Pt_Y_Delta0", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYDelta0->Sumw2();
  fOutputContainer->Add(fHistPtYDelta0);

  fHistPtYLambda              		= new TH2F("Pt_Y_Lambda","Pt_Y_Lambda", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYLambda->Sumw2();
  fOutputContainer->Add(fHistPtYLambda);

  fHistPtYPi0FromEta          		= new TH2F("Pt_Y_Pi0FromEta","Pt_Y_Pi0FromEta", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPi0FromEta->Sumw2();
  fOutputContainer->Add(fHistPtYPi0FromEta);

  fHistPtYPi0FromLambda       		= new TH2F("Pt_Y_Pi0FromLambda","Pt_Y_Pi0FromLambda", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPi0FromLambda->Sumw2();
  fOutputContainer->Add(fHistPtYPi0FromLambda);

  fHistPtYPi0FromK            		= new TH2F("Pt_Y_Pi0FromK","Pt_Y_Pi0FromK", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPi0FromK->Sumw2();
  fOutputContainer->Add(fHistPtYPi0FromK);

  fHistPtYPiPlFromK           		= new TH2F("Pt_Y_PiPlFromK","Pt_Y_PiPlFromK", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPiPlFromK->Sumw2();
  fOutputContainer->Add(fHistPtYPiPlFromK);

  fHistPtYPiMiFromK           		= new TH2F("Pt_Y_PiMiFromK","Pt_Y_PiMiFromK", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPiMiFromK->Sumw2();
  fOutputContainer->Add(fHistPtYPiMiFromK);


  fHistPtYPi0GG               		= new TH2F("Pt_Y_Pi0GG","Pt_Y_Pi0GG", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPi0GG->Sumw2();
  fOutputContainer->Add(fHistPtYPi0GG);
  fHistPtYPi0GGPCMAcc         		= new TH2F("Pt_Y_Pi0GGPCMAcc","Pt_Y_Pi0GGPCMAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPi0GGPCMAcc->Sumw2();
  fOutputContainer->Add(fHistPtYPi0GGPCMAcc);
  fHistPtYPi0GGEMCAcc         		= new TH2F("Pt_Y_Pi0GGEMCAcc","Pt_Y_Pi0GGEMCAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPi0GGEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtYPi0GGEMCAcc);
  fHistPtYPi0GGPHOAcc         		= new TH2F("Pt_Y_Pi0GGPHOAcc","Pt_Y_Pi0GGPHOAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPi0GGPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtYPi0GGPHOAcc);
  fHistPtYPi0GGPCMEMCAcc      		= new TH2F("Pt_Y_Pi0GGPCMEMCAcc","Pt_Y_Pi0GGPCMEMCAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPi0GGPCMEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtYPi0GGPCMEMCAcc);
  fHistPtYPi0GGPCMPHOAcc      		= new TH2F("Pt_Y_Pi0GGPCMPHOAcc","Pt_Y_Pi0GGPCMPHOAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYPi0GGPCMPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtYPi0GGPCMPHOAcc);

  fHistPtYEtaGG               		= new TH2F("Pt_Y_EtaGG","Pt_Y_EtaGG", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaGG->Sumw2();
  fOutputContainer->Add(fHistPtYEtaGG);
  fHistPtYEtaGGPCMAcc         		= new TH2F("Pt_Y_EtaGGPCMAcc","Pt_Y_EtaGGPCMAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaGGPCMAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaGGPCMAcc);
  fHistPtYEtaGGEMCAcc         		= new TH2F("Pt_Y_EtaGGEMCAcc","Pt_Y_EtaGGEMCAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaGGEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaGGEMCAcc);
  fHistPtYEtaGGPHOAcc         		= new TH2F("Pt_Y_EtaGGPHOAcc","Pt_Y_EtaGGPHOAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaGGPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaGGPHOAcc);
  fHistPtYEtaGGPCMEMCAcc      		= new TH2F("Pt_Y_EtaGGPCMEMCAcc","Pt_Y_EtaGGPCMEMCAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaGGPCMEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaGGPCMEMCAcc);
  fHistPtYEtaGGPCMPHOAcc      		= new TH2F("Pt_Y_EtaGGPCMPHOAcc","Pt_Y_EtaGGPCMPHOAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaGGPCMPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaGGPCMPHOAcc);

  fHistPtYEtaPrimGG           		= new TH2F("Pt_Y_EtaPrimGG","Pt_Y_EtaPrimGG", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaPrimGG->Sumw2();
  fOutputContainer->Add(fHistPtYEtaPrimGG);
  fHistPtYEtaPrimGGPCMAcc     		= new TH2F("Pt_Y_EtaPrimGGPCMAcc","Pt_Y_EtaPrimGGPCMAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaPrimGGPCMAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaPrimGGPCMAcc);
  fHistPtYEtaPrimGGEMCAcc     		= new TH2F("Pt_Y_EtaPrimGGEMCAcc","Pt_Y_EtaPrimGGEMCAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaPrimGGEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaPrimGGEMCAcc);
  fHistPtYEtaPrimGGPHOAcc     		= new TH2F("Pt_Y_EtaPrimGGPHOAcc","Pt_Y_EtaPrimGGPHOAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaPrimGGPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaPrimGGPHOAcc);
  fHistPtYEtaPrimGGPCMEMCAcc  		= new TH2F("Pt_Y_EtaPrimGGPCMEMCAcc","Pt_Y_EtaPrimGGPCMEMCAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaPrimGGPCMEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaPrimGGPCMEMCAcc);
  fHistPtYEtaPrimGGPCMPHOAcc  		= new TH2F("Pt_Y_EtaPrimGGPCMPHOAcc","Pt_Y_EtaPrimGGPCMPHOAcc", 1000,0, 100, 200, -1.0, 1.0);
  fHistPtYEtaPrimGGPCMPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtYEtaPrimGGPCMPHOAcc);

  fHistPtAlphaPi0GGPCMAcc     		= new TH2F("Pt_Alpha_Pi0GGPCMAcc","Pt_Alpha_Pi0GGPCMAcc", 500,0.1, 50, 100, 0., 1.);
  SetLogBinningXTH2(fHistPtAlphaPi0GGPCMAcc);
  fHistPtAlphaPi0GGPCMAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaPi0GGPCMAcc);
  fHistPtAlphaPi0GGEMCAcc     		= new TH2F("Pt_Alpha_Pi0GGEMCAcc","Pt_Alpha_Pi0GGEMCAcc", 500,0.1, 50, 100, 0., 1.);
  SetLogBinningXTH2(fHistPtAlphaPi0GGEMCAcc);
  fHistPtAlphaPi0GGEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaPi0GGEMCAcc);
  fHistPtAlphaPi0GGPHOAcc     		= new TH2F("Pt_Alpha_Pi0GGPHOAcc","Pt_Alpha_Pi0GGPHOAcc", 500,0.1, 50, 100, 0., 1.);
  SetLogBinningXTH2(fHistPtAlphaPi0GGPHOAcc);
  fHistPtAlphaPi0GGPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaPi0GGPHOAcc);
  fHistPtAlphaPi0GGPCMEMCAcc  		= new TH2F("Pt_Alpha_Pi0GGPCMEMCAcc","Pt_Alpha_Pi0GGPCMEMCAcc", 500,0.1, 50, 200, -1., 1.);
  SetLogBinningXTH2(fHistPtAlphaPi0GGPCMEMCAcc);
  fHistPtAlphaPi0GGPCMEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaPi0GGPCMEMCAcc);
  fHistPtAlphaPi0GGPCMPHOAcc  		= new TH2F("Pt_Alpha_Pi0GGPCMPHOAcc","Pt_Alpha_Pi0GGPCMPHOAcc", 500,0.1, 50, 200, -1., 1.);
  SetLogBinningXTH2(fHistPtAlphaPi0GGPCMPHOAcc);
  fHistPtAlphaPi0GGPCMPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaPi0GGPCMPHOAcc);

  fHistPtAlphaEtaGGPCMAcc     		= new TH2F("Pt_Alpha_EtaGGPCMAcc","Pt_Alpha_EtaGGPCMAcc", 500,0.1, 50, 100, 0., 1.);
  SetLogBinningXTH2(fHistPtAlphaEtaGGPCMAcc);
  fHistPtAlphaEtaGGPCMAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaEtaGGPCMAcc);
  fHistPtAlphaEtaGGEMCAcc     		= new TH2F("Pt_Alpha_EtaGGEMCAcc","Pt_Alpha_EtaGGEMCAcc", 500,0.1, 50, 100, 0., 1.);
  SetLogBinningXTH2(fHistPtAlphaEtaGGEMCAcc);
  fHistPtAlphaEtaGGEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaEtaGGEMCAcc);
  fHistPtAlphaEtaGGPHOAcc     		= new TH2F("Pt_Alpha_EtaGGPHOAcc","Pt_Alpha_EtaGGPHOAcc", 500,0.1, 50, 100, 0., 1.);
  SetLogBinningXTH2(fHistPtAlphaEtaGGPHOAcc);
  fHistPtAlphaEtaGGPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaEtaGGPHOAcc);
  fHistPtAlphaEtaGGPCMEMCAcc  		= new TH2F("Pt_Alpha_EtaGGPCMEMCAcc","Pt_Alpha_EtaGGPCMEMCAcc", 500,0.1, 50, 200, -1., 1.);
  SetLogBinningXTH2(fHistPtAlphaEtaGGPCMEMCAcc);
  fHistPtAlphaEtaGGPCMEMCAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaEtaGGPCMEMCAcc);
  fHistPtAlphaEtaGGPCMPHOAcc  		= new TH2F("Pt_Alpha_EtaGGPCMPHOAcc","Pt_Alpha_EtaGGPCMPHOAcc", 500,0.1, 50, 200, -1., 1.);
  SetLogBinningXTH2(fHistPtAlphaEtaGGPCMPHOAcc);
  fHistPtAlphaEtaGGPCMPHOAcc->Sumw2();
  fOutputContainer->Add(fHistPtAlphaEtaGGPCMPHOAcc);

  if (fIsK0 == 1){
        fHistPtYPi0FromKGG          		= new TH2F("Pt_Y_Pi0FromKGG","Pt_Y_Pi0FromKGG", 1000,0, 100, 200, -1.0, 1.0);
        fHistPtYPi0FromKGG->Sumw2();
        fOutputContainer->Add(fHistPtYPi0FromKGG);
        fHistPtYPi0FromKGGPCMAcc    		= new TH2F("Pt_Y_Pi0FromKGGPCMAcc","Pt_Y_Pi0FromKGGPCMAcc", 1000,0, 100, 200, -1.0, 1.0);
        fHistPtYPi0FromKGGPCMAcc->Sumw2();
        fOutputContainer->Add(fHistPtYPi0FromKGGPCMAcc);
        fHistPtYPi0FromKGGEMCAcc    		= new TH2F("Pt_Y_Pi0FromKGGEMCAcc","Pt_Y_Pi0FromKGGEMCAcc", 1000,0, 100, 200, -1.0, 1.0);
        fHistPtYPi0FromKGGEMCAcc->Sumw2();
        fOutputContainer->Add(fHistPtYPi0FromKGGEMCAcc);
        fHistPtYPi0FromKGGPCMEMCAcc 		= new TH2F("Pt_Y_Pi0FromKGGPCMEMCAcc","Pt_Y_Pi0FromKGGPCMEMCAcc", 1000,0, 100, 200, -1.0, 1.0);
        fHistPtYPi0FromKGGPCMEMCAcc->Sumw2();
        fOutputContainer->Add(fHistPtYPi0FromKGGPCMEMCAcc);
        fHistPtYPi0FromKGGEMCPCMAcc 		= new TH2F("Pt_Y_Pi0FromKGGEMCPCMAcc","Pt_Y_Pi0FromKGGEMCPCMAcc", 1000,0, 100, 200, -1.0, 1.0);
        fHistPtYPi0FromKGGEMCPCMAcc->Sumw2();
        fOutputContainer->Add(fHistPtYPi0FromKGGEMCPCMAcc);
        fHistPtYPi0FromKGGEMCAccSamePi0    = new TH2F("Pt_Y_Pi0FromKGGEMCAccSamePi0","Pt_Y_Pi0FromKGGEMCAccSamePi0", 1000,0, 100, 200, -1.0, 1.0);
        fHistPtYPi0FromKGGEMCAccSamePi0->Sumw2();
        fOutputContainer->Add(fHistPtYPi0FromKGGEMCAccSamePi0);
        fHistPtYPi0FromKGGEMCAccDiffPi0  	= new TH2F("Pt_Y_Pi0FromKGGEMCAccDiffPi0","Pt_Y_Pi0FromKGGEMCAccDiffPi0", 1000,0, 100, 200, -1.0, 1.0);
        fHistPtYPi0FromKGGEMCAccDiffPi0->Sumw2();
        fOutputContainer->Add(fHistPtYPi0FromKGGEMCAccDiffPi0);
        fHistPtAlphaPi0FromKGG              = new TH2F("Pt_Alpha_Pi0FromKGG","Pt_Alpha_Pi0FromKGG", 500,0.1, 50, 200, -1., 1.);
    	SetLogBinningXTH2(fHistPtAlphaPi0FromKGG);
        fHistPtAlphaPi0FromKGG->Sumw2();
        fOutputContainer->Add(fHistPtAlphaPi0FromKGG);
        fHistPtAlphaPi0FromKGGPCMAcc        = new TH2F("Pt_Alpha_Pi0FromKGGPCMAcc","Pt_Alpha_Pi0FromKGGPCMAcc", 500,0.1, 50, 200, -1., 1.);
    	SetLogBinningXTH2(fHistPtAlphaPi0FromKGGPCMAcc);
        fHistPtAlphaPi0FromKGGPCMAcc->Sumw2();
        fOutputContainer->Add(fHistPtAlphaPi0FromKGGPCMAcc);
        fHistPtAlphaPi0FromKGGEMCAcc      = new TH2F("Pt_Alpha_Pi0FromKGGEMCAcc","Pt_Alpha_Pi0FromKGGEMCAcc", 500,0.1, 50, 200, -1., 1.);
    	SetLogBinningXTH2(fHistPtAlphaPi0FromKGGEMCAcc);
        fHistPtAlphaPi0FromKGGEMCAcc->Sumw2();
        fOutputContainer->Add(fHistPtAlphaPi0FromKGGEMCAcc);
        fHistPtAlphaPi0FromKGGPCMEMCAcc      = new TH2F("Pt_Alpha_Pi0FromKGGPCMEMCAcc","Pt_Alpha_Pi0FromKGGPCMEMCAcc", 500,0.1, 50, 200, -1., 1.);
    	SetLogBinningXTH2(fHistPtAlphaPi0FromKGGPCMEMCAcc);
        fHistPtAlphaPi0FromKGGPCMEMCAcc->Sumw2();
        fOutputContainer->Add(fHistPtAlphaPi0FromKGGPCMEMCAcc);
        fHistPtAlphaPi0FromKGGEMCPCMAcc      = new TH2F("Pt_Alpha_Pi0FromKGGEMCPCMAcc","Pt_Alpha_Pi0FromKGGEMCPCMAcc", 500,0.1, 50, 200, -1., 1.);
    	SetLogBinningXTH2(fHistPtAlphaPi0FromKGGEMCPCMAcc);
        fHistPtAlphaPi0FromKGGEMCPCMAcc->Sumw2();
        fOutputContainer->Add(fHistPtAlphaPi0FromKGGEMCPCMAcc);
        fHistPtAlphaPi0FromKGGEMCAccSamePi0      = new TH2F("Pt_Alpha_Pi0FromKGGEMCAccSamePi0","Pt_Alpha_Pi0FromKGGEMCAccSamePi0", 500,0.1, 50, 200, -1., 1.);
    	SetLogBinningXTH2(fHistPtAlphaPi0FromKGGEMCAccSamePi0);
        fHistPtAlphaPi0FromKGGEMCAccSamePi0->Sumw2();
        fOutputContainer->Add(fHistPtAlphaPi0FromKGGEMCAccSamePi0);
        fHistPtAlphaPi0FromKGGEMCAccDiffPi0      = new TH2F("Pt_Alpha_Pi0FromKGGEMCAccDiffPi0","Pt_Alpha_Pi0FromKGGEMCAccDiffPi0", 500,0.1, 50, 200, -1., 1.);
    	SetLogBinningXTH2(fHistPtAlphaPi0FromKGGEMCAccDiffPi0);
        fHistPtAlphaPi0FromKGGEMCAccDiffPi0->Sumw2();
        fOutputContainer->Add(fHistPtAlphaPi0FromKGGEMCAccDiffPi0);
  }



  PostData(1, fOutputContainer);
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaPureMC::UserExec(Option_t *)
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
void AliAnalysisTaskGammaPureMC::ProcessMCParticles()
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

    const std::array<int, 18> kAcceptPdgCodes = {kPdgPi0, kPdgEta, kPdgEtaPrime, kPdgOmega, kPdgPiPlus, kPdgRho0, kPdgPhi, kPdgJPsi, kPdgSigma0, kPdgK0Short, kPdgDeltaPlus, kPdgDeltaPlusPlus, kPdgDeltaMinus, kPdgDelta0, kPdgRhoPlus, kPdgKStar, kPdgK0Long, kPdgLambda};
    if(std::find(kAcceptPdgCodes.begin(), kAcceptPdgCodes.end(), TMath::Abs(particle->GetPdgCode())) ==  kAcceptPdgCodes.end()) continue;  // species not supported

    if (!(TMath::Abs(particle->Energy()-particle->Pz())>0.)) continue;
    Double_t yPre = (particle->Energy()+particle->Pz())/(particle->Energy()-particle->Pz());
    //     cout << i << "\t"<< particle->GetPdgCode() << "\t"<< particle->Pz() << "\t" << particle->Energy()<< "\t" << particle->Energy()-particle->Pz() << "\t"<< yPre << endl;
    if (yPre == 0.) continue;

    Double_t y = 0.5*TMath::Log(yPre);


    if (y > 1.000) continue;
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
      fHistPtYEtaPrim->Fill(particle->Pt(), particle->Y());
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
    }

    // from here on, we are only intested in particles considered primaries in ALICE
    if ((particle->GetPdgCode()== kPdgPi0 || particle->GetPdgCode()== kPdgEta) && hasMother){
      if (TMath::Abs(motherParticle->GetPdgCode()) == kPdgK0Short ||
          TMath::Abs(motherParticle->GetPdgCode()) == kPdgK0Long  ||
          TMath::Abs(motherParticle->GetPdgCode()) == kPdgKPlus  ||
          TMath::Abs(motherParticle->GetPdgCode()) == kPdgLambda ||
          TMath::Abs(motherParticle->GetPdgCode()) == kPdgSigma0 ||
          TMath::Abs(motherParticle->GetPdgCode()) == kPdgSigmaPlus ||
          TMath::Abs(motherParticle->GetPdgCode()) == kPdgSigmaMinus ||
          TMath::Abs(motherParticle->GetPdgCode()) == kPdgXi0 ||
          TMath::Abs(motherParticle->GetPdgCode()) == kPdgXiMinus
      )
        continue;
    }

    // just looking at pi0, etas, etaprims
    if (particle->GetPdgCode()==kPdgPi0 || particle->GetPdgCode()==kPdgEta || particle->GetPdgCode() == kPdgEtaPrime){
      if (particle->GetNDaughters() != 2) continue;   // only the two particle decays
      UChar_t acceptanceGamma[2] = {0,0};
      Double_t energyGamma[2] = {0,0};
      Bool_t allOK[2] = {kFALSE,kFALSE};


      for(Int_t j=0;j<2;j++){
        TParticle *daughter=fMCEvent->Particle(particle->GetDaughter(j));
        if (!daughter) continue;

        // Is Daughter a Photon? 
        if(daughter->GetPdgCode() == 22) allOK[j] =kTRUE;
        if(IsInPCMAcceptance(daughter))  SETBIT(acceptanceGamma[j], kPCMAcceptance);
        if(IsInPHOSAcceptance(daughter)) SETBIT(acceptanceGamma[j], kPHOSAcceptance);
        if(IsInEMCalAcceptance(daughter)) SETBIT(acceptanceGamma[j], kEMCALAcceptance);
        energyGamma[j] = daughter->Energy();


      }

      if (!(allOK[0] && allOK[1])) continue;

      Double_t alpha = (energyGamma[0]-energyGamma[1])/(energyGamma[0]+energyGamma[1]);

      if (particle->GetPdgCode()==kPdgPi0){
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
      if (particle->GetPdgCode()==kPdgEta){
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
      if (particle->GetPdgCode()==kPdgEtaPrime){
        fHistPtYEtaPrimGG->Fill(particle->Pt(), particle->Y());
        if (TESTBIT(acceptanceGamma[0], kPCMAcceptance) && TESTBIT(acceptanceGamma[1], kPCMAcceptance))
          fHistPtYEtaPrimGGPCMAcc->Fill(particle->Pt(), particle->Y());
        if (TESTBIT(acceptanceGamma[0], kEMCALAcceptance) && TESTBIT(acceptanceGamma[1], kEMCALAcceptance))
          fHistPtYEtaPrimGGEMCAcc->Fill(particle->Pt(), particle->Y());
        if (TESTBIT(acceptanceGamma[0], kPHOSAcceptance) && TESTBIT(acceptanceGamma[1], kPHOSAcceptance))
          fHistPtYEtaPrimGGPHOAcc->Fill(particle->Pt(), particle->Y());
        if ( (TESTBIT(acceptanceGamma[0], kPCMAcceptance)  && TESTBIT(acceptanceGamma[1], kEMCALAcceptance)) ||
            (TESTBIT(acceptanceGamma[1], kPCMAcceptance) && TESTBIT(acceptanceGamma[0], kEMCALAcceptance))
        )
          fHistPtYEtaPrimGGPCMEMCAcc->Fill(particle->Pt(), particle->Y());
        if ( (TESTBIT(acceptanceGamma[0], kPCMAcceptance) && TESTBIT(acceptanceGamma[1], kPHOSAcceptance)) ||
            (TESTBIT(acceptanceGamma[1], kPCMAcceptance) && TESTBIT(acceptanceGamma[0], kPHOSAcceptance))
        )
          fHistPtYEtaPrimGGPCMPHOAcc->Fill(particle->Pt(), particle->Y());
      }
    }


    if(fIsK0 == 0) continue;
    if( particle->GetPdgCode() == kPdgK0Short){
      if (particle->GetNDaughters() != 2) continue;   
      UChar_t acceptanceGamma[2] = {0,0};
      Double_t energyPi0[2] = {0,0};
      Bool_t allOK[2] = {kFALSE,kFALSE};
      UChar_t gdAcceptanceGamma[4] = {0,0,0,0};
      Double_t gdEnergyGamma[4] = {0,0,0,0};
      Bool_t allGDOK[4] = {kFALSE, kFALSE, kFALSE,kFALSE};
      for(Int_t k=0;k<2;k++){
        TParticle *daughter=fMCEvent->Particle(particle->GetDaughter(k));
        if (!daughter) continue;
        
        // Is Daughter a pi0?
        if (daughter->GetPdgCode() == kPdgPi0){
          allOK[k] = kTRUE;
          if(daughter->GetNDaughters() != 2) continue;
          energyPi0[k] = daughter->Energy();
          for(Int_t h=0;h<2;h++){
            TParticle *granddaughter = fMCEvent->Particle(daughter->GetDaughter(k));
            if(granddaughter->GetPdgCode() == 22) allGDOK[2*k + h] = kTRUE;
            if(IsInPCMAcceptance(granddaughter))  SETBIT(gdAcceptanceGamma[2*k+h], kPCMAcceptance);
            if(IsInEMCalAcceptance(granddaughter)) SETBIT(gdAcceptanceGamma[2*k+h], kEMCALAcceptance);
            gdEnergyGamma[2*k+h] = granddaughter->Energy();
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
bool AliAnalysisTaskGammaPureMC::IsInPCMAcceptance(TParticle* part) const {
  const Double_t kBoundaryEta = 0.900001;
  if (part->Pt() > 0.050 && TMath::Abs(part->Eta()) < kBoundaryEta) return true;

  return false;
}

//________________________________________________________________________
bool AliAnalysisTaskGammaPureMC::IsInPHOSAcceptance(TParticle* part) const {
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
bool AliAnalysisTaskGammaPureMC::IsInEMCalAcceptance(TParticle* part) const {
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
