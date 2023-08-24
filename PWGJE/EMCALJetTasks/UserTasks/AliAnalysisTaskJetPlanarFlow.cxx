/*************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// C++
#include <array>
#include <sstream>

// Root
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TGrid.h>
#include <THashList.h>
#include <THnSparse.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <TParticlePDG.h>
#include <TRandom3.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TKey.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TTree.h>
#include <TChain.h>

// Aliroot general
#include "AliAnalysisDataSlot.h"
#include "AliEMCALGeometry.h"
#include "AliVEventHandler.h"


// Aliroot EMCal jet framework
#include "AliEmcalJetTask.h"
#include "AliFJWrapper.h"
#include "AliEmcalJetFinder.h"
#include "AliPicoTrack.h"
#include "FJ_includes.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliInputEventHandler.h"
#include "AliEmcalJet.h"
#include "AliEmcalParticle.h"
#include "AliEmcalPythiaInfo.h"
#include "AliGenPythiaEventHeader.h"
#include "AliJetContainer.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliParticleContainer.h"
#include "AliRhoParameter.h"
#include "AliVTrack.h"
#include "AliVCluster.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDSymEigen.h"
#include <AliAnalysisDataContainer.h>
#include "AliAODEvent.h"

#include "AliAnalysisTaskJetPlanarFlow.h"


using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskJetPlanarFlow)

    //________________________________________________________________________
    AliAnalysisTaskJetPlanarFlow::AliAnalysisTaskJetPlanarFlow()
    : AliAnalysisTaskEmcalJet("AliAnalysisTaskJetPlanarFlow", kTRUE),
      fJetAnalysisType(kData),
      fJetSubType(kNoSub),
      fJetRadius(0.4),
      fMinJetPt(10.0),
      fMaxJetPt(100.0),
      fMinJetConstiteuntPAPt(0.15),
      fMaxJetConstiteuntPAPt(100.0),
      fSharedFractionPtMin(0.5),
      fTrackingEfficiency(1.0),
      fCentralitySelection(kFALSE),
      fCentralityMin(0.0),
      fCentralityMax(10.0),
      fFillNsubjettiness(kTRUE),
      fFillDeltaR(kTRUE),
      fDoRotations(kTRUE),
      fRandom(0),
      fJetConstituentLabels(0),
      fShapesVar_Particles_E(0),
      fShapesVar_Particles_E_Truth(0),
      fShapesVar_Particles_pT(0),
      fShapesVar_Particles_pT_Truth(0),
      fShapesVar_Particles_Phi(0),
      fShapesVar_Particles_Phi_Truth(0),
      fShapesVar_Particles_Theta(0),
      fShapesVar_Particles_Theta_Truth(0),
      fShapesVar_Particles_InJet(0),
      fShapesVar_Particles_InJet_Truth(0),
      fShapesVar_Particles_DeltaR(0),
      fShapesVar_Particles_DeltaR_Truth(0),
      fShapesVar_Particles_NRPhi(0),
      fShapesVar_Particles_NRPhi_Truth(0),
      fShapesVar_Particles_Eta(0),
      fShapesVar_Particles_Eta_Truth(0),
      fTreeJet(0),
      fTreeParticles(0),
      fhEvent(0x0)

{
  for (Int_t i = 0; i < nVar; i++) {
    fShapesVar[i] = 0;
  }
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskJetPlanarFlow::AliAnalysisTaskJetPlanarFlow(const char* name)
    : AliAnalysisTaskEmcalJet(name, kTRUE),
      fJetAnalysisType(kData),
      fJetSubType(kNoSub),
      fJetRadius(0.4),
      fMinJetPt(10.0),
      fMaxJetPt(100.0),
      fMinJetConstiteuntPAPt(0.15),
      fMaxJetConstiteuntPAPt(100.0),
      fSharedFractionPtMin(0.5),
      fTrackingEfficiency(1.0),
      fCentralitySelection(kFALSE),
      fCentralityMin(0.0),
      fCentralityMax(10.0),
      fFillNsubjettiness(kTRUE),
      fFillDeltaR(kTRUE),
      fDoRotations(kTRUE),
      fRandom(0),
      fJetConstituentLabels(0),
      fShapesVar_Particles_E(0),
      fShapesVar_Particles_E_Truth(0),
      fShapesVar_Particles_pT(0),
      fShapesVar_Particles_pT_Truth(0),
      fShapesVar_Particles_Phi(0),
      fShapesVar_Particles_Phi_Truth(0),
      fShapesVar_Particles_Theta(0),
      fShapesVar_Particles_Theta_Truth(0),
      fShapesVar_Particles_InJet(0),
      fShapesVar_Particles_InJet_Truth(0),
      fShapesVar_Particles_DeltaR(0),
      fShapesVar_Particles_DeltaR_Truth(0),
      fShapesVar_Particles_NRPhi(0),
      fShapesVar_Particles_NRPhi_Truth(0),
      fShapesVar_Particles_Eta(0),
      fShapesVar_Particles_Eta_Truth(0),
      fTreeJet(0),
      fTreeParticles(0),
      fhEvent(0x0) {
  // Standard constructor.
  for (Int_t i = 0; i < nVar; i++) {
    fShapesVar[i] = 0;
  }
  SetMakeGeneralHistograms(kTRUE);
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskJetPlanarFlow::~AliAnalysisTaskJetPlanarFlow() {
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskJetPlanarFlow::UserCreateOutputObjects() {
  // Create user output.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  TH1::AddDirectory(oldStatus);
  // create a tree used for the MC data and making a 4D response matrix

  const char *nameoutput = GetOutputSlot(2)->GetContainer()->GetName();
  fTreeJet = new TTree(nameoutput, nameoutput);
  TString *fShapesVarNames = new TString[nVar];
  fShapesVarNames[0] = "pT_Jet";
  fShapesVarNames[1] = "pT_Jet_Truth";
  fShapesVarNames[2] = "Mass_Jet";
  fShapesVarNames[3] = "Mass_Jet_Truth";
  fShapesVarNames[4] = "Eta_Jet";
  fShapesVarNames[5] = "Eta_Jet_Truth";
  fShapesVarNames[6] = "DeltaR_Jet";
  fShapesVarNames[7] = "DeltaR_Jet_Truth";
  fShapesVarNames[8] = "Tau2to1_Jet";
  fShapesVarNames[9] = "Tau2to1_Jet_Truth";
  fShapesVarNames[10] = "PlanarFlow_Jet";
  fShapesVarNames[11] = "PlanarFlow_Jet_Truth";

  for (Int_t ivar = 0; ivar < nVar; ivar++) {
    cout << "looping over variables" << endl;
    fTreeJet->Branch(fShapesVarNames[ivar].Data(), &fShapesVar[ivar], Form("%s/D", fShapesVarNames[ivar].Data()));
  }

  const char *nameoutput_Particles = GetOutputSlot(3)->GetContainer()->GetName();
  fTreeParticles = new TTree(nameoutput_Particles, nameoutput_Particles);
  TString *fShapesVarNames_Particles = new TString[nVar_Particles];
  fShapesVarNames_Particles[0] = "E";
  fShapesVarNames_Particles[1] = "E_Truth";
  fShapesVarNames_Particles[2] = "pT";
  fShapesVarNames_Particles[3] = "pT_Truth";
  fShapesVarNames_Particles[4] = "Phi";
  fShapesVarNames_Particles[5] = "Phi_Truth";
  fShapesVarNames_Particles[6] = "Theta";
  fShapesVarNames_Particles[7] = "Theta_Truth";
  fShapesVarNames_Particles[8] = "InJet";
  fShapesVarNames_Particles[9] = "InJet_Truth";
  fShapesVarNames_Particles[10] = "DeltaR";
  fShapesVarNames_Particles[11] = "DeltaR_Truth";
  fShapesVarNames_Particles[12] = "NRPhi";
  fShapesVarNames_Particles[13] = "NRPhi_Truth";
  fShapesVarNames_Particles[14] = "Eta";
  fShapesVarNames_Particles[15] = "Eta_Truth";
  fTreeParticles->Branch(fShapesVarNames_Particles[0].Data(), &fShapesVar_Particles_E, 0, 1);
  fTreeParticles->Branch(fShapesVarNames_Particles[1].Data(), &fShapesVar_Particles_E_Truth, 0, 1);
  fTreeParticles->Branch(fShapesVarNames_Particles[2].Data(), &fShapesVar_Particles_pT, 0, 1);
  fTreeParticles->Branch(fShapesVarNames_Particles[3].Data(), &fShapesVar_Particles_pT_Truth, 0, 1);
  fTreeParticles->Branch(fShapesVarNames_Particles[4].Data(), &fShapesVar_Particles_Phi, 0, 1);
  fTreeParticles->Branch(fShapesVarNames_Particles[5].Data(), &fShapesVar_Particles_Phi_Truth, 0, 1);
  fTreeParticles->Branch(fShapesVarNames_Particles[6].Data(), &fShapesVar_Particles_Theta, 0, 1);
  fTreeParticles->Branch(fShapesVarNames_Particles[7].Data(), &fShapesVar_Particles_Theta_Truth, 0, 1);
  fTreeParticles->Branch(fShapesVarNames_Particles[8].Data(), &fShapesVar_Particles_InJet, 0, 1);
  fTreeParticles->Branch(fShapesVarNames_Particles[9].Data(), &fShapesVar_Particles_InJet_Truth, 0, 1);
  fTreeParticles->Branch(fShapesVarNames_Particles[10].Data(), &fShapesVar_Particles_DeltaR, 0, 1);
  fTreeParticles->Branch(fShapesVarNames_Particles[11].Data(), &fShapesVar_Particles_DeltaR_Truth, 0, 1);
  fTreeParticles->Branch(fShapesVarNames_Particles[10].Data(), &fShapesVar_Particles_NRPhi, 0, 1);
  fTreeParticles->Branch(fShapesVarNames_Particles[11].Data(), &fShapesVar_Particles_NRPhi_Truth, 0, 1);
  fTreeParticles->Branch(fShapesVarNames_Particles[10].Data(), &fShapesVar_Particles_Eta, 0, 1);
  fTreeParticles->Branch(fShapesVarNames_Particles[11].Data(), &fShapesVar_Particles_Eta_Truth, 0, 1);

  fhEvent = new TH1D("fhEvent", "fhEvent", 40, -0.5, 39.5);
  fOutput->Add(fhEvent);

  PostData(1, fOutput);
  PostData(2, fTreeJet);
  PostData(3, fTreeParticles);
  // delete [] fShapesVarNames;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetPlanarFlow::Run() {
  // Run analysis code here, if needed. It will be executed before FillHistograms().

  return kTRUE;
}

Float_t RelativePhi(Float_t mphi, Float_t vphi) {
  if (mphi < -TMath::Pi()){
    mphi += TMath::TwoPi();
  }
  else if (mphi > TMath::Pi()){
    mphi -= TMath::TwoPi();
  }
  if (vphi < -TMath::Pi()){
    vphi += TMath::TwoPi();
  }
  else if (vphi > TMath::Pi()){
    vphi -= TMath::TwoPi();
  }
  Float_t dphi = mphi - vphi;
  if (dphi < -TMath::Pi()){
    dphi += TMath::TwoPi();
  }
  else if (dphi > TMath::Pi()){
    dphi -= TMath::TwoPi();
  }
  return dphi;  // returns them in the range [-pi,pi]
}


void AliAnalysisTaskJetPlanarFlow::SetTree(AliEmcalJet *jet, AliJetContainer *jetContainer, AliTrackContainer *trackContainer, Float_t jetPt, Int_t level) {
  fShapesVar_Particles_E.clear();
  fShapesVar_Particles_pT.clear();
  fShapesVar_Particles_Phi.clear();
  fShapesVar_Particles_Theta.clear();
  fShapesVar_Particles_InJet.clear();
  fShapesVar_Particles_DeltaR.clear();
  fShapesVar_Particles_NRPhi.clear();
  fShapesVar_Particles_Eta.clear();

  fShapesVar_Particles_E_Truth.clear();
  fShapesVar_Particles_pT_Truth.clear();
  fShapesVar_Particles_Phi_Truth.clear();
  fShapesVar_Particles_Theta_Truth.clear();
  fShapesVar_Particles_InJet_Truth.clear();
  fShapesVar_Particles_DeltaR_Truth.clear();
  fShapesVar_Particles_NRPhi_Truth.clear();
  fShapesVar_Particles_Eta_Truth.clear();

      
  AliVParticle *jetConstituent = NULL;
  fJetConstituentLabels.clear();
  Float_t jetMass = TMath::Sqrt((jet->E() * jet->E()) - (jet->Pt() * jet->Pt()) - (jet->Pz() * jet->Pz()));

  Float_t Result_NSub1 = 10.0;
  Float_t Result_NSub2 = -100.0;
  Float_t deltaR = -10.0;
  if (fFillNsubjettiness) {
    AliEmcalJetFinder JetFinderNSub1("NSubjettiness1");
    JetFinderNSub1.SetJetMaxEta(0.9 - fJetRadius);
    JetFinderNSub1.SetRadius(fJetRadius);
    JetFinderNSub1.SetJetAlgorithm(0);
    JetFinderNSub1.SetRecombSheme(0);
    JetFinderNSub1.SetJetMinPt(jet->Pt());

    Double_t dVtxNSub1[3] = {1, 1, 1};
    Result_NSub1 = JetFinderNSub1.Nsubjettiness(jet, jetContainer, dVtxNSub1, 1, 0, fJetRadius, 0.0, 0, 0, 0.0, 0.1, 1);

    AliEmcalJetFinder JetFinderNSub2("NSubjettiness2");
    JetFinderNSub2.SetJetMaxEta(0.9 - fJetRadius);
    JetFinderNSub2.SetRadius(fJetRadius);
    JetFinderNSub2.SetJetAlgorithm(0);
    JetFinderNSub2.SetRecombSheme(0);
    JetFinderNSub2.SetJetMinPt(jet->Pt());

    Double_t dVtxNSub2[3] = {1, 1, 1};
    Result_NSub2 = JetFinderNSub2.Nsubjettiness(jet, jetContainer, dVtxNSub2, 2, 0, fJetRadius, 0.0, 0, 0, 0.0, 0.1, 1);
  }
  if (fFillDeltaR) {
    AliEmcalJetFinder JetFinderNSubdR("NSubjettinessdR");
    JetFinderNSubdR.SetJetMaxEta(0.9 - fJetRadius);
    JetFinderNSubdR.SetRadius(fJetRadius);
    JetFinderNSubdR.SetJetAlgorithm(0);
    JetFinderNSubdR.SetRecombSheme(0);
    JetFinderNSubdR.SetJetMinPt(jet->Pt());

    Double_t dVtxNSubdR[3] = {1, 1, 1};
    deltaR = JetFinderNSubdR.Nsubjettiness(jet, jetContainer, dVtxNSubdR, 2, 0, fJetRadius, 0.0, 2, 0, 0.0, 0.1, 1);
  }

  Float_t tau2to1 = Result_NSub2 / Result_NSub1;

  fShapesVar[0 + level] = jetPt;
  fShapesVar[2 + level] = jetMass;
  fShapesVar[4 + level] = jet->Eta();
  fShapesVar[6 + level] = deltaR;
  fShapesVar[8 + level] = tau2to1;

  Float_t thetaTrack = -1.0;
  Float_t phiTrack = -1.0;
  Float_t rotationMatrix[3][3];
  Float_t rotationMatrix2D[2][2];

  if (fDoRotations) {
    Float_t jetUnitVector[3] = {Float_t(TMath::Cos(jet->Phi()) / TMath::CosH(jet->Eta())), Float_t(TMath::Sin(jet->Phi()) / TMath::CosH(jet->Eta())), Float_t(TMath::SinH(jet->Eta()) / TMath::CosH(jet->Eta()))};
    Float_t magPt = TMath::Sqrt((jetUnitVector[0] * jetUnitVector[0]) + (jetUnitVector[1] * jetUnitVector[1]));
    Float_t cosTheta = jetUnitVector[2];
    Float_t sinTheta = magPt;
    Float_t cosPhi = TMath::Cos(jet->Phi());
    Float_t sinPhi = TMath::Sin(jet->Phi());

    rotationMatrix[0][0] = -1.0 * cosTheta * cosPhi;
    rotationMatrix[0][1] = -1.0 * cosTheta * sinPhi;
    rotationMatrix[0][2] = sinTheta;
    rotationMatrix[1][0] = sinPhi;
    rotationMatrix[1][1] = -1.0 * cosPhi;
    rotationMatrix[1][2] = 0.;
    rotationMatrix[2][0] = sinTheta * cosPhi;
    rotationMatrix[2][1] = sinTheta * sinPhi;
    rotationMatrix[2][2] = cosTheta;

    Float_t principleMatrix[2][2];
    for (int i = 0; i < 2; i++) {
      for (int j = 0; j < 2; j++) {
        principleMatrix[i][j] = 0.0;
      }
    }

    for (Int_t iConstituent = 0; iConstituent < jet->GetNumberOfTracks(); iConstituent++) {
      jetConstituent = static_cast<AliVParticle *>(jet->TrackAt(iConstituent, jetContainer->GetParticleContainer()->GetArray()));
      fJetConstituentLabels.push_back(jetConstituent->GetLabel());
      if (jetConstituent->Pt() < fMinJetConstiteuntPAPt || jetConstituent->Pt() >= fMaxJetConstiteuntPAPt) continue;
      if (fTrackingEfficiency != 1.0) {
        if (fRandom.Rndm() > fTrackingEfficiency) continue;
      }

      double normalisation_factor = (1.0 / (jetConstituent->E() * jetMass));
      Float_t pxRotated = (rotationMatrix[0][0] * jetConstituent->Px()) + (rotationMatrix[0][1] * jetConstituent->Py()) + (rotationMatrix[0][2] * jetConstituent->Pz());
      Float_t pyRotated = (rotationMatrix[1][0] * jetConstituent->Px()) + (rotationMatrix[1][1] * jetConstituent->Py()) + (rotationMatrix[1][2] * jetConstituent->Pz());
      Float_t pzRotated = (rotationMatrix[2][0] * jetConstituent->Px()) + (rotationMatrix[2][1] * jetConstituent->Py()) + (rotationMatrix[2][2] * jetConstituent->Pz());

      principleMatrix[0][0] += normalisation_factor * pxRotated * pxRotated;
      principleMatrix[0][1] += normalisation_factor * pxRotated * pyRotated;
      principleMatrix[1][0] += normalisation_factor * pyRotated * pxRotated;
      principleMatrix[1][1] += normalisation_factor * pyRotated * pyRotated;
    }

    Float_t principleMatrixTrace = principleMatrix[0][0] + principleMatrix[1][1];
    Float_t PrinciplMatrixDeterminant = (principleMatrix[0][0] * principleMatrix[1][1]) - (principleMatrix[0][1] * principleMatrix[1][0]);
    Float_t eigenValue1 = 0.5 * (principleMatrixTrace + TMath::Sqrt(principleMatrixTrace * principleMatrixTrace - 4 * PrinciplMatrixDeterminant));
    Float_t eigenValue2 = 0.5 * (principleMatrixTrace - TMath::Sqrt(principleMatrixTrace * principleMatrixTrace - 4 * PrinciplMatrixDeterminant));

    fShapesVar[10 + level] = (4.0 * PrinciplMatrixDeterminant) / (principleMatrixTrace * principleMatrixTrace);

    Float_t eigenVector1[2];
    Float_t eigenVector2[2];
    if (principleMatrix[1][0] == 0.0 || principleMatrix[0][1] == 0.0) {
      eigenVector1[0] = principleMatrix[0][0];
      eigenVector1[1] = principleMatrix[1][1];
      eigenVector2[0] = principleMatrix[0][0];
      eigenVector2[1] = principleMatrix[1][1];
    }
    else {
      eigenVector1[0] = eigenValue1 - principleMatrix[1][1];
      eigenVector1[1] = principleMatrix[1][0];
      eigenVector2[0] = principleMatrix[0][1];
      eigenVector2[1] = eigenValue2 - principleMatrix[0][0];
    }
    if (eigenValue1 < eigenValue2) {
      eigenVector1[0] = eigenVector2[0];
      eigenVector1[1] = eigenVector2[1];
    }

    Float_t theta = TMath::ATan(eigenVector1[1] / eigenVector1[0]);
    if (theta < 0) theta += TMath::Pi();
    rotationMatrix2D[0][0] = TMath::Cos(theta);
    rotationMatrix2D[0][1] = TMath::Sin(theta);
    rotationMatrix2D[1][0] = -TMath::Sin(theta);
    rotationMatrix2D[1][1] = TMath::Cos(theta);
  }

  int trackLabel = -1.0;
  for (Int_t iTrack = 0; iTrack < trackContainer->GetNTracks(); iTrack++) {
    AliAODTrack *track = static_cast<AliAODTrack *>(trackContainer->GetAcceptParticle(iTrack));
    if (!track) continue;
    if (TMath::Abs(track->Eta()) > 0.9) continue;
    trackLabel = track->GetLabel();

    Float_t isInJet = 0.0;
    for (Int_t iConstituent = 0; iConstituent < fJetConstituentLabels.size(); iConstituent++) {
      if (fJetConstituentLabels[iConstituent] == trackLabel) {  // Is this correct?
        isInJet = 1.0;
        break;
      }
    }
    Float_t pxRotatedPrincipleAxis = 0.0;
    Float_t pyRotatedPrincipleAxis = 0.0;
    if (fDoRotations) {
      Float_t pxRotated = (rotationMatrix[0][0] * track->Px()) + (rotationMatrix[0][1] * track->Py()) + (rotationMatrix[0][2] * track->Pz());
      Float_t pyRotated = (rotationMatrix[1][0] * track->Px()) + (rotationMatrix[1][1] * track->Py()) + (rotationMatrix[1][2] * track->Pz());
      Float_t pzRotated = (rotationMatrix[2][0] * track->Px()) + (rotationMatrix[2][1] * track->Py()) + (rotationMatrix[2][2] * track->Pz());

      pxRotatedPrincipleAxis = (rotationMatrix2D[0][0] * pxRotated) + (rotationMatrix2D[0][1] * pyRotated);
      pyRotatedPrincipleAxis = (rotationMatrix2D[1][0] * pxRotated) + (rotationMatrix2D[1][1] * pyRotated);
      thetaTrack = TMath::ACos(pzRotated / TMath::Sqrt((pxRotated * pxRotated) + (pyRotated * pyRotated) + (pzRotated * pzRotated)));
      phiTrack = TMath::ATan2(pyRotatedPrincipleAxis, pxRotatedPrincipleAxis);
    }
    else {
      thetaTrack = track->Eta();
      phiTrack = track->Phi();
    }

    Float_t particleDeltaR = TMath::Sqrt(TMath::Power(RelativePhi(track->Phi(),jet->Phi()),2)+TMath::Power((track->Eta()-jet->Eta()),2));

    if (level == 0) {
      fShapesVar_Particles_E.push_back(track->E());
      fShapesVar_Particles_pT.push_back(track->Pt());
      fShapesVar_Particles_Phi.push_back(phiTrack);
      fShapesVar_Particles_Theta.push_back(thetaTrack);
      fShapesVar_Particles_InJet.push_back(isInJet);
      fShapesVar_Particles_DeltaR.push_back(particleDeltaR);
      fShapesVar_Particles_NRPhi.push_back(track->Phi());
      fShapesVar_Particles_Eta.push_back(track->Eta());
    }

    if (level == 1) {
      fShapesVar_Particles_E_Truth.push_back(track->E());
      fShapesVar_Particles_pT_Truth.push_back(track->Pt());
      fShapesVar_Particles_Phi_Truth.push_back(phiTrack);
      fShapesVar_Particles_Theta_Truth.push_back(thetaTrack);
      fShapesVar_Particles_InJet_Truth.push_back(isInJet);
      fShapesVar_Particles_DeltaR_Truth.push_back(particleDeltaR);
      fShapesVar_Particles_NRPhi_Truth.push_back(track->Phi());
      fShapesVar_Particles_Eta_Truth.push_back(track->Eta());
    }
  }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetPlanarFlow::FillHistograms() {
  AliTrackContainer *trackContainer = GetTrackContainer(0);
  AliTrackContainer *trackContainerDet = NULL;
  AliTrackContainer *trackContainerTrue = NULL;

  AliJetContainer *jetContainer = GetJetContainer(0);
  AliJetContainer *jetContainerDet = NULL;
  AliJetContainer *jetContainerTrue = NULL;

  AliEmcalJet *jet = NULL;
  AliEmcalJet *jetDet = NULL;
  AliEmcalJet *jetTrue = NULL;

  if (fJetAnalysisType == kTrueDet) {
    jetContainerTrue = GetJetContainer(1);
    trackContainerTrue = GetTrackContainer(1);
  }
  if (fJetAnalysisType == kEmb) {
    jetContainerDet = GetJetContainer(1);
    trackContainerDet = GetTrackContainer(1);
    jetContainerTrue = GetJetContainer(2);
    trackContainerTrue = GetTrackContainer(2);
  }

  if (fCentralitySelection) {
    if ((fCent > fCentralityMax) || (fCent < fCentralityMin)) return kFALSE;
  }

  if (jetContainer) {
    jetContainer->ResetCurrentID();  //??
    // Jet Loop
    Float_t jetpT = 0.0;
    while ((jet = jetContainer->GetNextAcceptJet())) {
      if (!jet) continue;
      if (fJetSubType == kAreaSub)
        jetpT = jet->Pt() - GetRhoVal(0) * jet->Area();
      else
        jetpT = jet->Pt();
      if (jetpT < fMinJetPt || jetpT >= fMaxJetPt) continue;
      if (fJetAnalysisType == kTrueDet) {
        jetTrue = jet->ClosestJet();
        if (!jetTrue) continue;
        SetTree(jetTrue, jetContainerTrue, trackContainerTrue, jetTrue->Pt(), 1);
      }
      if (fJetAnalysisType == kEmb) {
        if (jetContainer->AliJetContainer::GetFractionSharedPt(jet) < fSharedFractionPtMin) continue;  // how does this work? should we not supply the other container?
        jetDet = jet->ClosestJet();
        if (!jetDet) continue;
        jetTrue = jetDet->ClosestJet();  // do labels need to be checked?
        if (!jetTrue) continue;
        SetTree(jetTrue, jetContainerTrue, trackContainerTrue, jetTrue->Pt(), 1);
      }
      SetTree(jet, jetContainer, trackContainer, jetpT, 0);
      fTreeJet->Fill();
      fTreeParticles->Fill();
    }
  }

  delete trackContainer;
  delete trackContainerTrue;
  delete trackContainerDet;
  delete jetContainer;
  delete jetContainerTrue;
  delete jetContainerDet;
  delete jet;
  delete jetTrue;
  delete jetDet;

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetPlanarFlow::RetrieveEventObjects() {
  //
  // retrieve event objects
  //
  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects()) return kFALSE;

  return kTRUE;
}

//_______________________________________________________________________
void AliAnalysisTaskJetPlanarFlow::Terminate(Option_t*) {
  // Called once at the end of the analysis.
}