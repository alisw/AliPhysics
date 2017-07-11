/**************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                   *
 * All rights reserved.                                                               *
 *                                                                                    *
 * Redistribution and use in source and binary forms, with or without                 *
 * modification, are permitted provided that the following conditions are met:        *
 *     * Redistributions of source code must retain the above copyright               *
 *       notice, this list of conditions and the following disclaimer.                *
 *     * Redistributions in binary form must reproduce the above copyright            *
 *       notice, this list of conditions and the following disclaimer in the          *
 *       documentation and/or other materials provided with the distribution.         *
 *     * Neither the name of the <organization> nor the                               *
 *       names of its contributors may be used to endorse or promote products         *
 *       derived from this software without specific prior written permission.        *
 *                                                                                    *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND    *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED      *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE             *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY                *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES         *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;       *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND        *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT         *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS      *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                       *
 **************************************************************************************/
#include <array>

#include <THistManager.h>
#include <TLinearBinning.h>
#include <TList.h>
#include "AliAnalysisTaskK0toPi0Pi0.h"
#include "AliAODConversionMother.h"
#include "AliAODConversionPhoton.h"
#include "AliCaloPhotonCuts.h"
#include "AliClusterContainer.h"
#include "AliConversionPhotonCuts.h"
#include "AliConversionMesonCuts.h"
#include "AliConvEventCuts.h"
#include "AliLog.h"
#include "AliVCluster.h"
#include "AliV0ReaderV1.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskK0toPi0Pi0)
/// \endcond

AliAnalysisTaskK0toPi0Pi0::AliAnalysisTaskK0toPi0Pi0():
  AliAnalysisTaskSE(),
  fLocalInitialized(kFALSE),
  fCurrentRun(-1),
  fNewFile(kFALSE),
  fV0Reader(nullptr),
  fV0ReaderName("V0ReaderV1"),
  fClusterContainer(nullptr),
  fIsMC(false),
  fWeightJetJetMC(1.),
  fEventPlaneAngle(0.),
  fEventCuts(nullptr),
  fConvPhotonCuts(nullptr),
  fCaloPhotonCuts(nullptr),
  fPi0Cuts(nullptr),
  fPi0CutsCaloCalo(nullptr),
  fK0Cuts(nullptr),
  fHistos(nullptr),
  fOutput(nullptr)
{

}

AliAnalysisTaskK0toPi0Pi0::AliAnalysisTaskK0toPi0Pi0(const char *name):
  AliAnalysisTaskSE(name),
  fLocalInitialized(kFALSE),
  fCurrentRun(-1),
  fNewFile(kFALSE),
  fV0Reader(nullptr),
  fV0ReaderName("V0ReaderV1"),
  fClusterContainer(nullptr),
  fIsMC(false),
  fWeightJetJetMC(1.),
  fEventPlaneAngle(0),
  fEventCuts(nullptr),
  fConvPhotonCuts(nullptr),
  fCaloPhotonCuts(nullptr),
  fPi0Cuts(nullptr),
  fPi0CutsCaloCalo(nullptr),
  fK0Cuts(nullptr),
  fHistos(nullptr),
  fOutput(nullptr)
{
  DefineOutput(1, TList::Class());
}



AliAnalysisTaskK0toPi0Pi0::~AliAnalysisTaskK0toPi0Pi0() {
  if(fHistos) delete fHistos;
}

void AliAnalysisTaskK0toPi0Pi0::UserCreateOutputObjects(){
  fOutput        = new TList();
  fOutput->SetOwner(kTRUE);
  
  
  // Connecting input V0 reader
  fV0Reader=dynamic_cast<AliV0ReaderV1*>(AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data()));
  if(!fV0Reader){
    AliFatal("Error: No V0 Reader");
  }// GetV0Reader

 // fEventCuts = fV0Reader->GetEventCuts();
 
  
  
  // Define histograms
  fHistos = new THistManager("K0stoPi0Pi0");
  
  AliDebug(2, "************* Defining Event Counter Histograms ********************");
  
  
  // Event counter histograms
  fHistos->CreateTH1("hEventQualityBefore", "Event Quality (0 = good)", 13, -0.5, 12.5);
  fHistos->CreateTH1("hEventQualityAfter", "Event Quality (0 = good)", 13, -0.5, 12.5);
  fHistos->CreateTH1("hVertexZ", "z-component of the primary vertex; z (cm); Number of events", 1000, -40., 40.);
  fHistos->CreateTH1("hCaloPhotonsBefore", "Number of Events", 13, -0.5, 12.5);
  fHistos->CreateTH1("hCaloPhotonsAfter", "Number of Events", 13, -0.5, 12.5);
  fHistos->CreateTH1("hConvPhotonsBefore", "Number of Events", 13, -0.5, 12.5);
  fHistos->CreateTH1("hConvPhotonsAfter", "Number of Events", 13, -0.5, 12.5);
  
  AliDebug(2, "************* Defining Photon QA Histograms ********************");
  
  // Photon QA
  fHistos->CreateTH1("hCaloPhotonPt", "p_{t}-distribution of the conversion photons; p_{t} (GeV); Yield", 300, 0., 30.);
  fHistos->CreateTH1("hConvPhotonPt", "p_{t}-distribution of the conversion photons; p_{t} (GeV); Yield", 300, 0., 30.);
  fHistos->CreateTH2("hConvPhotonEtaR", "#eta vs conversion radius of conversion photons; #eta; R (cm)", 200, -1.5, 1.5, 300, 0., 300);

  AliDebug(2, "************* Defining Pi0 Histograms ********************");
  
  // Pi0 invariant mass, alpha and opening angle distributions
  const std::array<TString, 3> pi0rec = {"ConvConv", "ConvCalo", "CaloCalo"}; // aka PCM, EMCAL, PCM-EMCAL
  for(const auto &reccase : pi0rec){
    // for candidates in a wide mass region
    fHistos->CreateTH2("hMassvsPtPi0" + reccase + "All", "inv. mass vs. p_{t} for all #pi^{0} (" + reccase + ") candidates; inv. mass (GeV/c^{2}); p_{t} (GeV/c)", 500, 0., 0.5, 300, 0.3, 30.);
    // only for candidates in the pi0 mass region
    fHistos->CreateTH2("hMassvsPtPi0" + reccase + "Sel", "inv. mass vs. p_{t} for selected #pi^{0} (" + reccase + ") candidates; inv. mass (GeV/c^{2}); p_{t} (GeV/c)", 500, 0., 0.5, 300, 0.3, 30.);
    fHistos->CreateTH2("hAlphavsPtPi0" + reccase, "#alpha vs p_{t} for selected #p^i{0} (" + reccase + ") candidates; #alpha; p_{t}", 200, -1., 1., 300, 0., 30.);
    fHistos->CreateTH2("hOpeningAnglevsPtPi0" + reccase, "Opening angle vs. p_{t} for selected #pi^{0} (" + reccase + ") candidates; opening angle; p_{t} (GeV/c)", 100, 0., 1., 300., 0.3, 30.);
  }

  AliDebug(2, "************* Defining K0 Histograms ********************");
  
  // K0short invariant mass distribution and opening angle distributions
  const std::array<TString, 6> k0Shortrec = {"AllConv", "AllCalo", "DiffMixed", "SameMixed", "ConvoCalo","CaloConvo" }; // aka 6 cases 
  for(const auto &reccase1 : k0Shortrec){
    fHistos->CreateTH2("hMassvsPtK0Short" + reccase1, "inv. mass vs. p_{t} for #k^{0}s (" + reccase1 + ") candidates; inv. mass (GeV/c^{2}); p_{t} (GeV/c)",  500, 0.3, 0.6, 300, 0.3, 30.); 
    fHistos->CreateTH2("hOpeningAnglevsPtK0Short"+ reccase1, "Opening angle vs. p_{t} for  k0Short (" + reccase1 + ") candidates; opening angle; p_{t} (GeV/c)",  100, 0., 1., 300., 0.3, 30.); 
  }
  
  for(auto hist : *(fHistos->GetListOfHistograms())) fOutput->Add(hist);
  
  // Adding cut QA
  
  TList *qaV0reader = new TList;
  qaV0reader->SetName("QA_V0reader");
  qaV0reader->SetOwner(kTRUE);
  qaV0reader->Add(fV0Reader->GetEventCutHistograms());
  qaV0reader->Add(fV0Reader->GetCutHistograms());
  fOutput->Add(qaV0reader);
  
  fOutput->Add(fEventCuts->GetCutHistograms());
  fOutput->Add(fConvPhotonCuts->GetCutHistograms());
  fOutput->Add(fCaloPhotonCuts->GetCutHistograms());
  
  PostData(1, fOutput);
  
}

void AliAnalysisTaskK0toPi0Pi0::ExecOnce() {
  if(!fClusterContainer) fClusterContainer = new AliClusterContainer(fInputEvent->IsA() == AliESDEvent::Class() ? "CaloClusters" : "caloClusters");
  fClusterContainer->SetArray(fInputEvent);
   
  if(fConvPhotonCuts){
    fConvPhotonCuts->InitPIDResponse();   
  }
   
}

void AliAnalysisTaskK0toPi0Pi0::UserExec(Option_t *){
  if(!fLocalInitialized) {
    ExecOnce();
    fLocalInitialized = kTRUE;
  }
  
  if(fCurrentRun != fInputEvent->GetRunNumber()) {
    RunChanged();
    fCurrentRun = fInputEvent->GetRunNumber();
  }
  
  
  // do event selection
  // Use the same event selection as for the v0 reader
  // Good events defined as events with event quality 0
  Int_t eventQuality = fEventCuts->GetEventQuality();
  fHistos->FillTH1("hEventQualityBefore", eventQuality);
  
  if(fEventCuts->IsEventAcceptedByCut(fV0Reader->GetEventCuts(), fInputEvent, fMCEvent, false, false)) return;
  fHistos->FillTH1("hEventQualityAfter", eventQuality);
  fHistos->FillTH1("hVertexZ", fInputEvent->GetPrimaryVertex()->GetZ());

  
  // get Photon candidates
  fHistos->FillTH1("hConvPhotonsBefore",fV0Reader->GetNReconstructedGammas());
  std::vector<AliAODConversionPhoton> conversionPhotons = MakeConversionPhotonCandidates(*fV0Reader, *fConvPhotonCuts);
  MakePhotonQAConv(conversionPhotons, *fEventCuts);
  Int_t numPhotons = conversionPhotons.size();
  fHistos->FillTH1("hConvPhotonsAfter", numPhotons); 
  
  fHistos->FillTH1("hCaloPhotonsBefore",fClusterContainer->GetNEntries());
  std::vector<AliAODConversionPhoton> caloPhotons = MakeCaloPhotonCandidates(*fClusterContainer, *fCaloPhotonCuts);
  MakePhotonQACalo(caloPhotons, *fEventCuts);
  Int_t numCaloPhotons = caloPhotons.size();
  fHistos->FillTH1("hCaloPhotonsAfter", numCaloPhotons);
  

  
  // get Pi0 candidates
  std::vector<AliAODConversionMother> samePi0PCM = MakePi0Candidates(&conversionPhotons, nullptr, *fPi0Cuts),
                                      samePi0EMCAL = MakePi0Candidates(&caloPhotons, nullptr, *fPi0CutsCaloCalo),
                                      mixedPi0 = MakePi0Candidates(&conversionPhotons, &caloPhotons, *fPi0Cuts);
                                      
  MakePi0QA(samePi0PCM, "ConvConv");
  MakePi0QA(samePi0EMCAL, "CaloCalo");
  MakePi0QA(mixedPi0, "ConvCalo");
  
  
                                      
  // get K0Short candidates
  
  std::vector<AliAODConversionMother> allPCM = MakeK0ShortCandidates(&samePi0PCM, nullptr, *fK0Cuts); 
  std::vector<AliAODConversionMother> allEMC = MakeK0ShortCandidates(&samePi0EMCAL, nullptr, *fK0Cuts);
  std::vector<AliAODConversionMother> PCMEMC = MakeK0ShortCandidates(&samePi0PCM, &mixedPi0, *fK0Cuts);
  std::vector<AliAODConversionMother> EMCPCM = MakeK0ShortCandidates(&samePi0EMCAL, &mixedPi0, *fK0Cuts);
  std::vector<AliAODConversionMother> mixedSame= MakeK0ShortCandidates(&samePi0EMCAL, &samePi0PCM, *fK0Cuts);
  std::vector<AliAODConversionMother> mixedDiff = MakeK0ShortCandidates(&mixedPi0, nullptr, *fK0Cuts);
  
  
  MakeK0ShortQA(allPCM, "AllConv");
  MakeK0ShortQA(allEMC,"AllCalo");
  MakeK0ShortQA(PCMEMC,"ConvoCalo");
  MakeK0ShortQA(EMCPCM,"CaloConvo");
  MakeK0ShortQA(mixedSame, "SameMixed");
  MakeK0ShortQA(mixedDiff, "DiffMixed");
  
  
   
  
  PostData(1, fOutput);
}

std::vector<AliAODConversionPhoton> AliAnalysisTaskK0toPi0Pi0::MakeCaloPhotonCandidates(const AliClusterContainer &inputcont, AliCaloPhotonCuts &cuts){
  std::vector<AliAODConversionPhoton> candidates;
  cuts.FillHistogramsExtendedQA(fInputEvent, fIsMC);

  // vertex
  Double_t vertex[3] = {0,0,0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  int clusterindex = 0;
  for(auto c : inputcont.all()) {
    clusterindex++;
    if(!cuts.ClusterIsSelected(c, fInputEvent, fMCEvent, fIsMC, 1., c->GetID())) continue;   // weight to be added later

    // TLorentzvector with cluster
    TLorentzVector clusterVector;
    c->GetMomentum(clusterVector,vertex);

    // convert to AODConversionPhoton
    AliAODConversionPhoton photonCandidate(&clusterVector);

    // Flag Photon as CaloPhoton
    photonCandidate.SetIsCaloPhoton();
    photonCandidate.SetCaloClusterRef(clusterindex);

    // get MC label
    if(fIsMC>0){
      photonCandidate.SetNCaloPhotonMCLabels(c->GetNLabels());
      for (UInt_t k = 0; k < c->GetNLabels(); k++){
        if(k < 50) photonCandidate.SetCaloPhotonMCLabel(k,c->GetLabels()[k]);
      }
    }

    candidates.push_back(photonCandidate);
  }

  return candidates;
}

std::vector<AliAODConversionPhoton> AliAnalysisTaskK0toPi0Pi0::MakeConversionPhotonCandidates(const AliV0ReaderV1 &reader, AliConversionPhotonCuts &cuts) {
  std::vector<AliAODConversionPhoton> candidates;
  Int_t nV0 = 0;
  std::vector<AliAODConversionPhoton *> GammaCandidatesStepOne;
  TList GammaCandidatesStepTwo;
  
  // Loop over Photon Candidates allocated by ReaderV1  
  for(auto photon : reader){
    if(!cuts.PhotonIsSelected(photon ,fInputEvent)) continue;
    if(!cuts.InPlaneOutOfPlaneCut(photon->GetPhotonPhi(), fEventPlaneAngle)) continue;
    if(!cuts.UseElecSharingCut() && !cuts.UseToCloseV0sCut()){
      candidates.push_back(*(static_cast<AliAODConversionPhoton *>(photon))); // if no second loop is required add to events good gammas
    }else if(cuts.UseElecSharingCut()){ // if Shared Electron cut is enabled, Fill array, add to step one
      cuts.FillElectonLabelArray(static_cast<AliAODConversionPhoton *>(photon), nV0);
      nV0++;
      GammaCandidatesStepOne.push_back(static_cast<AliAODConversionPhoton *>(photon));
    }else if(!cuts.UseElecSharingCut() && cuts.UseToCloseV0sCut()){ // shared electron is disabled, step one not needed -> step two
      GammaCandidatesStepTwo.Add(static_cast<AliAODConversionPhoton *>(photon));
    }
  }

  if(cuts.UseElecSharingCut()){
    int iV0 = 0;
    for(auto photon : GammaCandidatesStepOne){
      if(!cuts.RejectSharedElectronV0s(photon,iV0,GammaCandidatesStepOne.size())) continue;
      if(!cuts.UseToCloseV0sCut()){ // To Close v0s cut diabled, step two not needed
        candidates.push_back(*photon);
      } else GammaCandidatesStepTwo.Add(photon); // Close v0s cut enabled -> add to list two
      iV0++;
    }
  }

  if(cuts.UseToCloseV0sCut()){
    for(int i = 0; i < GammaCandidatesStepTwo.GetEntries(); i++){
      AliAODConversionPhoton *photon = static_cast<AliAODConversionPhoton *>(GammaCandidatesStepTwo.At(i));
      if(!cuts.RejectToCloseV0s(photon, &GammaCandidatesStepTwo,i)) continue;
      candidates.push_back(*photon); // Add gamma to current cut TList
    }
  }
  return candidates;
}

std::vector<AliAODConversionMother> AliAnalysisTaskK0toPi0Pi0::MakePi0Candidates(const std::vector<AliAODConversionPhoton> *primaryLeg,
                                                                                 const std::vector<AliAODConversionPhoton> *secondaryLeg,
                                                                                 AliConversionMesonCuts &cuts){
  // secondary leg: optional argument, if different methods for photon identification (i.e. PCM-EMCAL) is used
  std::vector<AliAODConversionMother> candidates;
  if(secondaryLeg){
    // Different methods for photon identification identification
    for(auto primphoton : *primaryLeg){
      for(auto secphoton : *secondaryLeg) {
        AliAODConversionMother candidate(&primphoton, &secphoton);
        if(!cuts.CheckWhetherInMassRange(candidate.M()))continue;
        // Do Pi0 selection
        if(!cuts.MesonIsSelected(&candidate, kTRUE, 0)) continue;      // Rapidity shift needed when going to asymmetric systems
        candidates.push_back(candidate);
      }
    }
  } else {
    // Same method for photon identification
    for(auto primiter = primaryLeg->begin(); primiter != primaryLeg->end(); ++primiter){
      for(auto seciter = primiter + 1; seciter != primaryLeg->end(); ++seciter){
        AliAODConversionMother candidate(&(*primiter), &(*seciter));
        if(!cuts.CheckWhetherInMassRange(candidate.M()))continue;
        // Do pi0 selection
        if(cuts.MesonIsSelected(&candidate, kTRUE, 0)) continue;      // Rapidity shift needed when going to asymmetric systems
        candidates.push_back(candidate);
      }
    }
  }
  return candidates;
}

std::vector<AliAODConversionMother> AliAnalysisTaskK0toPi0Pi0::MakeK0ShortCandidates(const std::vector<AliAODConversionMother> *primaryLeg,
                                                                                     const std::vector<AliAODConversionMother> *secondaryLeg,
                                                                                    AliConversionMesonCuts &cuts){
  std::vector<AliAODConversionMother> candidates;
  if(secondaryLeg) {
    // Different methods for Pi0 identification (one same, one mixed)
    for(auto primpi0 : *primaryLeg) {
      for(auto secpi0 : *secondaryLeg) {
        AliAODConversionMother candidate(&primpi0, &secpi0);
        if(!cuts.MesonIsSelected(&candidate, kTRUE, 0)) continue;
        candidates.push_back(candidate);
      }
    }
  } else {
    // Same methods for Pi0 identification (both same or both mixed)
    for(auto primpi0 = primaryLeg->begin(); primpi0 != primaryLeg->end(); ++primpi0) {
      for(auto secpi0 = primpi0 + 1; secpi0 != primaryLeg->end(); ++secpi0) {
        AliAODConversionMother candidate(&(*primpi0), &(*secpi0));
        if(!cuts.MesonIsSelected(&candidate, kTRUE, 0)) continue;
        candidates.push_back(candidate);
      }
    }
  }
  return candidates;
}


void AliAnalysisTaskK0toPi0Pi0::MakePhotonQACalo(const std::vector<AliAODConversionPhoton> &photons, AliConvEventCuts &cuts) {
  for(auto photon : photons) {
    fHistos->FillTH1("hCaloPhotonPt", photon.Pt());
  }
}

void AliAnalysisTaskK0toPi0Pi0::MakePhotonQAConv(const std::vector<AliAODConversionPhoton> &photons, AliConvEventCuts &cuts) {
  for(auto photon : photons) {
    fHistos->FillTH1("hConvPhotonPt", photon.Pt());
    fHistos->FillTH2("hConvPhotonEtaR", photon.Eta(), photon.GetConversionRadius());
  }
}

void AliAnalysisTaskK0toPi0Pi0::MakePi0QA(const std::vector<AliAODConversionMother> &pi0s, const char *reccase){
  TString reccaseString = reccase;
  for(auto pi0 : pi0s) {
  	fHistos->FillTH2("hMassvsPtPi0" + reccaseString + "All", pi0.M(), pi0.Pt());
  	// if in the pi0 mass region  
  	if((0.08 <=pi0.M()) && (pi0.M()<= 0.145)){
  		 fHistos->FillTH2("hMassvsPtPi0"+ reccaseString + "Sel", pi0.M(),pi0.Pt());
  		 fHistos->FillTH2("hAlphavsPtPi0" + reccaseString, pi0.GetAlpha(), pi0.Pt());
  		 fHistos->FillTH2("hOpeningAnglevsPtPi0" + reccaseString, pi0.GetOpeningAngle(), pi0.Pt());
  	}
  }
}

void AliAnalysisTaskK0toPi0Pi0::MakeK0ShortQA(const std::vector<AliAODConversionMother> &k0s,const char *reccase){
  TString reccaseString = reccase;
  for(auto k0: k0s) {
  	fHistos->FillTH2("hMassvsPtK0Short" + reccaseString, k0.M(), k0.Pt());
  	fHistos->FillTH2("hOpeningAnglevsPtK0Short" + reccaseString, k0.GetOpeningAngle(), k0.Pt());
  }
}

AliClusterContainer *AliAnalysisTaskK0toPi0Pi0::AddClusterContainer(const char *name) {
  fClusterContainer = new AliClusterContainer(name);
  return fClusterContainer;
}



