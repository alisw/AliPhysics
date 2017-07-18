/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/contrib/nsubjettiness.hh>
#include <fastjet/contrib/softdrop.hh>

#include <TMath.h>
#include <TString.h>

#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEmcalJetSubstructureTree.h"
#include "AliClusterContainer.h"
#include "AliJetContainer.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliEmcalJet.h"
#include "AliEmcalList.h"
#include "AliLog.h"
#include "AliParticleContainer.h"
#include "AliTrackContainer.h"
#include "AliVCluster.h"
#include "AliVParticle.h"

/// \cond CLASSIMP
ClassImp(EmcalTriggerJets::AliAnalysisTaskEmcalJetSubstructureTree)
/// \endcond

namespace EmcalTriggerJets {

AliAnalysisTaskEmcalJetSubstructureTree::AliAnalysisTaskEmcalJetSubstructureTree() :
    AliAnalysisTaskEmcalJet(),
    fJetSubstructureTree(nullptr),
    fJetSubstructureInfo(),
    fSDZCut(0.1),
    fSDBetaCut(0),
    fReclusterizer(kCAAlgo),
    fTriggerSelectionBits(AliVEvent::kAny),
    fTriggerSelectionString("")
{
}

AliAnalysisTaskEmcalJetSubstructureTree::AliAnalysisTaskEmcalJetSubstructureTree(const char *name) :
    AliAnalysisTaskEmcalJet(name, kTRUE),
    fJetSubstructureTree(nullptr),
    fJetSubstructureInfo(),
    fSDZCut(0.1),
    fSDBetaCut(0),
    fReclusterizer(kCAAlgo),
    fTriggerSelectionBits(AliVEvent::kAny),
    fTriggerSelectionString("")
{

}

AliAnalysisTaskEmcalJetSubstructureTree::~AliAnalysisTaskEmcalJetSubstructureTree() {

}

void AliAnalysisTaskEmcalJetSubstructureTree::UserCreateOutputObjects() {
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  fJetSubstructureTree = new TTree("jetSubstructure", "Tree with jet substructure information");
  std::stringstream leaflist;
  leaflist  << "fR/D:"
            << "fEventWeight/D:"
            << "fPtJetRec/D:"
            << "fPtJetSim/D:"
            << "fNCharged/I:"
            << "fNNeutral/I:"
            << "fNTrueConst/I:"
            << "fAreaRec/D:"
            << "fAreaSim/D:"
            << "fNEFRec/D:"
            << "fNEFSim/D:"
            << "fZgMeasured/D:"
            << "fZgTrue/D:"
            << "fRgMeasured/D:"
            << "fRgTrue/D:"
            << "fMgMeasured/D:"
            << "fMgTrue/D:"
            << "fPtgMeasured/D:"
            << "fPtgTrue/D:"
            << "fNDroppedMeasured/I:"
            << "fNDroppedTrue/I:"
            << "fOneSubjettinessMeasured/D:"
            << "fOneSubjettinessTrue/D:"
            << "fTwoSubjettinessMeasured/D:"
            << "fTwoSubjettinessTrue/D";
  fJetSubstructureTree->Branch("JetInfo", &fJetSubstructureInfo, leaflist.str().c_str());
  fOutput->Add(fJetSubstructureTree);
  PostData(1, fOutput);
}

bool AliAnalysisTaskEmcalJetSubstructureTree::Run(){
  AliClusterContainer *clusters = GetClusterContainer("caloClusters");
  AliTrackContainer *tracks = GetTrackContainer("tracks");
  AliParticleContainer *particles = GetParticleContainer("mcparticles");

  AliJetContainer *mcjets = GetJetContainer("mcjets");
  AliJetContainer *datajets = GetJetContainer("datajets");

  // Run trigger selection (not on pure MCgen train)
  if(datajets){
    if(!(fInputHandler->IsEventSelected() & fTriggerSelectionBits)) return false;
    if(fTriggerSelectionString.Length()) {
      if(!fInputEvent->GetFiredTriggerClasses().Contains(fTriggerSelectionString)) return false;
    }
  }

  Double_t weight = 1;

  AliSoftdropDefinition softdropSettings;
  softdropSettings.fBeta = fSDBetaCut;
  softdropSettings.fZ = fSDZCut;
  switch(fReclusterizer) {
  case kCAAlgo: softdropSettings.fRecluserAlgo = fastjet::cambridge_aachen_algorithm; break;
  case kKTAlgo: softdropSettings.fRecluserAlgo = fastjet::kt_algorithm; break;
  case kAKTAlgo: softdropSettings.fRecluserAlgo = fastjet::antikt_algorithm; break;
  };

  AliNSubjettinessDefinition nsubjettinessSettings;
  nsubjettinessSettings.fBeta = 1.;
  nsubjettinessSettings.fRadius = 0.4;

  if(mcjets && !datajets) {
    // pure MC (gen) train - run over MC jets
    for(auto jet : mcjets->accepted()) {
      try {
        AliJetSubstructureData structure = MakeJetSubstructure(*jet, mcjets->GetJetRadius() * 2., particles, nullptr,{softdropSettings, nsubjettinessSettings});
        FillTree(mcjets->GetJetRadius(), weight, nullptr, jet, nullptr, &(structure.fSoftDrop), nullptr, &(structure.fNsubjettiness));
      } catch (ReclusterizerException &e) {
        AliErrorStream() << "Error in reclusterization - skipping jet" << std::endl;
      }
    }
  }

  if(datajets) {
    for(auto jet : datajets->accepted()) {
      AliEmcalJet *associatedJet = jet->MatchedJet();
      if(mcjets) {
        if(!associatedJet) continue;
        try {
          AliJetSubstructureData structureData =  MakeJetSubstructure(*jet, mcjets->GetJetRadius() * 2., particles, nullptr, {softdropSettings, nsubjettinessSettings}),
                                 structureMC = MakeJetSubstructure(*associatedJet, mcjets->GetJetRadius() * 2, particles, nullptr, {softdropSettings, nsubjettinessSettings});
          FillTree(datajets->GetJetRadius(), weight, jet, associatedJet, &(structureData.fSoftDrop), &(structureMC.fSoftDrop), &(structureData.fNsubjettiness), &(structureMC.fNsubjettiness));
        } catch(ReclusterizerException &e) {
          AliErrorStream() << "Error in reclusterization - skipping jet" << std::endl;
        }
      } else {
        try {
          AliJetSubstructureData structure = MakeJetSubstructure(*jet, 0.4, tracks, clusters, {softdropSettings, nsubjettinessSettings});
          FillTree(datajets->GetJetRadius(), weight, jet, nullptr, &(structure.fSoftDrop), nullptr, &(structure.fNsubjettiness), nullptr);
        } catch(ReclusterizerException &e) {
          AliErrorStream() << "Error in reclusterization - skipping jet" << std::endl;
        }
      }
    }
  }

  return true;
}

void AliAnalysisTaskEmcalJetSubstructureTree::FillTree(double r, double weight,
                                                       const AliEmcalJet *datajet, const AliEmcalJet *mcjet,
                                                       AliSoftDropParameters *dataSoftdrop, AliSoftDropParameters *mcSoftdrop,
                                                       AliNSubjettinessParameters *dataSubjettiness, AliNSubjettinessParameters *mcSubjettiness){
  fJetSubstructureInfo.fR = r;
  fJetSubstructureInfo.fEventWeight = weight;
  if(datajet) {
    fJetSubstructureInfo.fPtJetRec = TMath::Abs(datajet->Pt());
    fJetSubstructureInfo.fNCharged = datajet->GetNumberOfTracks();
    fJetSubstructureInfo.fNNeutral = datajet->GetNumberOfClusters();
    fJetSubstructureInfo.fAreaRec = datajet->Area();
    fJetSubstructureInfo.fNEFRec = datajet->NEF();
  } else {
    fJetSubstructureInfo.fPtJetRec = 0.;
    fJetSubstructureInfo.fNCharged = 0;
    fJetSubstructureInfo.fNNeutral = 0;
    fJetSubstructureInfo.fAreaRec = 0.;
    fJetSubstructureInfo.fNEFRec = 0.;
  }

  if(mcjet) {
    fJetSubstructureInfo.fPtJetSim = TMath::Abs(mcjet->Pt());
    fJetSubstructureInfo.fNTrueConst = mcjet->GetNumberOfConstituents();
    fJetSubstructureInfo.fAreaSim = mcjet->Area();
    fJetSubstructureInfo.fNEFSim = mcjet->NEF();
  } else {
    fJetSubstructureInfo.fPtJetSim = 0.;
    fJetSubstructureInfo.fNTrueConst = 0;
    fJetSubstructureInfo.fAreaSim = 0.;
    fJetSubstructureInfo.fNEFSim = 0.;
  }

  if(dataSoftdrop) {
    fJetSubstructureInfo.fZgMeasured = dataSoftdrop->fZg;
    fJetSubstructureInfo.fRgMeasured = dataSoftdrop->fRg;
    fJetSubstructureInfo.fMgMeasured = dataSoftdrop->fMg;
    fJetSubstructureInfo.fPtgMeasured = dataSoftdrop->fPtg;
    fJetSubstructureInfo.fNDroppedMeasured = dataSoftdrop->fNDropped;
  } else {
    fJetSubstructureInfo.fZgMeasured = 0.;
    fJetSubstructureInfo.fRgMeasured = 0.;
    fJetSubstructureInfo.fMgMeasured = 0.;
    fJetSubstructureInfo.fPtgMeasured = 0.;
    fJetSubstructureInfo.fNDroppedMeasured = 0;
  }

  if(mcSoftdrop) {
    fJetSubstructureInfo.fZgTrue = mcSoftdrop->fZg;
    fJetSubstructureInfo.fRgTrue = mcSoftdrop->fRg;
    fJetSubstructureInfo.fMgTrue = mcSoftdrop->fMg;
    fJetSubstructureInfo.fPtgTrue = mcSoftdrop->fPtg;
    fJetSubstructureInfo.fNDroppedTrue = mcSoftdrop->fNDropped;
  } else {
    fJetSubstructureInfo.fZgTrue = 0.;
    fJetSubstructureInfo.fRgTrue = 0.;
    fJetSubstructureInfo.fMgTrue = 0.;
    fJetSubstructureInfo.fPtgTrue = 0.;
    fJetSubstructureInfo.fNDroppedTrue = 0;
  }

  if(dataSubjettiness) {
    fJetSubstructureInfo.fOneSubjettinessMeasured = dataSubjettiness->fOneSubjettiness;
    fJetSubstructureInfo.fTwoSubjettinessMeasured = dataSubjettiness->fTwoSubjettiness;
  } else {
    fJetSubstructureInfo.fOneSubjettinessMeasured = 0.;
    fJetSubstructureInfo.fTwoSubjettinessMeasured = 0.;
  }

  if(mcSubjettiness) {
    fJetSubstructureInfo.fOneSubjettinessTrue = mcSubjettiness->fOneSubjettiness;
    fJetSubstructureInfo.fTwoSubjettinessTrue = mcSubjettiness->fTwoSubjettiness;
  } else {
    fJetSubstructureInfo.fOneSubjettinessTrue = 0.;
    fJetSubstructureInfo.fTwoSubjettinessTrue = 0.;
  }

  fJetSubstructureTree->Fill();
}


AliJetSubstructureData AliAnalysisTaskEmcalJetSubstructureTree::MakeJetSubstructure(const AliEmcalJet &jet, double jetradius, const AliParticleContainer *tracks, const AliClusterContainer *clusters, const AliJetSubstructureSettings &settings) const {
  const int kClusterOffset = 30000; // In order to handle tracks and clusters in the same index space the cluster index needs and offset, large enough so that there is no overlap with track indices
  std::vector<fastjet::PseudoJet> constituents;
  for(int itrk = 0; itrk < jet.GetNumberOfTracks(); itrk++){
    AliVTrack *track = static_cast<AliVTrack *>(jet.TrackAt(itrk, tracks->GetArray()));
    fastjet::PseudoJet constituentTrack(track->Px(), track->Py(), track->Pz(), track->E());
    constituentTrack.set_user_index(jet.TrackAt(itrk));
    constituents.push_back(constituentTrack);
  }

  for(int icl = 0; icl < jet.GetNumberOfClusters(); icl++) {
    AliVCluster *cluster = jet.ClusterAt(icl, clusters->GetArray());
    TLorentzVector clustervec;
    cluster->GetMomentum(clustervec, fVertex);
    fastjet::PseudoJet constituentCluster(clustervec.Px(), clustervec.Py(), clustervec.Pz(), cluster->GetHadCorrEnergy());
    constituentCluster.set_user_index(jet.ClusterAt(icl) + kClusterOffset);
    constituents.push_back(constituentCluster);
  }

  // Redo jet finding on constituents with a
  fastjet::JetDefinition jetdef(fastjet::antikt_algorithm, jetradius*2, static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30 );
  std::vector<fastjet::PseudoJet> outputjets;
  try {
    fastjet::ClusterSequence jetfinder(constituents, jetdef);
    outputjets = jetfinder.inclusive_jets(0);
    AliJetSubstructureData result({MakeSoftDropParameters(outputjets[0], settings.fSoftdropSettings), MakeNsubjettinessParameters(outputjets[0], settings.fSubjettinessSettings)});
    return result;
  } catch (fastjet::Error &e) {
    AliErrorStream() << " FJ Exception caught: " << e.message() << std::endl;
    throw ReclusterizerException();
  }
}

AliSoftDropParameters AliAnalysisTaskEmcalJetSubstructureTree::MakeSoftDropParameters(const fastjet::PseudoJet &jet, const AliSoftdropDefinition &cutparameters) const {
  fastjet::contrib::SoftDrop softdropAlgorithm(cutparameters.fBeta, cutparameters.fZ);
  softdropAlgorithm.set_verbose_structure(kTRUE);
  std::unique_ptr<fastjet::contrib::Recluster> reclusterizer(new fastjet::contrib::Recluster(cutparameters.fRecluserAlgo, 1, true));
  softdropAlgorithm.set_reclustering(kTRUE, reclusterizer.get());
  fastjet::PseudoJet groomed = softdropAlgorithm(jet);

  AliSoftDropParameters result({groomed.structure_of<fastjet::contrib::SoftDrop>().symmetry(),
                                groomed.structure_of<fastjet::contrib::SoftDrop>().delta_R(),
                                groomed.structure_of<fastjet::contrib::SoftDrop>().mu(),
                                groomed.perp(),
                                groomed.structure_of<fastjet::contrib::SoftDrop>().dropped_count()});
  return result;
}

AliNSubjettinessParameters AliAnalysisTaskEmcalJetSubstructureTree::MakeNsubjettinessParameters(const fastjet::PseudoJet &jet, const AliNSubjettinessDefinition &cut) const {
  AliNSubjettinessParameters result({
    fastjet::contrib::Nsubjettiness (1,fastjet::contrib::KT_Axes(),fastjet::contrib::NormalizedMeasure(cut.fBeta,cut.fRadius)).result(jet),
    fastjet::contrib::Nsubjettiness (2,fastjet::contrib::KT_Axes(),fastjet::contrib::NormalizedMeasure(cut.fBeta,cut.fRadius)).result(jet)
  });
  return result;
}

AliAnalysisTaskEmcalJetSubstructureTree *AliAnalysisTaskEmcalJetSubstructureTree::AddEmcalJetSubstructureTreeMaker(Bool_t isMC, Bool_t isData, Double_t jetradius, const char *trigger){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  Bool_t isAOD(kFALSE);
  AliInputEventHandler *inputhandler = static_cast<AliInputEventHandler *>(mgr->GetInputEventHandler());
  if(inputhandler) {
    if(inputhandler->IsA() == AliAODInputHandler::Class()) isAOD = kTRUE;
  }

  AliAnalysisTaskEmcalJetSubstructureTree *treemaker = new AliAnalysisTaskEmcalJetSubstructureTree("JetSubstructureTreemaker_" + TString::Format("R%02d_", int(jetradius * 10.)) + trigger);
  mgr->AddTask(treemaker);
  treemaker->SetMakeGeneralHistograms(kTRUE);
  treemaker->SetVzRange(-10., 10);

  // Adding containers
  if(isMC) {
    AliParticleContainer *particles = treemaker->AddMCParticleContainer("mcparticles");

    AliJetContainer *mcjets = treemaker->AddJetContainer(
                              AliJetContainer::kFullJet,
                              AliJetContainer::antikt_algorithm,
                              AliJetContainer::pt_scheme,
                              jetradius,
                              AliEmcalJet::kEMCALfid,
                              particles, nullptr);
    mcjets->SetName("mcjets");
    mcjets->SetMinPt(20.);
  }

  if(isData) {
    AliTrackContainer *tracks = treemaker->AddTrackContainer(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::TrackContainerNameFactory(isAOD));
    tracks->SetMinPt(0.15);
    AliClusterContainer *clusters = treemaker->AddClusterContainer(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::ClusterContainerNameFactory(isAOD));
    clusters->SetMinE(0.3); // 300 MeV E-cut

    AliJetContainer *datajets = treemaker->AddJetContainer(
                              AliJetContainer::kFullJet,
                              AliJetContainer::antikt_algorithm,
                              AliJetContainer::pt_scheme,
                              jetradius,
                              AliEmcalJet::kEMCALfid,
                              tracks, clusters);
    datajets->SetName("datajets");
    datajets->SetMinPt(20.);

    treemaker->SetUseAliAnaUtils(true, true);

    // configure trigger selection
    TString triggerstring(trigger);
    if(triggerstring.Contains("INT7")) {
      treemaker->SetTriggerBits(AliVEvent::kINT7);
    } else if(triggerstring.Contains("EJ1")) {
      treemaker->SetTriggerBits(AliVEvent::kEMCEJE);
      treemaker->SetTriggerString("EJ1");
    } else if(triggerstring.Contains("EJ2")) {
      treemaker->SetTriggerBits(AliVEvent::kEMCEJE);
      treemaker->SetTriggerString("EJ2");
    }
  }

  // Connecting containers
  TString outputfile = mgr->GetCommonFileName();
  outputfile += TString::Format(":JetSubstructure_R%02d_%s", int(jetradius * 10.), trigger);
  mgr->ConnectInput(treemaker, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(treemaker, 1, mgr->CreateContainer("JetSubstructure_" + TString::Format("R%0d_", int(jetradius * 10.)) + trigger, AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, outputfile));

  return treemaker;
}

} /* namespace EmcalTriggerJets */
