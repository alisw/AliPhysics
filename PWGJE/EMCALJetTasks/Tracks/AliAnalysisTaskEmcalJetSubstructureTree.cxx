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
#include <fastjet/contrib/Nsubjettiness.hh>
#include <fastjet/contrib/SoftDrop.hh>

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
    fSDZCut(0.1),
    fSDBetaCut(0),
    fReclusterizer(kCAAlgo),
    fTriggerSelectionBits(AliVEvent::kAny),
    fTriggerSelectionString("")
{
  memset(fJetTreeData, 0, sizeof(Double_t) * kTNVar);
}

AliAnalysisTaskEmcalJetSubstructureTree::AliAnalysisTaskEmcalJetSubstructureTree(const char *name) :
    AliAnalysisTaskEmcalJet(name, kTRUE),
    fJetSubstructureTree(nullptr),
    fSDZCut(0.1),
    fSDBetaCut(0),
    fReclusterizer(kCAAlgo),
    fTriggerSelectionBits(AliVEvent::kAny),
    fTriggerSelectionString("")
{
  memset(fJetTreeData, 0, sizeof(Double_t) * kTNVar);
  DefineOutput(2, TTree::Class());
}

AliAnalysisTaskEmcalJetSubstructureTree::~AliAnalysisTaskEmcalJetSubstructureTree() {

}

void AliAnalysisTaskEmcalJetSubstructureTree::UserCreateOutputObjects() {
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  OpenFile(2);
  fJetSubstructureTree = new TTree("jetSubstructure", "Tree with jet substructure information");
  TString varnames[kTNVar];
  varnames[0] = "Radius";
  varnames[1] = "EventWeight";
  varnames[2] = "PtJetRec";
  varnames[3] = "PtJetSim";
  varnames[4] = "AreaRec";
  varnames[5] = "AreaSim";
  varnames[6] = "NEFRec";
  varnames[7] = "NEFSim";
  varnames[8] = "MassRec";
  varnames[9] = "MassSim";
  varnames[10] = "ZgMeasured";
  varnames[11] = "ZgTrue";
  varnames[12] = "RgMeasured";
  varnames[13] = "RgTrue";
  varnames[14] = "MgMeasured";
  varnames[15] = "MgTrue";
  varnames[16] = "PtgMeasured";
  varnames[17] = "PtgTrue";
  varnames[18] = "OneSubjettinessMeasured";
  varnames[19] = "OneSubjettinessTrue";
  varnames[20] = "TwoSubjettinessMeasured";
  varnames[21] = "TwoSubjettinessTrue";
  varnames[22] = "NCharged";
  varnames[23] = "NNeutral";
  varnames[24] = "NConstTrue";
  varnames[25] = "NDroppedMeasured";
  varnames[26] = "NDroppedTrue";

  for(int ib = 0; ib < kTNVar; ib++){
    fJetSubstructureTree->Branch(varnames[ib], fJetTreeData + ib, Form("%s/D", varnames[ib].Data()));
  }
  PostData(1, fOutput);
  PostData(2, fJetSubstructureTree);
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
  fJetTreeData[kTRadius] = r;
  fJetTreeData[kTWeight] = weight;
  if(datajet) {
    fJetTreeData[kTPtJetRec] = TMath::Abs(datajet->Pt());
    fJetTreeData[kTNCharged] = datajet->GetNumberOfTracks();
    fJetTreeData[kTNNeutral] = datajet->GetNumberOfClusters();
    fJetTreeData[kTAreaRec] = datajet->Area();
    fJetTreeData[kTNEFRec] = datajet->NEF();
    fJetTreeData[kTMassRec] = datajet->M();
  } else {
    fJetTreeData[kTPtJetRec] = 0.;
    fJetTreeData[kTNCharged] = 0;
    fJetTreeData[kTNNeutral] = 0;
    fJetTreeData[kTAreaRec] = 0.;
    fJetTreeData[kTNEFRec] = 0.;
    fJetTreeData[kTMassRec] = 0.;
  }

  if(mcjet) {
    fJetTreeData[kTPtJetSim] = TMath::Abs(mcjet->Pt());
    fJetTreeData[kTNConstTrue] = mcjet->GetNumberOfConstituents();
    fJetTreeData[kTAreaSim] = mcjet->Area();
    fJetTreeData[kTNEFSim] = mcjet->NEF();
    fJetTreeData[kTMassSim] = mcjet->M();
  } else {
    fJetTreeData[kTPtJetSim] = 0.;
    fJetTreeData[kTNConstTrue] = 0;
    fJetTreeData[kTAreaSim] = 0.;
    fJetTreeData[kTNEFSim] = 0.;
    fJetTreeData[kTMassSim] = 0;
  }

  if(dataSoftdrop) {
    fJetTreeData[kTZgMeasured] = dataSoftdrop->fZg;
    fJetTreeData[kTRgMeasured] = dataSoftdrop->fRg;
    fJetTreeData[kTMgMeasured] = dataSoftdrop->fMg;
    fJetTreeData[kTPtgMeasured] = dataSoftdrop->fPtg;
    fJetTreeData[kTNDroppedMeasured] = dataSoftdrop->fNDropped;
  } else {
    fJetTreeData[kTZgMeasured] = 0.;
    fJetTreeData[kTRgMeasured] = 0.;
    fJetTreeData[kTMgMeasured] = 0.;
    fJetTreeData[kTPtgMeasured] = 0.;
    fJetTreeData[kTNDroppedMeasured] = 0;
  }

  if(mcSoftdrop) {
    fJetTreeData[kTZgTrue] = mcSoftdrop->fZg;
    fJetTreeData[kTRgTrue] = mcSoftdrop->fRg;
    fJetTreeData[kTMgTrue] = mcSoftdrop->fMg;
    fJetTreeData[kTPtgTrue] = mcSoftdrop->fPtg;
    fJetTreeData[kTNDroppedTrue] = mcSoftdrop->fNDropped;
  } else {
    fJetTreeData[kTZgTrue] = 0.;
    fJetTreeData[kTRgTrue] = 0.;
    fJetTreeData[kTMgTrue] = 0.;
    fJetTreeData[kTPtgTrue] = 0.;
    fJetTreeData[kTNDroppedTrue] = 0;
  }

  if(dataSubjettiness) {
    fJetTreeData[kTOneNSubjettinessMeasured] = dataSubjettiness->fOneSubjettiness;
    fJetTreeData[kTTwoNSubjettinessMeasured] = dataSubjettiness->fTwoSubjettiness;
  } else {
    fJetTreeData[kTOneNSubjettinessMeasured] = 0.;
    fJetTreeData[kTTwoNSubjettinessMeasured] = 0.;
  }

  if(mcSubjettiness) {
    fJetTreeData[kTOneNSubjettinessTrue] = mcSubjettiness->fOneSubjettiness;
    fJetTreeData[kTTwoNSubjettinessTrue] = mcSubjettiness->fTwoSubjettiness;
  } else {
    fJetTreeData[kTOneNSubjettinessTrue] = 0.;
    fJetTreeData[kTTwoNSubjettinessTrue] = 0.;
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
    mcjets->SetJetPtCut(20.);
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
    datajets->SetJetPtCut(20.);

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
  mgr->ConnectOutput(treemaker, 1, mgr->CreateContainer("JetSubstructureHistos_" + TString::Format("R%0d_", int(jetradius * 10.)) + trigger, AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, outputfile));
  mgr->ConnectOutput(treemaker, 2, mgr->CreateContainer("JetSubstuctureTree_" + TString::Format("R%0d_", int(jetradius * 10.)) + trigger, TTree::Class(), AliAnalysisManager::kOutputContainer, Form("JetSubstructureTree_%s.root", trigger)));

  return treemaker;
}

} /* namespace EmcalTriggerJets */
