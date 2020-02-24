/************************************************************************************
 * Copyright (C) 2019, Copyright Holders of the ALICE Collaboration                 *
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
#include <memory>
#include <sstream>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/contrib/SoftDrop.hh>

#include <TArray.h>
#include <TCustomBinning.h>
#include <THistManager.h>
#include <TLinearBinning.h>

#include "AliAnalysisTaskEmcalSoftDropData.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliClusterContainer.h"
#include "AliEmcalDownscaleFactorsOCDB.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliEmcalJet.h"
#include "AliInputEventHandler.h"
#include "AliJetContainer.h"
#include "AliLog.h"
#include "AliTrackContainer.h"
#include "AliVCluster.h"
#include "AliVTrack.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalSoftDropData)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalSoftDropData::AliAnalysisTaskEmcalSoftDropData() : 
  AliAnalysisTaskEmcalJet(),
  fTriggerBits(AliVEvent::kAny),
  fTriggerString(""),
  fUseDownscaleWeight(false),
  fBeta(0.),
  fZcut(0.1),
  fReclusterizer(kCAAlgo),
  fUseChargedConstituents(kTRUE),
  fUseNeutralConstituents(kTRUE),
  fHistos(nullptr),
  fPtBinning(nullptr)
{

}

AliAnalysisTaskEmcalSoftDropData::AliAnalysisTaskEmcalSoftDropData(EMCAL_STRINGVIEW name) : 
  AliAnalysisTaskEmcalJet(name.data(), kTRUE),
  fTriggerBits(AliVEvent::kAny),
  fTriggerString(""),
  fUseDownscaleWeight(false),
  fBeta(0.),
  fZcut(0.1),
  fReclusterizer(kCAAlgo),
  fUseChargedConstituents(kTRUE),
  fUseNeutralConstituents(kTRUE),
  fHistos(nullptr),
  fPtBinning(nullptr)
{

}

AliAnalysisTaskEmcalSoftDropData::~AliAnalysisTaskEmcalSoftDropData() {
  if(fPtBinning) delete fPtBinning;
  if(fHistos) delete fHistos;
}

void AliAnalysisTaskEmcalSoftDropData::UserCreateOutputObjects() {
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  double R = GetJetContainer("datajets")->GetJetRadius();
  if(!fPtBinning) fPtBinning = new TLinearBinning(300, 0., 300.);  // Use fine binning for data, rebin offline
  std::unique_ptr<TBinning> zgBinning(GetZgBinning()),
                            rgBinning(GetRgBinning(R)),
                            nsdBinning(new TLinearBinning(22, -1.5, 20.5)),     // Negative bins for untagged jets
                            thetagBinning(new TLinearBinning(11, -0.1, 1.)); 
  TArrayD edgesPt;
  fPtBinning->CreateBinEdges(edgesPt);

  fHistos = new THistManager("histosSoftdrop");
  fHistos->CreateTH1("hEventCounter", "EventCounter", 1, 0.5, 1.5);
  fHistos->CreateTH1("hJetPtRaw", "raw jet pt", 300, 0., 300.);
  fHistos->CreateTH2("hZgVsPt", "zg vs pt", *zgBinning, *fPtBinning, "s");
  fHistos->CreateTH2("hRgVsPt", "rg vs pt", *rgBinning,  *fPtBinning, "s");
  fHistos->CreateTH2("hNsdVsPt", "nsd vs pt", *nsdBinning, *fPtBinning, "s");
  fHistos->CreateTH2("hThetagVsPt", "thetag vs pt", *thetagBinning,  *fPtBinning, "s");
  fHistos->CreateTH1("hSkippedJets", "Number of skipped jets", *fPtBinning);
  if(fUseDownscaleWeight){
    fHistos->CreateTH2("hZgVsPtWeighted", "zg vs pt (weighted)", *zgBinning, *fPtBinning, "s");
    fHistos->CreateTH2("hRgVsPtWeighted", "rg vs pt (weighted)", *rgBinning,  *fPtBinning, "s");
    fHistos->CreateTH2("hNsdVsPtWeighted", "nsd vs pt (weighted)", *nsdBinning, *fPtBinning, "s");
    fHistos->CreateTH2("hThetagVsPtWeighted", "thetag vs pt (weighted)", *thetagBinning,  *fPtBinning, "s");
    fHistos->CreateTH1("hEventCounterWeighted", "Event counter, weighted", 1., 0.5, 1.5);
    fHistos->CreateTH1("hJetPtRawWeighted", "raw jet pt", 300, 0., 300., "s");
    fHistos->CreateTH1("hSkippedJetsWeighted", "Number of skipped jets (weighted)", *fPtBinning);
  }

  // A bit of QA stuff
  fHistos->CreateTH2("hQANEFPt", "Neutral energy fraction; p_{t} (GeV/c); NEF", 350, 0., 350., 100, 0., 1.);
  fHistos->CreateTH2("hQAZchPt", "z_{ch,max}; p_{t} (GeV/c); z_{ch,max}", 350, 0., 350., 100, 0., 1.);
  fHistos->CreateTH2("hQAZnePt", "z_{ne,max}; p_{t} (GeV/c); z_{ne,max}", 350, 0., 350., 100, 0., 1.);
  fHistos->CreateTH2("hQANChPt", "Number of charged constituents; p_{t} (GeV/c); N_{ch}", 350, 0., 350., 100, 0., 100.);
  fHistos->CreateTH2("hQANnePt", "Number of neutral constituents; p_{t} (GeV/c); N_{ne}", 350, 0., 350., 100, 0., 100.);
  fHistos->CreateTH2("hSDUsedChargedPtjvPtc", "p_{t,j} vs. p_{t,const} for tracks used in SD", 350., 0., 350., 350, 0., 350.);
  fHistos->CreateTH2("hSDUsedNeutralPtjvPtc", "p_{t,j} vs. p_{t,const} for clusters used in SD", 350., 0., 350., 350, 0., 350.);

  for(auto h : *fHistos->GetListOfHistograms()) fOutput->Add(h);
  PostData(1, fOutput);
}

Bool_t AliAnalysisTaskEmcalSoftDropData::IsTriggerSelected(){
  if(!(fInputHandler->IsEventSelected() & fTriggerBits)) return false;
  if(fTriggerString.length()) {
    if(!fInputEvent->GetFiredTriggerClasses().Contains(fTriggerString)) return false;
  }
  return true;
}

void AliAnalysisTaskEmcalSoftDropData::RunChanged(Int_t newrun){
  if(fUseDownscaleWeight) 
    PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance()->SetRun(newrun);
}

Bool_t AliAnalysisTaskEmcalSoftDropData::Run() {
  auto jets = GetJetContainer("datajets");
  if(!jets) {
    AliErrorStream() << "Jet container not found" << std::endl;
    return false;
  } 
  auto clusters = GetClusterContainer(0);
  if(fUseNeutralConstituents && !clusters) {
    AliErrorStream() << "Cluster container not found, but neutral constituents requested" << std::endl; 
  }
  auto tracks  = GetTrackContainer(0);
  if(fUseChargedConstituents &&  !tracks) {
    AliErrorStream() << "Track container not found, but charged constituent requested." << std::endl;
    return false;
  }

  Double_t weight = fUseDownscaleWeight ? 1./GetDownscaleWeight() : 1.;
  Double_t Rjet = jets->GetJetRadius();
  fHistos->FillTH1("hEventCounter", 1., weight);
  if(fUseDownscaleWeight) fHistos->FillTH1("hEventCounterWeighted", 1., weight);

  for(auto jet : jets->accepted()){
    AliDebugStream(2) << "Next accepted jet with pt " << jet->Pt() << std::endl;
    fHistos->FillTH1("hJetPtRaw", jet->Pt());
    if(fUseDownscaleWeight) fHistos->FillTH1("hJetPtRawWeighted", jet->Pt(), weight);
    try {
      auto zgparams = MakeSoftdrop(*jet, jets->GetJetRadius(), tracks, clusters);
      bool untagged = zgparams[0] < fZcut;
      AliDebugStream(2) << "Found jet with pt " << jet->Pt() << " and zg " << zgparams[0] << std::endl;
      fHistos->FillTH2("hZgVsPt", zgparams[0], jet->Pt());
      fHistos->FillTH2("hRgVsPt", untagged ? -0.01 : zgparams[2], jet->Pt());
      fHistos->FillTH2("hNsdVsPt", untagged ? -1. : zgparams[5], jet->Pt());
      fHistos->FillTH2("hThetagVsPt", untagged ? -0.05 : zgparams[2]/Rjet, jet->Pt());
      if(fUseDownscaleWeight) {
        fHistos->FillTH2("hZgVsPtWeighted", zgparams[0], jet->Pt(), weight);
        fHistos->FillTH2("hRgVsPtWeighted", untagged ? -0.01 : zgparams[2], jet->Pt(), weight);
        fHistos->FillTH2("hNsdVsPtWeighted", untagged ? -1. : zgparams[5], jet->Pt(), weight);
        fHistos->FillTH2("hThetagVsPtWeighted", untagged ? -0.05 : zgparams[2]/Rjet, jet->Pt(), weight);
      } 
    } catch (int e) {
      if(fUseChargedConstituents && fUseNeutralConstituents) AliErrorStream() << "Softdrop error " << e << ": Having 0 constituents for reclustering" << std::endl;
      fHistos->FillTH1("hSkippedJets", jet->Pt());
      if(fUseDownscaleWeight) 
        fHistos->FillTH1("hSkippedJetsWeighted", jet->Pt(), weight);
    }

    // Fill QA plots - trigger cluster independent
    // Those plots have been in before (as part of the Tree) but were 
    // removed in order to reduce the disk space consumption.
    fHistos->FillTH2("hQANEFPt", jet->Pt(), jet->NEF(), weight);
    if(clusters){
      auto leadcluster = jet->GetLeadingCluster(clusters->GetArray());
      if(leadcluster){
        TLorentzVector ptvec;
        leadcluster->GetMomentum(ptvec, fVertex, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy());
        fHistos->FillTH2("hQAZnePt", jet->Pt(), jet->GetZ(ptvec.Px(), ptvec.Py(), ptvec.Pz()), weight);
      }
    }
    if(tracks){
      auto leadingtrack = jet->GetLeadingTrack(tracks->GetArray());
      if(leadingtrack) fHistos->FillTH2("hQAZchPt", jet->Pt(), jet->GetZ(leadingtrack->Px(), leadingtrack->Py(), leadingtrack->Pz()), weight);
    }
    fHistos->FillTH2("hQANChPt", jet->Pt(), jet->GetNumberOfTracks(), weight);
    fHistos->FillTH2("hQANnePt", jet->Pt(), jet->GetNumberOfClusters(), weight);
  }
  return true;
}

Double_t AliAnalysisTaskEmcalSoftDropData::GetDownscaleWeight() const {
  Double_t weight = 1.;
  TString triggerclass;
  if(fTriggerString == "INT7") triggerclass = "CINT7-B-NOPF-CENT";
  else if(fTriggerString == "EJ1") triggerclass = "CEMC7EJ1-B-NOPF-CENTNOTRD";
  else if(fTriggerString == "EJ2") triggerclass = "CEMC7EJ2-B-NOPF-CENT";
  if(triggerclass.Length()) weight = PWG::EMCAL::AliEmcalDownscaleFactorsOCDB::Instance()->GetDownscaleFactorForTriggerClass(triggerclass);
  return weight;
}

TBinning *AliAnalysisTaskEmcalSoftDropData::GetZgBinning() const {
  auto binning = new TCustomBinning;
  binning->SetMinimum(0.);
  binning->AddStep(fZcut, fZcut);
  binning->AddStep(0.5, 0.05);
  return binning;
}

TBinning *AliAnalysisTaskEmcalSoftDropData::GetRgBinning(double R) const {
  auto binning = new TCustomBinning;
  binning->SetMinimum(-0.05); // Negative bin for untagged jets
  binning->AddStep(R, 0.05);
  return binning;
}

std::vector<double> AliAnalysisTaskEmcalSoftDropData::MakeSoftdrop(const AliEmcalJet &jet, double jetradius, const AliParticleContainer *tracks, const AliClusterContainer *clusters) {
  const int kClusterOffset = 30000; // In order to handle tracks and clusters in the same index space the cluster index needs and offset, large enough so that there is no overlap with track indices
  std::vector<fastjet::PseudoJet> constituents;
  bool isMC = dynamic_cast<const AliMCParticleContainer *>(tracks);
  AliDebugStream(2) << "Make new jet substrucutre for " << (isMC ? "MC" : "data") << " jet: Number of tracks " << jet.GetNumberOfTracks() << ", clusters " << jet.GetNumberOfClusters() << std::endl;
  if(tracks && (fUseChargedConstituents || isMC)){                    // Neutral particles part of particle container in case of MC
    AliDebugStream(1) << "Jet substructure: Using charged constituents" << std::endl;
    for(int itrk = 0; itrk < jet.GetNumberOfTracks(); itrk++){
      auto track = jet.TrackAt(itrk, tracks->GetArray());
      if(!track->Charge() && !fUseNeutralConstituents) continue;      // Reject neutral constituents in case of using only charged consituents
      if(track->Charge() && !fUseChargedConstituents) continue;       // Reject charged constituents in case of using only neutral consituents
      fastjet::PseudoJet constituentTrack(track->Px(), track->Py(), track->Pz(), track->E());
      constituentTrack.set_user_index(jet.TrackAt(itrk));
      constituents.push_back(constituentTrack);
      fHistos->FillTH2("hSDUsedChargedPtjvPtc", jet.Pt(), constituentTrack.pt());
    }
  }

  if(clusters && fUseNeutralConstituents){
    AliDebugStream(1) << "Jet substructure: Using neutral constituents" << std::endl;
    for(int icl = 0; icl < jet.GetNumberOfClusters(); icl++) {
      auto cluster = jet.ClusterAt(icl, clusters->GetArray());
      TLorentzVector clustervec;
      cluster->GetMomentum(clustervec, fVertex, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy());
      fastjet::PseudoJet constituentCluster(clustervec.Px(), clustervec.Py(), clustervec.Pz(), cluster->GetHadCorrEnergy());
      constituentCluster.set_user_index(jet.ClusterAt(icl) + kClusterOffset);
      constituents.push_back(constituentCluster);
      fHistos->FillTH2("hSDUsedNeutralPtjvPtc", jet.Pt(), constituentCluster.pt());
    }
  }

  AliDebugStream(3) << "Found " << constituents.size() << " constituents for jet with pt=" << jet.Pt() << " GeV/c" << std::endl;
  if(!constituents.size()){
    if(fUseChargedConstituents && fUseNeutralConstituents) AliErrorStream() << "Jet has 0 constituents." << std::endl;
    throw 1;
  }
  // Redo jet finding on constituents with a
  fastjet::JetDefinition jetdef(fastjet::antikt_algorithm, jetradius*2, fastjet::E_scheme, fastjet::BestFJ30 );
  fastjet::ClusterSequence jetfinder(constituents, jetdef);
  std::vector<fastjet::PseudoJet> outputjets = jetfinder.inclusive_jets(0);
  auto sdjet = outputjets[0];
  fastjet::contrib::SoftDrop softdropAlgorithm(fBeta, fZcut);
  softdropAlgorithm.set_verbose_structure(kTRUE);
  fastjet::JetAlgorithm reclusterizingAlgorithm;
  switch(fReclusterizer) {
    case kCAAlgo: reclusterizingAlgorithm = fastjet::cambridge_aachen_algorithm; break;
    case kKTAlgo: reclusterizingAlgorithm = fastjet::kt_algorithm; break;
    case kAKTAlgo: reclusterizingAlgorithm = fastjet::antikt_algorithm; break;
  };
#if FASTJET_VERSION_NUMBER >= 30302
  fastjet::Recluster reclusterizer(reclusterizingAlgorithm, 1, fastjet::Recluster::keep_only_hardest);
#else
  fastjet::contrib::Recluster reclusterizer(reclusterizingAlgorithm, 1, true);
#endif
  softdropAlgorithm.set_reclustering(kTRUE, &reclusterizer);
  AliDebugStream(4) << "Jet has " << sdjet.constituents().size() << " constituents" << std::endl;
  auto groomed = softdropAlgorithm(sdjet);
  auto softdropstruct = groomed.structure_of<fastjet::contrib::SoftDrop>();

  std::vector<double> result ={softdropstruct.symmetry(),
                                  groomed.m(),
                                  softdropstruct.delta_R(),
                                  groomed.perp(),
                                  softdropstruct.mu(),
                                  static_cast<double>(softdropstruct.dropped_count())};
  return result;
}

AliAnalysisTaskEmcalSoftDropData *AliAnalysisTaskEmcalSoftDropData::AddTaskEmcalSoftDropData(Double_t jetradius, AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recombinationScheme, EMCAL_STRINGVIEW trigger) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  Bool_t isAOD(kFALSE);
  AliInputEventHandler *inputhandler = static_cast<AliInputEventHandler *>(mgr->GetInputEventHandler());
  if(inputhandler) {
    if(inputhandler->IsA() == AliAODInputHandler::Class()){
      std::cout << "Analysing AOD events\n";
      isAOD = kTRUE;
    } else {
      std::cout << "Analysing ESD events\n";
    }
  }

  std::stringstream taskname;
  taskname << "SoftdropDataMaker_R" << std::setw(2) << std::setfill('0') << int(jetradius*10) << trigger;  
  AliAnalysisTaskEmcalSoftDropData *datamaker = new AliAnalysisTaskEmcalSoftDropData(taskname.str().data());
  datamaker->SelectCollisionCandidates(AliVEvent::kAny);
  mgr->AddTask(datamaker);

  AliTrackContainer *tracks(nullptr);
  if((jettype == AliJetContainer::kChargedJet) || (jettype == AliJetContainer::kFullJet)){
      tracks = datamaker->AddTrackContainer(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::TrackContainerNameFactory(isAOD));
      std::cout << "Track container name: " << tracks->GetName() << std::endl;
      tracks->SetMinPt(0.15);
  }
  AliClusterContainer *clusters(nullptr);
  if((jettype == AliJetContainer::kFullJet) || (jettype == AliJetContainer::kNeutralJet)){
    std::cout << "Using full or neutral jets ..." << std::endl;
    clusters = datamaker->AddClusterContainer(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::ClusterContainerNameFactory(isAOD));
    std::cout << "Cluster container name: " << clusters->GetName() << std::endl;
    clusters->SetClusHadCorrEnergyCut(0.3); // 300 MeV E-cut
    clusters->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  } else {
    std::cout << "Using charged jets ... " << std::endl;
  }

  AliJetContainer *datajets = datamaker->AddJetContainer(
                              jettype,
                              AliJetContainer::antikt_algorithm,
                              recombinationScheme,
                              jetradius,
                              ((jettype == AliJetContainer::kFullJet) || (jettype == AliJetContainer::kNeutralJet)) ? AliEmcalJet::kEMCALfid : AliEmcalJet::kTPCfid,
                              tracks, clusters);
  datajets->SetName("datajets");
  datajets->SetJetPtCut(0.);
  datajets->SetMaxTrackPt(1000.);

  std::string jettypestring;
  switch(jettype) {
    case AliJetContainer::kFullJet: jettypestring = "FullJets"; break;
    case AliJetContainer::kChargedJet: jettypestring = "ChargedJets"; break;
    case AliJetContainer::kNeutralJet: jettypestring = "NeutralJets"; break;
    default: jettypestring = "Undef";
  };

  ULong_t triggerbits(AliVEvent::kINT7);
  std::string triggerstring(trigger);
  std::cout << "Found trigger " << triggerstring << std::endl;
  if(triggerstring == "EJ1") {
    std::cout << "Setting binning mode for EJ1" << std::endl;
    triggerbits = AliVEvent::kEMCEJE;
  } else if(triggerstring == "EJ2") {
    std::cout << "Setting binning mode for EJ2" << std::endl;
    triggerbits = AliVEvent::kEMCEJE;
  }
  datamaker->SetSelectTrigger(triggerbits, trigger.data());

  // Connecting containers
  std::stringstream outputfile, histname;
  outputfile << mgr->GetCommonFileName() << ":SoftDropResponse_" << jettypestring << "_R" << std::setw(2) << std::setfill('0') << int(jetradius * 10.) << "_" << trigger;
  histname << "SoftDropResponseHistos_" << jettypestring << "_R" << std::setw(2) << std::setfill('0') << int(jetradius * 10.) << "_" << trigger;
  mgr->ConnectInput(datamaker, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(datamaker, 1, mgr->CreateContainer(histname.str().data(), AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, outputfile.str().data()));

  return datamaker;
}