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
#include <array>
#include <algorithm>
#include <iostream>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/contrib/SoftDrop.hh>

#include <TLinearBinning.h>
#include <TCustomBinning.h>

#include "AliAnalysisEmcalSoftdropHelper.h"
#include "AliEmcalJet.h"
#include "AliLog.h"
#include "AliVParticle.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisEmcalSoftdropHelperImpl)
ClassImp(PWGJE::EMCALJetTasks::AliAnalysisEmcalSoftdropHelper)

using namespace PWGJE::EMCALJetTasks;

TBinning *AliAnalysisEmcalSoftdropHelperImpl::GetDefaultPartLevelPtBinning(EBinningMode_t binmode) const
{
  auto binning = new TCustomBinning;
  binning->SetMinimum(0);
  switch (binmode)
  {
  case kSDModeINT7:
  {
    binning->AddStep(20., 20.);
    binning->AddStep(40., 10.);
    binning->AddStep(80., 20.);
    binning->AddStep(120., 40.);
    binning->AddStep(240., 120.);
    break;
  }
  case kSDModeEJ1:
  {
    binning->AddStep(80., 80.);
    binning->AddStep(140., 10.);
    binning->AddStep(200., 20.);
    binning->AddStep(240., 40.);
    binning->AddStep(400., 160.);
    break;
  }
  case kSDModeEJ2:
  {
    binning->AddStep(70., 70.);
    binning->AddStep(100., 10.);
    binning->AddStep(140., 20.);
    binning->AddStep(400., 260.);
    break;
  }
  };
  return binning;
}

TBinning *AliAnalysisEmcalSoftdropHelperImpl::GetDefaultDetLevelPtBinning(EBinningMode_t binmode) const
{
  auto binning = new TCustomBinning;
  switch (binmode)
  {
  case kSDModeINT7:
  {
    binning->SetMinimum(20);
    binning->AddStep(40., 5.);
    binning->AddStep(60., 10.);
    binning->AddStep(80., 20.);
    binning->AddStep(120., 40.);
    break;
  }
  case kSDModeEJ1:
  {
    binning->SetMinimum(80.);
    binning->AddStep(120., 5.);
    binning->AddStep(160., 10.);
    binning->AddStep(200., 20.);
    break;
  }
  case kSDModeEJ2:
  {
    binning->SetMinimum(70.);
    binning->AddStep(100., 5.);
    binning->AddStep(120., 10.);
    binning->AddStep(140., 20.);
    break;
  }
  };
  return binning;
}

TBinning *AliAnalysisEmcalSoftdropHelperImpl::GetZgBinning(double zcut) const
{
  auto binning = new TCustomBinning;
  binning->SetMinimum(0.);
  binning->AddStep(zcut, zcut);
  binning->AddStep(0.5, 0.05);
  return binning;
}

TBinning *AliAnalysisEmcalSoftdropHelperImpl::GetRgBinning(double R) const {
  auto binning = new TCustomBinning;
  binning->SetMinimum(-0.05);    // Negative bins are for untagged jets
  binning->AddStep(R + 0.05, 0.05); // Adding outer margin of 0.05 above jet R
  return binning;
}

AliAnalysisEmcalSoftdropHelperImpl::SoftdropResults AliAnalysisEmcalSoftdropHelperImpl::MakeSoftdrop(const AliEmcalJet &jet, double jetradius, bool isPartLevel, SoftdropParams sdparams, AliVCluster::VCluUserDefEnergy_t energydef, double *vertex, bool dropMass0Jets) {
  const int kClusterOffset = 30000; // In order to handle tracks and clusters in the same index space the cluster index needs and offset, large enough so that there is no overlap with track indices
  std::vector<fastjet::PseudoJet> constituents;
  fastjet::PseudoJet inputjet(jet.Px(), jet.Py(), jet.Pz(), jet.E());
  AliDebugGeneralStream("MakeSoftdrop", 2) <<  "Make new jet substrucutre for " << (isPartLevel ? "part" : "det") << " jet: Number of tracks " << jet.GetNumberOfTracks() << ", clusters " << jet.GetNumberOfClusters() << std::endl;
  if(sdparams.fUseChargedConstituents || isPartLevel){                    // Neutral particles part of particle container in case of MC
    AliDebugGeneralStream("MakeSoftdrop",1) << "Jet substructure: Using charged constituents" << std::endl;
    for(int itrk = 0; itrk < jet.GetNumberOfTracks(); itrk++){
      auto track = jet.Track(itrk);
      if(!track->Charge() && !sdparams.fUseNeutralConstituents) continue;      // Reject neutral constituents in case of using only charged consituents
      if(track->Charge() && !sdparams.fUseChargedConstituents) continue;       // Reject charged constituents in case of using only neutral consituents
      fastjet::PseudoJet constituentTrack(track->Px(), track->Py(), track->Pz(), track->E());
      constituentTrack.set_user_index(jet.TrackAt(itrk));
      constituents.push_back(constituentTrack);
    }
  }

  if(sdparams.fUseNeutralConstituents){
    AliDebugGeneralStream("MakeSoftDrop",1) << "Jet substructure: Using neutral constituents" << std::endl;
    for(int icl = 0; icl < jet.GetNumberOfClusters(); icl++) {
      auto cluster = jet.Cluster(icl);
      TLorentzVector clustervec;
      cluster->GetMomentum(clustervec, vertex, energydef);
      fastjet::PseudoJet constituentCluster(clustervec.Px(), clustervec.Py(), clustervec.Pz(), cluster->GetHadCorrEnergy());
      constituentCluster.set_user_index(jet.ClusterAt(icl) + kClusterOffset);
      constituents.push_back(constituentCluster);
    }
  }

  AliDebugGeneralStream("MakeSoftDrop",3) << "Found " << constituents.size() << " constituents for jet with pt=" << jet.Pt() << " GeV/c" << std::endl;
  if(!constituents.size()){
    if(sdparams.fUseChargedConstituents && sdparams.fUseNeutralConstituents) AliErrorGeneralStream("MakeSoft") << "Jet has 0 constituents." << std::endl;
    throw 1;
  }
  // Redo jet finding on constituents with a
  fastjet::JetDefinition jetdef(fastjet::antikt_algorithm, jetradius*2, fastjet::E_scheme, fastjet::BestFJ30 );
  fastjet::ClusterSequence jetfinder(constituents, jetdef);
  std::vector<fastjet::PseudoJet> outputjets = jetfinder.inclusive_jets(0);
  auto sdjet = outputjets[0];
  if(dropMass0Jets && (TMath::Abs(sdjet.m2()) < 1e-5)) {
    AliErrorGeneralStream("MakeSoft") << "Detected jet with mass 0, not possible to evaluate SoftDrop" << std::endl;
    throw  2;
  }
  fastjet::contrib::SoftDrop softdropAlgorithm(sdparams.fBeta, sdparams.fZcut, jetradius);
  softdropAlgorithm.set_verbose_structure(kTRUE);
  fastjet::JetAlgorithm reclusterizingAlgorithm;
  switch(sdparams.fReclusterizer) {
    case kCAAlgo: reclusterizingAlgorithm = fastjet::cambridge_aachen_algorithm; break;
    case kKTAlgo: reclusterizingAlgorithm = fastjet::kt_algorithm; break;
    case kAKTAlgo: reclusterizingAlgorithm = fastjet::antikt_algorithm; break;
    default:
      AliErrorGeneralStream("MakeSoftdrop") << "Non-supported reclusterizer type" << std::endl;
      throw 3;
  };
#if FASTJET_VERSION_NUMBER >= 30302
  fastjet::Recluster reclusterizer(reclusterizingAlgorithm, 1, fastjet::Recluster::keep_only_hardest);
#else
  fastjet::contrib::Recluster reclusterizer(reclusterizingAlgorithm, 1, true);
#endif
  softdropAlgorithm.set_reclustering(kTRUE, &reclusterizer);
  AliDebugGeneralStream("MakeSoftdrop", 4) << "Jet has " << sdjet.constituents().size() << " constituents" << std::endl;
  auto groomed = softdropAlgorithm(sdjet);
  auto softdropstruct = groomed.structure_of<fastjet::contrib::SoftDrop>();

  // fjcontrib 1.024 uses 0. for untagged jets, fjcontrib 1.041 and newer -1. 
  // implementation is supposed to be consistent with putting 0. for untagged 
  // jets in SD, as histogram ranges for the untagged bin is usually adapted
  // to this
  return {softdropstruct.symmetry() < 0. ? 0. : softdropstruct.symmetry(),
          groomed.m(),
          softdropstruct.delta_R() < 0. ? 0. : softdropstruct.delta_R(),
          groomed.perp(),
          softdropstruct.mu() < 0. ? 0. : softdropstruct.mu(),
          softdropstruct.dropped_count()};
}

AliAnalysisEmcalSoftdropHelperImpl::SoftdropResults AliAnalysisEmcalSoftdropHelperImpl::MakeSoftdropStandAlone(const AliEmcalJet &jet, double jetradius, bool isPartLevel, SoftdropParams sdparams, AliVCluster::VCluUserDefEnergy_t energydef, double *vertex, bool dropMass0Jets) {
  const int kClusterOffset = 30000; // In order to handle tracks and clusters in the same index space the cluster index needs and offset, large enough so that there is no overlap with track indices
  std::vector<fastjet::PseudoJet> constituents;
  fastjet::PseudoJet inputjet(jet.Px(), jet.Py(), jet.Pz(), jet.E());
  AliDebugGeneralStream("MakeSoftdrop", 2) <<  "Make new jet substrucutre for " << (isPartLevel ? "part" : "det") << " jet: Number of tracks " << jet.GetNumberOfTracks() << ", clusters " << jet.GetNumberOfClusters() << std::endl;
  if(sdparams.fUseChargedConstituents || isPartLevel){                    // Neutral particles part of particle container in case of MC
    AliDebugGeneralStream("MakeSoftdrop",1) << "Jet substructure: Using charged constituents" << std::endl;
    for(int itrk = 0; itrk < jet.GetNumberOfTracks(); itrk++){
      auto track = jet.Track(itrk);
      if(!track->Charge() && !sdparams.fUseNeutralConstituents) continue;      // Reject neutral constituents in case of using only charged consituents
      if(track->Charge() && !sdparams.fUseChargedConstituents) continue;       // Reject charged constituents in case of using only neutral consituents
      fastjet::PseudoJet constituentTrack(track->Px(), track->Py(), track->Pz(), track->E());
      constituentTrack.set_user_index(jet.TrackAt(itrk));
      constituents.push_back(constituentTrack);
    }
  }

  if(sdparams.fUseNeutralConstituents){
    AliDebugGeneralStream("MakeSoftDrop",1) << "Jet substructure: Using neutral constituents" << std::endl;
    for(int icl = 0; icl < jet.GetNumberOfClusters(); icl++) {
      auto cluster = jet.Cluster(icl);
      TLorentzVector clustervec;
      cluster->GetMomentum(clustervec, vertex, energydef);
      fastjet::PseudoJet constituentCluster(clustervec.Px(), clustervec.Py(), clustervec.Pz(), cluster->GetHadCorrEnergy());
      constituentCluster.set_user_index(jet.ClusterAt(icl) + kClusterOffset);
      constituents.push_back(constituentCluster);
    }
  }

  AliDebugGeneralStream("MakeSoftDrop",3) << "Found " << constituents.size() << " constituents for jet with pt=" << jet.Pt() << " GeV/c" << std::endl;
  if(!constituents.size()){
    if(sdparams.fUseChargedConstituents && sdparams.fUseNeutralConstituents) AliErrorGeneralStream("MakeSoft") << "Jet has 0 constituents." << std::endl;
    throw 1;
  }
  // Redo jet finding on constituents with a
  fastjet::JetDefinition jetdef(fastjet::antikt_algorithm, jetradius*2, fastjet::E_scheme, fastjet::BestFJ30 );
  fastjet::ClusterSequence jetfinder(constituents, jetdef);
  std::vector<fastjet::PseudoJet> outputjets = jetfinder.inclusive_jets(0);
  auto sdjet = outputjets[0];
  if(dropMass0Jets && (TMath::Abs(sdjet.m2()) < 1e-5)) {
    AliErrorGeneralStream("MakeSoft") << "Detected jet with mass 0, not possible to evaluate SoftDrop" << std::endl;
    throw  2;
  }

  fastjet::JetAlgorithm jetalgo;
  switch (sdparams.fReclusterizer)
  {
  case EReclusterizer_t::kCAAlgo : 
    jetalgo = fastjet::cambridge_algorithm;
    break;
  case EReclusterizer_t::kAKTAlgo : 
    jetalgo = fastjet::antikt_algorithm;
    break;
  case EReclusterizer_t::kKTAlgo : 
    jetalgo = fastjet::kt_algorithm;
    break;
  default:
    AliErrorGeneralStream("IterativeDecluster") << "Non-supported reclusterizer type" << std::endl;
    throw 3;
  };

  try {
    fastjet::JetDefinition fJetDef(jetalgo, 1., static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30 ); 
    fastjet::ClusterSequence recluster(constituents, fJetDef);
    auto outputJets = recluster.inclusive_jets(0);
    fastjet::PseudoJet harder, softer, splitting = outputJets[0];
    if(dropMass0Jets && (TMath::Abs(splitting.m2()) < 1e-5)) {
      AliErrorGeneralStream("MakeSoft") << "Detected jet with mass 0, not possible to evaluate SoftDrop" << std::endl;
      throw  2;
    } 
    int drop_count = 0;
    while(splitting.has_parents(harder,softer)){
      if(harder.perp() < softer.perp()) std::swap(harder,softer);
      drop_count += 1;
      if(splitting.m2() < 1e-5 || splitting.pt() < 1e-5) {
        // treat as untagged jet
        // use 0. for symmetry parameter corresponding to implemetation in fjcontrib 1.042
        return SoftdropResults {
          0.,
          harder.m(),
          0.,
          harder.perp(),
          0.,
          drop_count
        };
      }
      auto sym = softer.perp()/(harder.perp()+softer.perp()),
           geoterm = sdparams.fBeta > 0 ? std::pow(harder.delta_R(softer) / jetradius, sdparams.fBeta) : 1.,
           zcut = sdparams.fZcut * geoterm; 
      if(sym > zcut) {
        // accept splitting
        double mu2 = std::max(harder.m2(), softer.m2())/splitting.m2();
        SoftdropResults acceptedSplitting {
          sym,
          harder.m(),
          harder.delta_R(softer),
          harder.perp(),
          mu2 >= 0 ? std::sqrt(mu2) : -std::sqrt(-mu2),
          drop_count 
        };
        return acceptedSplitting;
      }
      splitting = harder;
    }
    // Untagged jet
    return SoftdropResults {
      0.,
      splitting.m(),
      0.,
      splitting.perp(),
      0.,
      0
    };
  } catch(...) {
    AliErrorGeneralStream("IterativeDecluster") << "Fastjet error " << std::endl;
    throw 4;
  }
}

std::vector<AliAnalysisEmcalSoftdropHelperImpl::SoftdropResults> AliAnalysisEmcalSoftdropHelperImpl::IterativeDecluster(const AliEmcalJet &jet, double jetradius, bool isPartLevel, SoftdropParams sdparams, AliVCluster::VCluUserDefEnergy_t energydef, double *vertex, bool dropMass0Jets) {
  const int kClusterOffset = 30000; // In order to handle tracks and clusters in the same index space the cluster index needs and offset, large enough so that there is no overlap with track indices
  std::vector<fastjet::PseudoJet> constituents;
  fastjet::PseudoJet inputjet(jet.Px(), jet.Py(), jet.Pz(), jet.E());
  AliDebugGeneralStream("MakeSoftdrop", 2) <<  "Make new jet substrucutre for " << (isPartLevel ? "part" : "det") << " jet: Number of tracks " << jet.GetNumberOfTracks() << ", clusters " << jet.GetNumberOfClusters() << std::endl;
  if(sdparams.fUseChargedConstituents || isPartLevel){                    // Neutral particles part of particle container in case of MC
    AliDebugGeneralStream("MakeSoftdrop",1) << "Jet substructure: Using charged constituents" << std::endl;
    for(int itrk = 0; itrk < jet.GetNumberOfTracks(); itrk++){
      auto track = jet.Track(itrk);
      if(!track->Charge() && !sdparams.fUseNeutralConstituents) continue;      // Reject neutral constituents in case of using only charged consituents
      if(track->Charge() && !sdparams.fUseChargedConstituents) continue;       // Reject charged constituents in case of using only neutral consituents
      fastjet::PseudoJet constituentTrack(track->Px(), track->Py(), track->Pz(), track->E());
      constituentTrack.set_user_index(jet.TrackAt(itrk));
      constituents.push_back(constituentTrack);
    }
  }

  if(sdparams.fUseNeutralConstituents){
    AliDebugGeneralStream("MakeSoftDrop",1) << "Jet substructure: Using neutral constituents" << std::endl;
    for(int icl = 0; icl < jet.GetNumberOfClusters(); icl++) {
      auto cluster = jet.Cluster(icl);
      TLorentzVector clustervec;
      cluster->GetMomentum(clustervec, vertex, energydef);
      fastjet::PseudoJet constituentCluster(clustervec.Px(), clustervec.Py(), clustervec.Pz(), cluster->GetHadCorrEnergy());
      constituentCluster.set_user_index(jet.ClusterAt(icl) + kClusterOffset);
      constituents.push_back(constituentCluster);
    }
  }

  AliDebugGeneralStream("IterativeDecluster",3) << "Found " << constituents.size() << " constituents for jet with pt=" << jet.Pt() << " GeV/c" << std::endl;
  if(!constituents.size()){
    if(sdparams.fUseChargedConstituents && sdparams.fUseNeutralConstituents) AliErrorGeneralStream("MakeSoft") << "Jet has 0 constituents." << std::endl;
    throw 1;
  } 

  std::vector<SoftdropResults> result;

  fastjet::JetAlgorithm jetalgo;
  switch (sdparams.fReclusterizer)
  {
  case EReclusterizer_t::kCAAlgo : 
    jetalgo = fastjet::cambridge_algorithm;
    break;
  case EReclusterizer_t::kAKTAlgo : 
    jetalgo = fastjet::antikt_algorithm;
    break;
  case EReclusterizer_t::kKTAlgo : 
    jetalgo = fastjet::kt_algorithm;
    break;
  default:
    AliErrorGeneralStream("IterativeDecluster") << "Non-supported reclusterizer type" << std::endl;
    throw 3;
  };
  

  try {
    fastjet::JetDefinition fJetDef(jetalgo, 1., static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30 ); 
    fastjet::ClusterSequence recluster(constituents, fJetDef);
    auto outputJets = recluster.inclusive_jets(0);
    fastjet::PseudoJet harder, softer, splitting = outputJets[0];
    if(dropMass0Jets && (TMath::Abs(splitting.m2()) < 1e-5)) {
      AliErrorGeneralStream("MakeSoft") << "Detected jet with mass 0, not possible to evaluate SoftDrop" << std::endl;
      throw  2;
    } 
    int drop_count = 0;
    while(splitting.has_parents(harder,softer)){
      if(harder.perp() < softer.perp()) std::swap(harder,softer);
      drop_count += 1;
      auto sym = softer.perp()/(harder.perp()+softer.perp()),
           geoterm = sdparams.fBeta > 0 ? std::pow(harder.delta_R(softer) / jetradius, sdparams.fBeta) : 1.,
           zcut = sdparams.fZcut * geoterm; 
      if(sym > zcut) {
        // accept splitting
        double mu2 = TMath::Abs(splitting.m2()) < 1e-5 ? 100000. : std::max(harder.m2(), softer.m2())/splitting.m2();
        SoftdropResults acceptedSplitting{
          sym,
          harder.m(),
          harder.delta_R(softer),
          harder.perp(),
          mu2 >= 0 ? std::sqrt(mu2) : -std::sqrt(-mu2),
          drop_count 
        };
        result.push_back(acceptedSplitting);
      }
      splitting = harder;
    }
  } catch(...) {
    AliErrorGeneralStream("IterativeDecluster") << "Fastjet error " << std::endl;
    throw 4;
  }

  return result;
}
