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
#include <algorithm>
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/PseudoJet.hh"

#include "AliClusterContainer.h"
#include "AliEmcalJet.h"
#include "AliParticleContainer.h"
#include "AliVParticle.h"

#include "AliLundPlaneHelper.h"

ClassImp(PWGJE::EMCALJetTasks::AliLundPlaneHelper)

using namespace PWGJE::EMCALJetTasks;

AliLundPlaneHelper::AliLundPlaneHelper(): 
  TObject(),
  fHardCutoff(0)
{

}

AliLundPlaneData AliLundPlaneHelper::Evaluate(const AliEmcalJet &jet, const AliParticleContainer *tracks, const AliClusterContainer *clusters, Double_t *vertexpos) {
  std::vector<fastjet::PseudoJet>  fInputVectors;
  fInputVectors.clear();
  fastjet::PseudoJet  PseudoTracks;
  Int_t nall=0;
  Double_t z = 0;

  if (tracks){
    for (Int_t i=0; i < jet.GetNumberOfTracks(); i++) {
      AliVParticle *fTrk = jet.TrackAt(i, tracks->GetArray());
      if (!fTrk) continue;
      PseudoTracks.reset(fTrk->Px(), fTrk->Py(), fTrk->Pz(),fTrk->E());
      PseudoTracks.set_user_index(jet.TrackAt(i)+100);
      fInputVectors.push_back(PseudoTracks);
    }
  }

  if(clusters){
    if(!vertexpos) throw AliLundPlaneException(AliLundPlaneException::kVertexNotSet, "vertex position not defined");
    for(int icl = 0; icl < jet.GetNumberOfClusters(); icl++) {
      auto cluster = jet.ClusterAt(icl, clusters->GetArray());
      TLorentzVector clustervec;
      cluster->GetMomentum(clustervec, vertexpos, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy());
      fastjet::PseudoJet constituentCluster(clustervec.Px(), clustervec.Py(), clustervec.Pz(), cluster->GetHadCorrEnergy());
      constituentCluster.set_user_index(jet.ClusterAt(icl) + 1000);
      fInputVectors.push_back(constituentCluster);
    }
  }

  fastjet::JetAlgorithm jetalgo(fastjet::cambridge_algorithm);
  fastjet::GhostedAreaSpec ghost_spec(1, 1, 0.05);
   
  fastjet::JetDefinition fJetDef(jetalgo, 1., static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30 ); 
  fastjet::AreaDefinition fAreaDef(fastjet::passive_area,ghost_spec); 

  AliLundPlaneData result;
  try {
    fastjet::ClusterSequenceArea fClustSeqSA(fInputVectors, fJetDef, fAreaDef);
    std::vector<fastjet::PseudoJet> fOutputJets;
    fOutputJets.clear();
    fOutputJets=fClustSeqSA.inclusive_jets(0);
  
    fastjet::PseudoJet j1, j2, jj = fOutputJets[0];
    while(jj.has_parents(j1,j2)){
      nall++;
      if(j1.perp() < j2.perp()) std::swap(j1,j2);
      z=j2.perp()/(j1.perp()+j2.perp());
      double delta_R=j1.delta_R(j2);
      double lndeltaR =log(1.0/delta_R);
      double lnpt_rel=log(j2.perp()*delta_R);
      AliLundPlaneParameters currentsplitting(lndeltaR, lnpt_rel, j2.perp(), nall);
      result.InsertSplitting(currentsplitting);
      jj=j1;
    }
  } catch (fastjet::Error &e) {
    throw AliLundPlaneException(AliLundPlaneException::kFastjetError, e.message().data());
  }

  return result;
}