/**************************************************************************
 * Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved. *
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

#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TClonesArray.h>
#include <TGrid.h>

#include <AliLog.h>
#include <AliVEventHandler.h>
#include <AliAnalysisManager.h>

#include "AliRhoParameter.h"
#include "AliEmcalJet.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliJetContainer.h"

#include "AliAnalysisTaskJetUE.h"

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskJetUE);
/// \endcond

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskJetUE::AliAnalysisTaskJetUE() :
  AliAnalysisTaskEmcalJetLight(),
  fBackToBackJetPtFraction(0.6),
  fMaxMomentumThridJet(12),
  fPtBinWidth(0.5),
  fMaxPt(250),
  fNtracks(0),
  fNclusters(0),
  fNjets(),
  fTotJetArea(),
  fLeadingParticle(0),
  fLeadingCluster(0),
  fLeadingJet(),
  fSortedJets()
{
}

/**
 * Standard constructor. Should be used by the user.
 *
 * @param[in] name  Name of the task
 * @param[in] histo If kTRUE, the task will also produce QA histograms
 */
AliAnalysisTaskJetUE::AliAnalysisTaskJetUE(const char *name, Bool_t histo) :
  AliAnalysisTaskEmcalJetLight(name, histo),
  fBackToBackJetPtFraction(0.6),
  fMaxMomentumThridJet(12),
  fPtBinWidth(0.5),
  fMaxPt(250),
  fNtracks(0),
  fNclusters(0),
  fNjets(),
  fTotJetArea(),
  fLeadingParticle(0),
  fLeadingCluster(0),
  fLeadingJet(),
  fSortedJets()
{
  SetMakeGeneralHistograms(histo);
}


/**
 * Calculates event properties such as number of tracks, cluster, jets, leading track, cluster, jets etc.
 */
void AliAnalysisTaskJetUE::CalculateEventProperties()
{
  fNtracks = 0;
  fNclusters = 0;
  fLeadingParticle = nullptr;
  fLeadingCluster = nullptr;

  // Loop over all possible containers
  for (auto partCont : fParticleCollArray) {
    for (auto track : partCont.second->accepted()) {
      if (!fLeadingParticle || track->Pt() > fLeadingParticle->Pt()) fLeadingParticle = track;
      fNtracks++;
    }
  }

  // Loop over all possible containers
  for (auto clusCont : fClusterCollArray) {
    for (auto clus : clusCont.second->accepted()) {
      if (!fLeadingCluster || clus->E() > fLeadingCluster->E()) fLeadingCluster = clus;
      fNclusters++;
    }
  }

  SortJets();
}

/**
 * Sort jets according to their transverse momentum
 */
void AliAnalysisTaskJetUE::SortJets()
{
  for (auto jetCont : fJetCollArray) {
    fLeadingJet[jetCont.first] = nullptr;
    fNjets[jetCont.first] = 0;
    fTotJetArea[jetCont.first] = 0;
    auto itSortedJets = fSortedJets.find(jetCont.first);
    if (itSortedJets == fSortedJets.end()) {
      auto res = fSortedJets.emplace(jetCont.first, std::list<AliEmcalJet*>());
      if (res.second) {
        itSortedJets = res.first;
      }
      else {
        AliError(Form("Error while trying to insert in the map fSortedJets for jet collection %s", jetCont.first.c_str()));
        continue;
      }
    }
    else {
      itSortedJets->second.clear();
    }
    for (auto jet1 : jetCont.second->accepted()) {
      if (!jet1->IsGhost()) {
        fNjets[jetCont.first]++;
        fTotJetArea[jetCont.first] += jet1->Area();
      }
      std::list<AliEmcalJet*>::iterator itJet2 = itSortedJets->second.begin();
      while (itJet2 != itSortedJets->second.end()) {
        AliEmcalJet* jet2 = *itJet2;
        if (jet1->Pt() > jet2->Pt()) break;
        itJet2++;
      }
      itSortedJets->second.insert(itJet2, jet1);
    }
    if (!itSortedJets->second.empty()) fLeadingJet[jetCont.first] = *(itSortedJets->second.begin());
  }
}

/**
 * Determines whether the current event is a "back-to-back" event, using the following definition.
 * An event is back-to-back if:
 * 1) The event contains a jet whose azimuthal angle difference with the leading jet is > 5/6 pi (the "back-to-back" jet), and
 * 2) the back-to-back jet carries a minimum fraction (fBackToBackJetPtFraction) of the leading jet momentum, and
 * 3) there are no jets with momentum greater than a given threshold (fMaxMomentumThridJet) that have an azimuthal angle difference smaller than 5/6 pi with the leading jet.
 * @param jetCollName name of the jet collection used to establish the back-to-back nature of the event
 * @return A boolean indicating whether the event is classified as back-to-back (true) or not (false).
 */
Bool_t AliAnalysisTaskJetUE::IsB2BEvent(std::string jetCollName)
{
  static Float_t minPhi = (5.0/6.0) * TMath::Pi();

  Bool_t b2bJet = kFALSE;
  Bool_t thirdJetOverThreshold = kFALSE;

  auto itJet = fSortedJets[jetCollName].begin();
  if (itJet == fSortedJets[jetCollName].end()) return kFALSE;
  AliEmcalJet* leadingJet = *itJet;
  Double_t minB2Bpt = leadingJet->Pt() * fBackToBackJetPtFraction;
  itJet++;
  while (itJet != fSortedJets[jetCollName].end()) {
    auto jet = *itJet;
    Double_t phidiff = TMath::Abs(AliEmcalContainer::RelativePhi(jet->Phi(), leadingJet->Phi()));
    if (phidiff > minPhi) {
      if (jet->Pt() > minB2Bpt) b2bJet = kTRUE;
    }
    else if (jet->Pt() > fMaxMomentumThridJet) {
      thirdJetOverThreshold = kTRUE;
      break;
    }
    if (jet->Pt() < fMaxMomentumThridJet && (b2bJet || jet->Pt() < minB2Bpt)) break;
    itJet++;
  }
  return b2bJet && !thirdJetOverThreshold;
}

/**
 * Verify whether two jets have any track in common.
 * @param jet1 First jet
 * @param jet2 Second jet
 * @return kTRUE if the two jets have at least a track in common, kFALSE otherwise
 */
Bool_t AliAnalysisTaskJetUE::AreJetsOverlapping(AliEmcalJet* jet1, AliEmcalJet* jet2)
{
  for (Int_t i = 0; i < jet1->GetNumberOfTracks(); ++i) {
    Int_t jet1Track = jet1->TrackAt(i);
    for (Int_t j = 0; j < jet2->GetNumberOfTracks(); ++j) {
      Int_t jet2Track = jet2->TrackAt(j);
      if (jet1Track == jet2Track) return kTRUE;
    }
  }
  return kFALSE;
}
