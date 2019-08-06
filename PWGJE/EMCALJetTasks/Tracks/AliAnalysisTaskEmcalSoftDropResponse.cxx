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
#include <iostream>
#include <vector>

#include <TArrayD.h>
#include <TBinning.h>
#include <TCustomBinning.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TRandom.h>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/contrib/SoftDrop.hh>
#include <fastjet/config.h>
#if FASJET_VERSION_NUMBER >= 30302
#include <fastjet/tools/Recluster.hh>
#else
#include <fastjet/contrib/Recluster.hh>
#endif

#include <RooUnfoldResponse.h>

#include <AliAODEvent.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisManager.h>
#include "AliAnalysisTaskEmcalSoftDropResponse.h"
#include <AliClusterContainer.h>
#include <AliEmcalJet.h>
#include <AliEmcalAnalysisFactory.h>
#include <AliJetContainer.h>
#include <AliMCParticleContainer.h>
#include <AliTrackContainer.h>
#include <AliVCluster.h>
#include <AliVParticle.h>
#include <AliVTrack.h>

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalSoftDropResponse)

    using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalSoftDropResponse::AliAnalysisTaskEmcalSoftDropResponse() : AliAnalysisTaskEmcalJet(),
                                                                               fBinningMode(kSDModeINT7),
                                                                               fFractionResponseClosure(0.5),
                                                                               fZcut(0.1),
                                                                               fBeta(0.),
                                                                               fReclusterizer(kCAAlgo),
                                                                               fSampleFraction(1.),
                                                                               fMinFractionShared(0),
                                                                               fUseChargedConstituents(true),
                                                                               fUseNeutralConstituents(true),
                                                                               fNameMCParticles("mcparticles"),
                                                                               fSampleSplitter(nullptr),
                                                                               fSampleTrimmer(nullptr),
                                                                               fPartLevelPtBinning(nullptr),
                                                                               fDetLevelPtBinning(nullptr),
                                                                               fNamePartLevelJetContainer(""),
                                                                               fNameDetLevelJetContainer(""),
                                                                               fNameUnSubLevelJetContainer(""),
                                                                               fIsEmbeddedEvent(false),
                                                                               fHistManager("AliAnalysisTaskSoftDropResponse")
{
}

AliAnalysisTaskEmcalSoftDropResponse::AliAnalysisTaskEmcalSoftDropResponse(const char *name) : AliAnalysisTaskEmcalJet(name, kTRUE),
                                                                                               fBinningMode(kSDModeINT7),
                                                                                               fFractionResponseClosure(0.5),
                                                                                               fZcut(0.1),
                                                                                               fBeta(0.),
                                                                                               fReclusterizer(kCAAlgo),
                                                                                               fSampleFraction(1.),
                                                                                               fMinFractionShared(0),
                                                                                               fUseChargedConstituents(true),
                                                                                               fUseNeutralConstituents(true),
                                                                                               fNameMCParticles("mcparticles"),
                                                                                               fSampleSplitter(nullptr),
                                                                                               fSampleTrimmer(nullptr),
                                                                                               fPartLevelPtBinning(nullptr),
                                                                                               fDetLevelPtBinning(nullptr),
                                                                                               fNamePartLevelJetContainer(""),
                                                                                               fNameDetLevelJetContainer(""),
                                                                                               fNameUnSubLevelJetContainer(""),
                                                                                               fIsEmbeddedEvent(false),
                                                                                               fHistManager(name)
{
  SetMakeGeneralHistograms(true);
}

AliAnalysisTaskEmcalSoftDropResponse::~AliAnalysisTaskEmcalSoftDropResponse()
{
  if (fPartLevelPtBinning)
    delete fPartLevelPtBinning;
  if (fDetLevelPtBinning)
    delete fDetLevelPtBinning;
  if (fSampleSplitter)
    delete fSampleSplitter;
  if (fSampleTrimmer)
    delete fSampleTrimmer;
}

void AliAnalysisTaskEmcalSoftDropResponse::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  fSampleSplitter = new TRandom;
  if (fSampleFraction < 1.)
    fSampleTrimmer = new TRandom;

  if (!fPartLevelPtBinning)
    fPartLevelPtBinning = GetDefaultPartLevelPtBinning();
  if (!fDetLevelPtBinning)
    fDetLevelPtBinning = GetDefaultDetLevelPtBinning();
  auto zgbinning = GetZgBinning();
  TArrayD binEdgesZg, binEdgesPtPart, binEdgesPtDet;
  zgbinning->CreateBinEdges(binEdgesZg);
  fPartLevelPtBinning->CreateBinEdges(binEdgesPtPart);
  fDetLevelPtBinning->CreateBinEdges(binEdgesPtDet);

  //Need to do centrality bins in the histograms if it is not pp
  if (fForceBeamType != kpp)
  {
    TString histname;
    TString histtitle;
    for (Int_t cent = 0; cent < fNcentBins; cent++)
    {
      histname = TString::Format("hZgDetLevel_%d", cent);
      histtitle = TString::Format("Zg response at detector level, %d centrality bin", cent);
      fHistManager.CreateTH2(histname, histtitle, binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
      histname = TString::Format("hZgPartLevel_%d", cent);
      histtitle = TString::Format("Zg response at particle level, %d centrality bin", cent);
      fHistManager.CreateTH2(histname, histtitle, binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
      histname = TString::Format("hZgPartLevelTruncated_%d", cent);
      histtitle = TString::Format("Zg response at particle level (truncated), %d centrality bin", cent);
      fHistManager.CreateTH2(histname, histtitle, binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());

      // For closure test
      histname = TString::Format("hZgPartLevelClosureNoResp_%d", cent);
      histtitle = TString::Format("Zg response at particle level (closure test, jets not used for the response matrix), %d centrality bin", cent);
      fHistManager.CreateTH2(histname, histtitle, binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
      histname = TString::Format("hZgDetLevelClosureNoResp_%d", cent);
      histtitle = TString::Format("Zg response at detector level (closure test, jets not used for the response matrix), %d centrality bin", cent);
      fHistManager.CreateTH2(histname, histtitle, binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
      histname = TString::Format("hZgPartLevelClosureResp_%d", cent);
      histtitle = TString::Format("Zg response at particle level (closure test, jets used for the response matrix), %d centrality bin", cent);
      fHistManager.CreateTH2(histname, histtitle, binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
      histname = TString::Format("hZgDetLevelClosureResp_%d", cent);
      histtitle = TString::Format("Zg response at detector level (closure test, jets used for the response matrix), %d centrality bin", cent);
      fHistManager.CreateTH2(histname, histtitle, binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());

      histname = TString::Format("hZgResponse_%d", cent);
      histtitle = TString::Format("z_{g} response matrix, %d centrality bin", cent);
      RooUnfoldResponse *r = new RooUnfoldResponse(histname, histtitle);
      TString fZgDetLevel = TString::Format("hZgDetLevel_%d", cent);
      TString fZgPartLevel = TString::Format("hZgPartLevel_%d", cent);
      r->Setup((TH1 *)fHistManager.FindObject(fZgDetLevel), (TH1 *)fHistManager.FindObject(fZgPartLevel));
      fZgResponse.push_back(r);
      histname = TString::Format("hZgResponseClosure_%d", cent);
      histtitle = TString::Format("z_{g} response matrix for the closure test, %d centrality bin", cent);
      RooUnfoldResponse *r_closure = new RooUnfoldResponse(histname, histtitle);
      r_closure->Setup((TH1 *)fHistManager.FindObject(fZgDetLevel), (TH1 *)fHistManager.FindObject(fZgPartLevel));
      fZgResponseClosure.push_back(r_closure);
    }
  }
  else
  {
    fHistManager.CreateTH2("hZgDetLevel", "Zg response at detector level", binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
    fHistManager.CreateTH2("hZgPartLevel", "Zg response at particle level", binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
    fHistManager.CreateTH2("hZgPartLevelTruncated", "Zg response at particle level (truncated)", binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());

    // For closure test
    fHistManager.CreateTH2("hZgPartLevelClosureNoResp", "Zg response at particle level (closure test, jets not used for the response matrix)", binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
    fHistManager.CreateTH2("hZgDetLevelClosureNoResp", "Zg response at detector level (closure test, jets not used for the response matrix)", binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
    fHistManager.CreateTH2("hZgPartLevelClosureResp", "Zg response at particle level (closure test, jets used for the response matrix)", binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
    fHistManager.CreateTH2("hZgDetLevelClosureResp", "Zg response at detector level (closure test, jets used for the response matrix)", binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());

    RooUnfoldResponse *r = new RooUnfoldResponse("hZgResponse", "z_{g} response matrix");
    r->Setup((TH1 *)fHistManager.FindObject("hZgDetLevel"), (TH1 *)fHistManager.FindObject("hZgPartLevel"));
    fZgResponse.push_back(r);
    RooUnfoldResponse *r_closure = new RooUnfoldResponse("hZgResponseClosure", "z_{g} response matrix for the closure test");
    r_closure->Setup((TH1 *)fHistManager.FindObject("hZgDetLevel"), (TH1 *)fHistManager.FindObject("hZgPartLevel"));
    fZgResponseClosure.push_back(r_closure);
  }

  TIter next(fHistManager.GetListOfHistograms());
  TObject *obj = 0;
  while ((obj = next()))
  {
    fOutput->Add(obj);
  }
  for (int i = 0; i < fZgResponse.size(); i++)
  {
    fOutput->Add(fZgResponse.at(i));
  }
  for (int i = 0; i < fZgResponseClosure.size(); i++)
  {
    fOutput->Add(fZgResponseClosure.at(i));
  }

  PostData(1, fOutput);
}

Bool_t AliAnalysisTaskEmcalSoftDropResponse::CheckMCOutliers()
{
  if (!fMCRejectFilter)
    return true;
  if (!(fIsPythia || fIsHerwig))
    return true; // Only relevant for pt-hard production
  AliDebugStream(1) << "Using custom MC outlier rejection" << std::endl;
  auto partjets = GetJetContainer(fNamePartLevelJetContainer);
  if (!partjets)
    return true;

  // Check whether there is at least one particle level jet with pt above n * event pt-hard
  auto jetiter = partjets->accepted();
  auto max = std::max_element(jetiter.begin(), jetiter.end(), [](const AliEmcalJet *lhs, const AliEmcalJet *rhs) { return lhs->Pt() < rhs->Pt(); });
  if (max != jetiter.end())
  {
    // At least one jet found with pt > n * pt-hard
    AliDebugStream(1) << "Found max jet with pt " << (*max)->Pt() << " GeV/c" << std::endl;
    if ((*max)->Pt() > fPtHardAndJetPtFactor * fPtHard)
      return false;
  }
  return true;
}

bool AliAnalysisTaskEmcalSoftDropResponse::Run()
{
  AliJetContainer *partLevelJets = this->GetJetContainer(fNamePartLevelJetContainer),
                  *detLevelJets = GetJetContainer(fNameDetLevelJetContainer);
  AliClusterContainer *clusters = GetClusterContainer(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::ClusterContainerNameFactory(fInputEvent->IsA() == AliAODEvent::Class()));
  AliTrackContainer *tracks = GetTrackContainer(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::TrackContainerNameFactory(fInputEvent->IsA() == AliAODEvent::Class()));
  AliParticleContainer *particles = GetParticleContainer(fNameMCParticles.Data());
  if (!(partLevelJets || detLevelJets))
  {
    AliErrorStream() << "Either of the jet containers not found" << std::endl;
    return kFALSE;
  }

  if (fSampleFraction < 1.)
  {
    if (fSampleTrimmer->Uniform() > fSampleFraction)
      return false;
  }

  // get truncations at detector level
  TString histname;
  if (fForceBeamType != kpp)
  {
    histname = TString::Format("hZgDetLevel_%d", fCentBin);
  }
  else
  {
    histname = "hZgDetLevel";
  }
  TH2D *fZgDetLevel = (TH2D *)fHistManager.FindObject(histname);
  auto ptmindet = fZgDetLevel->GetYaxis()->GetBinLowEdge(1),
       ptmaxdet = fZgDetLevel->GetYaxis()->GetBinUpEdge(fZgDetLevel->GetYaxis()->GetNbins());

  //when  embedding and doing the constituent subtraction there is an additional step to the detector (or hybrid) to particle level matching because the detector jet (or hybrid) is the constituent subtracted jet which is not matched so we need to find the unsubtracted jet that it corresponds to and get the matched jets from there
  AliJetContainer *jetContUS = nullptr;
  if (fIsEmbeddedEvent)
    jetContUS = GetJetContainer(fNameUnSubLevelJetContainer);

  for (auto detjet : detLevelJets->accepted())
  {
    AliEmcalJet *partjet = nullptr;
    //variables for embedded pbpb data
    AliEmcalJet *jetUS = nullptr;
    Int_t ilab = -1;
    //for embedding, find the unsubtracted jet and get it's matched detector level jet
    if (fIsEmbeddedEvent)
    {
      for (Int_t i = 0; i < jetContUS->GetNJets(); i++)
      {
        jetUS = jetContUS->GetJet(i);
        if (jetUS->GetLabel() == detjet->GetLabel())
        {
          ilab = i;
          break;
        }
      }
      if (ilab == -1)
        continue;
      jetUS = jetContUS->GetJet(ilab);
      partjet = jetUS->ClosestJet();
    }
    //if we aren't embedding then just find the matched jet
    else
      partjet = detjet->ClosestJet();
    if (!partjet)
      continue;
    //one extra level of matching needed for embedding to go from detector to particle level
    if (fIsEmbeddedEvent)
    {
      partjet = partjet->ClosestJet();
      if (!partjet)
        continue;
    }

    //cut on the shared pt fraction, when embedding the unsubtracted jet should be used
    Double_t fraction = 0;
    if (fIsEmbeddedEvent)
      fraction = jetContUS->GetFractionSharedPt(jetUS);
    else
      fraction = detLevelJets->GetFractionSharedPt(detjet);
    if (fraction < fMinFractionShared)
      continue;

    // sample splitting (for closure test)
    bool closureUseResponse = (fSampleSplitter->Uniform() < fFractionResponseClosure);

    // Get the softdrop response
    std::vector<double> softdropDet, softdropPart;

    try
    {
      softdropDet = MakeSoftdrop(*detjet, detLevelJets->GetJetRadius(), tracks, clusters);
      softdropPart = MakeSoftdrop(*partjet, partLevelJets->GetJetRadius(), particles, nullptr);
      TString histname;
      if (fForceBeamType != kpp)
      {
        histname = TString::Format("hZgPartLevel_%d", fCentBin);
        fHistManager.FillTH1(histname, softdropPart[0], partjet->Pt());
      }
      else
      {
        fHistManager.FillTH1("hZgPartLevel", softdropPart[0], partjet->Pt());
      }
      if (detjet->Pt() >= ptmindet && detjet->Pt() <= ptmaxdet)
      {
        if (fForceBeamType != kpp)
        {
          histname = TString::Format("hZgPartLevelTruncated_%d", fCentBin);
          fHistManager.FillTH2(histname, softdropPart[0], partjet->Pt());
          histname = TString::Format("hZgDetLevel_%d", fCentBin);
          fHistManager.FillTH2(histname, softdropDet[0], detjet->Pt());
          fZgResponse.at(fCentBin)->Fill(softdropDet[0], detjet->Pt(), softdropPart[0], partjet->Pt());
        }
        else
        {
          fHistManager.FillTH2("hZgPartLevelTruncated", softdropPart[0], partjet->Pt());
          fHistManager.FillTH2("hZgDetLevel", softdropDet[0], detjet->Pt());
          fZgResponse.at(0)->Fill(softdropDet[0], detjet->Pt(), softdropPart[0], partjet->Pt());
        }
        if (closureUseResponse)
        {
          if (fForceBeamType != kpp)
          {
            fZgResponseClosure.at(fCentBin)->Fill(softdropDet[0], detjet->Pt(), softdropPart[0], partjet->Pt());
            histname = TString::Format("hZgDetLevelClosureResp_%d", fCentBin);
            fHistManager.FillTH2(histname, softdropDet[0], detjet->Pt());
            histname = TString::Format("hZgPartLevelClosureResp_%d", fCentBin);
            fHistManager.FillTH2(histname, softdropPart[0], partjet->Pt());
          }
          else
          {
            fZgResponseClosure.at(0)->Fill(softdropDet[0], detjet->Pt(), softdropPart[0], partjet->Pt());
            fHistManager.FillTH2("hZgDetLevelClosureResp", softdropDet[0], detjet->Pt());
            fHistManager.FillTH2("hZgPartLevelClosureResp", softdropPart[0], partjet->Pt());
          }
        }
        else
        {
          if (fForceBeamType != kpp)
          {
            histname = TString::Format("hZgPartLevelClosureNoResp_%d", fCentBin);
            fHistManager.FillTH2(histname, softdropPart[0], partjet->Pt());
            histname = TString::Format("hZgDetLevelClosureNoResp_%d", fCentBin);
            fHistManager.FillTH2(histname, softdropDet[0], detjet->Pt());
          }
          else
          {
            fHistManager.FillTH2("hZgDetLevelClosureNoResp", softdropDet[0], detjet->Pt());
            fHistManager.FillTH2("hZgPartLevelClosureNoResp", softdropPart[0], partjet->Pt());
          }
        }
      }
    }
    catch (...)
    {
      AliErrorStream() << "Error in softdrop evaluation - jet will be ignored" << std::endl;
      continue;
    }
  }
  return kTRUE;
}

std::vector<double> AliAnalysisTaskEmcalSoftDropResponse::MakeSoftdrop(const AliEmcalJet &jet, double jetradius, const AliParticleContainer *tracks, const AliClusterContainer *clusters) const
{
  const int kClusterOffset = 30000; // In order to handle tracks and clusters in the same index space the cluster index needs and offset, large enough so that there is no overlap with track indices
  std::vector<fastjet::PseudoJet> constituents;
  bool isMC = dynamic_cast<const AliMCParticleContainer *>(tracks);
  AliDebugStream(2) << "Make new jet substrucutre for " << (isMC ? "MC" : "data") << " jet: Number of tracks " << jet.GetNumberOfTracks() << ", clusters " << jet.GetNumberOfClusters() << std::endl;
  if (tracks && (fUseChargedConstituents || isMC))
  { // Neutral particles part of particle container in case of MC
    AliDebugStream(1) << "Jet substructure: Using charged constituents" << std::endl;
    for (int itrk = 0; itrk < jet.GetNumberOfTracks(); itrk++)
    {
      auto track = jet.Track(itrk);
      if (!track->Charge() && !fUseNeutralConstituents)
        continue; // Reject neutral constituents in case of using only charged consituents
      if (track->Charge() && !fUseChargedConstituents)
        continue; // Reject charged constituents in case of using only neutral consituents
      fastjet::PseudoJet constituentTrack(track->Px(), track->Py(), track->Pz(), track->E());
      constituentTrack.set_user_index(jet.TrackAt(itrk));
      constituents.push_back(constituentTrack);
    }
  }

  if (clusters && fUseNeutralConstituents)
  {
    AliDebugStream(1) << "Jet substructure: Using neutral constituents" << std::endl;
    for (int icl = 0; icl < jet.GetNumberOfClusters(); icl++)
    {
      auto cluster = jet.ClusterAt(icl, clusters->GetArray());
      TLorentzVector clustervec;
      cluster->GetMomentum(clustervec, fVertex, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy());
      fastjet::PseudoJet constituentCluster(clustervec.Px(), clustervec.Py(), clustervec.Pz(), cluster->GetHadCorrEnergy());
      constituentCluster.set_user_index(jet.ClusterAt(icl) + kClusterOffset);
      constituents.push_back(constituentCluster);
    }
  }

  AliDebugStream(3) << "Found " << constituents.size() << " constituents for jet with pt=" << jet.Pt() << " GeV/c" << std::endl;
  if (!constituents.size())
  {
    AliErrorStream() << "Jet has 0 constituents." << std::endl;
    throw 1;
  }
  // Redo jet finding on constituents with a
  fastjet::JetDefinition jetdef(fastjet::antikt_algorithm, jetradius * 2, static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30);
  fastjet::ClusterSequence jetfinder(constituents, jetdef);
  std::vector<fastjet::PseudoJet> outputjets = jetfinder.inclusive_jets(0);
  auto sdjet = outputjets[0];
  fastjet::contrib::SoftDrop softdropAlgorithm(fBeta, fZcut);
  softdropAlgorithm.set_verbose_structure(kTRUE);
  fastjet::JetAlgorithm reclusterizingAlgorithm;
  switch (fReclusterizer)
  {
  case kCAAlgo:
    reclusterizingAlgorithm = fastjet::cambridge_aachen_algorithm;
    break;
  case kKTAlgo:
    reclusterizingAlgorithm = fastjet::kt_algorithm;
    break;
  case kAKTAlgo:
    reclusterizingAlgorithm = fastjet::antikt_algorithm;
    break;
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

  std::vector<double> result = {softdropstruct.symmetry(),
                                groomed.m(),
                                softdropstruct.delta_R(),
                                groomed.perp(),
                                softdropstruct.mu(),
                                static_cast<double>(softdropstruct.dropped_count())};
  return result;
}

TBinning *AliAnalysisTaskEmcalSoftDropResponse::GetDefaultPartLevelPtBinning() const
{
  auto binning = new TCustomBinning;
  binning->SetMinimum(0);
  switch (fBinningMode)
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

TBinning *AliAnalysisTaskEmcalSoftDropResponse::GetDefaultDetLevelPtBinning() const
{
  auto binning = new TCustomBinning;
  switch (fBinningMode)
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

TBinning *AliAnalysisTaskEmcalSoftDropResponse::GetZgBinning() const
{
  auto binning = new TCustomBinning;
  binning->SetMinimum(0.);
  binning->AddStep(0.1, 0.1);
  binning->AddStep(0.5, 0.05);
  return binning;
}

AliAnalysisTaskEmcalSoftDropResponse *AliAnalysisTaskEmcalSoftDropResponse::AddTaskEmcalSoftDropResponse(Double_t jetradius, AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recombinationScheme, bool ifembed, const char *namepartcont, const char *trigger)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  Bool_t isAOD(kFALSE);
  AliInputEventHandler *inputhandler = static_cast<AliInputEventHandler *>(mgr->GetInputEventHandler());
  if (inputhandler)
  {
    if (inputhandler->IsA() == AliAODInputHandler::Class())
    {
      std::cout << "Analysing AOD events\n";
      isAOD = kTRUE;
    }
    else
    {
      std::cout << "Analysing ESD events\n";
    }
  }

  std::stringstream taskname;
  taskname << "SoftdropResponsemaker_R" << std::setw(2) << std::setfill('0') << int(jetradius * 10) << trigger;
  AliAnalysisTaskEmcalSoftDropResponse *responsemaker = new AliAnalysisTaskEmcalSoftDropResponse(taskname.str().data());
  responsemaker->SelectCollisionCandidates(AliVEvent::kINT7);
  responsemaker->SetIsEmbeddedEvent(ifembed);
  mgr->AddTask(responsemaker);

  TString partcontname(namepartcont);
  if (partcontname == "usedefault")
    partcontname = "mcparticles";
  AliParticleContainer *particles = responsemaker->AddMCParticleContainer(partcontname.Data());
  //if embedding need to specify that the particles are embedded
  if (ifembed)
    particles->SetIsEmbedding(true);
  particles->SetMinPt(0.);
  responsemaker->SetNameMCParticleContainer(partcontname.Data());

  TString partLevelTag("Jet");
  //if embedding need to specify the tag of the particle level jet from the embedding framework
  if (ifembed)
    partLevelTag = "partLevelJets";
  AliJetContainer *mcjets = responsemaker->AddJetContainer(
      jettype,
      AliJetContainer::antikt_algorithm,
      recombinationScheme,
      jetradius,
      ((jettype == AliJetContainer::kFullJet) || (jettype == AliJetContainer::kNeutralJet)) ? AliEmcalJet::kEMCALfid : AliEmcalJet::kTPC,
      particles, nullptr, partLevelTag);
  mcjets->SetJetPtCut(0.);
  mcjets->SetMaxTrackPt(1000.);
  responsemaker->SetNamePartLevelJetContainer(mcjets->GetName());

  AliTrackContainer *tracks(nullptr);
  if ((jettype == AliJetContainer::kChargedJet) || (jettype == AliJetContainer::kFullJet))
  {
    tracks = responsemaker->AddTrackContainer(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::TrackContainerNameFactory(isAOD));
    std::cout << "Track container name: " << tracks->GetName() << std::endl;
    //if embedding need to specify that the tracks are embedded
    if (ifembed)
      tracks->SetIsEmbedding(true);
    tracks->SetMinPt(0.15);
  }
  AliClusterContainer *clusters(nullptr);
  if ((jettype == AliJetContainer::kFullJet) || (jettype == AliJetContainer::kNeutralJet))
  {
    std::cout << "Using full or neutral jets ..." << std::endl;
    clusters = responsemaker->AddClusterContainer(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::ClusterContainerNameFactory(isAOD));
    std::cout << "Cluster container name: " << clusters->GetName() << std::endl;
    clusters->SetClusHadCorrEnergyCut(0.3); // 300 MeV E-cut
    clusters->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }
  else
  {
    std::cout << "Using charged jets ... " << std::endl;
  }

  TString detLevelTag("Jet");
  //if embedding need to specify the tag of the detector (or hybrid) level jet from the embedding framework
  if (ifembed)
    detLevelTag = "hybridLevelJets";
  AliJetContainer *datajets = responsemaker->AddJetContainer(
      jettype,
      AliJetContainer::antikt_algorithm,
      recombinationScheme,
      jetradius,
      ((jettype == AliJetContainer::kFullJet) || (jettype == AliJetContainer::kNeutralJet)) ? AliEmcalJet::kEMCALfid : AliEmcalJet::kTPCfid,
      tracks, clusters, detLevelTag);
  datajets->SetJetPtCut(0.);
  datajets->SetMaxTrackPt(1000.);
  //if embedding then this jet is the unsubtracted jet and the subtracted jet is added in the run macro, if not then it is the detector level jet
  if (!ifembed)
    responsemaker->SetNameDetLevelJetContainer(datajets->GetName());
  else
    responsemaker->SetNameUnSubLevelJetContainer(datajets->GetName());

  std::string jettypestring;
  switch (jettype)
  {
  case AliJetContainer::kFullJet:
    jettypestring = "FullJets";
    break;
  case AliJetContainer::kChargedJet:
    jettypestring = "ChargedJets";
    break;
  case AliJetContainer::kNeutralJet:
    jettypestring = "NeutralJets";
    break;
  default:
    jettypestring = "Undef";
  };

  EBinningMode_t binmode(kSDModeINT7);
  std::string triggerstring(trigger);
  if (triggerstring == "EJ1")
    binmode = kSDModeEJ1;
  else if (triggerstring == "EJ2")
    binmode = kSDModeEJ2;
  responsemaker->SetBinningMode(binmode);

  // Connecting containers
  std::stringstream outputfile, histname;
  outputfile << mgr->GetCommonFileName() << ":SoftDropResponse_" << jettypestring << "_R" << std::setw(2) << std::setfill('0') << int(jetradius * 10.) << "_" << trigger;
  histname << "SoftDropResponseHistos_" << jettypestring << "_R" << std::setw(2) << std::setfill('0') << int(jetradius * 10.) << "_" << trigger;
  mgr->ConnectInput(responsemaker, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(responsemaker, 1, mgr->CreateContainer(histname.str().data(), AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, outputfile.str().data()));

  return responsemaker;
}
