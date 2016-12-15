/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: R. Haake.                                                      *
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

#include <algorithm>
#include <vector>
#include <TClonesArray.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THn.h>
#include <TTree.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliEmcalPythiaInfo.h"
#include "AliMCEvent.h"
#include "AliPythia.h"
#include "AliStack.h"

#include "AliVTrack.h"
#include "AliVHeader.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliTrackContainer.h"
#include "AliAODTrack.h"
#include "AliPicoTrack.h"
#include "AliVParticle.h"
#include "TRandom3.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliAnalysisTaskJetExtractor.h"

#include "AliAnalysisTaskChargedJetsHadronCF.h"

//________________________________________________________________________
AliChargedJetsHadronCFCuts::~AliChargedJetsHadronCFCuts() 
{
// dummy destructor
}

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskChargedJetsHadronCF)
/// \endcond
//________________________________________________________________________
AliAnalysisTaskChargedJetsHadronCF::AliAnalysisTaskChargedJetsHadronCF() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskChargedJetsHadronCF", kTRUE),
  fJetsCont(0),
  fTracksCont(0),
  fJetsTree(0),
  fJetsTreeBuffer(0),
  fExtractionPercentage(0),
  fExtractionMinPt(0),
  fExtractionMaxPt(0),
  fEventExtractionPercentage(0),
  fEventExtractionMinJetPt(0),
  fEventExtractionMaxJetPt(0),
  fHadronMatchingRadius(0.5),
  fConstPtFilterBit(1024),
  fNumberOfCentralityBins(10),
  fJetsOutput(),
  fTracksOutput(),
  fJetParticleArrayName("JetsDPhiBasicParticles"),
  fTrackParticleArrayName(""),
  fJetEmbeddingArray(),
  fJetEmbeddingArrayName(""),
  fJetEmbeddingTrackArrayName(""),
  fJetEmbeddingMaxDistance(0.3),
  fJetEmbeddingNumMatchedJets(2),
  fJetEmbeddingUsePerTrackMCPercentage(kTRUE),
  fJetEmbeddingUseBgrdForMCPercentage(kFALSE),
  fJetEmbeddingCreatePtPlotPerCut(kFALSE),
  fJetEmbeddingCuts(),
  fJetVetoArray(),
  fJetVetoArrayName(""),
  fJetVetoJetByJet(1),
  fMatchedJets(),
  fRandom(0),
  fJetOutputMode(0),
  fPythiaExtractionMode(0),
  fPythiaExtractionUseHadronMatching(kFALSE),
  fLeadingJet(0),
  fSubleadingJet(0),
  fInitialPartonMatchedJet1(0),
  fInitialPartonMatchedJet2(0),
  fAcceptedJets(0),
  fAcceptedTracks(0)
{
  // Default constructor.
  SetMakeGeneralHistograms(kTRUE);
  fRandom = new TRandom3(0);
}


//________________________________________________________________________
AliAnalysisTaskChargedJetsHadronCF::AliAnalysisTaskChargedJetsHadronCF(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fJetsCont(0),
  fTracksCont(0),
  fJetsTree(0),
  fJetsTreeBuffer(0),
  fExtractionPercentage(0),
  fExtractionMinPt(0),
  fExtractionMaxPt(0),
  fEventExtractionPercentage(0),
  fEventExtractionMinJetPt(0),
  fEventExtractionMaxJetPt(0),
  fHadronMatchingRadius(0.5),
  fConstPtFilterBit(1024),
  fNumberOfCentralityBins(10),
  fJetsOutput(),
  fTracksOutput(),
  fJetParticleArrayName("JetsDPhiBasicParticles"),
  fTrackParticleArrayName(""),
  fJetEmbeddingArray(),
  fJetEmbeddingArrayName(""),
  fJetEmbeddingTrackArrayName(""),
  fJetEmbeddingMaxDistance(0.3),
  fJetEmbeddingNumMatchedJets(2),
  fJetEmbeddingUsePerTrackMCPercentage(kTRUE),
  fJetEmbeddingUseBgrdForMCPercentage(kFALSE),
  fJetEmbeddingCreatePtPlotPerCut(kFALSE),
  fJetEmbeddingCuts(),
  fJetVetoArray(),
  fJetVetoArrayName(""),
  fJetVetoJetByJet(1),
  fMatchedJets(),
  fRandom(0),
  fJetOutputMode(0),
  fPythiaExtractionMode(0),
  fPythiaExtractionUseHadronMatching(kFALSE),
  fLeadingJet(0),
  fSubleadingJet(0),
  fInitialPartonMatchedJet1(0),
  fInitialPartonMatchedJet2(0),
  fAcceptedJets(0),
  fAcceptedTracks(0)
{
  // Constructor
  SetMakeGeneralHistograms(kTRUE);
  fRandom = new TRandom3(0);
}

//________________________________________________________________________
AliAnalysisTaskChargedJetsHadronCF::~AliAnalysisTaskChargedJetsHadronCF()
{
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronCF::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  // ### Basic container settings
  fJetsCont           = GetJetContainer(0);
  if(fJetsCont) { //get particles connected to jets
    fJetsCont->PrintCuts();
    fTracksCont       = static_cast<AliTrackContainer*>(fJetsCont->GetParticleContainer());
  } else {        //no jets, just analysis tracks
    fTracksCont       = static_cast<AliTrackContainer*>(GetParticleContainer(0));
  }
  if(fTracksCont) fTracksCont->SetClassName("AliVTrack");

  // ### Create all histograms
  fHistEventRejection->GetXaxis()->SetBinLabel(15,"JetVeto");

  // Track QA plots
  AddHistogram2D<TH2D>("hTrackPt", "Tracks p_{T} distribution", "", 300, 0., 300., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)", "Centrality", "dN^{Tracks}/dp_{T}");
  AddHistogram2D<TH2D>("hTrackPhi", "Track angular distribution in #phi", "LEGO2", 180, 0., 2*TMath::Pi(), fNumberOfCentralityBins, 0, 100, "#phi", "Centrality", "dN^{Tracks}/(d#phi)");
  AddHistogram2D<TH2D>("hTrackEta", "Track angular distribution in #eta", "LEGO2", 100, -2.5, 2.5, fNumberOfCentralityBins, 0, 100, "#eta", "Centrality", "dN^{Tracks}/(d#eta)");
  AddHistogram2D<TH2D>("hTrackPhiEta", "Track angular distribution #phi/#eta", "COLZ", 180, 0., 2*TMath::Pi(), 100, -2.5, 2.5, "#phi", "#eta", "dN^{Tracks}/d#phi d#eta");

  AddHistogram2D<TH2D>("hLeadingTrackPt", "Leading tracks p_{T} distribution", "", 300, 0., 300., fNumberOfCentralityBins, 0, 100, "p_{T} (GeV/c)", "Centrality", "dN^{Tracks}/dp_{T}");
  AddHistogram2D<TH2D>("hLeadingTrackPhi", "Leading tracks angular distribution in #phi", "LEGO2", 180, 0., 2*TMath::Pi(), fNumberOfCentralityBins, 0, 100, "#phi", "Centrality", "dN^{Tracks}/(d#phi)");
  AddHistogram2D<TH2D>("hLeadingTrackEta", "Leading tracks angular distribution in #eta", "LEGO2", 100, -2.5, 2.5, fNumberOfCentralityBins, 0, 100, "#eta", "Centrality", "dN^{Tracks}/(d#eta)");
  AddHistogram2D<TH2D>("hLeadingTrackPhiEta", "Track angular distribution #phi/#eta", "COLZ", 180, 0., 2*TMath::Pi(), 100, -2.5, 2.5, "#phi", "#eta", "dN^{Tracks}/d#phi d#eta");

  AddHistogram2D<TH2D>("hTrackEtaPt", "Track angular distribution in #eta vs. p_{T}", "LEGO2", 100, -2.5, 2.5, 300, 0., 300., "#eta", "p_{T} (GeV/c)", "dN^{Tracks}/(d#eta dp_{T})");
  AddHistogram2D<TH2D>("hTrackPhiPt", "Track angular distribution in #phi vs. p_{T}", "LEGO2", 180, 0, 2*TMath::Pi(), 300, 0., 300., "#phi", "p_{T} (GeV/c)", "dN^{Tracks}/(d#phi dp_{T})");

  // Create plots for each embedding cut
  for(Int_t i = -1; i<static_cast<Int_t>(fJetEmbeddingCuts.size()); i++)
  {
    const char* appendix = "";
    if(i>-1)
    {
      AliChargedJetsHadronCFCuts currentCut = fJetEmbeddingCuts.at(i);
      appendix = Form("_%s", currentCut.fCutName.Data());

      // Don't double-add cuts
      if( static_cast<TH1*>(fOutput->FindObject(Form("hJetPtRaw%s", appendix))) )
        continue;
    }
    // Jet QA plots
    AddHistogram2D<TH2D>(Form("hJetPtRaw%s", appendix), "Jets p_{T} distribution (no bgrd. corr.)", "", 300, 0., 300., fNumberOfCentralityBins, 0, 100, "p_{T, jet} (GeV/c)", "Centrality", "dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>(Form("hJetPt%s", appendix), "Jets p_{T} distribution (background subtracted)", "", 400, -100., 300., fNumberOfCentralityBins, 0, 100, "p_{T, jet} (GeV/c)", "Centrality", "dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>(Form("hJetPhi%s", appendix), "Jet angular distribution #phi", "LEGO2", 180, 0., 2*TMath::Pi(), fNumberOfCentralityBins, 0, 100, "#phi", "Centrality", "dN^{Jets}/d#phi");
    AddHistogram2D<TH2D>(Form("hJetEta%s", appendix), "Jet angular distribution #eta", "LEGO2", 100, -2.5, 2.5, fNumberOfCentralityBins, 0, 100, "#eta","Centrality","dN^{Jets}/d#eta");
    AddHistogram2D<TH2D>(Form("hJetPhiPt%s", appendix), "Jet angular distribution #phi vs. p_{T}", "LEGO2", 180, 0., 2*TMath::Pi(), 400, -100., 300., "#phi", "p_{T, jet} (GeV/c)", "dN^{Jets}/d#phi dp_{T}");
    AddHistogram2D<TH2D>(Form("hJetEtaPt%s", appendix), "Jet angular distribution #eta  vs. p_{T}", "LEGO2", 100, -2.5, 2.5, 400, -100., 300., "#eta","p_{T, jet} (GeV/c)","dN^{Jets}/d#eta dp_{T}");
    AddHistogram2D<TH2D>(Form("hJetPhiEta%s", appendix), "Jet angular distribution #phi/#eta", "COLZ", 180, 0., 2*TMath::Pi(), 100, -2.5, 2.5, "#phi", "#eta", "dN^{Jets}/d#phi d#eta");
    AddHistogram2D<TH2D>(Form("hJetArea%s", appendix), "Jet area", "LEGO2", 200, 0., 2., fNumberOfCentralityBins, 0, 100, "Jet A", "Centrality", "dN^{Jets}/dA");
    AddHistogram2D<TH2D>(Form("hJetAreaPt%s", appendix), "Jet area vs. p_{T}", "LEGO2", 200, 0., 2., 400, -100., 300., "Jet A", "p_{T, jet} (GeV/c)", "dN^{Jets}/dA dp_{T}");
    AddHistogram2D<TH2D>(Form("hJetPtLeadingHadron%s", appendix), "Jet leading hadron p_{T} distribution vs. jet p_{T}", "", 300, 0., 300., 300, 0., 300., "p_{T, jet} (GeV/c)", "p_{T,lead had} (GeV/c)", "dN^{Jets}/dp_{T}dp_{T,had}");

    // Leading/subleading ...
    AddHistogram2D<TH2D>(Form("hLeadingJetPtRaw%s", appendix), "Jets p_{T} distribution (no bgrd. corr.)", "", 300, 0., 300., fNumberOfCentralityBins, 0, 100, "p_{T, jet} (GeV/c)", "Centrality", "dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>(Form("hLeadingJetPt%s", appendix), "Jets p_{T} distribution (background subtracted)", "", 400, -100., 300., fNumberOfCentralityBins, 0, 100, "p_{T, jet} (GeV/c)", "Centrality", "dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>(Form("hLeadingJetPhi%s", appendix), "Jet angular distribution #phi", "LEGO2", 180, 0., 2*TMath::Pi(), fNumberOfCentralityBins, 0, 100, "#phi", "Centrality", "dN^{Jets}/d#phi");
    AddHistogram2D<TH2D>(Form("hLeadingJetEta%s", appendix), "Jet angular distribution #eta", "LEGO2", 100, -2.5, 2.5, fNumberOfCentralityBins, 0, 100, "#eta","Centrality","dN^{Jets}/d#eta");
    AddHistogram2D<TH2D>(Form("hLeadingJetPhiPt%s", appendix), "Jet angular distribution #phi vs. p_{T}", "LEGO2", 180, 0., 2*TMath::Pi(), 400, -100., 300., "#phi", "p_{T, jet} (GeV/c)", "dN^{Jets}/d#phi dp_{T}");
    AddHistogram2D<TH2D>(Form("hLeadingJetEtaPt%s", appendix), "Jet angular distribution #eta  vs. p_{T}", "LEGO2", 100, -2.5, 2.5, 400, -100., 300., "#eta","p_{T, jet} (GeV/c)","dN^{Jets}/d#eta dp_{T}");
    AddHistogram2D<TH2D>(Form("hLeadingJetPhiEta%s", appendix), "Jet angular distribution #phi/#eta", "COLZ", 180, 0., 2*TMath::Pi(), 100, -2.5, 2.5, "#phi", "#eta", "dN^{Jets}/d#phi d#eta");
    AddHistogram2D<TH2D>(Form("hLeadingJetArea%s", appendix), "Jet area", "LEGO2", 200, 0., 2., fNumberOfCentralityBins, 0, 100, "Jet A", "Centrality", "dN^{Jets}/dA");
    AddHistogram2D<TH2D>(Form("hLeadingJetAreaPt%s", appendix), "Jet area vs. p_{T}", "LEGO2", 200, 0., 2., 400, -100., 300., "Jet A", "p_{T, jet} (GeV/c)", "dN^{Jets}/dA dp_{T}");
    AddHistogram2D<TH2D>(Form("hLeadingJetPtLeadingHadron%s", appendix), "Jet leading hadron p_{T} distribution vs. jet p_{T}", "", 300, 0., 300., 300, 0., 300., "p_{T, jet} (GeV/c)", "p_{T,lead had} (GeV/c)", "dN^{Jets}/dp_{T}dp_{T,had}");

    AddHistogram2D<TH2D>(Form("hSubleadingJetPtRaw%s", appendix), "Jets p_{T} distribution (no bgrd. corr.)", "", 300, 0., 300., fNumberOfCentralityBins, 0, 100, "p_{T, jet} (GeV/c)", "Centrality", "dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>(Form("hSubleadingJetPt%s", appendix), "Jets p_{T} distribution (background subtracted)", "", 400, -100., 300., fNumberOfCentralityBins, 0, 100, "p_{T, jet} (GeV/c)", "Centrality", "dN^{Jets}/dp_{T}");
    AddHistogram2D<TH2D>(Form("hSubleadingJetPhi%s", appendix), "Jet angular distribution #phi", "LEGO2", 180, 0., 2*TMath::Pi(), fNumberOfCentralityBins, 0, 100, "#phi", "Centrality", "dN^{Jets}/d#phi");
    AddHistogram2D<TH2D>(Form("hSubleadingJetEta%s", appendix), "Jet angular distribution #eta", "LEGO2", 100, -2.5, 2.5, fNumberOfCentralityBins, 0, 100, "#eta","Centrality","dN^{Jets}/d#eta");
    AddHistogram2D<TH2D>(Form("hSubleadingJetPhiPt%s", appendix), "Jet angular distribution #phi vs. p_{T}", "LEGO2", 180, 0., 2*TMath::Pi(), 400, -100., 300., "#phi", "p_{T, jet} (GeV/c)", "dN^{Jets}/d#phi dp_{T}");
    AddHistogram2D<TH2D>(Form("hSubleadingJetEtaPt%s", appendix), "Jet angular distribution #eta  vs. p_{T}", "LEGO2", 100, -2.5, 2.5, 400, -100., 300., "#eta","p_{T, jet} (GeV/c)","dN^{Jets}/d#eta dp_{T}");
    AddHistogram2D<TH2D>(Form("hSubleadingJetPhiEta%s", appendix), "Jet angular distribution #phi/#eta", "COLZ", 180, 0., 2*TMath::Pi(), 100, -2.5, 2.5, "#phi", "#eta", "dN^{Jets}/d#phi d#eta");
    AddHistogram2D<TH2D>(Form("hSubleadingJetArea%s", appendix), "Jet area", "LEGO2", 200, 0., 2., fNumberOfCentralityBins, 0, 100, "Jet A", "Centrality", "dN^{Jets}/dA");
    AddHistogram2D<TH2D>(Form("hSubleadingJetAreaPt%s", appendix), "Jet area vs. p_{T}", "LEGO2", 200, 0., 2., 400, -100., 300., "Jet A", "p_{T, jet} (GeV/c)", "dN^{Jets}/dA dp_{T}");
    AddHistogram2D<TH2D>(Form("hSubleadingJetPtLeadingHadron%s", appendix), "Jet leading hadron p_{T} distribution vs. jet p_{T}", "", 300, 0., 300., 300, 0., 300., "p_{T, jet} (GeV/c)", "p_{T,lead had} (GeV/c)", "dN^{Jets}/dp_{T}dp_{T,had}");

    AddHistogram2D<TH2D>(Form("hJetConstituentPt_Cent0_100%s", appendix), "Jet constituent p_{T} distribution vs. jet p_T (background subtracted)", "", 400, -100., 300., 300, 0., 300., "p_{T, jet} (GeV/c)", "p_{T, track} (GeV/c)", "dN^{Tracks}/d^{2}p_{T}");
    AddHistogram2D<TH2D>(Form("hJetConstituentPt_Cent0_10%s", appendix), "Jet constituent p_{T} distribution vs. jet p_T (background subtracted), 0-10 centrality", "", 400, -100., 300., 300, 0., 300., "p_{T, jet} (GeV/c)", "p_{T, track} (GeV/c)", "dN^{Tracks}/d^{2}p_{T}");
    AddHistogram2D<TH2D>(Form("hJetConstituentPt_Cent10_30%s", appendix), "Jet constituent p_{T} distribution vs. jet p_T (background subtracted), 10-30 centrality", "", 400, -100., 300., 300, 0., 300., "p_{T, jet} (GeV/c)", "p_{T, track} (GeV/c)", "dN^{Tracks}/d^{2}p_{T}");
    AddHistogram2D<TH2D>(Form("hJetConstituentPt_Cent30_50%s", appendix), "Jet constituent p_{T} distribution vs. jet p_T (background subtracted), 30-50 centrality", "", 400, -100., 300., 300, 0., 300., "p_{T, jet} (GeV/c)", "p_{T, track} (GeV/c)", "dN^{Tracks}/d^{2}p_{T}");
    AddHistogram2D<TH2D>(Form("hJetConstituentPt_Cent50_90%s", appendix), "Jet constituent p_{T} distribution vs. jet p_T (background subtracted), 50-90 centrality", "", 400, -100., 300., 300, 0., 300., "p_{T, jet} (GeV/c)", "p_{T, track} (GeV/c)", "dN^{Tracks}/d^{2}p_{T}");

    AddHistogram2D<TH2D>(Form("hJetConstituentPt_Cent0_100_FilterBit%i%s", fConstPtFilterBit, appendix), "Jet constituent p_{T} distribution vs. jet p_T (background subtracted)", "", 400, -100., 300., 300, 0., 300., "p_{T, jet} (GeV/c)", "p_{T, track} (GeV/c)", "dN^{Tracks}/d^{2}p_{T}");
    AddHistogram2D<TH2D>(Form("hJetConstituentPt_Cent0_10_FilterBit%i%s", fConstPtFilterBit, appendix), "Jet constituent p_{T} distribution vs. jet p_T (background subtracted), 0-10 centrality", "", 400, -100., 300., 300, 0., 300., "p_{T, jet} (GeV/c)", "p_{T, track} (GeV/c)", "dN^{Tracks}/d^{2}p_{T}");
    AddHistogram2D<TH2D>(Form("hJetConstituentPt_Cent10_30_FilterBit%i%s", fConstPtFilterBit, appendix), "Jet constituent p_{T} distribution vs. jet p_T (background subtracted), 10-30 centrality", "", 400, -100., 300., 300, 0., 300., "p_{T, jet} (GeV/c)", "p_{T, track} (GeV/c)", "dN^{Tracks}/d^{2}p_{T}");
    AddHistogram2D<TH2D>(Form("hJetConstituentPt_Cent30_50_FilterBit%i%s", fConstPtFilterBit, appendix), "Jet constituent p_{T} distribution vs. jet p_T (background subtracted), 30-50 centrality", "", 400, -100., 300., 300, 0., 300., "p_{T, jet} (GeV/c)", "p_{T, track} (GeV/c)", "dN^{Tracks}/d^{2}p_{T}");
    AddHistogram2D<TH2D>(Form("hJetConstituentPt_Cent50_90_FilterBit%i%s", fConstPtFilterBit, appendix), "Jet constituent p_{T} distribution vs. jet p_T (background subtracted), 50-90 centrality", "", 400, -100., 300., 300, 0., 300., "p_{T, jet} (GeV/c)", "p_{T, track} (GeV/c)", "dN^{Tracks}/d^{2}p_{T}");

    AddHistogram2D<TH2D>(Form("hJetConstituentCount_Cent0_100%s", appendix), "Jet constituent count vs. jet p_T (background subtracted)", "", 400, -100., 300., 200, 0., 200., "p_{T, jet} (GeV/c)", "Count", "dN^{Jets}/dNdp_{T}");
    AddHistogram2D<TH2D>(Form("hJetConstituentCount_Cent0_10%s", appendix), "Jet constituent count vs. jet p_T (background subtracted), 0-10 centrality", "", 400, -100., 300., 200, 0., 200., "p_{T, jet} (GeV/c)", "Count", "dN^{Jets}/dNdp_{T}");
  }

  // Embedding plots
  if(fJetOutputMode == 4  || fJetOutputMode == 5)
  {
    Double_t maxRatio = 1.;
    if(fJetEmbeddingUseBgrdForMCPercentage)
      maxRatio = 2.;

    for(Int_t i = -1; i<static_cast<Int_t>(fJetEmbeddingCuts.size()); i++)
    {
      const char* appendix = "";
      if(i>-1)
      {
        AliChargedJetsHadronCFCuts currentCut = fJetEmbeddingCuts.at(i);
        appendix = Form("_%s", currentCut.fCutName.Data());

        // Don't double-add cuts
        if( static_cast<TH1*>(fOutput->FindObject(Form("hEmbeddingDeltaR%s", appendix))) )
          continue;
      }
      AddHistogram2D<TH2D>(Form("hEmbeddingDeltaR%s", appendix), "Matched jet #Delta R distribution", "", 200, -50., 150., 100, 0, 1.0, "p_{T, jet} (GeV/c)", "#Delta R", "dN^{Matched}/dp_{T}dR");
      AddHistogram2D<TH2D>(Form("hEmbeddingDeltaEta%s", appendix), "Matched jet #Delta #eta distribution", "", 200, -50., 150., 100, -1.0, 1.0, "p_{T, jet} (GeV/c)", "#Delta #eta", "dN^{Matched}/dp_{T}d#eta");
      AddHistogram2D<TH2D>(Form("hEmbeddingDeltaPhi%s", appendix), "Matched jet #Delta #phi distribution", "", 200, -50., 150., 100, -1.0, 1.0, "p_{T, jet} (GeV/c)", "#Delta #phi", "dN^{Matched}/dp_{T}d#phi");
      AddHistogram1D<TH1D>(Form("hEmbeddingJetPt%s", appendix), "Embedded jets p_{T} distribution", "", 200, -50., 150., "p_{T, jet} (GeV/c)", "dN/dp_{T}");
      AddHistogram2D<TH2D>(Form("hEmbeddingJetPhiEta%s", appendix), "Embedded jet angular distribution #phi/#eta", "COLZ", 180, 0., 2*TMath::Pi(), 100, -2.5, 2.5, "#phi", "#eta", "dN^{Jets}/d#phi d#eta");

      if(fJetEmbeddingCreatePtPlotPerCut)
      {
        AddHistogram3D<TH3D>(Form("hEmbeddingPtCorr010%s", appendix), "Matched jet p_{T} distributions (0-10% centrality)", "", 150, 0., 150., 150, 0., 150., 100, 0., maxRatio, "p_{T, MC jet} (GeV/c)", "p_{T, emb} (GeV/c)", "% MC");
        AddHistogram3D<TH3D>(Form("hEmbeddingPtCorr1030%s", appendix), "Matched jet p_{T} distributions (10-30% centrality)", "", 150, 0., 150., 150, 0., 150., 100, 0., maxRatio, "p_{T, MC jet} (GeV/c)", "p_{T, emb} (GeV/c)", "% MC");
        AddHistogram3D<TH3D>(Form("hEmbeddingPtCorr3050%s", appendix), "Matched jet p_{T} distributions (30-50% centrality)", "", 150, 0., 150., 150, 0., 150., 100, 0., maxRatio, "p_{T, MC jet} (GeV/c)", "p_{T, emb} (GeV/c)", "% MC");
        AddHistogram3D<TH3D>(Form("hEmbeddingPtCorr5090%s", appendix), "Matched jet p_{T} distributions (50-90% centrality)", "", 150, 0., 150., 150, 0., 150., 100, 0., maxRatio, "p_{T, MC jet} (GeV/c)", "p_{T, emb} (GeV/c)", "% MC");
      }
      else
      {
        AddHistogram2D<TH2D>(Form("hEmbeddingPtCorr010_Above20%s", appendix), "Matched jet p_{T} distributions, MC ratio > 20% (0-10% centrality)", "", 150, 0., 150., 150, 0., 150.,"p_{T, MC jet} (GeV/c)", "p_{T, emb} (GeV/c)", "% MC");
        AddHistogram2D<TH2D>(Form("hEmbeddingPtCorr1030_Above20%s", appendix),"Matched jet p_{T} distributions, MC ratio > 20% (10-30% centrality)","", 150, 0., 150., 150, 0., 150.,"p_{T, MC jet} (GeV/c)", "p_{T, emb} (GeV/c)", "% MC");
        AddHistogram2D<TH2D>(Form("hEmbeddingPtCorr3050_Above20%s", appendix),"Matched jet p_{T} distributions, MC ratio > 20% (30-50% centrality)","", 150, 0., 150., 150, 0., 150.,"p_{T, MC jet} (GeV/c)", "p_{T, emb} (GeV/c)", "% MC");
        AddHistogram2D<TH2D>(Form("hEmbeddingPtCorr5090_Above20%s", appendix),"Matched jet p_{T} distributions, MC ratio > 20% (50-90% centrality)","", 150, 0., 150., 150, 0., 150.,"p_{T, MC jet} (GeV/c)", "p_{T, emb} (GeV/c)", "% MC");
      }
    }

    if(!fJetEmbeddingCreatePtPlotPerCut)
    {
      AddHistogram3D<TH3D>("hEmbeddingPtCorr010", "Matched jet p_{T} distributions (0-10% centrality)", "", 150, 0., 150., 150, 0., 150., 100, 0., maxRatio, "p_{T, MC jet} (GeV/c)", "p_{T, emb} (GeV/c)", "% MC");
      AddHistogram3D<TH3D>("hEmbeddingPtCorr1030", "Matched jet p_{T} distributions (10-30% centrality)", "", 150, 0., 150., 150, 0., 150., 100, 0., maxRatio, "p_{T, MC jet} (GeV/c)", "p_{T, emb} (GeV/c)", "% MC");
      AddHistogram3D<TH3D>("hEmbeddingPtCorr3050", "Matched jet p_{T} distributions (30-50% centrality)", "", 150, 0., 150., 150, 0., 150., 100, 0., maxRatio, "p_{T, MC jet} (GeV/c)", "p_{T, emb} (GeV/c)", "% MC");
      AddHistogram3D<TH3D>("hEmbeddingPtCorr5090", "Matched jet p_{T} distributions (50-90% centrality)", "", 150, 0., 150., 150, 0., 150., 100, 0., maxRatio, "p_{T, MC jet} (GeV/c)", "p_{T, emb} (GeV/c)", "% MC");
    }
  }

  // Random cone plots, background, ...
  AddHistogram2D<TH2D>("hRandomConePt", "Random cone p_{T} distribution", "", 400, -100., 300., fNumberOfCentralityBins, 0, 100, "p_{T, cone} (GeV/c)", "Centrality", "dN^{Tracks}/dp_{T}");
  AddHistogram2D<TH2D>("hRandomConePtCut3GeV", "Random cone p_{T} distribution, cut p_{T} > 3 GeV/c", "", 400, -100., 300., fNumberOfCentralityBins, 0, 100, "p_{T, cone} (GeV/c)", "Centrality", "dN^{Tracks}/dp_{T}");
  AddHistogram2D<TH2D>("hRandomConeRawPt", "Random cone p_{T} distribution (no bgrd. correction)", "", 300, 0., 300., fNumberOfCentralityBins, 0, 100, "p_{T, cone} (GeV/c)", "Centrality", "dN^{Tracks}/dp_{T}");
  AddHistogram2D<TH2D>("hRandomConeRawPtCut3GeV", "Random cone p_{T} distribution (no bgrd. correction), cut p_{T} > 3 GeV/c", "", 300, 0., 300., fNumberOfCentralityBins, 0, 100, "p_{T, cone} (GeV/c)", "Centrality", "dN^{Tracks}/dp_{T}");

  AddHistogram2D<TH2D>("hTrackCount", "Number of tracks in acceptance vs. centrality", "LEGO2", 500, 0., 5000., fNumberOfCentralityBins, 0, 100, "N tracks","Centrality", "dN^{Events}/dN^{Tracks}");
  AddHistogram2D<TH2D>("hJetCount", "Number of jets in acceptance vs. centrality", "LEGO2", 100, 0., 100., fNumberOfCentralityBins, 0, 100, "N Jets","Centrality", "dN^{Events}/dN^{Jets}");
  AddHistogram2D<TH2D>("hBackgroundPt", "Background p_{T} distribution", "", 150, 0., 150., fNumberOfCentralityBins, 0, 100, "Background p_{T} (GeV/c)", "Centrality", "dN^{Events}/dp_{T}");


  for(Int_t i = -1; i<static_cast<Int_t>(fJetEmbeddingCuts.size()); i++)
  {
    const char* appendix = "";
    if(i>-1)
    {
      AliChargedJetsHadronCFCuts currentCut = fJetEmbeddingCuts.at(i);
      appendix = Form("_%s", currentCut.fCutName.Data());

      // Don't double-add cuts
      if( static_cast<TH1*>(fOutput->FindObject(Form("hBackgroundPtJetPt_Cent0_100%s", appendix))) )
        continue;
    }
    AddHistogram2D<TH2D>(Form("hBackgroundPtJetPt_Cent0_100%s", appendix), "Background p_{T} distribution vs. jet p_{T}", "", 150, 0., 150.,  400, -100., 300., "Background p_{T} (GeV/c)", "Jet p_{T} (GeV/c)", "dN^{Events}/dp_{T}");
    AddHistogram2D<TH2D>(Form("hBackgroundPtJetPt_Cent0_10%s", appendix), "Background p_{T} distribution vs. jet p_{T}", "", 150, 0., 150.,  400, -100., 300., "Background p_{T} (GeV/c)", "Jet p_{T} (GeV/c)", "dN^{Events}/dp_{T}");
    AddHistogram2D<TH2D>(Form("hBackgroundPtConstCount_Cent0_100%s", appendix), "Background p_{T} distribution vs. const. count", "", 150, 0., 150., 200, 0., 200., "Background p_{T} (GeV/c)", "Count", "dN^{Events}/dp_{T}");
    AddHistogram2D<TH2D>(Form("hBackgroundPtConstCount_Cent0_10%s", appendix), "Background p_{T} distribution vs. const. count", "", 150, 0., 150., 200, 0., 200., "Background p_{T} (GeV/c)", "Count", "dN^{Events}/dp_{T}");
  }

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}


//________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronCF::ExecOnce() {

  AliAnalysisTaskEmcalJet::ExecOnce();


  // ### Add the jets as basic correlation particles to the event
  // This output object carries all accepted jets
  if (!(fInputEvent->FindListObject(Form("%s", fJetParticleArrayName.Data()))))
  {
    fJetsOutput.push_back(new TClonesArray("AliPicoTrack"));
    fJetsOutput.at(0)->SetName(fJetParticleArrayName.Data());
    fInputEvent->AddObject(fJetsOutput.at(0));
  }
  else
    AliFatal(Form("%s: Object with name %s already in event!", GetName(), Form("%s", fJetParticleArrayName.Data())));

  // These output objects carry all jets that pass certain cuts
  if( (fJetOutputMode==4 || fJetOutputMode==5) )
  {
    // before, check if all given output names are OK
    for(Int_t i = 0; i<fJetEmbeddingCuts.size(); i++)
      if (fInputEvent->FindListObject(Form("%s", fJetEmbeddingCuts.at(i).fOutputName.Data())))
        AliFatal(Form("%s: Object with name %s already in event!", GetName(), Form("%s", fJetEmbeddingCuts.at(i).fOutputName.Data())));

    for(Int_t i = 0; i<fJetEmbeddingCuts.size(); i++)
    {
      // If the cut demands a new output stream, add it
      if (!fInputEvent->FindListObject(Form("%s", fJetEmbeddingCuts.at(i).fOutputName.Data())))
      {
        fJetsOutput.push_back(new TClonesArray("AliPicoTrack"));
        fJetsOutput.at(fJetsOutput.size()-1)->SetName(fJetEmbeddingCuts.at(i).fOutputName.Data());
        fInputEvent->AddObject(fJetsOutput.at(fJetsOutput.size()-1));
        fJetEmbeddingCuts.at(i).fArrayIndex = fJetsOutput.size()-1;

        // Set the array indices for all cuts that use this output stream
        for(Int_t j = 0; j<fJetEmbeddingCuts.size(); j++)
        {
          if(fJetEmbeddingCuts.at(j).fArrayIndex != -1)
            continue;
          if(fJetEmbeddingCuts.at(j).fOutputName != fJetEmbeddingCuts.at(i).fOutputName)
            continue;
          fJetEmbeddingCuts.at(j).fArrayIndex = fJetEmbeddingCuts.at(i).fArrayIndex;
        }
      }
    }
  }

  // ##############################################
  // ##############################################

  // ### Add the tracks as basic correlation particles to the event (optional)
  if(fTrackParticleArrayName != "")
  {
    if (!(fInputEvent->FindListObject(Form("%s", fTrackParticleArrayName.Data()))))
    {
      fTracksOutput = new TClonesArray("AliPicoTrack");
      fTracksOutput->SetName(fTrackParticleArrayName.Data());
      fInputEvent->AddObject(fTracksOutput);
    }
    else
      AliFatal(Form("%s: Object with name %s already in event!", GetName(), Form("%s", fTrackParticleArrayName.Data())));
  }

  // ### Import jets for embedding (optional)
  if(fJetEmbeddingArrayName != "")
  {
    fJetEmbeddingArray = static_cast<TClonesArray*>(InputEvent()->FindListObject(Form("%s", fJetEmbeddingArrayName.Data())));
    if(!fJetEmbeddingArray)
      AliFatal(Form("Importing jets for embedding failed! Array '%s' not found!", fJetEmbeddingArrayName.Data()));
  }
  else if(fJetOutputMode==4 || fJetOutputMode==5)
    AliFatal(Form("fJetEmbeddingArrayName must be set in jet output mode 4 or 5."));

  // ### Import veto jets for matching (optional)
  if(fJetVetoArrayName != "")
  {
    fJetVetoArray = static_cast<TClonesArray*>(InputEvent()->FindListObject(Form("%s", fJetVetoArrayName.Data())));
    if(!fJetVetoArray)
      AliFatal(Form("Importing jets for veto failed! Array '%s' not found!", fJetVetoArrayName.Data()));
  }

  // ### Jets tree (optional)
  if(fExtractionPercentage)
  {
    fJetsTree = new TTree("ExtractedJets", "ExtractedJets");
    fJetsTree->Branch("Jets", "AliBasicJet", &fJetsTreeBuffer, 1000);
    fOutput->Add(fJetsTree);
  }


}

//________________________________________________________________________
Bool_t AliAnalysisTaskChargedJetsHadronCF::IsJetSelected(AliEmcalJet* jet)
{
  if(fJetOutputMode==3) // output leading&subleading jet
  {
    if((jet!=fLeadingJet) && (jet!=fSubleadingJet))
      return kFALSE;
  }
  else if(fJetOutputMode==1) // output the leading jet
  {
    if(jet!=fLeadingJet)
      return kFALSE;
  }
  else if(fJetOutputMode==2) // output the subleading jet
  {
    if(jet!=fSubleadingJet)
      return kFALSE;
  }
  else if(fJetOutputMode==6)
  {
    if(jet!=fInitialPartonMatchedJet1 && jet!=fInitialPartonMatchedJet2)
      return kFALSE;
  }

  if(fJetOutputMode==4) // matching jets only
    return (std::find(fMatchedJets.begin(), fMatchedJets.end(), jet) != fMatchedJets.end());
  else if(fJetOutputMode==5) // non-matching jets only
    return (std::find(fMatchedJets.begin(), fMatchedJets.end(), jet) == fMatchedJets.end());

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronCF::FillHistogramsJets(AliEmcalJet* jet, const char* cutName)
{
  TString appendix("");
  if(cutName)
    appendix = Form("_%s", cutName);

  // All jets
  FillHistogram(Form("hJetPtRaw%s", appendix.Data()), jet->Pt(), fCent); 
  FillHistogram(Form("hJetPt%s", appendix.Data()), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), fCent); 
  FillHistogram(Form("hJetPhi%s", appendix.Data()), jet->Phi(), fCent); 
  FillHistogram(Form("hJetEta%s", appendix.Data()), jet->Eta(), fCent); 
  FillHistogram(Form("hJetEtaPt%s", appendix.Data()), jet->Eta(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area()); 
  FillHistogram(Form("hJetPhiPt%s", appendix.Data()), jet->Phi(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area()); 
  FillHistogram(Form("hJetPhiEta%s", appendix.Data()), jet->Phi(), jet->Eta()); 
  FillHistogram(Form("hJetArea%s", appendix.Data()), jet->Area(), fCent); 
  FillHistogram(Form("hJetAreaPt%s", appendix.Data()), jet->Area(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area()); 
  FillHistogram(Form("hJetPtLeadingHadron%s", appendix.Data()), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), fJetsCont->GetLeadingHadronPt(jet));

  FillHistogram(Form("hBackgroundPtJetPt_Cent0_100%s", appendix.Data()),  fJetsCont->GetRhoVal(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area()); 
  if( (fCent >= 0) && (fCent < 10) )
    FillHistogram(Form("hBackgroundPtJetPt_Cent0_10%s", appendix.Data()),  fJetsCont->GetRhoVal(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area()); 

  // Leading jet plots
  if(jet==fLeadingJet)
  {
    FillHistogram(Form("hLeadingJetPtRaw%s", appendix.Data()), jet->Pt(), fCent); 
    FillHistogram(Form("hLeadingJetPt%s", appendix.Data()), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), fCent); 
    FillHistogram(Form("hLeadingJetPhi%s", appendix.Data()), jet->Phi(), fCent); 
    FillHistogram(Form("hLeadingJetEta%s", appendix.Data()), jet->Eta(), fCent); 
    FillHistogram(Form("hLeadingJetEtaPt%s", appendix.Data()), jet->Eta(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area()); 
    FillHistogram(Form("hLeadingJetPhiPt%s", appendix.Data()), jet->Phi(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area()); 
    FillHistogram(Form("hLeadingJetPhiEta%s", appendix.Data()), jet->Phi(), jet->Eta()); 
    FillHistogram(Form("hLeadingJetArea%s", appendix.Data()), jet->Area(), fCent); 
    FillHistogram(Form("hLeadingJetAreaPt%s", appendix.Data()), jet->Area(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area()); 
    FillHistogram(Form("hLeadingJetPtLeadingHadron%s", appendix.Data()), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), fJetsCont->GetLeadingHadronPt(jet));
  }

  // Subleading jet plot
  else if(jet==fSubleadingJet)
  {
    FillHistogram(Form("hSubleadingJetPtRaw%s", appendix.Data()), jet->Pt(), fCent); 
    FillHistogram(Form("hSubleadingJetPt%s", appendix.Data()), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), fCent); 
    FillHistogram(Form("hSubleadingJetPhi%s", appendix.Data()), jet->Phi(), fCent); 
    FillHistogram(Form("hSubleadingJetEta%s", appendix.Data()), jet->Eta(), fCent); 
    FillHistogram(Form("hSubleadingJetEtaPt%s", appendix.Data()), jet->Eta(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area()); 
    FillHistogram(Form("hSubleadingJetPhiPt%s", appendix.Data()), jet->Phi(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area()); 
    FillHistogram(Form("hSubleadingJetPhiEta%s", appendix.Data()), jet->Phi(), jet->Eta());
    FillHistogram(Form("hSubleadingJetArea%s", appendix.Data()), jet->Area(), fCent); 
    FillHistogram(Form("hSubleadingJetAreaPt%s", appendix.Data()), jet->Area(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area()); 
    FillHistogram(Form("hSubleadingJetPtLeadingHadron%s", appendix.Data()), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), fJetsCont->GetLeadingHadronPt(jet));
  }

  // ####### Jet constituent plots
  for(Int_t i = 0; i < jet->GetNumberOfTracks(); i++)
  {
    AliVParticle* constituent = static_cast<AliVParticle*>(jet->TrackAt(i, fTracksCont->GetArray()));
    if(!constituent) 
      continue;

    Bool_t filterConditionFulfilled = kFALSE;
    AliAODTrack* aodTrack = static_cast<AliAODTrack*>(constituent);
    if (aodTrack)
      filterConditionFulfilled = aodTrack->TestFilterBit(fConstPtFilterBit);

    // Fill jet constituent plots
    FillHistogram(Form("hJetConstituentPt_Cent0_100%s", appendix.Data()), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), constituent->Pt());
    if(filterConditionFulfilled)
      FillHistogram(Form("hJetConstituentPt_Cent0_100_FilterBit%i%s", fConstPtFilterBit, appendix.Data()), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), constituent->Pt());
    if( (fCent >= 0) && (fCent < 10) )
    {
      FillHistogram(Form("hJetConstituentPt_Cent0_10%s", appendix.Data()), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), constituent->Pt());
      if(filterConditionFulfilled)
        FillHistogram(Form("hJetConstituentPt_Cent0_10_FilterBit%i%s", fConstPtFilterBit, appendix.Data()), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), constituent->Pt());
    }
    else if( (fCent >= 10) && (fCent < 30) )
    {
      FillHistogram(Form("hJetConstituentPt_Cent10_30%s", appendix.Data()), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), constituent->Pt());
      if(filterConditionFulfilled)
        FillHistogram(Form("hJetConstituentPt_Cent10_30_FilterBit%i%s", fConstPtFilterBit, appendix.Data()), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), constituent->Pt());
    }
    else if( (fCent >= 30) && (fCent < 50) )
    {
      FillHistogram(Form("hJetConstituentPt_Cent30_50%s", appendix.Data()), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), constituent->Pt());
      if(filterConditionFulfilled)
        FillHistogram(Form("hJetConstituentPt_Cent30_50_FilterBit%i%s", fConstPtFilterBit, appendix.Data()), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), constituent->Pt());
    }
    else if( (fCent >= 50) && (fCent < 90) )
    {
      FillHistogram(Form("hJetConstituentPt_Cent50_90%s", appendix.Data()), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), constituent->Pt());
      if(filterConditionFulfilled)
        FillHistogram(Form("hJetConstituentPt_Cent50_90_FilterBit%i%s", fConstPtFilterBit, appendix.Data()), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), constituent->Pt());
    }
  }

  FillHistogram(Form("hJetConstituentCount_Cent0_100%s", appendix.Data()), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), jet->GetNumberOfTracks()); 
  if( (fCent >= 0) && (fCent < 10) )
    FillHistogram(Form("hJetConstituentCount_Cent0_10%s", appendix.Data()), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), jet->GetNumberOfTracks()); 

  FillHistogram(Form("hBackgroundPtConstCount_Cent0_100%s", appendix.Data()), fJetsCont->GetRhoVal(), jet->GetNumberOfTracks()); 
  if( (fCent >= 0) && (fCent < 10) )
    FillHistogram(Form("hBackgroundPtConstCount_Cent0_10%s", appendix.Data()), fJetsCont->GetRhoVal(), jet->GetNumberOfTracks()); 

  // ####### Embedding plots
  if( (fJetOutputMode == 4) || (fJetOutputMode == 5))
  {
    AliEmcalJet* refJet = GetReferenceJet(jet);
    Double_t deltaEta = (jet->Eta()-refJet->Eta());
    Double_t deltaPhi = TMath::Min(TMath::Abs(jet->Phi()-refJet->Phi()),TMath::TwoPi() - TMath::Abs(jet->Phi()-refJet->Phi()));
    if(jet->Phi() < refJet->Phi())
      deltaPhi = -deltaPhi;

    Double_t deltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
    FillHistogram(Form("hEmbeddingDeltaR%s", appendix.Data()), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), deltaR);
    FillHistogram(Form("hEmbeddingDeltaEta%s", appendix.Data()), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), deltaPhi);
    FillHistogram(Form("hEmbeddingDeltaPhi%s", appendix.Data()), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), deltaEta);
    FillHistogram(Form("hEmbeddingJetPt%s", appendix.Data()), refJet->Pt());
    FillHistogram(Form("hEmbeddingJetPhiEta%s", appendix.Data()), refJet->Phi(), refJet->Eta()); 

    // Only create 3D plots for each cut on demand
    Double_t trackRatio = 0.;
    Double_t ptRatio = 0.;
    GetTrackMCRatios(jet, refJet, trackRatio, ptRatio);

    if(fCent >= 0 && fCent < 10)
    {
      if((appendix == "") || fJetEmbeddingCreatePtPlotPerCut)
        FillHistogram3D(Form("hEmbeddingPtCorr010%s", appendix.Data()), refJet->Pt(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), ptRatio);
      if(ptRatio >= 0.2)
        FillHistogram(Form("hEmbeddingPtCorr010_Above20%s", appendix.Data()), refJet->Pt(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area());
    }
    else if (fCent >= 10 && fCent < 30)
    {
      if((appendix == "") || fJetEmbeddingCreatePtPlotPerCut)
        FillHistogram3D(Form("hEmbeddingPtCorr1030%s", appendix.Data()), refJet->Pt(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), ptRatio);
      if(ptRatio >= 0.2)
        FillHistogram(Form("hEmbeddingPtCorr1030_Above20%s", appendix.Data()), refJet->Pt(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area());
    }
    else if (fCent >= 30 && fCent < 50)
    {
      if((appendix == "") || fJetEmbeddingCreatePtPlotPerCut)
        FillHistogram3D(Form("hEmbeddingPtCorr3050%s", appendix.Data()), refJet->Pt(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), ptRatio);
      if(ptRatio >= 0.2)
        FillHistogram(Form("hEmbeddingPtCorr3050_Above20%s", appendix.Data()), refJet->Pt(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area());
    }
    else if (fCent >= 50 && fCent < 90)
    {
      if((appendix == "") || fJetEmbeddingCreatePtPlotPerCut)
        FillHistogram3D(Form("hEmbeddingPtCorr5090%s", appendix.Data()), refJet->Pt(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), ptRatio);
      if(ptRatio >= 0.2)
        FillHistogram(Form("hEmbeddingPtCorr5090_Above20%s", appendix.Data()), refJet->Pt(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area());
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronCF::FillHistogramsTracks(AliVTrack* track)
{
  FillHistogram("hTrackPt", track->Pt(), fCent); 
  FillHistogram("hTrackPhi", track->Phi(), fCent); 
  FillHistogram("hTrackEta", track->Eta(), fCent); 
  FillHistogram("hTrackEtaPt", track->Eta(), track->Pt()); 
  FillHistogram("hTrackPhiPt", track->Phi(), track->Pt()); 
  FillHistogram("hTrackPhiEta", track->Phi(), track->Eta()); 
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronCF::AddJetToOutputArray(AliEmcalJet* jet, Int_t arrayIndex, Int_t& jetsAlreadyInArray)
{
  Double_t tmpPt = jet->Pt() - fJetsCont->GetRhoVal()*jet->Area();
  new ((*fJetsOutput.at(arrayIndex))[jetsAlreadyInArray]) AliPicoTrack(tmpPt, jet->Eta(), jet->Phi(), jet->Charge(), 0, 0);
  jetsAlreadyInArray++;
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronCF::AddTrackToOutputArray(AliVTrack* track)
{
  if(fTrackParticleArrayName != "")
  {
    new ((*fTracksOutput)[fAcceptedTracks]) AliPicoTrack(track->Pt(), track->Eta(), track->Phi(), track->Charge(), 0, 0); // only Pt,Eta,Phi are interesting for correlations;
    fAcceptedTracks++;
  }
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronCF::AddJetToTree(AliEmcalJet* jet)
{
  // Check pT threshold
  if( ((jet->Pt()-jet->Area()*fJetsCont->GetRhoVal()) < fExtractionMinPt) || ((jet->Pt()-jet->Area()*fJetsCont->GetRhoVal()) >= fExtractionMaxPt) )
    return;

  AliVHeader* eventIDHeader = InputEvent()->GetHeader();
  Long64_t eventID = 0;
  if(eventIDHeader)
    eventID = eventIDHeader->GetEventIdAsLong();

  // if only the two initial collision partons will be added, get info on them
  Int_t matchedIC = 0;
  Int_t matchedHadron = 0;
  if(fJetOutputMode==6)
  {
    // Get type of jet
    CalculateJetType(jet, matchedIC, matchedHadron);
    Int_t partid = matchedIC;
    if (fPythiaExtractionUseHadronMatching)
      partid = matchedHadron;
    // If fPythiaExtractionMode is set, only extract certain jets
    if( (fPythiaExtractionMode==1) && not (partid>=1 && partid<=6)) // all quark-jet extraction
      return;
    else if( (fPythiaExtractionMode==2) && not (partid==21)) // gluon-jet extraction
      return;
    else if( (fPythiaExtractionMode==3) && not (partid==0)) // extract only those w/o hadron matching
      return;
    else if( (fPythiaExtractionMode<0) && (fPythiaExtractionMode!=-partid) ) // custom type jet extraction by given a negative number
      return;
  }

  // Load vertex if possible
  Double_t vtxX = 0;
  Double_t vtxY = 0;
  Double_t vtxZ = 0;
  const AliVVertex* myVertex = InputEvent()->GetPrimaryVertex();
  if(!myVertex && MCEvent())
    myVertex = MCEvent()->GetPrimaryVertex();
  if(myVertex)
  {
    vtxX = myVertex->GetX();
    vtxY = myVertex->GetY();
    vtxZ = myVertex->GetZ();
  }

  // Discard jets statistically
  if(fRandom->Rndm() >= fExtractionPercentage)
    return;

  AliBasicJet basicJet(jet->Eta(), jet->Phi(), jet->Pt(), jet->Charge(), fJetsCont->GetJetRadius(), jet->Area(), matchedIC, matchedHadron, fJetsCont->GetRhoVal(), InputEvent()->GetMagneticField(), vtxX, vtxY, vtxZ, eventID, fCent);

  // Add constituents
  for(Int_t i = 0; i < jet->GetNumberOfTracks(); i++)
  {
    AliVParticle* particle = static_cast<AliVParticle*>(jet->TrackAt(i, fTracksCont->GetArray()));
    if(!particle) continue;

    AliAODTrack*  aodtrack = static_cast<AliAODTrack*>(jet->TrackAt(i, fTracksCont->GetArray()));
    Int_t constid = 9; // 9 mean unknown
    if(fJetOutputMode==6) // MC
    {
      // Use same convention as PID in AODs
      if(TMath::Abs(particle->PdgCode()) == 2212) // proton
        constid = 4;
      else if (TMath::Abs(particle->PdgCode()) == 211) // pion
        constid = 2;
      else if (TMath::Abs(particle->PdgCode()) == 321) // kaon
        constid = 3;
      else if (TMath::Abs(particle->PdgCode()) == 11) // electron
        constid = 0;
      else if (TMath::Abs(particle->PdgCode()) == 13) // muon
        constid = 1;
    }
    else if (aodtrack) // data
      constid = aodtrack->GetMostProbablePID();

    basicJet.AddJetConstituent(particle->Eta(), particle->Phi(), particle->Pt(), particle->Charge(), constid, particle->Xv(), particle->Yv(), particle->Zv(), 0);
  }
  if(std::find(fMatchedJets.begin(), fMatchedJets.end(), jet) != fMatchedJets.end()) // set the true pT from the matched jets (only possible in modes 4 & 7)
    basicJet.SetTruePt(fMatchedJetsReference[std::find(fMatchedJets.begin(), fMatchedJets.end(), jet)-fMatchedJets.begin()]->Pt());

  fJetsTreeBuffer = &basicJet;
  fJetsTree->Fill();
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronCF::AddEventToTree()
{
  // Check jet pT threshold
  if(fLeadingJet && ( ((fLeadingJet->Pt() - fLeadingJet->Area()*fJetsCont->GetRhoVal()) < fEventExtractionMinJetPt) || ((fLeadingJet->Pt() - fLeadingJet->Area()*fJetsCont->GetRhoVal()) >= fEventExtractionMaxJetPt)))
    return;

  // Discard jets statistically
  if(fRandom->Rndm() >= fEventExtractionPercentage)
    return;

  static Int_t numSavedEvents = 0;
  numSavedEvents++;


  AddHistogram2D<TH2D>(Form("Event%i", numSavedEvents), "Event display", "COLZ", 180, 0., 2*TMath::Pi(), 100, -2.5, 2.5, "#phi", "#eta", "dN^{Tracks}/d#phi d#eta");
  fTracksCont->ResetCurrentID();
  while(AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle()))
    FillHistogram(Form("Event%i", numSavedEvents), track->Phi(), track->Eta(), track->Pt());

}

//________________________________________________________________________
Bool_t AliAnalysisTaskChargedJetsHadronCF::Run()
{
  CalculateEventProperties();

  // Jet veto:
  // 1. jet-by-jet mode: veto is active if jets overlaps with a suitable jet
  // 2. other mode: veto is active if suitable jet is in sample
  AliEmcalJet* vetoJet = 0;
  if(!fJetVetoJetByJet)
    vetoJet = GetLeadingVetoJet();

  // ####### Jet loop
  fAcceptedJets = 0;
  for(Int_t i = 0; i<fJetEmbeddingCuts.size(); i++)
    fJetEmbeddingCuts.at(i).fAcceptedJets = 0;
  fJetsCont->ResetCurrentID();
  while(AliEmcalJet *jet = fJetsCont->GetNextAcceptJet())
  {
    if(!IsJetSelected(jet))
      continue;

    // Plots + output jets regardless of embedding cut
    FillHistogramsJets(jet, 0);
    AddJetToOutputArray(jet, 0, fAcceptedJets);

    // Plots + output jets for each embedding cut
    if(fJetVetoJetByJet)
      vetoJet = GetVetoJet(jet);
    for(Int_t i = 0; i<fJetEmbeddingCuts.size(); i++)
    {
      AliChargedJetsHadronCFCuts currentCut = fJetEmbeddingCuts.at(i);
      AliEmcalJet* refJet = GetReferenceJet(jet);
      Double_t trackRatio = 0.;
      Double_t ptRatio = 0.;
      GetTrackMCRatios(jet, refJet, trackRatio, ptRatio);

      Double_t vetoJetPt = 0.;
      if(vetoJet)
        vetoJetPt = vetoJet->Pt();

      if(!currentCut.IsCutFulfilled(jet->Pt(), refJet->Pt(), fCent, ptRatio, vetoJetPt))
        continue;

      FillHistogramsJets(jet, currentCut.fCutName.Data());
      AddJetToOutputArray(jet, currentCut.fArrayIndex, currentCut.fAcceptedJets);
    }

    // Add jet to output array
    if(fExtractionPercentage)
      AddJetToTree(jet);
  }

  // ####### Particle loop
  // Throw random cone
  Double_t tmpRandConeEta = fJetsCont->GetJetEtaMin() + fRandom->Rndm()*TMath::Abs(fJetsCont->GetJetEtaMax()-fJetsCont->GetJetEtaMin());
  Double_t tmpRandConePhi = fRandom->Rndm()*TMath::TwoPi();
  Double_t tmpRandConePt  = 0; // to be determined
  Double_t tmpRandConePt3GeV = 0; // to be determined

  fAcceptedTracks = 0;
  fTracksCont->ResetCurrentID();
  Int_t trackcount = 0;
  while(AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle()))
  {
    // Track plots
    FillHistogramsTracks(track);

    if(IsTrackInCone(track, tmpRandConeEta, tmpRandConePhi, fJetsCont->GetJetRadius()))
    {
      tmpRandConePt += track->Pt();
      if (track->Pt() > 3.0)
        tmpRandConePt3GeV += track->Pt();
    }

    // Add track to output array
    trackcount++;
    AddTrackToOutputArray(track);
  }

  // Add event to output tree
  if(fEventExtractionPercentage)
    AddEventToTree();

  // ####### Event properties
  FillHistogram("hRandomConePt", tmpRandConePt - fJetsCont->GetRhoVal()*fJetsCont->GetJetRadius()*fJetsCont->GetJetRadius()*TMath::Pi(), fCent);
  FillHistogram("hRandomConePtCut3GeV", tmpRandConePt3GeV - fJetsCont->GetRhoVal()*fJetsCont->GetJetRadius()*fJetsCont->GetJetRadius()*TMath::Pi(), fCent);
  FillHistogram("hRandomConeRawPt", tmpRandConePt, fCent); 
  FillHistogram("hRandomConeRawPtCut3GeV", tmpRandConePt3GeV, fCent);

  FillHistogram("hBackgroundPt", fJetsCont->GetRhoVal(), fCent);
  FillHistogram("hJetCount", fAcceptedJets, fCent);
  FillHistogram("hTrackCount", trackcount, fCent);
  // NOTE: It is possible to use fTracksCont->GetLeadingParticle() since we do not apply additional track cuts
  AliVTrack* leadTrack = static_cast<AliVTrack*>(fTracksCont->GetLeadingParticle());
  if(leadTrack)
  {
    FillHistogram("hLeadingTrackPt", leadTrack->Pt(), fCent); 
    FillHistogram("hLeadingTrackPhi", leadTrack->Phi(), fCent); 
    FillHistogram("hLeadingTrackEta", leadTrack->Eta(), fCent); 
    FillHistogram("hLeadingTrackPhiEta", leadTrack->Phi(), leadTrack->Eta()); 
  }

  return kTRUE;
}

//########################################################################
// HELPERS
//########################################################################

//________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronCF::GetInitialCollisionJets()
{
  if(!fPythiaInfo)
    AliError("fPythiaInfo object not available. Is it activated with SetGeneratePythiaInfoObject()?");

  Double_t bestMatchDeltaR1 = 999.;
  Double_t bestMatchDeltaR2 = 999.;

  fJetsCont->ResetCurrentID();
  while(AliEmcalJet *jet = fJetsCont->GetNextAcceptJet())
  {
    // Check via geometrical matching if jet is connected to the initial collision
    Double_t deltaEta1 = TMath::Abs(jet->Eta()-fPythiaInfo->GetPartonEta6());
    Double_t deltaEta2 = TMath::Abs(jet->Eta()-fPythiaInfo->GetPartonEta7());
    Double_t deltaPhi1 = TMath::Min(TMath::Abs(jet->Phi()-fPythiaInfo->GetPartonPhi6()),TMath::TwoPi() - TMath::Abs(jet->Phi()-fPythiaInfo->GetPartonPhi6()));
    Double_t deltaPhi2 = TMath::Min(TMath::Abs(jet->Phi()-fPythiaInfo->GetPartonPhi7()),TMath::TwoPi() - TMath::Abs(jet->Phi()-fPythiaInfo->GetPartonPhi7()));

    Double_t deltaR1 = TMath::Sqrt(deltaEta1*deltaEta1 + deltaPhi1*deltaPhi1);
    Double_t deltaR2 = TMath::Sqrt(deltaEta2*deltaEta2 + deltaPhi2*deltaPhi2);

    if(deltaR1 < bestMatchDeltaR1)
    {
      bestMatchDeltaR1 = deltaR1;
      fInitialPartonMatchedJet1 = jet;
    }
    if(deltaR2 < bestMatchDeltaR2)
    {
      bestMatchDeltaR2 = deltaR2;
      fInitialPartonMatchedJet2 = jet;
    }
  }

  if(bestMatchDeltaR1 > 0.3)
    fInitialPartonMatchedJet1 = 0;
  if(bestMatchDeltaR2 > 0.3)
    fInitialPartonMatchedJet2 = 0;
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronCF::GetMatchingJets()
{
  fMatchedJets.clear();
  fMatchedJetsReference.clear();

  AliEmcalJet* jetLeading = 0;
  AliEmcalJet* jetSubLeading = 0;
  GetLeadingJetsInArray(fJetEmbeddingArray, "rho", jetLeading, jetSubLeading);

  // ############ Search for all matches with leading/subleading jets
  for(Int_t i=0; i<fJetEmbeddingNumMatchedJets; i++)
  {
    AliEmcalJet* probeJet = 0;
    if(i==0)
      probeJet = jetLeading;
    else if(i==1)
      probeJet = jetSubLeading;

    if(!probeJet)
      continue;

    AliEmcalJet* matchedJet = 0;
    AliEmcalJet* matchedJetReference = 0;
    Double_t     bestMatchDeltaR = 999.;
    fJetsCont->ResetCurrentID();
    // Loop over all embedded jets to find the best match
    while(AliEmcalJet* embeddedJet = fJetsCont->GetNextAcceptJet())
    {
      Double_t deltaEta = (embeddedJet->Eta()-probeJet->Eta());
      Double_t deltaPhi = TMath::Min(TMath::Abs(embeddedJet->Phi()-probeJet->Phi()),TMath::TwoPi() - TMath::Abs(embeddedJet->Phi()-probeJet->Phi()));
      Double_t deltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);

      // Cut jets too far away
      if (deltaR > fJetEmbeddingMaxDistance)
        continue;

      // Search for the best match
      if(deltaR < bestMatchDeltaR)
      {
        bestMatchDeltaR = deltaR;
        matchedJet = embeddedJet;
        matchedJetReference = probeJet;
      }
    }
    // Put matched jet to a list
    if(matchedJet && matchedJetReference)
    {
      fMatchedJets.push_back(matchedJet);
      fMatchedJetsReference.push_back(matchedJetReference);
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronCF::CalculateJetType(AliEmcalJet* jet, Int_t& typeIC, Int_t& typeHM)
{
  typeIC = 0;
  if(!fPythiaInfo)
    AliError("fPythiaInfo object not available. Is it activated with SetGeneratePythiaInfoObject()?");
  else if(jet==fInitialPartonMatchedJet1)
    typeIC = fPythiaInfo->GetPartonFlag6();
  else if (jet==fInitialPartonMatchedJet2)
    typeIC = fPythiaInfo->GetPartonFlag7();


  typeHM = 0;
  AliStack* stack = MCEvent()->Stack();
  // Go through the whole particle stack
  for(Int_t i=0; i<stack->GetNtrack(); i++)
  {
    TParticle *part = stack->Particle(i);
    if(!part) continue;

    // Check if particle is in a radius around the jet
    Double_t rsquared = (part->Eta() - jet->Eta())*(part->Eta() - jet->Eta()) + (part->Phi() - jet->Phi())*(part->Phi() - jet->Phi());
    if(rsquared >= fHadronMatchingRadius*fHadronMatchingRadius)
    {
      continue;
    }

    // Check if the particle has beauty, charm or strangeness
    // If it has beauty, we are done (exclusive definition)
    Int_t absPDG = TMath::Abs(part->GetPdgCode());
    // Particle has beauty
    if ((absPDG > 500 && absPDG < 600) || (absPDG > 5000 && absPDG < 6000))
    {
      typeHM = 5; // beauty
      break;
    }
    // Particle has charm
    else if ((absPDG > 400 && absPDG < 500) || (absPDG > 4000 && absPDG < 5000))
      typeHM = 4; // charm
    // Particle has strangeness: Only search for strangeness, if charm was not already found
    else if (typeHM != 4 && (absPDG > 300 && absPDG < 400) || (absPDG > 3000 && absPDG < 4000))
      typeHM = 3; // strange
  }
}


//________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronCF::GetTrackMCRatios(AliEmcalJet* jet, AliEmcalJet* mcJet, Double_t& trackRatio, Double_t& ptRatio)
{
  Int_t tracksFromMC = 0;
  Int_t tracksTotal  = 0;
  Double_t ptFromMC  = 0;
  Double_t ptTotal   = 0;

  if(fJetEmbeddingUsePerTrackMCPercentage) // Calculate MC track percentage from tracks of the matched jet
  {
    for(Int_t j = 0; j < jet->GetNumberOfTracks(); j++)
    {
      AliVParticle* constituent = static_cast<AliVParticle*>(jet->TrackAt(j, fTracksCont->GetArray()));
      if(!constituent) 
        continue;

      TClonesArray* mcArray = static_cast<TClonesArray*>(InputEvent()->FindListObject(Form("%s", fJetEmbeddingTrackArrayName.Data())));
      Bool_t foundInMC = kFALSE;
      for(Int_t k=0; k<mcJet->GetNumberOfTracks(); k++)
      {
        AliVParticle* mcConstituent = static_cast<AliVParticle*>(mcJet->TrackAt(k, mcArray));
        if(!mcConstituent) 
          continue;
        if(mcConstituent->GetLabel() == constituent->GetLabel())
        {
          foundInMC = kTRUE;
          break;
        }
      }

      if(foundInMC)
      {
        tracksFromMC++;
        ptFromMC += constituent->Pt();
      }
      tracksTotal++;
      ptTotal += constituent->Pt();

    }
  }
  else // Calculate MC track percentage from all tracks in MC
  {
    for(Int_t j = 0; j < jet->GetNumberOfTracks(); j++)
    {
      AliVParticle* constituent = static_cast<AliVParticle*>(jet->TrackAt(j, fTracksCont->GetArray()));
      if(!constituent) 
        continue;

      // Plots on MC percentage in jet
      if(constituent->GetLabel() > 10000)
      {
        tracksFromMC++;
        ptFromMC += constituent->Pt();
      }
      tracksTotal++;
      ptTotal += constituent->Pt();
    }
  }


  trackRatio = 0.;
  if(tracksTotal)
    trackRatio = tracksFromMC/((Double_t)tracksTotal);

  if(fJetEmbeddingUseBgrdForMCPercentage)
    ptTotal = jet->Pt() - fJetsCont->GetRhoVal()*jet->Area();

  ptRatio = 0.;
  if(ptTotal)
    ptRatio = ptFromMC/ptTotal;

}

//________________________________________________________________________
AliEmcalJet* AliAnalysisTaskChargedJetsHadronCF::GetVetoJet(AliEmcalJet* jet)
{
  Double_t     leadingVetoJetPt   = -999.;
  AliEmcalJet* leadingVetoJet     = 0;
  // Search for the 'leading' overlapping jet in the veto sample
  if(fJetVetoArray && jet)
  {
    for(Int_t j=0; j<fJetVetoArray->GetEntries(); j++)
    {
      UInt_t dummy = 0;
      AliEmcalJet* vetoJet = static_cast<AliEmcalJet*>(fJetVetoArray->At(j));
      // Check if veto jet is accepted
      if(!fJetsCont->AcceptJet(vetoJet , dummy))
        continue;

      // Check matching distance
      Double_t vetoPt  = vetoJet->Pt() - vetoJet->Area()*fJetsCont->GetRhoVal();
      Double_t deltaEta = (jet->Eta()-vetoJet->Eta());
      Double_t deltaPhi = TMath::Min(TMath::Abs(jet->Phi()-vetoJet->Phi()),TMath::TwoPi() - TMath::Abs(jet->Phi()-vetoJet->Phi()));
      Double_t deltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);

      if ((vetoPt > leadingVetoJetPt) && (deltaR <= fJetEmbeddingMaxDistance))
      {
        leadingVetoJetPt = vetoPt;
        leadingVetoJet = vetoJet;
      }
    }
  }

  return leadingVetoJet;
}

//________________________________________________________________________
AliEmcalJet* AliAnalysisTaskChargedJetsHadronCF::GetLeadingVetoJet()
{
  Double_t     leadingVetoJetPt   = -999.;
  AliEmcalJet* leadingVetoJet     = 0;
  // Search for the 'leading' jet in the veto sample
  if(fJetVetoArray)
  {
    for(Int_t j=0; j<fJetVetoArray->GetEntries(); j++)
    {
      UInt_t dummy = 0;
      AliEmcalJet* vetoJet = static_cast<AliEmcalJet*>(fJetVetoArray->At(j));
      // Check if veto jet is accepted
      if(!fJetsCont->AcceptJet(vetoJet , dummy))
        continue;

      // Check matching distance
      Double_t vetoPt  = vetoJet->Pt() - vetoJet->Area()*fJetsCont->GetRhoVal();
      if (vetoPt > leadingVetoJetPt)
      {
        leadingVetoJetPt = vetoPt;
        leadingVetoJet = vetoJet;
      }
    }
  }
  return leadingVetoJet;
}

//________________________________________________________________________
inline Bool_t AliAnalysisTaskChargedJetsHadronCF::IsTrackInCone(AliVParticle* track, Double_t eta, Double_t phi, Double_t radius)
{
  // This is to use a full cone in phi even at the edges of phi (2pi -> 0) (0 -> 2pi)
  Double_t trackPhi = 0.0;
  if (track->Phi() > (TMath::TwoPi() - (radius-phi)))
    trackPhi = track->Phi() - TMath::TwoPi();
  else if (track->Phi() < (phi+radius - TMath::TwoPi()))
    trackPhi = track->Phi() + TMath::TwoPi();
  else
    trackPhi = track->Phi();
  
  if ( TMath::Abs(trackPhi-phi)*TMath::Abs(trackPhi-phi) + TMath::Abs(track->Eta()-eta)*TMath::Abs(track->Eta()-eta) <= radius*radius)
    return kTRUE;
  
  return kFALSE;
}


//________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronCF::CalculateEventProperties()
{
  // Calculate leading + subleading jet
  GetLeadingJets("rho", fLeadingJet, fSubleadingJet);
  if(fJetOutputMode==6)
    GetInitialCollisionJets();
  else if(fJetOutputMode==4 || fJetOutputMode==5)
    GetMatchingJets();
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronCF::GetLeadingJets(const char* opt, AliEmcalJet*& jetLeading, AliEmcalJet*& jetSubLeading)
{
  // Customized from AliJetContainer::GetLeadingJet()
  // Get the leading+subleading jet; if opt contains "rho" the sorting is according to pt-A*rho

  TString option(opt);
  option.ToLower();

  jetLeading = 0;
  jetSubLeading = 0;

  fJetsCont->ResetCurrentID();
  Double_t     tmpLeadingPt = 0;
  Double_t     tmpSubleadingPt = 0;

  if (option.Contains("rho")) {
    while (AliEmcalJet* jet = fJetsCont->GetNextAcceptJet()) {
      if      ( (jet->Pt()-jet->Area()*fJetsCont->GetRhoVal()) > tmpLeadingPt )
      {
        jetSubLeading = jetLeading;
        jetLeading = jet;
        tmpSubleadingPt = tmpLeadingPt;
        tmpLeadingPt = jet->Pt()-jet->Area()*fJetsCont->GetRhoVal();
      }
      else if ( (jet->Pt()-jet->Area()*fJetsCont->GetRhoVal()) > tmpSubleadingPt )
      {
        jetSubLeading = jet;
        tmpSubleadingPt = jet->Pt()-jet->Area()*fJetsCont->GetRhoVal();
      }
    }
  }
  else {
    while (AliEmcalJet* jet = fJetsCont->GetNextAcceptJet()) {
      if      ( (jet->Pt()) > tmpLeadingPt )
      {
        jetSubLeading = jetLeading;
        jetLeading = jet;
        tmpSubleadingPt = tmpLeadingPt;
        tmpLeadingPt = jet->Pt();
      }
      else if ( (jet->Pt()) > tmpSubleadingPt )
      {
        jetSubLeading = jet;
        tmpSubleadingPt = jet->Pt();
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronCF::GetLeadingJetsInArray(TClonesArray* arr, const char* opt, AliEmcalJet*& jetLeading, AliEmcalJet*& jetSubLeading)
{
  // Customized from AliJetContainer::GetLeadingJet()
  // Get the leading+subleading jet; if opt contains "rho" the sorting is according to pt-A*rho

  TString option(opt);
  option.ToLower();

  jetLeading = 0;
  jetSubLeading = 0;
  Double_t     tmpLeadingPt = -999.;
  Double_t     tmpSubleadingPt = -999.;

  for(Int_t i=0; i<arr->GetEntries(); i++)
  {
    AliEmcalJet* jet = static_cast<AliEmcalJet*>(arr->At(i));
    UInt_t   dummy = 0;
    if(!fJetsCont->AcceptJet(jet , dummy))
      continue;

    Double_t jetPt = jet->Pt();
    if (option.Contains("rho"))
      jetPt -= jet->Area()*fJetsCont->GetRhoVal();
    
    if      (  jetPt > tmpLeadingPt )
    {
      jetSubLeading = jetLeading;
      jetLeading = jet;
      tmpSubleadingPt = tmpLeadingPt;
      tmpLeadingPt = jetPt;
    }
    else if ( jetPt > tmpSubleadingPt )
    {
      jetSubLeading = jet;
      tmpSubleadingPt = jetPt;
    }

  }
}


//________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronCF::BinLogAxis(const THn *h, Int_t axisNumber)
{
  // Method for the correct logarithmic binning of histograms
  TAxis *axis = h->GetAxis(axisNumber);
  int bins = axis->GetNbins();

  Double_t from = axis->GetXmin();
  Double_t to = axis->GetXmax();
  Double_t *newBins = new Double_t[bins + 1];
   
  newBins[0] = from;
  Double_t factor = pow(to/from, 1./bins);
  
  for (int i = 1; i <= bins; i++) {
   newBins[i] = factor * newBins[i-1];
  }
  axis->Set(bins, newBins);
  delete [] newBins;
}

//________________________________________________________________________
AliEmcalJet* AliAnalysisTaskChargedJetsHadronCF::GetReferenceJet(AliEmcalJet* jet)
{
  std::vector<AliEmcalJet*>::iterator matchedJetFindResult = std::find(fMatchedJets.begin(), fMatchedJets.end(), jet);
  if(matchedJetFindResult == fMatchedJets.end())
  {
    AliError("Checked for a reference jet but it was not found. Check code.");
    return 0;
  }

  Int_t matchedJetIndex = (matchedJetFindResult - fMatchedJets.begin());
  AliEmcalJet* refJet = fMatchedJetsReference[matchedJetIndex];
  return refJet;
}

//________________________________________________________________________
inline void AliAnalysisTaskChargedJetsHadronCF::FillHistogram(const char * key, Double_t x)
{
  TH1* tmpHist = static_cast<TH1*>(fOutput->FindObject(key));
  if(!tmpHist)
  {
    AliError(Form("Cannot find histogram <%s> ",key)) ;
    return;
  }

  tmpHist->Fill(x);
}

//________________________________________________________________________
inline void AliAnalysisTaskChargedJetsHadronCF::FillHistogram(const char * key, Double_t x, Double_t y)
{
  TH1* tmpHist = static_cast<TH1*>(fOutput->FindObject(key));
  if(!tmpHist)
  {
    AliError(Form("Cannot find histogram <%s> ",key));
    return;
  }

  if (tmpHist->IsA()->GetBaseClass("TH1"))
    static_cast<TH1*>(tmpHist)->Fill(x,y); // Fill x with y
  else if (tmpHist->IsA()->GetBaseClass("TH2"))
    static_cast<TH2*>(tmpHist)->Fill(x,y); // Fill x,y with 1
}

//________________________________________________________________________
inline void AliAnalysisTaskChargedJetsHadronCF::FillHistogram(const char * key, Double_t x, Double_t y, Double_t add)
{
  TH2* tmpHist = static_cast<TH2*>(fOutput->FindObject(key));
  if(!tmpHist)
  {
    AliError(Form("Cannot find histogram <%s> ",key));
    return;
  }
  
  tmpHist->Fill(x,y,add);
}

//________________________________________________________________________
inline void AliAnalysisTaskChargedJetsHadronCF::FillHistogram3D(const char * key, Double_t x, Double_t y, Double_t z, Double_t add)
{
  TH3* tmpHist = static_cast<TH3*>(fOutput->FindObject(key));
  if(!tmpHist)
  {
    AliError(Form("Cannot find histogram <%s> ",key));
    return;
  }
  
  if(add)
    tmpHist->Fill(x,y,z,add);
  else
    tmpHist->Fill(x,y,z);
}


//________________________________________________________________________
template <class T> T* AliAnalysisTaskChargedJetsHadronCF::AddHistogram1D(const char* name, const char* title, const char* options, Int_t xBins, Double_t xMin, Double_t xMax, const char* xTitle, const char* yTitle)
{
  T* tmpHist = new T(name, title, xBins, xMin, xMax);

  tmpHist->GetXaxis()->SetTitle(xTitle);
  tmpHist->GetYaxis()->SetTitle(yTitle);
  tmpHist->SetOption(options);
  tmpHist->SetMarkerStyle(kFullCircle);
  tmpHist->Sumw2();

  fOutput->Add(tmpHist);

  return tmpHist;
}

//________________________________________________________________________
template <class T> T* AliAnalysisTaskChargedJetsHadronCF::AddHistogram2D(const char* name, const char* title, const char* options, Int_t xBins, Double_t xMin, Double_t xMax, Int_t yBins, Double_t yMin, Double_t yMax, const char* xTitle, const char* yTitle, const char* zTitle)
{
  T* tmpHist = new T(name, title, xBins, xMin, xMax, yBins, yMin, yMax);
  tmpHist->GetXaxis()->SetTitle(xTitle);
  tmpHist->GetYaxis()->SetTitle(yTitle);
  tmpHist->GetZaxis()->SetTitle(zTitle);
  tmpHist->SetOption(options);
  tmpHist->SetMarkerStyle(kFullCircle);
  tmpHist->Sumw2();

  fOutput->Add(tmpHist);

  return tmpHist;
}

//________________________________________________________________________
template <class T> T* AliAnalysisTaskChargedJetsHadronCF::AddHistogram3D(const char* name, const char* title, const char* options, Int_t xBins, Double_t xMin, Double_t xMax, Int_t yBins, Double_t yMin, Double_t yMax, Int_t zBins, Double_t zMin, Double_t zMax, const char* xTitle, const char* yTitle, const char* zTitle)
{
  T* tmpHist = new T(name, title, xBins, xMin, xMax, yBins, yMin, yMax, zBins, zMin, zMax);
  tmpHist->GetXaxis()->SetTitle(xTitle);
  tmpHist->GetYaxis()->SetTitle(yTitle);
  tmpHist->GetZaxis()->SetTitle(zTitle);
  tmpHist->SetOption(options);
  tmpHist->SetMarkerStyle(kFullCircle);
  tmpHist->Sumw2();

  fOutput->Add(tmpHist);

  return tmpHist;
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronCF::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

