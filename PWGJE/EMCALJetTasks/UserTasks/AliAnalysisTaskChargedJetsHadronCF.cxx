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

#include <TClonesArray.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THn.h>
#include <TTree.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliTrackContainer.h"
#include "AliPicoTrack.h"
#include "AliVParticle.h"
#include "TRandom3.h"
#include "AliAnalysisTaskEmcalJet.h"

#include "AliAnalysisTaskChargedJetsHadronCF.h"

/// \cond CLASSIMP
ClassImp(AliBasicJet)
/// \endcond
//________________________________________________________________________
AliBasicJet::~AliBasicJet() 
{
// dummy destructor
}

/// \cond CLASSIMP
ClassImp(AliBasicJetConstituent)
/// \endcond
//________________________________________________________________________
AliBasicJetConstituent::~AliBasicJetConstituent() 
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
  fNumberOfCentralityBins(10),
  fJetsOutput(),
  fTracksOutput(),
  fJetsInput(),
  fJetParticleArrayName("JetsDPhiBasicParticles"),
  fTrackParticleArrayName(""),
  fJetMatchingArrayName(""),
  fRandom(0),
  fRejectionFunction(0),
  fJetOutputMode(0),
  fMinFakeFactorPercentage(0),
  fMaxFakeFactorPercentage(0),
  fEventCriteriumMode(0),
  fEventCriteriumMinBackground(0),
  fEventCriteriumMaxBackground(0),
  fEventCriteriumMinLeadingJetPt(0),
  fEventCriteriumMinSubleadingJetPt(0),
  fLeadingJet(),
  fSubleadingJet(),
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
  fNumberOfCentralityBins(10),
  fJetsOutput(),
  fTracksOutput(),
  fJetsInput(),
  fJetParticleArrayName("JetsDPhiBasicParticles"),
  fTrackParticleArrayName(""),
  fJetMatchingArrayName(""),
  fRandom(0),
  fRejectionFunction(0),
  fJetOutputMode(0),
  fMinFakeFactorPercentage(0),
  fMaxFakeFactorPercentage(0),
  fEventCriteriumMode(0),
  fEventCriteriumMinBackground(0),
  fEventCriteriumMaxBackground(0),
  fEventCriteriumMinLeadingJetPt(0),
  fEventCriteriumMinSubleadingJetPt(0),
  fLeadingJet(),
  fSubleadingJet(),
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

  // Change the event rejection histogram -> Add a custom value
  fHistEventRejection->GetXaxis()->SetBinLabel(14,"JetCrit");

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


  // Jet QA plots
  AddHistogram2D<TH2D>("hJetPtRaw", "Jets p_{T} distribution (no bgrd. corr.)", "", 300, 0., 300., fNumberOfCentralityBins, 0, 100, "p_{T, jet} (GeV/c)", "Centrality", "dN^{Jets}/dp_{T}");
  AddHistogram2D<TH2D>("hJetPt", "Jets p_{T} distribution (background subtracted)", "", 400, -100., 300., fNumberOfCentralityBins, 0, 100, "p_{T, jet} (GeV/c)", "Centrality", "dN^{Jets}/dp_{T}");
  AddHistogram2D<TH2D>("hJetPhi", "Jet angular distribution #phi", "LEGO2", 180, 0., 2*TMath::Pi(), fNumberOfCentralityBins, 0, 100, "#phi", "Centrality", "dN^{Jets}/d#phi");
  AddHistogram2D<TH2D>("hJetEta", "Jet angular distribution #eta", "LEGO2", 100, -2.5, 2.5, fNumberOfCentralityBins, 0, 100, "#eta","Centrality","dN^{Jets}/d#eta");
  AddHistogram2D<TH2D>("hJetPhiPt", "Jet angular distribution #phi vs. p_{T}", "LEGO2", 180, 0., 2*TMath::Pi(), 400, -100., 300., "#phi", "p_{T, jet} (GeV/c)", "dN^{Jets}/d#phi dp_{T}");
  AddHistogram2D<TH2D>("hJetEtaPt", "Jet angular distribution #eta  vs. p_{T}", "LEGO2", 100, -2.5, 2.5, 400, -100., 300., "#eta","p_{T, jet} (GeV/c)","dN^{Jets}/d#eta dp_{T}");
  AddHistogram2D<TH2D>("hJetPhiEta", "Jet angular distribution #phi/#eta", "COLZ", 180, 0., 2*TMath::Pi(), 100, -2.5, 2.5, "#phi", "#eta", "dN^{Jets}/d#phi d#eta");
  AddHistogram2D<TH2D>("hJetArea", "Jet area", "LEGO2", 200, 0., 2., fNumberOfCentralityBins, 0, 100, "Jet A", "Centrality", "dN^{Jets}/dA");
  AddHistogram2D<TH2D>("hJetAreaPt", "Jet area vs. p_{T}", "LEGO2", 200, 0., 2., 400, -100., 300., "Jet A", "p_{T, jet} (GeV/c)", "dN^{Jets}/dA dp_{T}");
  AddHistogram2D<TH2D>("hJetPtLeadingHadron", "Jet leading hadron p_{T} distribution vs. jet p_{T}", "", 300, 0., 300., 300, 0., 300., "p_{T, jet} (GeV/c)", "p_{T,lead had} (GeV/c)", "dN^{Jets}/dp_{T}dp_{T,had}");

  AddHistogram2D<TH2D>("hJetConstituentPt_Cent0_100", "Jet constituent p_{T} distribution vs. jet p_T (background subtracted)", "", 400, -100., 300., 300, 0., 300., "p_{T, jet} (GeV/c)", "p_{T, track} (GeV/c)", "dN^{Tracks}/d^{2}p_{T}");
  AddHistogram2D<TH2D>("hJetConstituentPt_Cent0_10", "Jet constituent p_{T} distribution vs. jet p_T (background subtracted), 0-10 centrality", "", 400, -100., 300., 300, 0., 300., "p_{T, jet} (GeV/c)", "p_{T, track} (GeV/c)", "dN^{Tracks}/d^{2}p_{T}");

  AddHistogram2D<TH2D>("hJetConstituentCount_Cent0_100", "Jet constituent count vs. jet p_T (background subtracted)", "", 400, -100., 300., 200, 0., 200., "p_{T, jet} (GeV/c)", "Count", "dN^{Jets}/dNdp_{T}");
  AddHistogram2D<TH2D>("hJetConstituentCount_Cent0_10", "Jet constituent count vs. jet p_T (background subtracted), 0-10 centrality", "", 400, -100., 300., 200, 0., 200., "p_{T, jet} (GeV/c)", "Count", "dN^{Jets}/dNdp_{T}");

  AddHistogram2D<TH2D>("hLeadingJetPtRaw", "Jets p_{T} distribution (no bgrd. corr.)", "", 300, 0., 300., fNumberOfCentralityBins, 0, 100, "p_{T, jet} (GeV/c)", "Centrality", "dN^{Jets}/dp_{T}");
  AddHistogram2D<TH2D>("hLeadingJetPt", "Jets p_{T} distribution (background subtracted)", "", 400, -100., 300., fNumberOfCentralityBins, 0, 100, "p_{T, jet} (GeV/c)", "Centrality", "dN^{Jets}/dp_{T}");
  AddHistogram2D<TH2D>("hLeadingJetPhi", "Jet angular distribution #phi", "LEGO2", 180, 0., 2*TMath::Pi(), fNumberOfCentralityBins, 0, 100, "#phi", "Centrality", "dN^{Jets}/d#phi");
  AddHistogram2D<TH2D>("hLeadingJetEta", "Jet angular distribution #eta", "LEGO2", 100, -2.5, 2.5, fNumberOfCentralityBins, 0, 100, "#eta","Centrality","dN^{Jets}/d#eta");
  AddHistogram2D<TH2D>("hLeadingJetPhiPt", "Jet angular distribution #phi vs. p_{T}", "LEGO2", 180, 0., 2*TMath::Pi(), 400, -100., 300., "#phi", "p_{T, jet} (GeV/c)", "dN^{Jets}/d#phi dp_{T}");
  AddHistogram2D<TH2D>("hLeadingJetEtaPt", "Jet angular distribution #eta  vs. p_{T}", "LEGO2", 100, -2.5, 2.5, 400, -100., 300., "#eta","p_{T, jet} (GeV/c)","dN^{Jets}/d#eta dp_{T}");
  AddHistogram2D<TH2D>("hLeadingJetPhiEta", "Jet angular distribution #phi/#eta", "COLZ", 180, 0., 2*TMath::Pi(), 100, -2.5, 2.5, "#phi", "#eta", "dN^{Jets}/d#phi d#eta");
  AddHistogram2D<TH2D>("hLeadingJetArea", "Jet area", "LEGO2", 200, 0., 2., fNumberOfCentralityBins, 0, 100, "Jet A", "Centrality", "dN^{Jets}/dA");
  AddHistogram2D<TH2D>("hLeadingJetAreaPt", "Jet area vs. p_{T}", "LEGO2", 200, 0., 2., 400, -100., 300., "Jet A", "p_{T, jet} (GeV/c)", "dN^{Jets}/dA dp_{T}");
  AddHistogram2D<TH2D>("hLeadingJetPtLeadingHadron", "Jet leading hadron p_{T} distribution vs. jet p_{T}", "", 300, 0., 300., 300, 0., 300., "p_{T, jet} (GeV/c)", "p_{T,lead had} (GeV/c)", "dN^{Jets}/dp_{T}dp_{T,had}");

  AddHistogram2D<TH2D>("hSubleadingJetPtRaw", "Jets p_{T} distribution (no bgrd. corr.)", "", 300, 0., 300., fNumberOfCentralityBins, 0, 100, "p_{T, jet} (GeV/c)", "Centrality", "dN^{Jets}/dp_{T}");
  AddHistogram2D<TH2D>("hSubleadingJetPt", "Jets p_{T} distribution (background subtracted)", "", 400, -100., 300., fNumberOfCentralityBins, 0, 100, "p_{T, jet} (GeV/c)", "Centrality", "dN^{Jets}/dp_{T}");
  AddHistogram2D<TH2D>("hSubleadingJetPhi", "Jet angular distribution #phi", "LEGO2", 180, 0., 2*TMath::Pi(), fNumberOfCentralityBins, 0, 100, "#phi", "Centrality", "dN^{Jets}/d#phi");
  AddHistogram2D<TH2D>("hSubleadingJetEta", "Jet angular distribution #eta", "LEGO2", 100, -2.5, 2.5, fNumberOfCentralityBins, 0, 100, "#eta","Centrality","dN^{Jets}/d#eta");
  AddHistogram2D<TH2D>("hSubleadingJetPhiPt", "Jet angular distribution #phi vs. p_{T}", "LEGO2", 180, 0., 2*TMath::Pi(), 400, -100., 300., "#phi", "p_{T, jet} (GeV/c)", "dN^{Jets}/d#phi dp_{T}");
  AddHistogram2D<TH2D>("hSubleadingJetEtaPt", "Jet angular distribution #eta  vs. p_{T}", "LEGO2", 100, -2.5, 2.5, 400, -100., 300., "#eta","p_{T, jet} (GeV/c)","dN^{Jets}/d#eta dp_{T}");
  AddHistogram2D<TH2D>("hSubleadingJetPhiEta", "Jet angular distribution #phi/#eta", "COLZ", 180, 0., 2*TMath::Pi(), 100, -2.5, 2.5, "#phi", "#eta", "dN^{Jets}/d#phi d#eta");
  AddHistogram2D<TH2D>("hSubleadingJetArea", "Jet area", "LEGO2", 200, 0., 2., fNumberOfCentralityBins, 0, 100, "Jet A", "Centrality", "dN^{Jets}/dA");
  AddHistogram2D<TH2D>("hSubleadingJetAreaPt", "Jet area vs. p_{T}", "LEGO2", 200, 0., 2., 400, -100., 300., "Jet A", "p_{T, jet} (GeV/c)", "dN^{Jets}/dA dp_{T}");
  AddHistogram2D<TH2D>("hSubleadingJetPtLeadingHadron", "Jet leading hadron p_{T} distribution vs. jet p_{T}", "", 300, 0., 300., 300, 0., 300., "p_{T, jet} (GeV/c)", "p_{T,lead had} (GeV/c)", "dN^{Jets}/dp_{T}dp_{T,had}");

  AddHistogram2D<TH2D>("hTrackCount", "Number of tracks in acceptance vs. centrality", "LEGO2", 500, 0., 5000., fNumberOfCentralityBins, 0, 100, "N tracks","Centrality", "dN^{Events}/dN^{Tracks}");
  AddHistogram2D<TH2D>("hJetCount", "Number of jets in acceptance vs. centrality", "LEGO2", 100, 0., 100., fNumberOfCentralityBins, 0, 100, "N Jets","Centrality", "dN^{Events}/dN^{Jets}");
  AddHistogram2D<TH2D>("hFakeFactor", "Fake factor distribution", "LEGO2", 1000, 0., 100., fNumberOfCentralityBins, 0, 100, "Fake factor","Centrality", "dN^{Jets}/df");
  AddHistogram2D<TH2D>("hBackgroundPt", "Background p_{T} distribution", "", 150, 0., 150., fNumberOfCentralityBins, 0, 100, "Background p_{T} (GeV/c)", "Centrality", "dN^{Events}/dp_{T}");
  AddHistogram2D<TH2D>("hBackgroundPtJetPt_Cent0_100", "Background p_{T} distribution vs. jet p_{T}", "", 150, 0., 150.,  400, -100., 300., "Background p_{T} (GeV/c)", "Jet p_{T} (GeV/c)", "dN^{Events}/dp_{T}");
  AddHistogram2D<TH2D>("hBackgroundPtJetPt_Cent0_10", "Background p_{T} distribution vs. jet p_{T}", "", 150, 0., 150.,  400, -100., 300., "Background p_{T} (GeV/c)", "Jet p_{T} (GeV/c)", "dN^{Events}/dp_{T}");
  AddHistogram2D<TH2D>("hBackgroundPtConstCount_Cent0_100", "Background p_{T} distribution vs. const. count", "", 150, 0., 150., 200, 0., 200., "Background p_{T} (GeV/c)", "Count", "dN^{Events}/dp_{T}");
  AddHistogram2D<TH2D>("hBackgroundPtConstCount_Cent0_10", "Background p_{T} distribution vs. const. count", "", 150, 0., 150., 200, 0., 200., "Background p_{T} (GeV/c)", "Count", "dN^{Events}/dp_{T}");

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}


//________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronCF::ExecOnce() {

  AliAnalysisTaskEmcalJet::ExecOnce();

  // ### Add the jets as basic correlation particles to the event
  if (!(fInputEvent->FindListObject(Form("%s", fJetParticleArrayName.Data()))))
  {
    fJetsOutput = new TClonesArray("AliPicoTrack");
    fJetsOutput->SetName(fJetParticleArrayName.Data());
    fInputEvent->AddObject(fJetsOutput);
  }
  else
    AliError(Form("%s: Object with name %s already in event!", GetName(), Form("%s", fJetParticleArrayName.Data())));

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
      AliError(Form("%s: Object with name %s already in event!", GetName(), Form("%s", fTrackParticleArrayName.Data())));
  }

  // ### Import generated jets from toymodel for matching (optional)
  if(fJetMatchingArrayName != "")
  {
    fJetsInput = static_cast<TClonesArray*>(InputEvent()->FindListObject(Form("%s", fJetMatchingArrayName.Data())));
    if(!fJetsInput)
      AliFatal(Form("Importing jets for matching failed! Array '%s' not found!", fJetMatchingArrayName.Data()));
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
Bool_t AliAnalysisTaskChargedJetsHadronCF::IsEventCriteriumFulfilled()
{

    // In case of special selection criteria, trigger on certain events
    if(fEventCriteriumMode==0) // "minimum bias"
    {
      // do nothing
    }
    else if(fEventCriteriumMode==1) // background constraints
    {
      if( (fJetsCont->GetRhoVal() < fEventCriteriumMinBackground) || (fJetsCont->GetRhoVal() > fEventCriteriumMaxBackground) )
      {
        fHistEventRejection->Fill("JetCrit", 1);
        return kFALSE;
      }
    }
    else if(fEventCriteriumMode==2) // Minimum leading jet pT
    {
      if(fLeadingJet)
      {
        if(fLeadingJet->Pt() - fJetsCont->GetRhoVal()*fLeadingJet->Area() <= fEventCriteriumMinLeadingJetPt)
        {
          fHistEventRejection->Fill("JetCrit", 1);
          return kFALSE;
        }
      }
    }
    else if(fEventCriteriumMode==3) // Simple dijet trigger
    {
      if(fLeadingJet && fSubleadingJet)
      {
        if((fLeadingJet->Pt() - fJetsCont->GetRhoVal()*fLeadingJet->Area() <= fEventCriteriumMinLeadingJetPt) || (fSubleadingJet->Pt() - fJetsCont->GetRhoVal()*fSubleadingJet->Area() <= fEventCriteriumMinSubleadingJetPt))
        {
          fHistEventRejection->Fill("JetCrit", 1);
          return kFALSE;
      }
      }
    }
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskChargedJetsHadronCF::IsJetSelected(AliEmcalJet* jet)
{
  if( (fJetOutputMode==1) || (fJetOutputMode==3) ) // output the leading jet
    if(jet!=fLeadingJet)
      return kFALSE;

  if( (fJetOutputMode==2) || (fJetOutputMode==3) ) // output the subleading jet
    if(jet!=fSubleadingJet)
      return kFALSE;

  // Fake jet rejection (0810.1219)
/*
  if(fFakeFactorCutProfile)
  {
    Double_t fakeFactor = CalculateFakeFactor(jet);
    FillHistogram("hFakeFactor", fakeFactor, fCent);
    if( (fakeFactor >= fMinFakeFactorPercentage*fFakeFactorCutProfile->GetBinContent(fFakeFactorCutProfile->GetXaxis()->FindBin(fCent))) && (fakeFactor < fMaxFakeFactorPercentage*fFakeFactorCutProfile->GetBinContent(fFakeFactorCutProfile->GetXaxis()->FindBin(fCent))) )
      return kFALSE;
  }
*/

  // Poor man's fake jet rejection (according to jet const.)
  if(fRejectionFunction)
  {
    if( jet->GetNumberOfTracks() <  fRejectionFunction->Eval(jet->Pt() - fJetsCont->GetRhoVal()*jet->Area()) )
      return kFALSE;
  }

  // Jet matching. Only done if SetJetMatchingArrayName() called
  Bool_t matchedFound = kFALSE;
  if(fJetsInput)
  {
    // Go through all jets and check if the matching condition is fulfiled by at least one jet
    Double_t bestMatchDeltaR = 999.;
    for(Int_t i=0; i<fJetsInput->GetEntries(); i++)
    {
      AliEmcalJet* matchJet = static_cast<AliEmcalJet*>(fJetsInput->At(i));
      Double_t deltaPhi = TMath::Min(TMath::Abs(jet->Phi()-matchJet->Phi()),TMath::TwoPi() - TMath::Abs(jet->Phi()-matchJet->Phi()));
      Double_t deltaEta = TMath::Abs(jet->Eta() - matchJet->Eta());
      Double_t deltaR   = TMath::Sqrt((deltaPhi*deltaPhi) + (deltaEta*deltaEta));

      if(deltaR < bestMatchDeltaR)
      {
        bestMatchDeltaR = deltaR;
      }
    }
    // Check if a matching jet is found.
    if(bestMatchDeltaR < 0.9*fJetsCont->GetJetRadius())
      matchedFound = kTRUE;
  }

  if(fJetOutputMode==4) // matching jets only
    return matchedFound;
  else if(fJetOutputMode==5) // non-matching jets only
    return !matchedFound;

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronCF::FillHistogramsJets(AliEmcalJet* jet)
{
    // All jets
    FillHistogram("hJetPtRaw", jet->Pt(), fCent); 
    FillHistogram("hJetPt", jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), fCent); 
    FillHistogram("hJetPhi", jet->Phi(), fCent); 
    FillHistogram("hJetEta", jet->Eta(), fCent); 
    FillHistogram("hJetEtaPt", jet->Eta(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area()); 
    FillHistogram("hJetPhiPt", jet->Phi(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area()); 
    FillHistogram("hJetPhiEta", jet->Phi(), jet->Eta()); 
    FillHistogram("hJetArea", jet->Area(), fCent); 
    FillHistogram("hJetAreaPt", jet->Area(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area()); 
    FillHistogram("hJetPtLeadingHadron", jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), fJetsCont->GetLeadingHadronPt(jet));

    FillHistogram("hBackgroundPtJetPt_Cent0_100",  fJetsCont->GetRhoVal(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area()); 
    if( (fCent >= 0) && (fCent < 10) )
      FillHistogram("hBackgroundPtJetPt_Cent0_10",  fJetsCont->GetRhoVal(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area()); 

    // Leading jet plots
    if(jet==fLeadingJet)
    {
      FillHistogram("hLeadingJetPtRaw", jet->Pt(), fCent); 
      FillHistogram("hLeadingJetPt", jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), fCent); 
      FillHistogram("hLeadingJetPhi", jet->Phi(), fCent); 
      FillHistogram("hLeadingJetEta", jet->Eta(), fCent); 
      FillHistogram("hLeadingJetEtaPt", jet->Eta(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area()); 
      FillHistogram("hLeadingJetPhiPt", jet->Phi(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area()); 
      FillHistogram("hLeadingJetPhiEta", jet->Phi(), jet->Eta()); 
      FillHistogram("hLeadingJetArea", jet->Area(), fCent); 
      FillHistogram("hLeadingJetAreaPt", jet->Area(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area()); 
      FillHistogram("hLeadingJetPtLeadingHadron", jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), fJetsCont->GetLeadingHadronPt(jet));
    }

    // Subleading jet plot
    else if(jet==fSubleadingJet)
    {
      FillHistogram("hSubleadingJetPtRaw", jet->Pt(), fCent); 
      FillHistogram("hSubleadingJetPt", jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), fCent); 
      FillHistogram("hSubleadingJetPhi", jet->Phi(), fCent); 
      FillHistogram("hSubleadingJetEta", jet->Eta(), fCent); 
      FillHistogram("hSubleadingJetEtaPt", jet->Eta(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area()); 
      FillHistogram("hSubleadingJetPhiPt", jet->Phi(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area()); 
      FillHistogram("hSubleadingJetPhiEta", jet->Phi(), jet->Eta());
      FillHistogram("hSubleadingJetArea", jet->Area(), fCent); 
      FillHistogram("hSubleadingJetAreaPt", jet->Area(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area()); 
      FillHistogram("hSubleadingJetPtLeadingHadron", jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), fJetsCont->GetLeadingHadronPt(jet));
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
void AliAnalysisTaskChargedJetsHadronCF::FillHistogramsJetConstituents(AliEmcalJet* jet)
{
  // Loop over all jet constituents
  for(Int_t i = 0; i < jet->GetNumberOfTracks(); i++)
  {
    AliVParticle* constituent = static_cast<AliVParticle*>(jet->TrackAt(i, fTracksCont->GetArray()));
    if(!constituent) 
      continue;

    // Fill jet constituent plots
    FillHistogram("hJetConstituentPt_Cent0_100", jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), constituent->Pt()); 
    if( (fCent >= 0) && (fCent < 10) )
      FillHistogram("hJetConstituentPt_Cent0_10", jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), constituent->Pt()); 
  }

  FillHistogram("hJetConstituentCount_Cent0_100", jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), jet->GetNumberOfTracks()); 
  if( (fCent >= 0) && (fCent < 10) )
    FillHistogram("hJetConstituentCount_Cent0_10", jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), jet->GetNumberOfTracks()); 

  FillHistogram("hBackgroundPtConstCount_Cent0_100", fJetsCont->GetRhoVal(), jet->GetNumberOfTracks()); 
  if( (fCent >= 0) && (fCent < 10) )
    FillHistogram("hBackgroundPtConstCount_Cent0_10", fJetsCont->GetRhoVal(), jet->GetNumberOfTracks()); 
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronCF::AddJetToOutputArray(AliEmcalJet* jet)
{
  new ((*fJetsOutput)[fAcceptedJets]) AliPicoTrack(jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), jet->Eta(), jet->Phi(), jet->Charge(), 0, 0);
  fAcceptedJets++;
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
  if( (jet->Pt()-jet->Area()*fJetsCont->GetRhoVal()) < fExtractionMinPt )
    return;

  // Discard jets statistically
  if(fRandom->Rndm() >= fExtractionPercentage)
    return;

  Long64_t eventID = InputEvent()->GetHeader()->GetEventIdAsLong();
  AliBasicJet basicJet(jet->Eta(), jet->Phi(), jet->Pt(), jet->Charge(), fJetsCont->GetJetRadius(), jet->Area(), fJetsCont->GetRhoVal(), eventID, fCent);
  // Add constituents
  for(Int_t i = 0; i < jet->GetNumberOfTracks(); i++)
  {
    AliVParticle* particle = static_cast<AliVParticle*>(jet->TrackAt(i, fTracksCont->GetArray()));
    if(!particle) continue;
    basicJet.AddJetConstituent(particle->Eta(), particle->Phi(), particle->Pt(), particle->Charge());
  }
  fJetsTreeBuffer = &basicJet;
  fJetsTree->Fill();
}

//________________________________________________________________________
Bool_t AliAnalysisTaskChargedJetsHadronCF::Run()
{
  CalculateEventProperties();

  if(!IsEventCriteriumFulfilled())
    return kFALSE;

  // ####### Jet loop
  fAcceptedJets = 0;
  fJetsCont->ResetCurrentID();
  while(AliEmcalJet *jet = fJetsCont->GetNextAcceptJet())
  {
    if(!IsJetSelected(jet))
      continue;

    // Jet plots
    FillHistogramsJets(jet);
    FillHistogramsJetConstituents(jet);

    // Add jet to output array
    if(fExtractionPercentage)
      AddJetToTree(jet);
    AddJetToOutputArray(jet);
  }


  // ####### Particle loop
  fAcceptedTracks = 0;
  fTracksCont->ResetCurrentID();
  Int_t trackcount = 0;
  while(AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle()))
  {
    // Track plots
    FillHistogramsTracks(track);

    // Add track to output array
    trackcount++;
    AddTrackToOutputArray(track);
  }

  // ####### Event properties
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
void AliAnalysisTaskChargedJetsHadronCF::CalculateEventProperties()
{
  // Calculate leading + subleading jet
  fLeadingJet    = fJetsCont->GetLeadingJet("rho");
  fSubleadingJet = GetSubleadingJet("rho");
}

//________________________________________________________________________
Double_t AliAnalysisTaskChargedJetsHadronCF::CalculateFakeFactor(AliEmcalJet* jet)
{
  Double_t fakeFactor = 0;

  // Loop over all jet constituents
  for(Int_t i = 0; i < jet->GetNumberOfTracks(); i++)
  {
    AliVParticle* constituent = static_cast<AliVParticle*>(jet->TrackAt(i, fTracksCont->GetArray()));

    Double_t deltaPhi = TMath::Min(TMath::Abs(jet->Phi()-constituent->Phi()),TMath::TwoPi() - TMath::Abs(jet->Phi()-constituent->Phi()));
    Double_t deltaR = TMath::Sqrt( (jet->Eta() - constituent->Eta())*(jet->Eta() - constituent->Eta()) + deltaPhi*deltaPhi );
    fakeFactor += constituent->Pt() * TMath::Sin(deltaR);
  }

  return fakeFactor;
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronCF::SetEventCriteriumSelection(Int_t type)
{
  fEventCriteriumMode = type;

  if(fEventCriteriumMode==0)
    AliWarning("Set event criterium to 'default'              -- no further selection criterium.");
  else if(fEventCriteriumMode==1)
    AliWarning("Set event criterium to 'background'           -- select events with certain backgrounds");
  else if(fEventCriteriumMode==2)
    AliWarning("Set event criterium to 'simple jet trigger'   -- select events with certain minimum leading jet pT (bgrd corr.)");
  else if(fEventCriteriumMode==3)
    AliWarning("Set event criterium to 'simple dijet trigger' -- select events with certain minimum leading + subleading jet pT (bgrd corr.)");
  else
  {
    AliFatal("Event criterium not valid.");
  }
}


//________________________________________________________________________
AliEmcalJet* AliAnalysisTaskChargedJetsHadronCF::GetSubleadingJet(const char* opt)
{
  // Customized from AliJetContainer::GetLeadingJet()
  // Get the subleading jet; if opt contains "rho" the sorting is according to pt-A*rho

  TString option(opt);
  option.ToLower();

  AliEmcalJet *jetLeading = fJetsCont->GetLeadingJet(opt);
  AliEmcalJet *jetSubLeading = 0;

  fJetsCont->ResetCurrentID();
  AliEmcalJet *jet = fJetsCont->GetNextAcceptJet();
  Double_t     tmpPt = 0;

  if (option.Contains("rho")) {
    while ((jet = fJetsCont->GetNextAcceptJet())) {
      if(jet == jetLeading)
        continue;
      else if ( (jet->Pt()-jet->Area()*fJetsCont->GetRhoVal()) > tmpPt )
      {
        jetSubLeading = jet;
        tmpPt = jet->Pt()-jet->Area()*fJetsCont->GetRhoVal();
      }

    }
  }
  else {
    while ((jet = fJetsCont->GetNextAcceptJet())) {
      if(jet == jetLeading)
        continue;
      else if ( jet->Pt() > tmpPt )
      {
        jetSubLeading = jet;
        tmpPt = jet->Pt();
      }
    }
  }

  return jetSubLeading;
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

