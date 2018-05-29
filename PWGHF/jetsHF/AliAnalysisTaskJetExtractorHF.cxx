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

#include <cstdlib>
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

#include "TObjArray.h"
#include "AliESDVertex.h"
#include "AliAODVertex.h"
#include "AliAODPid.h"

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
#include "AliESDtrackCuts.h"
#include "AliRhoParameter.h"

#include "AliHFJetsTaggingVertex.h"
#include "AliAnalysisTaskJetExtractorHF.h"

//________________________________________________________________________
AliBasicJet::~AliBasicJet() 
{
// dummy destructor
}

//________________________________________________________________________
AliBasicPID::~AliBasicPID() 
{
// dummy destructor
}

//________________________________________________________________________
AliBasicJetConstituent::~AliBasicJetConstituent() 
{
// dummy destructor
}

//________________________________________________________________________
AliBasicJetSecondaryVertex::~AliBasicJetSecondaryVertex() 
{
// dummy destructor
}

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskJetExtractorHF)
/// \endcond
//________________________________________________________________________
AliAnalysisTaskJetExtractorHF::AliAnalysisTaskJetExtractorHF() :
  AliAnalysisTaskEmcalJet("AliAnalysisTaskJetExtractorHF", kTRUE),
  fTruthParticleArray(0),
  fJetsCont(0),
  fTracksCont(0),
  fJetsTree(0),
  fJetsTreeBuffer(0),
  fRandom(0),
  fVtxTagger(0),
  fCurrentNJetsInEvents(0),
  fCurrentNJetsInEventsTotal(0),
  fCurrentJetTypeHM(0),
  fCurrentJetTypePM(0),
  fCurrentLeadingJet(0),
  fCurrentSubleadingJet(0),
  fCurrentJetHFHadron(0),
  fCurrentTrueJetPt(0),
  fUnderflowBinContents(),
  fExtractionCutMinCent(-1),
  fExtractionCutMaxCent(-1),
  fExtractionCutUsePM(kFALSE),
  fExtractionCutUseHM(kFALSE),
  fExtractionCutMinPt(0.),
  fExtractionCutMaxPt(200.),
  fExtractionPercentage(1.0),
  fExtractionListPIDsHM(),
  fSetEmcalJetFlavour(0),
  fUnderflowNumBins(1),
  fUnderflowCutOff(5.0),
  fUnderflowPercentage(0.5),
  fHadronMatchingRadius(0.4),
  fTruthJetsArrayName(""),
  fTruthJetsRhoName(""),
  fTruthParticleArrayName("mcparticles"),
  fSecondaryVertexMaxChi2(1e10),
  fSecondaryVertexMaxDispersion(0.05),
  fAddPIDSignal(kFALSE),
  fCalculateSecondaryVertices(kFALSE),
  fVertexerCuts(0)
{
  // Default constructor.
  SetMakeGeneralHistograms(kTRUE);
  fRandom = new TRandom3(0);
  for(size_t i=0;i<100;i++)
    fUnderflowBinContents[i] = 0;
}

//________________________________________________________________________
AliAnalysisTaskJetExtractorHF::AliAnalysisTaskJetExtractorHF(const char *name) :
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fTruthParticleArray(0),
  fJetsCont(0),
  fTracksCont(0),
  fJetsTree(0),
  fJetsTreeBuffer(0),
  fRandom(0),
  fVtxTagger(0),
  fCurrentNJetsInEvents(0),
  fCurrentNJetsInEventsTotal(0),
  fCurrentJetTypeHM(0),
  fCurrentJetTypePM(0),
  fCurrentLeadingJet(0),
  fCurrentSubleadingJet(0),
  fCurrentJetHFHadron(0),
  fCurrentTrueJetPt(0),
  fUnderflowBinContents(),
  fExtractionCutMinCent(-1),
  fExtractionCutMaxCent(-1),
  fExtractionCutUsePM(kFALSE),
  fExtractionCutUseHM(kFALSE),
  fExtractionCutMinPt(0.),
  fExtractionCutMaxPt(200.),
  fExtractionPercentage(1.0),
  fExtractionListPIDsHM(),
  fSetEmcalJetFlavour(0),
  fUnderflowNumBins(1),
  fUnderflowCutOff(5.0),
  fUnderflowPercentage(0.5),
  fHadronMatchingRadius(0.4),
  fTruthJetsArrayName(""),
  fTruthJetsRhoName(""),
  fTruthParticleArrayName("mcparticles"),
  fSecondaryVertexMaxChi2(1e10),
  fSecondaryVertexMaxDispersion(0.05),
  fAddPIDSignal(kFALSE),
  fCalculateSecondaryVertices(kFALSE),
  fVertexerCuts(0)
{
  // Default constructor.
  SetMakeGeneralHistograms(kTRUE);
  fRandom = new TRandom3(0);
  for(size_t i=0;i<100;i++)
    fUnderflowBinContents[i] = 0;
}

//________________________________________________________________________
AliAnalysisTaskJetExtractorHF::~AliAnalysisTaskJetExtractorHF()
{
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskJetExtractorHF::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  // ### Basic container settings
  fJetsCont           = GetJetContainer(0);
  if(!fJetsCont)
    AliFatal("Jet input container not found!");
  fJetsCont->PrintCuts();
  fTracksCont       = static_cast<AliTrackContainer*>(fJetsCont->GetParticleContainer());
  if(!fTracksCont)
    AliFatal("Particle input container not found attached to jets!");

  //if(fTracksCont) fTracksCont->SetClassName("AliVTrack");

  // ### Add control histograms (already some created in base task)
  AddHistogram2D<TH2D>("hJetCount", "Number of tagged jets in acceptance vs. centrality", "COLZ", 100, 0., 100., 100, 0, 100, "N Jets","Centrality", "dN^{Events}/dN^{Jets}");
  AddHistogram2D<TH2D>("hJetCountTotal", "Number of inclusive jets in acceptance vs. centrality", "COLZ", 100, 0., 100., 100, 0, 100, "N Jets","Centrality", "dN^{Events}/dN^{Jets}");
  AddHistogram2D<TH2D>("hTrackCount", "Number of tracks in acceptance vs. centrality", "COLZ", 500, 0., 5000., 100, 0, 100, "N tracks","Centrality", "dN^{Events}/dN^{Tracks}");
  AddHistogram2D<TH2D>("hBackgroundPt", "Background p_{T} distribution", "", 150, 0., 150., 100, 0, 100, "Background p_{T} (GeV/c)", "Centrality", "dN^{Events}/dp_{T}");

  AddHistogram2D<TH2D>("hJetPtRaw", "Jets p_{T} distribution (raw)", "COLZ", 300, 0., 300., 100, 0, 100, "p_{T, jet} (GeV/c)", "Centrality", "dN^{Jets}/dp_{T}");
  AddHistogram2D<TH2D>("hJetPt", "Jets p_{T} distribution (background subtracted)", "COLZ", 400, -100., 300., 100, 0, 100, "p_{T, jet} (GeV/c)", "Centrality", "dN^{Jets}/dp_{T}");
  AddHistogram2D<TH2D>("hJetPtAcceptance", "Statistically discarded jets p_{T} distribution (background subtracted)", "COLZ", 400, -100., 300., 100, 0, 100, "p_{T, jet} (GeV/c)", "Centrality", "dN^{Discarded jets}/dp_{T}");
  AddHistogram2D<TH2D>("hJetPhiEta", "Jet angular distribution #phi/#eta", "COLZ", 180, 0., 2*TMath::Pi(), 100, -2.5, 2.5, "#phi", "#eta", "dN^{Jets}/d#phi d#eta");
  AddHistogram2D<TH2D>("hJetArea", "Jet area", "COLZ", 200, 0., 2., 100, 0, 100, "Jet A", "Centrality", "dN^{Jets}/dA");

  TH1* tmpHisto = AddHistogram1D<TH1D>("hJetAcceptance", "Accepted jets", "", 5, 0, 5, "stage","N^{jets}/cut");
  tmpHisto->GetXaxis()->SetBinLabel(1, "Before cuts");
  tmpHisto->GetXaxis()->SetBinLabel(2, "Centrality");
  tmpHisto->GetXaxis()->SetBinLabel(3, "p_{T}");
  tmpHisto->GetXaxis()->SetBinLabel(4, "PID");
  tmpHisto->GetXaxis()->SetBinLabel(5, "Extraction percentage");

  tmpHisto = AddHistogram1D<TH1D>("hBJetFromInitialCollision", "B-jets directly from a b/anti b quark", "", 2, 0, 2, "IC","N^{jets}");
  tmpHisto->GetXaxis()->SetBinLabel(1, "From gluon splittings etc.");
  tmpHisto->GetXaxis()->SetBinLabel(2, "From initial collision");


  AddHistogram1D<TH1D>("hJetTypeHM", "Jet types from HM (before tagging)", "COLZ", 100, 0, 100, "Jet type from hadron matching", "dN^{Jets}/dType");
  AddHistogram1D<TH1D>("hJetTypePM", "Jet types from PM (before tagging)", "COLZ", 100, 0, 100, "Jet type from parton matching", "dN^{Jets}/dType");

  AddHistogram1D<TH1D>("hJetHFHadronHM_C", "Hadron content for HM c-jets (before tagging)", "COLZ", 30000, 0, 30000,  "PDG code", "dN^{Jets}/dType");
  AddHistogram1D<TH1D>("hJetHFHadronHM_B", "Hadron content for HM b-jets (before tagging)", "COLZ", 30000, 0, 30000,  "PDG code", "dN^{Jets}/dType");
  AddHistogram1D<TH1D>("hJetHFHadronPM_C", "Hadron content for PM c-jets (before tagging)", "COLZ", 30000, 0, 30000,  "PDG code", "dN^{Jets}/dType");
  AddHistogram1D<TH1D>("hJetHFHadronPM_B", "Hadron content for PM b-jets (before tagging)", "COLZ", 30000, 0, 30000,  "PDG code", "dN^{Jets}/dType");

  AddHistogram2D<TH2D>("hJetTypeCorr", "Jet type correlation", "COLZ", 7, 0, 7, 7, 0, 7, "Jet type from hadron matching", "Jet type from parton matching", "dN^{Jets}/dType");
  AddHistogram2D<TH2D>("hDeltaPt", "#delta p_{T} distribution", "", 400, -100., 300., 100, 0, 100, "p_{T, cone} (GeV/c)", "Centrality", "dN^{Tracks}/dp_{T}");

  AddHistogram2D<TH2D>("hConstituentPt", "Jet constituent p_{T} distribution", "COLZ", 400, 0., 300., 100, 0, 100, "p_{T, const} (GeV/c)", "Centrality", "dN^{Const}/dp_{T}");
  AddHistogram2D<TH2D>("hConstituentPhiEta", "Jet constituent relative #phi/#eta distribution", "COLZ", 120, -0.6, 0.6, 120, -0.6, 0.6, "#Delta#phi", "#Delta#eta", "dN^{Const}/d#phi d#eta");

  // Track QA plots
  AddHistogram2D<TH2D>("hTrackPt", "Tracks p_{T} distribution", "", 300, 0., 300., 100, 0, 100, "p_{T} (GeV/c)", "Centrality", "dN^{Tracks}/dp_{T}");
  AddHistogram2D<TH2D>("hTrackPhi", "Track angular distribution in #phi", "LEGO2", 180, 0., 2*TMath::Pi(), 100, 0, 100, "#phi", "Centrality", "dN^{Tracks}/(d#phi)");
  AddHistogram2D<TH2D>("hTrackEta", "Track angular distribution in #eta", "LEGO2", 100, -2.5, 2.5, 100, 0, 100, "#eta", "Centrality", "dN^{Tracks}/(d#eta)");
  AddHistogram2D<TH2D>("hTrackPhiEta", "Track angular distribution #phi/#eta", "COLZ", 180, 0., 2*TMath::Pi(), 100, -2.5, 2.5, "#phi", "#eta", "dN^{Tracks}/d#phi d#eta");

  AddHistogram2D<TH2D>("hTrackEtaPt", "Track angular distribution in #eta vs. p_{T}", "LEGO2", 100, -2.5, 2.5, 300, 0., 300., "#eta", "p_{T} (GeV/c)", "dN^{Tracks}/(d#eta dp_{T})");
  AddHistogram2D<TH2D>("hTrackPhiPt", "Track angular distribution in #phi vs. p_{T}", "LEGO2", 180, 0, 2*TMath::Pi(), 300, 0., 300., "#phi", "p_{T} (GeV/c)", "dN^{Tracks}/(d#phi dp_{T})");

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}


//________________________________________________________________________
void AliAnalysisTaskJetExtractorHF::ExecOnce() {

  AliAnalysisTaskEmcalJet::ExecOnce();

  // ### Prepare the jet tree
  fJetsTree = new TTree("ExtractedJets", "ExtractedJets");
  fJetsTree->Branch("Jets", "AliBasicJet", &fJetsTreeBuffer, 1000);
  fOutput->Add(fJetsTree);

  if (fTruthParticleArrayName != "")
    fTruthParticleArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTruthParticleArrayName.Data()));

  // ### Prepare vertexer
  if(fCalculateSecondaryVertices)
  {
    if(!fVertexerCuts)
      AliFatal("VertexerCuts not given but secondary vertex calculation turned on.");
    fVtxTagger = new AliHFJetsTaggingVertex();
    fVtxTagger->SetCuts(fVertexerCuts);
  }

  // Status message on low-pt extraction
  if(fUnderflowNumBins)
    AliWarning(Form("Low-pT extraction active: %2.0f%% of jets will be extracted in %d bin(s) from %2.2f to %2.2f GeV/c.", fUnderflowPercentage*100., fUnderflowNumBins, fUnderflowCutOff, fExtractionCutMinPt));
  else
    AliWarning(Form("Low-pT extraction not active: Jets below %2.2f GeV/c will be discarded completely.", fExtractionCutMinPt));

}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetExtractorHF::IsJetSelected(AliEmcalJet* jet)
{
  // ######## Check if cuts for extraction have been passed

  // Extraction pt, jet PID (according to initial parton, according to HM), centrality
  FillHistogram("hJetAcceptance", 0.5);

  // ### CENTRALITY
  if(fExtractionCutMinCent != -1 && (fCent < fExtractionCutMinCent || fCent >= fExtractionCutMaxCent))
    return kFALSE;
  FillHistogram("hJetAcceptance", 1.5);

  Double_t jetPt = jet->Pt()-jet->Area()*fJetsCont->GetRhoVal();

  // ### PT
  if(jetPt >= fExtractionCutMaxPt)
  {
    FillHistogram("hJetPtAcceptance", jetPt, fCent);
    return kFALSE;
  }

  if(jetPt < fExtractionCutMinPt && jetPt >= fUnderflowCutOff)
  {
    if(!fUnderflowNumBins) // If low-pT extraction is not active, discard jet
    {
      FillHistogram("hJetPtAcceptance", jetPt, fCent);
      return kFALSE;
    }
    else // If active, put in underflow bin
    {
      if(fRandom->Rndm() >= fUnderflowPercentage)
      {
        FillHistogram("hJetPtAcceptance", jetPt, fCent);
        return kFALSE;
      }
      Double_t extractionBinSize = (fExtractionCutMinPt-fUnderflowCutOff)/fUnderflowNumBins;
      // Low-pT extraction
      for(Int_t i=0; i<fUnderflowNumBins; i++)
        if( (jetPt >= i*extractionBinSize + fUnderflowCutOff) && (jetPt < (i+1)*extractionBinSize + fUnderflowCutOff))
        {
          if (fUnderflowBinContents[i+1] >= fUnderflowBinContents[0])
          {
            FillHistogram("hJetPtAcceptance", jetPt, fCent);
            return kFALSE;
          }
        }
    }
  }
  fCurrentNJetsInEventsTotal++; // number inclusive jets within kinematic extraction criterium
  FillHistogram("hJetAcceptance", 2.5);

  // Fill histograms to evaluate fractions of different jet types
  // This is after the kinematic extraction cuts, but before jet ID filtering
  FillHistogram("hJetTypeHM", fCurrentJetTypeHM);
  FillHistogram("hJetTypePM", fCurrentJetTypePM);

  // Jet hadron content of b- and c-jets
  if(fCurrentJetTypeHM == 4)
    FillHistogram("hJetHFHadronHM_C", fCurrentJetHFHadron);
  else if (fCurrentJetTypeHM == 5)
    FillHistogram("hJetHFHadronHM_B", fCurrentJetHFHadron);
  if(fCurrentJetTypePM == 4)
    FillHistogram("hJetHFHadronPM_C", fCurrentJetHFHadron);
  else if (fCurrentJetTypePM == 5)
    FillHistogram("hJetHFHadronPM_B", fCurrentJetHFHadron);

  Bool_t passedCutPID = kTRUE;
  // ### PID (from parton matching)
  if(fExtractionCutUsePM)
  {
    passedCutPID = kFALSE;
    for(size_t i=0; i<fExtractionListPIDsPM.size(); i++)
    {
      if (fExtractionListPIDsPM.at(i) == fCurrentJetTypePM)
        passedCutPID = kTRUE;
    }
  }
  // ### PID (from hadron matching)
  else if(fExtractionCutUseHM)
  {
    passedCutPID = kFALSE;
    for(size_t i=0; i<fExtractionListPIDsHM.size(); i++)
    {
      if (fExtractionListPIDsHM.at(i) == fCurrentJetTypeHM)
        passedCutPID = kTRUE;
    }
  }

  // Note: When either of the cuts is passed, accepted it
  if(!passedCutPID)
    return kFALSE;
  FillHistogram("hJetAcceptance", 3.5);
  fCurrentNJetsInEvents++; // number tagged jets within kinematic extraction criterium


  // Discard jets statistically
  if(fRandom->Rndm() >= fExtractionPercentage)
  {
    FillHistogram("hJetPtAcceptance", jetPt, fCent);
    return kFALSE;
  }
  FillHistogram("hJetAcceptance", 4.5);

  if(fUnderflowNumBins)
  {
    // Fill low-pt extraction bin contents
    // Those contents decide whether we extract further jets in an "underflow" bin
    Double_t extractionBinSize = (fExtractionCutMinPt-fUnderflowCutOff)/fUnderflowNumBins;
    for(Int_t i=0; i<fUnderflowNumBins; i++)
      if( (jetPt >= i*extractionBinSize + fUnderflowCutOff) && (jetPt < (i+1)*extractionBinSize + fUnderflowCutOff))
        fUnderflowBinContents[i+1]++;

    if( (jetPt >= fExtractionCutMinPt) && (jetPt < fExtractionCutMinPt+extractionBinSize) )
      fUnderflowBinContents[0]++;
  }

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskJetExtractorHF::AddJetToTree(AliEmcalJet* jet)
{
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

  // Get event ID from header
  AliVHeader* eventIDHeader = InputEvent()->GetHeader();
  Long64_t eventID = 0;
  if(eventIDHeader)
    eventID = eventIDHeader->GetEventIdAsLong();

  // ### Actually add the basic jet
  AliBasicJet basicJet(jet->Eta(), jet->Phi(), jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), jet->Charge(), fJetsCont->GetJetRadius(), jet->Area(), fCurrentJetTypePM, fCurrentJetTypeHM, fJetsCont->GetRhoVal(), InputEvent()->GetMagneticField(), vtxX, vtxY, vtxZ, eventID, fCent, jet->M(), fPtHard);

  // Add constituents
  for(Int_t i = 0; i < jet->GetNumberOfTracks(); i++)
  {
    AliVParticle* particle = static_cast<AliVParticle*>(jet->TrackAt(i, fTracksCont->GetArray()));
    if(!particle) continue;

    // Impact parameter analysis
    Double_t d0 = 0;
    Double_t z0 = 0;
    Double_t d0cov = 0;
    Double_t z0cov = 0;
    // Note: significant IPs are given by d0/sqrt(d0cov) & z0/sqrt(z0cov)
    GetTrackImpactParameters(myVertex, dynamic_cast<AliAODTrack*>(particle), d0, d0cov, z0, z0cov);

    basicJet.AddJetConstituent(particle->Eta(), particle->Phi(), particle->Pt(), particle->Charge(), particle->Xv(), particle->Yv(), particle->Zv(), d0, d0cov, z0, z0cov);
    AddPIDInformation(particle, *basicJet.GetJetConstituent(basicJet.GetNumbersOfConstituents()-1));
  }

  // Add secondary vertices to jet (from AOD.VertexingHF or calculated)
  AddSecondaryVertices(myVertex, jet, basicJet);

  // Add here pt from matched jet
  basicJet.SetTruePt(fCurrentTrueJetPt);

  fJetsTreeBuffer = &basicJet;
  fJetsTree->Fill();
}

//________________________________________________________________________
void AliAnalysisTaskJetExtractorHF::FillEventControlHistograms()
{
  // ### Event control plots
  FillHistogram("hTrackCount", fTracksCont->GetNAcceptedParticles(), fCent);
  FillHistogram("hJetCount", fCurrentNJetsInEvents, fCent);
  FillHistogram("hJetCountTotal", fCurrentNJetsInEventsTotal, fCent);
  FillHistogram("hBackgroundPt", fJetsCont->GetRhoVal(), fCent);
}

//________________________________________________________________________
void AliAnalysisTaskJetExtractorHF::FillJetControlHistograms(AliEmcalJet* jet)
{
  // ### Jet control plots
  FillHistogram("hJetPtRaw", jet->Pt(), fCent); 
  FillHistogram("hJetPt", jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), fCent);
  FillHistogram("hJetPhiEta", jet->Phi(), jet->Eta());
  FillHistogram("hJetArea", jet->Area(), fCent);
  FillHistogram("hJetTypeCorr", fCurrentJetTypeHM, fCurrentJetTypePM);

  // ### Jet constituent plots
  for(Int_t i = 0; i < jet->GetNumberOfTracks(); i++)
  {
    AliVParticle* jetConst = static_cast<AliVParticle*>(jet->TrackAt(i, fTracksCont->GetArray()));
    if(!jetConst) continue;

    // Constituent eta/phi (relative to jet)
    Double_t deltaEta = jet->Eta() - jetConst->Eta();
    Double_t deltaPhi = TMath::Min(TMath::Abs(jet->Phi() - jetConst->Phi()), TMath::TwoPi() - TMath::Abs(jet->Phi() - jetConst->Phi()));
    if(!((jet->Phi() - jetConst->Phi() < 0) && (jet->Phi() - jetConst->Phi() <= TMath::Pi())))
    deltaPhi = -deltaPhi;

    FillHistogram("hConstituentPt", jetConst->Pt(), fCent);
    FillHistogram("hConstituentPhiEta", deltaPhi, deltaEta);
  }

  // ### Random cone / delta pT plots
  const Int_t kNumRandomConesPerEvent = 4;
  for(Int_t iCone=0; iCone<kNumRandomConesPerEvent; iCone++)
  {
    // Throw random cone
    Double_t tmpRandConeEta = fJetsCont->GetJetEtaMin() + fRandom->Rndm()*TMath::Abs(fJetsCont->GetJetEtaMax()-fJetsCont->GetJetEtaMin());
    Double_t tmpRandConePhi = fRandom->Rndm()*TMath::TwoPi();
    Double_t tmpRandConePt  = 0;
    // Fill pT that is in cone
    fTracksCont->ResetCurrentID();
    while(AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle()))
      if(IsTrackInCone(track, tmpRandConeEta, tmpRandConePhi, fJetsCont->GetJetRadius()))
        tmpRandConePt += track->Pt();

    // Fill histograms
    FillHistogram("hDeltaPt", tmpRandConePt - fJetsCont->GetRhoVal()*fJetsCont->GetJetRadius()*fJetsCont->GetJetRadius()*TMath::Pi(), fCent);
  }

}

//________________________________________________________________________
void AliAnalysisTaskJetExtractorHF::FillTrackControlHistograms(AliVTrack* track)
{
  FillHistogram("hTrackPt", track->Pt(), fCent);
  FillHistogram("hTrackPhi", track->Phi(), fCent);
  FillHistogram("hTrackEta", track->Eta(), fCent);
  FillHistogram("hTrackEtaPt", track->Eta(), track->Pt());
  FillHistogram("hTrackPhiPt", track->Phi(), track->Pt());
  FillHistogram("hTrackPhiEta", track->Phi(), track->Eta());
}
//________________________________________________________________________
Bool_t AliAnalysisTaskJetExtractorHF::Run()
{
  CalculateEventProperties();

  // ### Jet loop
  fJetsCont->ResetCurrentID();
  while(AliEmcalJet *jet = fJetsCont->GetNextAcceptJet())
  {
    CalculateJetProperties(jet);

    if(!IsJetSelected(jet))
      continue;

    FillJetControlHistograms(jet);
    AddJetToTree(jet);
  }

  // ### Particle loop
  fTracksCont->ResetCurrentID();
  while(AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle()))
    FillTrackControlHistograms(track);

  FillEventControlHistograms();

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskJetExtractorHF::AddSecondaryVertices(const AliVVertex* primVtx, const AliEmcalJet* jet, AliBasicJet& basicJet)
{
  if(!primVtx)
    return;

  // Create ESD vertex from the existing AliVVertex
  Double_t vtxPos[3]   = {primVtx->GetX(), primVtx->GetY(), primVtx->GetZ()};
  Double_t covMatrix[6] = {0};
  primVtx->GetCovarianceMatrix(covMatrix);
  AliESDVertex* esdVtx = new AliESDVertex(vtxPos, covMatrix, primVtx->GetChi2(), primVtx->GetNContributors());

  TClonesArray* secVertexArr = 0;
  vctr_pair_dbl_int arrDispersion;
  arrDispersion.reserve(5);
  if(fCalculateSecondaryVertices)
  {
    //###########################################################################
    // ********* Calculate secondary vertices
    // Derived from AliAnalysisTaskEmcalJetBtagSV
    secVertexArr = new TClonesArray("AliAODVertex");
    Int_t nDauRejCount = 0;
    Int_t nVtx = fVtxTagger->FindVertices(jet,
                                         fTracksCont->GetArray(),
                                         (AliAODEvent*)InputEvent(),
                                         esdVtx,
                                         InputEvent()->GetMagneticField(),
                                         secVertexArr,
                                         0,
                                         arrDispersion,
                                         nDauRejCount);


    if(nVtx < 0)
    {
      secVertexArr->Clear();
      delete secVertexArr;
      return;
    }
    //###########################################################################
  }
  else // Load HF vertex branch from AOD event, if possible
  {
    secVertexArr = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("VerticesHF"));
    if(!secVertexArr)
      return;
  }

  // Loop over all potential secondary vertices
  for(Int_t i=0; i<secVertexArr->GetEntriesFast(); i++)
  {
    AliAODVertex* secVtx = (AliAODVertex*)(secVertexArr->UncheckedAt(i));
    if(!fCalculateSecondaryVertices)
      if((strcmp(secVtx->GetParent()->ClassName(), "AliAODRecoDecayHF3Prong")))
        continue;

    // Calculate vtx distance
    Double_t effX = secVtx->GetX() - esdVtx->GetX();
    Double_t effY = secVtx->GetY() - esdVtx->GetY();
    //Double_t effZ = secVtx->GetZ() - esdVtx->GetZ();

    // ##### Vertex properties
    // vertex dispersion
    Double_t dispersion = arrDispersion[i].first;

    // invariant mass
    Double_t mass = fVtxTagger->GetVertexInvariantMass(secVtx);

    // signed length
    Double_t Lxy  = TMath::Sqrt(effX*effX + effY*effY);
    Double_t jetP[3]; jet->PxPyPz(jetP);
    Double_t signLxy = effX * jetP[0] + effY * jetP[1];
    if (signLxy < 0.) Lxy *= -1.;

    Double_t sigmaLxy  = 0;
    AliAODVertex* aodVtx = (AliAODVertex*)(primVtx);
    if (aodVtx)
      sigmaLxy = aodVtx->ErrorDistanceXYToVertex(secVtx);

    // Add secondary vertices if they fulfill the conditions
    if( (dispersion > fSecondaryVertexMaxDispersion) || (TMath::Abs(secVtx->GetChi2perNDF()) > fSecondaryVertexMaxChi2) )
      continue;

    basicJet.AddSecondaryVertex(secVtx->GetX(), secVtx->GetY(), secVtx->GetZ(), TMath::Abs(secVtx->GetChi2perNDF()), dispersion, mass, Lxy, sigmaLxy);

    //cout << Form("(%3.3f, %3.3f, %3.3f), r=%3.3f, chi2=%3.3f, dispersion=%3.3f", effX, effY, effZ, vtxDistance, TMath::Abs(vtx->GetChi2perNDF()), dispersion) << endl;
  }

  if(fCalculateSecondaryVertices)
  {
    secVertexArr->Clear();
    delete secVertexArr;
  }
  delete esdVtx;
}

//________________________________________________________________________
void AliAnalysisTaskJetExtractorHF::GetTrackImpactParameters(const AliVVertex* vtx, AliAODTrack* track, Double_t& d0, Double_t& d0cov, Double_t& z0, Double_t& z0cov)
{
  if (track)
  {
    Double_t d0rphiz[2],covd0[3];
    Bool_t isDCA=track->PropagateToDCA(vtx,InputEvent()->GetMagneticField(),3.0,d0rphiz,covd0);
    if(isDCA)
    {
      d0 = d0rphiz[0];
      z0 = d0rphiz[1];
      d0cov = covd0[0];
      z0cov = covd0[2];
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskJetExtractorHF::CalculateJetProperties(AliEmcalJet* jet)
{
  CalculateJetType(jet, fCurrentJetTypeHM, fCurrentJetTypePM);

  // If fTruthJetsArrayName is set, the true pt field is calculated by searching the matching jet
  if(fTruthJetsArrayName != "")
    FindMatchingJet(jet);
  else // otherwise, the true jet pt is set as fraction of non-MC particles in jets
  {
    Double_t pt_nonMC = 0.;
    Double_t pt_all   = 0.;

    for(Int_t iConst = 0; iConst < jet->GetNumberOfTracks(); iConst++)
    {
      // Loop over all valid jet constituents
      AliVParticle* particle = static_cast<AliVParticle*>(jet->TrackAt(iConst, fTracksCont->GetArray()));
      if(!particle) continue;
      if(particle->Pt() < 1e-6) continue;

      // Particles marked w/ labels above 100000 are considered from toy
      if(particle->GetLabel() < 100000)
        pt_nonMC += particle->Pt();
      pt_all += particle->Pt();
    }
    if(pt_all)
      fCurrentTrueJetPt = (pt_nonMC/pt_all);
  }
}

//________________________________________________________________________
void AliAnalysisTaskJetExtractorHF::CalculateEventProperties()
{
  fCurrentNJetsInEvents = 0;
  fCurrentNJetsInEventsTotal = 0;
  GetLeadingJets("rho", fCurrentLeadingJet, fCurrentSubleadingJet);
  if(fCent==-1)
    fCent = 99;
}


//________________________________________________________________________
Bool_t AliAnalysisTaskJetExtractorHF::IsStrangeJet(AliEmcalJet* jet)
{
  // Do hadron matching jet type tagging using mcparticles
  // ... if not explicitly deactivated
  if (fTruthParticleArray)
  {
    for(Int_t i=0; i<fTruthParticleArray->GetEntries();i++)
    {
      AliAODMCParticle* part = (AliAODMCParticle*)fTruthParticleArray->At(i);
      if(!part) continue;

      // Check if the particle has strangeness
      Int_t absPDG = TMath::Abs(part->PdgCode());
      if ((absPDG > 300 && absPDG < 400) || (absPDG > 3000 && absPDG < 4000))
      {
        // Check if particle is in a radius around the jet
        Double_t rsquared = (part->Eta() - jet->Eta())*(part->Eta() - jet->Eta()) + (part->Phi() - jet->Phi())*(part->Phi() - jet->Phi());
        if(rsquared >= fHadronMatchingRadius*fHadronMatchingRadius)
          continue;
        else
          return kTRUE;
      }
    }
  }
  // As fallback, the MC stack will be tried
  else if(MCEvent() && (MCEvent()->Stack()))
  {
    AliStack* stack = MCEvent()->Stack();
    // Go through the whole particle stack
    for(Int_t i=0; i<stack->GetNtrack(); i++)
    {
      TParticle *part = stack->Particle(i);
      if(!part) continue;

      // Check if the particle has strangeness
      Int_t absPDG = TMath::Abs(part->GetPdgCode());
      if ((absPDG > 300 && absPDG < 400) || (absPDG > 3000 && absPDG < 4000))
      {
        // Check if particle is in a radius around the jet
        Double_t rsquared = (part->Eta() - jet->Eta())*(part->Eta() - jet->Eta()) + (part->Phi() - jet->Phi())*(part->Phi() - jet->Phi());
        if(rsquared >= fHadronMatchingRadius*fHadronMatchingRadius)
          continue;
        else
          return kTRUE;
      }
    }
  }
  return kFALSE;

}

//________________________________________________________________________
void AliAnalysisTaskJetExtractorHF::CalculateJetType(AliEmcalJet* jet, Int_t& typeHM, Int_t& typePM)
{
  Double_t radius = fHadronMatchingRadius;

  if(!fTruthParticleArray)
    return;

  typeHM = 0;
  typePM = 0;

  AliAODMCParticle* parton[2];
  parton[0] = (AliAODMCParticle*) fVtxTagger->IsMCJetParton(fTruthParticleArray, jet, radius);  // method 1 (parton matching)
  parton[1] = (AliAODMCParticle*) fVtxTagger->IsMCJetMeson(fTruthParticleArray, jet, radius);   // method 2 (hadron matching)

  if (parton[0]) {
    Int_t pdg = TMath::Abs(parton[0]->PdgCode());
    typePM = pdg;
  }

  if (!parton[1])
  {
    // No HF jet -- now also separate jets in udg (1) and s-jets (3) (on demand)
    if((std::find(fExtractionListPIDsHM.begin(), fExtractionListPIDsHM.end(), 3) != fExtractionListPIDsHM.end()) && IsStrangeJet(jet))
      typeHM = 3;
    else
      typeHM = 1;
  }
  else {
    Int_t pdg = TMath::Abs(parton[1]->PdgCode());
    if(fVtxTagger->IsDMeson(pdg)) typeHM = 4;
    else if (fVtxTagger->IsBMeson(pdg)) typeHM = 5;
    fCurrentJetHFHadron = pdg;

    if(typeHM == 5)
    {
      Int_t motherCode = parton[1]->GetMother();
      if(motherCode >= 0)
      {
        if(TMath::Abs( ((AliAODMCParticle*)fTruthParticleArray->At(motherCode))->PdgCode()) == 5)
          FillHistogram("hBJetFromInitialCollision", 1.5);
        else
          FillHistogram("hBJetFromInitialCollision", 0.5);
      }
      else
        FillHistogram("hBJetFromInitialCollision", 1.5);
    }
  }

  // Set flavour of AliEmcalJet object (set ith bit while i corresponds to type)
  if(fSetEmcalJetFlavour)
    jet->AddFlavourTag(static_cast<Int_t>(TMath::Power(2, typeHM)));
}


//________________________________________________________________________
void AliAnalysisTaskJetExtractorHF::FindMatchingJet(AliEmcalJet* jet)
{
  // "True" background
  AliRhoParameter* rho = static_cast<AliRhoParameter*>(InputEvent()->FindListObject(fTruthJetsRhoName.Data()));
  Double_t trueRho = rho->GetVal();

  TClonesArray* truthArray = static_cast<TClonesArray*>(InputEvent()->FindListObject(Form("%s", fTruthJetsArrayName.Data())));
  Double_t     bestMatchDeltaR = 999.;

  // Loop over all true jets to find the best match
  fCurrentTrueJetPt = 0;
  for(Int_t i=0; i<truthArray->GetEntries(); i++)
  {
    AliEmcalJet* truthJet = static_cast<AliEmcalJet*>(truthArray->At(i));
    if(truthJet->Pt() < 1.0)
      continue;

    Double_t deltaEta = (truthJet->Eta()-jet->Eta());
    Double_t deltaPhi = TMath::Min(TMath::Abs(truthJet->Phi()-jet->Phi()),TMath::TwoPi() - TMath::Abs(truthJet->Phi()-jet->Phi()));
    Double_t deltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);

    // Cut jets too far away
    if (deltaR > fJetsCont->GetJetRadius()*2)
      continue;

    // Search for the best match
    if(deltaR < bestMatchDeltaR)
    {
      bestMatchDeltaR = deltaR;
      fCurrentTrueJetPt = truthJet->Pt() - truthJet->Area()* trueRho;
    }
  }
}


//________________________________________________________________________
void AliAnalysisTaskJetExtractorHF::AddPIDInformation(AliVParticle* particle, AliBasicJetConstituent& constituent)
{
  // On demand and if we have AODs, add the particle PID signal object
  // Add the truth values if available
  AliAODTrack* aodtrack = dynamic_cast<AliAODTrack*>(particle);

  if(!fAddPIDSignal)
    return;
  if(!aodtrack)
    return;

  // Get AOD value from reco
  Int_t recoPID  = aodtrack->GetMostProbablePID();
  Int_t truthPID = 9;

  // Get truth values if we are on MC
  if(fTruthParticleArray)
  {
    for(Int_t i=0; i<fTruthParticleArray->GetEntries();i++)
    {
      AliAODMCParticle* mcParticle = dynamic_cast<AliAODMCParticle*>(fTruthParticleArray->At(i));
      if(!mcParticle) continue;

      if (mcParticle->GetLabel() == particle->GetLabel())
      {
        // Use same convention as PID in AODs
        if(TMath::Abs(mcParticle->PdgCode()) == 2212) // proton
          truthPID = 4;
        else if (TMath::Abs(mcParticle->PdgCode()) == 211) // pion
          truthPID = 2;
        else if (TMath::Abs(mcParticle->PdgCode()) == 321) // kaon
          truthPID = 3;
        else if (TMath::Abs(mcParticle->PdgCode()) == 11) // electron
          truthPID = 0;
        else if (TMath::Abs(mcParticle->PdgCode()) == 13) // muon
          truthPID = 1;
        else if (TMath::Abs(mcParticle->PdgCode()) == 700201) // deuteron
          truthPID = 5;
        else if (TMath::Abs(mcParticle->PdgCode()) == 700301) // triton
          truthPID = 6;
        else if (TMath::Abs(mcParticle->PdgCode()) == 700302) // He3
          truthPID = 7;
        else if (TMath::Abs(mcParticle->PdgCode()) == 700202) // alpha
          truthPID = 8;
        else
          truthPID = 9;

        break;
      }
    }
  }

  AliAODPid* pidObj = aodtrack->GetDetPid();
  constituent.SetPIDSignal(pidObj->GetITSsignal(), pidObj->GetTPCsignal(), pidObj->GetTOFsignal(), pidObj->GetTRDsignal(), truthPID, recoPID);
}

//________________________________________________________________________
inline Bool_t AliAnalysisTaskJetExtractorHF::IsTrackInCone(AliVParticle* track, Double_t eta, Double_t phi, Double_t radius)
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
void AliAnalysisTaskJetExtractorHF::GetLeadingJets(const char* opt, AliEmcalJet*& jetLeading, AliEmcalJet*& jetSubLeading)
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


//########################################################################
// HELPERS
//########################################################################


//________________________________________________________________________
inline void AliAnalysisTaskJetExtractorHF::FillHistogram(const char * key, Double_t x)
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
inline void AliAnalysisTaskJetExtractorHF::FillHistogram(const char * key, Double_t x, Double_t y)
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
inline void AliAnalysisTaskJetExtractorHF::FillHistogram(const char * key, Double_t x, Double_t y, Double_t add)
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
inline void AliAnalysisTaskJetExtractorHF::FillHistogram3D(const char * key, Double_t x, Double_t y, Double_t z, Double_t add)
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
template <class T> T* AliAnalysisTaskJetExtractorHF::AddHistogram1D(const char* name, const char* title, const char* options, Int_t xBins, Double_t xMin, Double_t xMax, const char* xTitle, const char* yTitle)
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
template <class T> T* AliAnalysisTaskJetExtractorHF::AddHistogram2D(const char* name, const char* title, const char* options, Int_t xBins, Double_t xMin, Double_t xMax, Int_t yBins, Double_t yMin, Double_t yMax, const char* xTitle, const char* yTitle, const char* zTitle)
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
template <class T> T* AliAnalysisTaskJetExtractorHF::AddHistogram3D(const char* name, const char* title, const char* options, Int_t xBins, Double_t xMin, Double_t xMax, Int_t yBins, Double_t yMin, Double_t yMax, Int_t zBins, Double_t zMin, Double_t zMax, const char* xTitle, const char* yTitle, const char* zTitle)
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
void AliAnalysisTaskJetExtractorHF::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

