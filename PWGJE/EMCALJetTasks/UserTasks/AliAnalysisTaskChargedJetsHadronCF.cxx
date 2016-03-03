// $Id$
//
// Jet+h correlation task
//
// Author: R. Haake

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliPicoTrack.h"

#include "AliAnalysisTaskChargedJetsHadronCF.h"

ClassImp(AliAnalysisTaskChargedJetsHadronCF)

//________________________________________________________________________
AliAnalysisTaskChargedJetsHadronCF::AliAnalysisTaskChargedJetsHadronCF() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskChargedJetsHadronCF", kTRUE),
  fJetsCont(0),
  fTracksCont(0),
  fNumberOfCentralityBins(10),
  fJetsOutput(),
  fJetParticleArrayName("JetsDPhiBasicParticles"),
  fEventCriteriumMode(0),
  fEventCriteriumMinBackground(0),
  fEventCriteriumMaxBackground(0),
  fEventCriteriumMinLeadingJetPt(0),
  fEventCriteriumMinSubleadingJetPt(0)
{
  // Default constructor.
  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskChargedJetsHadronCF::AliAnalysisTaskChargedJetsHadronCF(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fJetsCont(0),
  fTracksCont(0),
  fNumberOfCentralityBins(10),
  fJetsOutput(),
  fJetParticleArrayName("JetsDPhiBasicParticles"),
  fEventCriteriumMode(0),
  fEventCriteriumMinBackground(0),
  fEventCriteriumMaxBackground(0),
  fEventCriteriumMinLeadingJetPt(0),
  fEventCriteriumMinSubleadingJetPt(0)
{
  SetMakeGeneralHistograms(kTRUE);
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
  fJetsCont->PrintCuts();
  if(fJetsCont) { //get particles connected to jets
    fTracksCont       = fJetsCont->GetParticleContainer();
  } else {        //no jets, just analysis tracks
    fTracksCont       = GetParticleContainer(0);
  }
  if(fTracksCont) fTracksCont->SetClassName("AliAODTrack");

  // ### Create all histograms

  // Change the event rejection histogram -> Add a custom value
  fHistEventRejection->GetXaxis()->SetBinLabel(14,"JetCrit");

  // Track QA plots
  AddHistogram2D<TH2D>("hTrackCount", "Number of tracks in acceptance vs. centrality", "LEGO2", 500, 0., 5000., fNumberOfCentralityBins, 0, 100, "N tracks","Centrality", "dN^{Events}/dN^{Tracks}");
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
  AddHistogram2D<TH2D>("hJetCount", "Number of jets in acceptance vs. centrality", "LEGO2", 100, 0., 100., fNumberOfCentralityBins, 0, 100, "N Jets","Centrality", "dN^{Events}/dN^{Jets}");

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


  AddHistogram2D<TH2D>("hBackgroundPt", "Background p_{T} distribution", "", 1000, 0., 50., fNumberOfCentralityBins, 0, 100, "Background p_{T} (GeV/c)", "Centrality", "dN^{Events}/dp_{T}");


  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}

//________________________________________________________________________
Bool_t AliAnalysisTaskChargedJetsHadronCF::FillHistograms()
{
  // Fill histograms.

  // ################# TRACK PLOTS
  if (fTracksCont)
  {
    FillHistogram("hTrackCount", fTracksCont->GetNAcceptedParticles(), fCent); 
    fTracksCont->ResetCurrentID();
    AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());

    // All tracks plots
    while(track) {
      FillHistogram("hTrackPt", track->Pt(), fCent); 
      FillHistogram("hTrackPhi", track->Phi(), fCent); 
      FillHistogram("hTrackEta", track->Eta(), fCent); 
      FillHistogram("hTrackEtaPt", track->Eta(), track->Pt()); 
      FillHistogram("hTrackPhiPt", track->Phi(), track->Pt()); 
      FillHistogram("hTrackPhiEta", track->Phi(), track->Eta()); 

      track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle());
    }

    // Leading track plots
    track = static_cast<AliVTrack*>(fTracksCont->GetLeadingParticle());
    if(track)
    {
      FillHistogram("hLeadingTrackPt", track->Pt(), fCent); 
      FillHistogram("hLeadingTrackPhi", track->Phi(), fCent); 
      FillHistogram("hLeadingTrackEta", track->Eta(), fCent); 
      FillHistogram("hLeadingTrackPhiEta", track->Phi(), track->Eta()); 
    }
  }

  // ################# JET PLOTS

  if (fJetsCont)
  {
    FillHistogram("hBackgroundPt", fJetsCont->GetRhoVal(), fCent);
    Int_t count = 0;
    fJetsCont->ResetCurrentID();
    AliEmcalJet *jet = fJetsCont->GetNextAcceptJet(); 
    // All jets
    while(jet) {
      count++;
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

      jet = fJetsCont->GetNextAcceptJet(); 
    }

    // Leading jet plots
    jet = fJetsCont->GetLeadingJet("rho"); // Get the leading bgrd-corrected jet
    if(jet)
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
    jet = GetSubleadingJet("rho"); // Get the subleading bgrd-corrected jet
    if(jet)
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

    // Jet count
    FillHistogram("hJetCount", count, fCent);
  }

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskChargedJetsHadronCF::ExecOnce() {

  AliAnalysisTaskEmcalJet::ExecOnce();
  if (fJetsCont && fJetsCont->GetArray() == 0) fJetsCont = 0;
  if (fTracksCont && fTracksCont->GetArray() == 0) fTracksCont = 0;


  // ### Add the jets as basic correlation particles to the event
  if (!(fInputEvent->FindListObject(Form("%s_%s", GetName(), fJetParticleArrayName.Data()))))
  {
    fJetsOutput = new TClonesArray("AliPicoTrack");
    fJetsOutput->SetName(fJetParticleArrayName.Data());
    fInputEvent->AddObject(fJetsOutput);
  }
  else
    AliFatal(Form("%s: Object with name %s already in event!", GetName(), Form("%s_%s", GetName(), fJetParticleArrayName.Data())));

}

//________________________________________________________________________
Bool_t AliAnalysisTaskChargedJetsHadronCF::Run()
{
  // Note: Normal event, track, and jet selection is done in the framework and not here
  if(fJetsOutput && fJetsCont)
  {
    // Clear the array containing the correlation-like objects. Should be done in Event::Reset()
    fJetsOutput->Delete();

    // ##################
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
      AliEmcalJet* jet = fJetsCont->GetLeadingJet("rho"); // Get the leading jet (w/ bgrd-corrected)
      if(jet)
      {
        if(jet->Pt() - fJetsCont->GetRhoVal()*jet->Area() <= fEventCriteriumMinLeadingJetPt)
        {
          fHistEventRejection->Fill("JetCrit", 1);
          return kFALSE;
        }
      }
    }
    else if(fEventCriteriumMode==3) // Simple dijet trigger
    {
      AliEmcalJet* leadingJet    = fJetsCont->GetLeadingJet("rho");
      AliEmcalJet* subleadingJet = GetSubleadingJet("rho");

      if(leadingJet && subleadingJet)
      {
        if((leadingJet->Pt() - fJetsCont->GetRhoVal()*leadingJet->Area() <= fEventCriteriumMinLeadingJetPt) || (subleadingJet->Pt() - fJetsCont->GetRhoVal()*subleadingJet->Area() <= fEventCriteriumMinSubleadingJetPt))
        {
          fHistEventRejection->Fill("JetCrit", 1);
          return kFALSE;
      }
      }
    }

    // ##################
    // Create correlation-like objects and save them in a TClonesArray
    fJetsCont->ResetCurrentID();
    AliEmcalJet *jet = fJetsCont->GetNextAcceptJet();
    Int_t count = 0;
    while(jet) {

      // Add the particle object to the clones array
      new ((*fJetsOutput)[count]) AliPicoTrack(jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), jet->Eta(), jet->Phi(), jet->Charge(), 0, 0); // only Pt,Eta,Phi are interesting for correlations;

      jet = fJetsCont->GetNextAcceptJet(); 
      count++;
    }
    return kTRUE;
  }
  else
    return kFALSE;
}

//########################################################################
// HELPERS
//########################################################################

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
void AliAnalysisTaskChargedJetsHadronCF::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

