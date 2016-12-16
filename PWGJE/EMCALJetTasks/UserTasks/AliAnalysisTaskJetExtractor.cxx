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

//________________________________________________________________________
AliBasicJet::~AliBasicJet() 
{
// dummy destructor
}

//________________________________________________________________________
AliBasicJetConstituent::~AliBasicJetConstituent() 
{
// dummy destructor
}

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskJetExtractor)
/// \endcond
//________________________________________________________________________
AliAnalysisTaskJetExtractor::AliAnalysisTaskJetExtractor() :
  AliAnalysisTaskEmcalJet("AliAnalysisTaskJetExtractor", kTRUE),
  fJetsCont(0),
  fTracksCont(0),
  fJetsTree(0),
  fJetsTreeBuffer(0),
  fRandom(0),
  fCurrentNJetsInEvents(0),
  fCurrentJetTypeHM(0),
  fCurrentJetTypeIC(0),
  fCurrentLeadingJet(0),
  fCurrentSubleadingJet(0),
  fCurrentInitialParton1(0),
  fCurrentInitialParton2(0),
  fExtractionCutMinCent(-1),
  fExtractionCutMaxCent(-1),
  fExtractionCutUseIC(kFALSE),
  fExtractionCutUseHM(kFALSE),
  fExtractionCutMinPt(0.),
  fExtractionCutMaxPt(200.),
  fExtractionPercentage(1.0),
  fExtractionListPIDsHM(),
  fHadronMatchingRadius(0.5),
  fInitialCollisionMatchingRadius(0.3),
  fTruthParticleArrayName("mcparticles")
{
  // Default constructor.
  SetMakeGeneralHistograms(kTRUE);
  fRandom = new TRandom3(0);
}

//________________________________________________________________________
AliAnalysisTaskJetExtractor::AliAnalysisTaskJetExtractor(const char *name) :
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fJetsCont(0),
  fTracksCont(0),
  fJetsTree(0),
  fJetsTreeBuffer(0),
  fRandom(0),
  fCurrentNJetsInEvents(0),
  fCurrentJetTypeHM(0),
  fCurrentJetTypeIC(0),
  fCurrentLeadingJet(0),
  fCurrentSubleadingJet(0),
  fCurrentInitialParton1(0),
  fCurrentInitialParton2(0),
  fExtractionCutMinCent(-1),
  fExtractionCutMaxCent(-1),
  fExtractionCutUseIC(kFALSE),
  fExtractionCutUseHM(kFALSE),
  fExtractionCutMinPt(0.),
  fExtractionCutMaxPt(200.),
  fExtractionPercentage(1.0),
  fExtractionListPIDsHM(),
  fHadronMatchingRadius(0.5),
  fInitialCollisionMatchingRadius(0.3),
  fTruthParticleArrayName("mcparticles")
{
  // Default constructor.
  SetMakeGeneralHistograms(kTRUE);
  fRandom = new TRandom3(0);
}

//________________________________________________________________________
AliAnalysisTaskJetExtractor::~AliAnalysisTaskJetExtractor()
{
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskJetExtractor::UserCreateOutputObjects()
{
  std::cout << "########################################\n";
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

  // ### Add control histograms (already some create in base task)
  AddHistogram2D<TH2D>("hJetCount", "Number of jets in acceptance vs. centrality", "COLZ", 100, 0., 100., 100, 0, 100, "N Jets","Centrality", "dN^{Events}/dN^{Jets}");
  AddHistogram2D<TH2D>("hTrackCount", "Number of tracks in acceptance vs. centrality", "COLZ", 500, 0., 5000., 100, 0, 100, "N tracks","Centrality", "dN^{Events}/dN^{Tracks}");
  AddHistogram2D<TH2D>("hBackgroundPt", "Background p_{T} distribution", "", 150, 0., 150., 100, 0, 100, "Background p_{T} (GeV/c)", "Centrality", "dN^{Events}/dp_{T}");
  AddHistogram2D<TH2D>("hJetPtRaw", "Jets p_{T} distribution (raw)", "COLZ", 300, 0., 300., 100, 0, 100, "p_{T, jet} (GeV/c)", "Centrality", "dN^{Jets}/dp_{T}");
  AddHistogram2D<TH2D>("hJetPt", "Jets p_{T} distribution (background subtracted)", "COLZ", 400, -100., 300., 100, 0, 100, "p_{T, jet} (GeV/c)", "Centrality", "dN^{Jets}/dp_{T}");
  AddHistogram2D<TH2D>("hJetPhiEta", "Jet angular distribution #phi/#eta", "COLZ", 180, 0., 2*TMath::Pi(), 100, -2.5, 2.5, "#phi", "#eta", "dN^{Jets}/d#phi d#eta");
  AddHistogram2D<TH2D>("hJetArea", "Jet area", "COLZ", 200, 0., 2., 100, 0, 100, "Jet A", "Centrality", "dN^{Jets}/dA");

  TH1* tmpHisto = AddHistogram1D<TH1D>("hJetAcceptance", "Accepted jets", "", 5, 0, 5, "stage","N^{jets}/cut");
  tmpHisto->GetXaxis()->SetBinLabel(1, "Before cuts");
  tmpHisto->GetXaxis()->SetBinLabel(2, "Centrality");
  tmpHisto->GetXaxis()->SetBinLabel(3, "p_{T}");
  tmpHisto->GetXaxis()->SetBinLabel(4, "PID");
  tmpHisto->GetXaxis()->SetBinLabel(5, "Extraction percentage");

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}


//________________________________________________________________________
void AliAnalysisTaskJetExtractor::ExecOnce() {

  AliAnalysisTaskEmcalJet::ExecOnce();

  // ### Prepare the jet tree
  fJetsTree = new TTree("ExtractedJets", "ExtractedJets");
  fJetsTree->Branch("Jets", "AliBasicJet", &fJetsTreeBuffer, 1000);
  fOutput->Add(fJetsTree);
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetExtractor::IsJetSelected(AliEmcalJet* jet)
{
  // ######## Check if cuts for extraction have been passed

  // Extraction pt, jet PID (according to initial parton, according to HM), centrality
  FillHistogram("hJetAcceptance", 0.5);

  // ### CENTRALITY
  if(fExtractionCutMinCent != -1 && (fCent < fExtractionCutMinCent || fCent >= fExtractionCutMaxCent))
    return kFALSE;
  FillHistogram("hJetAcceptance", 1.5);

  // ### PT
  if( ((jet->Pt()-jet->Area()*fJetsCont->GetRhoVal()) < fExtractionCutMinPt) || ((jet->Pt()-jet->Area()*fJetsCont->GetRhoVal()) >= fExtractionCutMaxPt) )
    return kFALSE;
  FillHistogram("hJetAcceptance", 2.5);

  Bool_t passedCutPID = kTRUE;
  // ### PID (from initial collision)
  if(fExtractionCutUseIC)
  {
    passedCutPID = kFALSE;
    for(Int_t i=0; i<fExtractionListPIDsIC.size(); i++)
    {
      if (fExtractionListPIDsIC.at(i) == fCurrentJetTypeIC)
        passedCutPID = kTRUE;
    }
  }
  // ### PID (from hadron matching)
  else if(fExtractionCutUseHM)
  {
    passedCutPID = kFALSE;
    for(Int_t i=0; i<fExtractionListPIDsHM.size(); i++)
    {
      if (fExtractionListPIDsHM.at(i) == fCurrentJetTypeHM)
        passedCutPID = kTRUE;
    }
  }

  // Note: When either of the cuts is passed, accepted it
  if(!passedCutPID)
    return kFALSE;
  FillHistogram("hJetAcceptance", 3.5);

  // Discard jets statistically
  if(fRandom->Rndm() >= fExtractionPercentage)
    return kFALSE;
  FillHistogram("hJetAcceptance", 4.5);

  fCurrentNJetsInEvents++;
  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskJetExtractor::AddJetToTree(AliEmcalJet* jet)
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
  AliBasicJet basicJet(jet->Eta(), jet->Phi(), jet->Pt(), jet->Charge(), fJetsCont->GetJetRadius(), jet->Area(), fCurrentJetTypeIC, fCurrentJetTypeHM, fJetsCont->GetRhoVal(), InputEvent()->GetMagneticField(), vtxX, vtxY, vtxZ, eventID, fCent);

  // Add constituents
  for(Int_t i = 0; i < jet->GetNumberOfTracks(); i++)
  {
    AliVParticle* particle = static_cast<AliVParticle*>(jet->TrackAt(i, fTracksCont->GetArray()));
    if(!particle) continue;

    // Secondary vertex analysis
    Double_t z = GetTrackImpactParameter(myVertex, dynamic_cast<AliAODTrack*>(particle));

    Int_t constid = GetTrackPID(particle);
    basicJet.AddJetConstituent(particle->Eta(), particle->Phi(), particle->Pt(), particle->Charge(), constid, particle->Xv(), particle->Yv(), particle->Zv(), z);
  }

  // Field currently not used
  basicJet.SetTruePt(0);

  fJetsTreeBuffer = &basicJet;
  fJetsTree->Fill();
}

//________________________________________________________________________
void AliAnalysisTaskJetExtractor::FillEventControlHistograms()
{
  // ### Event control plots
  FillHistogram("hTrackCount", fTracksCont->GetNAcceptedParticles(), fCent);
  FillHistogram("hJetCount", fCurrentNJetsInEvents, fCent);
  FillHistogram("hBackgroundPt", fJetsCont->GetRhoVal(), fCent);
}

//________________________________________________________________________
void AliAnalysisTaskJetExtractor::FillJetControlHistograms(AliEmcalJet* jet)
{
  // ### Jet control plots
  FillHistogram("hJetPtRaw", jet->Pt(), fCent); 
  FillHistogram("hJetPt", jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), fCent);
  FillHistogram("hJetPhiEta", jet->Phi(), jet->Eta());
  FillHistogram("hJetArea", jet->Area(), fCent);
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetExtractor::Run()
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

  FillEventControlHistograms();

  return kTRUE;
}

//________________________________________________________________________
std::vector<std::vector<Float_t>> AliAnalysisTaskJetExtractor::GetSecondaryVertices(const AliVVertex* vtx, AliEmcalJet* jet)
{
  std::vector<std::vector<Float_t>> dummy;
  return dummy;
}


//________________________________________________________________________
Double_t AliAnalysisTaskJetExtractor::GetTrackImpactParameter(const AliVVertex* vtx, AliAODTrack* track)
{
  Double_t z = 0;
  // Impact parameter z
  if (track)
  {
    Double_t d0rphiz[2],covd0[3];
    Bool_t isDCA=track->PropagateToDCA(vtx,InputEvent()->GetMagneticField(),9999.,d0rphiz,covd0);
    if(isDCA)
      z = d0rphiz[1]; // impact parameter z
  }

  return z;
}

//________________________________________________________________________
void AliAnalysisTaskJetExtractor::CalculateJetProperties(AliEmcalJet* jet)
{
  CalculateJetType(jet, fCurrentJetTypeIC, fCurrentJetTypeHM);
}

//________________________________________________________________________
void AliAnalysisTaskJetExtractor::CalculateEventProperties()
{
  fCurrentNJetsInEvents = 0;
  GetLeadingJets("rho", fCurrentLeadingJet, fCurrentSubleadingJet);
  if(fPythiaInfo)
    CalculatePYTHIAInitialCollisionJets();
}

//________________________________________________________________________
void AliAnalysisTaskJetExtractor::CalculateJetType(AliEmcalJet* jet, Int_t& typeIC, Int_t& typeHM)
{
  typeIC = -1;
  typeHM = -1;

  // Only do this in case we have PYTHIA information available
  if(fPythiaInfo)
  {
    typeIC = 0;
    if(jet==fCurrentInitialParton1)
      typeIC = fPythiaInfo->GetPartonFlag6();
    else if (jet==fCurrentInitialParton2)
      typeIC = fPythiaInfo->GetPartonFlag7();
  }

  // Do hadron matching jet type tagging using mcparticles
  // ... if not explicitly deactivated
  if (fTruthParticleArrayName != "")
  {
    TClonesArray* fTruthParticleArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTruthParticleArrayName.Data()));
    if(!fTruthParticleArray)
      return;

    typeHM = 0;
    for(Int_t i=0; i<fTruthParticleArray->GetEntries();i++)
    {
      AliAODMCParticle* part = dynamic_cast<AliAODMCParticle*>(fTruthParticleArray->At(i));
      if(!part) continue;

      // Check if particle is in a radius around the jet
      Double_t rsquared = (part->Eta() - jet->Eta())*(part->Eta() - jet->Eta()) + (part->Phi() - jet->Phi())*(part->Phi() - jet->Phi());
      if(rsquared >= fHadronMatchingRadius*fHadronMatchingRadius)
        continue;

      // Check if the particle has beauty, charm or strangeness
      // If it has beauty, we are done (exclusive definition)
      Int_t absPDG = TMath::Abs(part->PdgCode());
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
  // As fallback, the MC stack will be tried
  else if(MCEvent() && (MCEvent()->Stack()))
  {
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
        continue;

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
}

//________________________________________________________________________
Int_t AliAnalysisTaskJetExtractor::GetTrackPID(AliVParticle* particle)
{
  AliAODTrack*  aodtrack = static_cast<AliAODTrack*>(particle);
  Int_t constid = 9; // 9 means unknown

  if(fPythiaInfo || MCEvent()) // in case of PYTHIA or MC, use PDG codes for PID
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
  return constid;
}

//________________________________________________________________________
void AliAnalysisTaskJetExtractor::GetLeadingJets(const char* opt, AliEmcalJet*& jetLeading, AliEmcalJet*& jetSubLeading)
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
void AliAnalysisTaskJetExtractor::CalculatePYTHIAInitialCollisionJets()
{
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
      fCurrentInitialParton1 = jet;
    }
    if(deltaR2 < bestMatchDeltaR2)
    {
      bestMatchDeltaR2 = deltaR2;
      fCurrentInitialParton2 = jet;
    }
  }

  if(bestMatchDeltaR1 > fInitialCollisionMatchingRadius)
    fCurrentInitialParton1 = 0;
  if(bestMatchDeltaR2 > fInitialCollisionMatchingRadius)
    fCurrentInitialParton2 = 0;
}



//########################################################################
// HELPERS
//########################################################################


//________________________________________________________________________
inline void AliAnalysisTaskJetExtractor::FillHistogram(const char * key, Double_t x)
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
inline void AliAnalysisTaskJetExtractor::FillHistogram(const char * key, Double_t x, Double_t y)
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
inline void AliAnalysisTaskJetExtractor::FillHistogram(const char * key, Double_t x, Double_t y, Double_t add)
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
inline void AliAnalysisTaskJetExtractor::FillHistogram3D(const char * key, Double_t x, Double_t y, Double_t z, Double_t add)
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
template <class T> T* AliAnalysisTaskJetExtractor::AddHistogram1D(const char* name, const char* title, const char* options, Int_t xBins, Double_t xMin, Double_t xMax, const char* xTitle, const char* yTitle)
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
template <class T> T* AliAnalysisTaskJetExtractor::AddHistogram2D(const char* name, const char* title, const char* options, Int_t xBins, Double_t xMin, Double_t xMax, Int_t yBins, Double_t yMin, Double_t yMax, const char* xTitle, const char* yTitle, const char* zTitle)
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
template <class T> T* AliAnalysisTaskJetExtractor::AddHistogram3D(const char* name, const char* title, const char* options, Int_t xBins, Double_t xMin, Double_t xMax, Int_t yBins, Double_t yMin, Double_t yMax, Int_t zBins, Double_t zMin, Double_t zMax, const char* xTitle, const char* yTitle, const char* zTitle)
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
void AliAnalysisTaskJetExtractor::Terminate(Option_t *) 
{
  // Called once at the end of the analysis.
}

