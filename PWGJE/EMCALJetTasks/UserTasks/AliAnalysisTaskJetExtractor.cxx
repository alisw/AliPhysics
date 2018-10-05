/**************************************************************************
 * Copyright(c) 1998-2018, ALICE Experiment at CERN, All rights reserved. *
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
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TTree.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TNamed.h>

#include "AliMCEvent.h"
#include "AliESDVertex.h"
#include "AliAODVertex.h"
#include "AliAODPid.h"

#include "AliRhoParameter.h"

#include "AliVTrack.h"
#include "AliVHeader.h"
#include "AliEmcalJet.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliTrackContainer.h"
#include "AliAODTrack.h"
#include "AliVParticle.h"
#include "TRandom3.h"
#include "AliEmcalPythiaInfo.h"
#include "AliAnalysisTaskEmcalJet.h"

#include "AliGenHijingEventHeader.h"
#include "AliHFJetsTaggingVertex.h"


#include "AliAnalysisTaskJetExtractor.h"

/// \cond CLASSIMP
ClassImp(AliEmcalJetTree)
/// \endcond

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskJetExtractor)
/// \endcond

//________________________________________________________________________
AliEmcalJetTree::AliEmcalJetTree() : TNamed("CustomTree", "CustomTree"), fJetTree(0), fInitialized(0), fSaveEventProperties(0), fSaveConstituents(0), fSaveConstituentsIP(0), fSaveConstituentPID(0), fSaveJetShapes(0), fSaveMCInformation(0), fSaveSecondaryVertices(0), fSaveTriggerTracks(0), fExtractionPercentages(), fExtractionPercentagePtBins(), fExtractionJetTypes_HM(), fExtractionJetTypes_PM()
{
  // For these arrays, we need to reserve memory
  fBuffer_Const_Pt         = new Float_t[kMaxNumConstituents];
  fBuffer_Const_Eta        = new Float_t[kMaxNumConstituents];
  fBuffer_Const_Phi        = new Float_t[kMaxNumConstituents];
  fBuffer_Const_Charge     = new Float_t[kMaxNumConstituents];
  fBuffer_Const_ProdVtx_X  = new Float_t[kMaxNumConstituents];
  fBuffer_Const_ProdVtx_Y  = new Float_t[kMaxNumConstituents];
  fBuffer_Const_ProdVtx_Z  = new Float_t[kMaxNumConstituents];

}

//________________________________________________________________________
AliEmcalJetTree::AliEmcalJetTree(const char* name) : TNamed(name, name), fJetTree(0), fInitialized(0), fSaveEventProperties(0), fSaveConstituents(0), fSaveConstituentsIP(0), fSaveConstituentPID(0), fSaveJetShapes(0), fSaveMCInformation(0), fSaveSecondaryVertices(0), fSaveTriggerTracks(0), fExtractionPercentages(), fExtractionPercentagePtBins(), fExtractionJetTypes_HM(), fExtractionJetTypes_PM()
{
  // For these arrays, we need to reserve memory
  fBuffer_Const_Pt         = new Float_t[kMaxNumConstituents];
  fBuffer_Const_Eta        = new Float_t[kMaxNumConstituents];
  fBuffer_Const_Phi        = new Float_t[kMaxNumConstituents];
  fBuffer_Const_Charge     = new Float_t[kMaxNumConstituents];
  fBuffer_Const_ProdVtx_X  = new Float_t[kMaxNumConstituents];
  fBuffer_Const_ProdVtx_Y  = new Float_t[kMaxNumConstituents];
  fBuffer_Const_ProdVtx_Z  = new Float_t[kMaxNumConstituents];
}


//________________________________________________________________________
Bool_t AliEmcalJetTree::AddJetToTree(AliEmcalJet* jet, Float_t bgrdDensity, Float_t vertexX, Float_t vertexY, Float_t vertexZ, Float_t centrality, Long64_t eventID, Float_t magField,
  AliTrackContainer* trackCont, Int_t motherParton, Int_t motherHadron, Int_t partonInitialCollision, Float_t matchedPt, Float_t truePtFraction, Float_t ptHard,
  Float_t* trackPID_ITS, Float_t* trackPID_TPC, Float_t* trackPID_TOF, Float_t* trackPID_TRD, Short_t* trackPID_Reco, Short_t* trackPID_Truth,
  Int_t numTriggerTracks, Float_t* triggerTrackPt, Float_t* triggerTrackDeltaEta, Float_t* triggerTrackDeltaPhi,
  Float_t* trackIP_d0, Float_t* trackIP_z0, Float_t* trackIP_d0cov, Float_t* trackIP_z0cov,
  Int_t numSecVertices, Float_t* secVtx_X, Float_t* secVtx_Y, Float_t* secVtx_Z, Float_t* secVtx_Mass, Float_t* secVtx_Lxy, Float_t* secVtx_SigmaLxy, Float_t* secVtx_Chi2, Float_t* secVtx_Dispersion)
{
  // ############ 
  // Approach: Fill buffer and then add to tree
  //
  if(!fInitialized)
    AliFatal("Tree is not initialized.");

  fBuffer_JetPt                                   = jet->Pt() - bgrdDensity*jet->Area();

  // Check if jet type is contained in extraction list
  if( (fExtractionJetTypes_PM.size() || fExtractionJetTypes_HM.size()) &&
      (std::find(fExtractionJetTypes_PM.begin(), fExtractionJetTypes_PM.end(), motherParton) == fExtractionJetTypes_PM.end()) &&
      (std::find(fExtractionJetTypes_HM.begin(), fExtractionJetTypes_HM.end(), motherHadron) == fExtractionJetTypes_HM.end()) )
    return kFALSE;

  // Check acceptance percentage for the given jet and discard statistically on demand
  Bool_t inPtRange = kFALSE;
  for(size_t i = 0; i<fExtractionPercentages.size(); i++)
  {
    if(fBuffer_JetPt>=fExtractionPercentagePtBins[i*2] && fBuffer_JetPt<fExtractionPercentagePtBins[i*2+1])
    {
      inPtRange = kTRUE;
      if(gRandom->Rndm() >= fExtractionPercentages[i])
        return kFALSE;
      break;
    }
  }
  if(!inPtRange && fExtractionPercentagePtBins.size()) // only discard jet if fExtractionPercentagePtBins was given
    return kFALSE;

  // Set basic properties
  fBuffer_JetEta                                  = jet->Eta();
  fBuffer_JetPhi                                  = jet->Phi();
  fBuffer_JetArea                                 = jet->Area();

  // Set event properties
  if(fSaveEventProperties)
  {
    fBuffer_Event_BackgroundDensity               = bgrdDensity;
    fBuffer_Event_Vertex_X                        = vertexX;
    fBuffer_Event_Vertex_Y                        = vertexY;
    fBuffer_Event_Vertex_Z                        = vertexZ;
    fBuffer_Event_Centrality                      = centrality;
    fBuffer_Event_ID                              = eventID;
    fBuffer_Event_MagneticField                   = magField;
    if(fSaveMCInformation)
      fBuffer_Event_PtHard                          = ptHard;
  }

  // Extract basic constituent properties directly from AliEmcalJet object
  fBuffer_NumConstituents = 0;
  if(trackCont && (fSaveConstituents || fSaveConstituentsIP))
    for(Int_t i = 0; i < jet->GetNumberOfTracks(); i++)
    {
      AliVParticle* particle = static_cast<AliVParticle*>(jet->TrackAt(i, trackCont->GetArray()));
      if(!particle) continue;

      if(fSaveConstituents)
      {
        fBuffer_Const_Pt[fBuffer_NumConstituents] = particle->Pt();
        fBuffer_Const_Eta[fBuffer_NumConstituents] = particle->Eta();
        fBuffer_Const_Phi[fBuffer_NumConstituents] = particle->Phi();
        fBuffer_Const_Charge[fBuffer_NumConstituents] = particle->Charge();
      }
      if(fSaveConstituentsIP)
      {
        fBuffer_Const_ProdVtx_X[fBuffer_NumConstituents] = particle->Xv();
        fBuffer_Const_ProdVtx_Y[fBuffer_NumConstituents] = particle->Yv();
        fBuffer_Const_ProdVtx_Z[fBuffer_NumConstituents] = particle->Zv();
      }
      fBuffer_NumConstituents++;
    }


  // Set constituent arrays for impact parameters
  if(fSaveConstituentsIP)
  {
    fJetTree->SetBranchAddress("Jet_Const_CovIPd", trackIP_d0cov);
    fJetTree->SetBranchAddress("Jet_Const_CovIPz", trackIP_z0cov);
    fJetTree->SetBranchAddress("Jet_Const_IPd", trackIP_d0);
    fJetTree->SetBranchAddress("Jet_Const_IPz", trackIP_z0);
  }
  // Set constituent arrays for PID parameters
  if(fSaveConstituentPID)
  {
    fJetTree->SetBranchAddress("Jet_Const_PID_ITS", trackPID_ITS);
    fJetTree->SetBranchAddress("Jet_Const_PID_TPC", trackPID_TPC);
    fJetTree->SetBranchAddress("Jet_Const_PID_TOF", trackPID_TOF);
    fJetTree->SetBranchAddress("Jet_Const_PID_TRD", trackPID_TRD);
    fJetTree->SetBranchAddress("Jet_Const_PID_Reconstructed", trackPID_Reco);
    if(fSaveMCInformation)
      fJetTree->SetBranchAddress("Jet_Const_PID_Truth", trackPID_Truth);
  }

  // Set secondary vertices arrays
  if(fSaveSecondaryVertices)
  {
    fBuffer_NumSecVertices = numSecVertices;
    fJetTree->SetBranchAddress("Jet_SecVtx_X", secVtx_X);
    fJetTree->SetBranchAddress("Jet_SecVtx_Y", secVtx_Y);
    fJetTree->SetBranchAddress("Jet_SecVtx_Z", secVtx_Z);
    fJetTree->SetBranchAddress("Jet_SecVtx_Mass", secVtx_Mass);
    fJetTree->SetBranchAddress("Jet_SecVtx_Lxy", secVtx_Lxy);
    fJetTree->SetBranchAddress("Jet_SecVtx_SigmaLxy", secVtx_SigmaLxy);
    fJetTree->SetBranchAddress("Jet_SecVtx_Chi2", secVtx_Chi2);
    fJetTree->SetBranchAddress("Jet_SecVtx_Dispersion", secVtx_Dispersion);
  }

  // Set jet shape observables
  fBuffer_Shape_Mass  = jet->M();

  // Set Monte Carlo information
  if(fSaveMCInformation)
  {
    fBuffer_Jet_MC_MotherParton = motherParton;
    fBuffer_Jet_MC_MotherHadron = motherHadron;
    fBuffer_Jet_MC_MotherIC = partonInitialCollision;
    fBuffer_Jet_MC_MatchedPt = matchedPt;
    fBuffer_Jet_MC_TruePtFraction = truePtFraction;
  }

  if(fSaveTriggerTracks)
  {
    fBuffer_NumTriggerTracks = numTriggerTracks;
    fJetTree->SetBranchAddress("Jet_TriggerTrack_Pt", triggerTrackPt);
    fJetTree->SetBranchAddress("Jet_TriggerTrack_dEta", triggerTrackDeltaEta);
    fJetTree->SetBranchAddress("Jet_TriggerTrack_dPhi", triggerTrackDeltaPhi);
  }

  // Add buffered jet to tree
  fJetTree->Fill();

  return kTRUE;
}


//________________________________________________________________________
void AliEmcalJetTree::InitializeTree()
{
  // Create the tree with active branches
  
  void* dummy = 0; // for some branches, do not need a buffer pointer yet
  

  if(fInitialized)
    AliFatal("Tree is already initialized.");

  // ### Prepare the jet tree
  fJetTree = new TTree(Form("JetTree_%s", GetName()), "");

  fJetTree->Branch("Jet_Pt",&fBuffer_JetPt,"Jet_Pt/F");
  fJetTree->Branch("Jet_Phi",&fBuffer_JetPhi,"Jet_Phi/F");
  fJetTree->Branch("Jet_Eta",&fBuffer_JetEta,"Jet_Eta/F");
  fJetTree->Branch("Jet_Area",&fBuffer_JetArea,"Jet_Area/F");
  fJetTree->Branch("Jet_NumConstituents",&fBuffer_NumConstituents,"Jet_NumConstituents/I");

  if(fSaveEventProperties)
  {
    fJetTree->Branch("Event_BackgroundDensity",&fBuffer_Event_BackgroundDensity,"Event_BackgroundDensity/F");
    fJetTree->Branch("Event_Vertex_X",&fBuffer_Event_Vertex_X,"Event_Vertex_X/F");
    fJetTree->Branch("Event_Vertex_Y",&fBuffer_Event_Vertex_Y,"Event_Vertex_Y/F");
    fJetTree->Branch("Event_Vertex_Z",&fBuffer_Event_Vertex_Z,"Event_Vertex_Z/F");
    fJetTree->Branch("Event_Centrality",&fBuffer_Event_Centrality,"Event_Centrality/F");
    fJetTree->Branch("Event_ID",&fBuffer_Event_ID,"Event_ID/L");
    fJetTree->Branch("Event_MagneticField",&fBuffer_Event_MagneticField,"Event_MagneticField/F");

    if(fSaveMCInformation)
      fJetTree->Branch("Event_PtHard",&fBuffer_Event_PtHard,"Event_PtHard/F");
  }

  if(fSaveConstituents)
  {
    fJetTree->Branch("Jet_Const_Pt",fBuffer_Const_Pt,"Jet_Const_Pt[Jet_NumConstituents]/F");
    fJetTree->Branch("Jet_Const_Phi",fBuffer_Const_Phi,"Jet_Const_Phi[Jet_NumConstituents]/F");
    fJetTree->Branch("Jet_Const_Eta",fBuffer_Const_Eta,"Jet_Const_Eta[Jet_NumConstituents]/F");
    fJetTree->Branch("Jet_Const_Charge",fBuffer_Const_Charge,"Jet_Const_Charge[Jet_NumConstituents]/F");
  }

  if(fSaveConstituentsIP)
  {
    fJetTree->Branch("Jet_Const_IPd",&dummy,"Jet_Const_IPd[Jet_NumConstituents]/F");
    fJetTree->Branch("Jet_Const_IPz",&dummy,"Jet_Const_IPz[Jet_NumConstituents]/F");
    fJetTree->Branch("Jet_Const_CovIPd",&dummy,"Jet_Const_CovIPd[Jet_NumConstituents]/F");
    fJetTree->Branch("Jet_Const_CovIPz",&dummy,"Jet_Const_CovIPz[Jet_NumConstituents]/F");

    fJetTree->Branch("Jet_Const_ProdVtx_X",fBuffer_Const_ProdVtx_X,"Jet_Const_ProdVtx_X[Jet_NumConstituents]/F");
    fJetTree->Branch("Jet_Const_ProdVtx_Y",fBuffer_Const_ProdVtx_Y,"Jet_Const_ProdVtx_Y[Jet_NumConstituents]/F");
    fJetTree->Branch("Jet_Const_ProdVtx_Z",fBuffer_Const_ProdVtx_Z,"Jet_Const_ProdVtx_Z[Jet_NumConstituents]/F");
  }

  if(fSaveConstituentPID)
  {
    fJetTree->Branch("Jet_Const_PID_ITS",&dummy,"Jet_Const_PID_ITS[Jet_NumConstituents]/F");
    fJetTree->Branch("Jet_Const_PID_TPC",&dummy,"Jet_Const_PID_TPC[Jet_NumConstituents]/F");
    fJetTree->Branch("Jet_Const_PID_TOF",&dummy,"Jet_Const_PID_TOF[Jet_NumConstituents]/F");
    fJetTree->Branch("Jet_Const_PID_TRD",&dummy,"Jet_Const_PID_TRD[Jet_NumConstituents]/F");

    fJetTree->Branch("Jet_Const_PID_Reconstructed",&dummy,"Jet_Const_PID_Reconstructed[Jet_NumConstituents]/S");
    if(fSaveMCInformation)
      fJetTree->Branch("Jet_Const_PID_Truth",&dummy,"Jet_Const_PID_Truth[Jet_NumConstituents]/S");
  }

  if(fSaveJetShapes)
  {
    fJetTree->Branch("Jet_Shape_Mass",&fBuffer_Shape_Mass,"Jet_Shape_Mass/F");
  }

  if(fSaveMCInformation)
  {
    fJetTree->Branch("Jet_MC_MotherParton",&fBuffer_Jet_MC_MotherParton,"Jet_MC_MotherParton/I");
    fJetTree->Branch("Jet_MC_MotherHadron",&fBuffer_Jet_MC_MotherHadron,"Jet_MC_MotherHadron/I");
    fJetTree->Branch("Jet_MC_MotherIC",&fBuffer_Jet_MC_MotherIC,"Jet_MC_MotherIC/I");
    fJetTree->Branch("Jet_MC_MatchedPt",&fBuffer_Jet_MC_MatchedPt,"Jet_MC_MatchedPt/F");
    fJetTree->Branch("Jet_MC_TruePtFraction",&fBuffer_Jet_MC_TruePtFraction,"Jet_MC_TruePtFraction/F");
  }

  if(fSaveSecondaryVertices)
  {
    fJetTree->Branch("Jet_NumSecVertices",&fBuffer_NumSecVertices,"Jet_NumSecVertices/I");

    fJetTree->Branch("Jet_SecVtx_X",&dummy,"Jet_SecVtx_X[Jet_NumSecVertices]/F");
    fJetTree->Branch("Jet_SecVtx_Y",&dummy,"Jet_SecVtx_Y[Jet_NumSecVertices]/F");
    fJetTree->Branch("Jet_SecVtx_Z",&dummy,"Jet_SecVtx_Z[Jet_NumSecVertices]/F");
    fJetTree->Branch("Jet_SecVtx_Mass",&dummy,"Jet_SecVtx_Mass[Jet_NumSecVertices]/F");
    fJetTree->Branch("Jet_SecVtx_Lxy",&dummy,"Jet_SecVtx_Lxy[Jet_NumSecVertices]/F");
    fJetTree->Branch("Jet_SecVtx_SigmaLxy",&dummy,"Jet_SecVtx_SigmaLxy[Jet_NumSecVertices]/F");
    fJetTree->Branch("Jet_SecVtx_Chi2",&dummy,"Jet_SecVtx_Chi2[Jet_NumSecVertices]/F");
    fJetTree->Branch("Jet_SecVtx_Dispersion",&dummy,"Jet_SecVtx_Dispersion[Jet_NumSecVertices]/F");
  }

  // Trigger track requirement active -> Save trigger track
  if(fSaveTriggerTracks)
  {
    fJetTree->Branch("Jet_NumTriggerTracks",&fBuffer_NumTriggerTracks,"Jet_NumTriggerTracks/I");
    fJetTree->Branch("Jet_TriggerTrack_Pt",&dummy,"Jet_TriggerTrack_Pt[Jet_NumTriggerTracks]/F");
    fJetTree->Branch("Jet_TriggerTrack_dEta",&dummy,"Jet_TriggerTrack_dEta[Jet_NumTriggerTracks]/F");
    fJetTree->Branch("Jet_TriggerTrack_dPhi",&dummy,"Jet_TriggerTrack_dPhi[Jet_NumTriggerTracks]/F");
  }

  fInitialized = kTRUE;
}

//________________________________________________________________________
AliAnalysisTaskJetExtractor::AliAnalysisTaskJetExtractor() :
  AliAnalysisTaskEmcalJet("AliAnalysisTaskJetExtractor", kTRUE),
  fJetTree(0),
  fVtxTagger(0),
  fHadronMatchingRadius(0.4),
  fTrueJetMatchingRadius(0.4),
  fSecondaryVertexMaxChi2(1e10),
  fSecondaryVertexMaxDispersion(0.05),
  fCalculateSecondaryVertices(kTRUE),
  fVertexerCuts(0),
  fSetEmcalJetFlavour(0),
  fEventPercentage(1.0),
  fTruthMinLabel(0),
  fTruthMaxLabel(100000),
  fEventCut_TriggerTrackMinPt(0),
  fEventCut_TriggerTrackMaxPt(0),
  fEventCut_TriggerTrackMinLabel(-9999999),
  fEventCut_TriggerTrackMaxLabel(+9999999),
  fJetsCont(0),
  fTracksCont(0),
  fTruthParticleArray(0),
  fTruthJetsArrayName(""),
  fTruthJetsRhoName(""),
  fTruthParticleArrayName("mcparticles"),
  fTriggerTracks_Pt(),
  fTriggerTracks_Eta(),
  fTriggerTracks_Phi()
{
  SetMakeGeneralHistograms(kTRUE);
  fJetTree = new AliEmcalJetTree(GetName());
  fJetTree->SetSaveEventProperties(kTRUE);
  fJetTree->SetSaveConstituents(kTRUE);
  fJetTree->SetSaveJetShapes(kTRUE);
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskJetExtractor::AliAnalysisTaskJetExtractor(const char *name) :
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fJetTree(0),
  fVtxTagger(0),
  fHadronMatchingRadius(0.4),
  fTrueJetMatchingRadius(0.4),
  fSecondaryVertexMaxChi2(1e10),
  fSecondaryVertexMaxDispersion(0.05),
  fCalculateSecondaryVertices(kTRUE),
  fVertexerCuts(0),
  fSetEmcalJetFlavour(0),
  fEventPercentage(1.0),
  fTruthMinLabel(0),
  fTruthMaxLabel(100000),
  fEventCut_TriggerTrackMinPt(0),
  fEventCut_TriggerTrackMaxPt(0),
  fEventCut_TriggerTrackMinLabel(-9999999),
  fEventCut_TriggerTrackMaxLabel(+9999999),
  fJetsCont(0),
  fTracksCont(0),
  fTruthParticleArray(0),
  fTruthJetsArrayName(""),
  fTruthJetsRhoName(""),
  fTruthParticleArrayName("mcparticles"),
  fTriggerTracks_Pt(),
  fTriggerTracks_Eta(),
  fTriggerTracks_Phi()
{
  SetMakeGeneralHistograms(kTRUE);
  fJetTree = new AliEmcalJetTree(GetName());
  fJetTree->SetSaveEventProperties(kTRUE);
  fJetTree->SetSaveConstituents(kTRUE);
  fJetTree->SetSaveJetShapes(kTRUE);
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskJetExtractor::~AliAnalysisTaskJetExtractor()
{
  // Destructor.
}

//________________________________________________________________________
void AliAnalysisTaskJetExtractor::UserCreateOutputObjects()
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

  // Activate saving of trigger tracks if this is demanded
  if(fEventCut_TriggerTrackMinPt || fEventCut_TriggerTrackMaxPt)
    fJetTree->SetSaveTriggerTracks(kTRUE);

  // ### Initialize the jet tree (settings must all be given at this stage)
  fJetTree->InitializeTree();
  OpenFile(2);
  PostData(2, fJetTree->GetTreePointer());

  PrintConfig();

  // ### Add control histograms (already some created in base task)
  AddHistogram2D<TH2D>("hTrackCount", "Number of tracks in acceptance vs. centrality", "COLZ", 500, 0., 5000., 100, 0, 100, "N tracks","Centrality", "dN^{Events}/dN^{Tracks}");
  AddHistogram2D<TH2D>("hBackgroundPt", "Background p_{T} distribution", "", 150, 0., 150., 100, 0, 100, "Background p_{T} (GeV/c)", "Centrality", "dN^{Events}/dp_{T}");

  AddHistogram2D<TH2D>("hJetPtRaw", "Jets p_{T} distribution (raw)", "COLZ", 300, 0., 300., 100, 0, 100, "p_{T, jet} (GeV/c)", "Centrality", "dN^{Jets}/dp_{T}");
  AddHistogram2D<TH2D>("hJetPt", "Jets p_{T} distribution (background subtracted)", "COLZ", 400, -100., 300., 100, 0, 100, "p_{T, jet} (GeV/c)", "Centrality", "dN^{Jets}/dp_{T}");
  AddHistogram2D<TH2D>("hJetPtExtracted", "Extracted jets p_{T} distribution (background subtracted)", "COLZ", 400, -100., 300., 100, 0, 100, "p_{T, jet} (GeV/c)", "Centrality", "dN^{Jets}/dp_{T}");
  AddHistogram2D<TH2D>("hJetPhiEta", "Jet angular distribution #phi/#eta", "COLZ", 180, 0., 2*TMath::Pi(), 100, -2.5, 2.5, "#phi", "#eta", "dN^{Jets}/d#phi d#eta");
  AddHistogram2D<TH2D>("hJetArea", "Jet area", "COLZ", 200, 0., 2., 100, 0, 100, "Jet A", "Centrality", "dN^{Jets}/dA");

  AddHistogram2D<TH2D>("hDeltaPt", "#delta p_{T} distribution", "", 400, -100., 300., 100, 0, 100, "p_{T, cone} (GeV/c)", "Centrality", "dN^{Tracks}/dp_{T}");

  AddHistogram2D<TH2D>("hConstituentPt", "Jet constituent p_{T} distribution", "COLZ", 400, 0., 300., 100, 0, 100, "p_{T, const} (GeV/c)", "Centrality", "dN^{Const}/dp_{T}");
  AddHistogram2D<TH2D>("hConstituentPhiEta", "Jet constituent relative #phi/#eta distribution", "COLZ", 120, -0.6, 0.6, 120, -0.6, 0.6, "#Delta#phi", "#Delta#eta", "dN^{Const}/d#phi d#eta");

  AddHistogram1D<TH1D>("hExtractionPercentage", "Extracted jets p_{T} distribution (background subtracted)", "e1p", 400, -100., 300., "p_{T, jet} (GeV/c)", "Extracted percentage");

  // Track QA plots
  AddHistogram2D<TH2D>("hTrackPt", "Tracks p_{T} distribution", "", 300, 0., 300., 100, 0, 100, "p_{T} (GeV/c)", "Centrality", "dN^{Tracks}/dp_{T}");
  AddHistogram2D<TH2D>("hTrackPhi", "Track angular distribution in #phi", "LEGO2", 180, 0., 2*TMath::Pi(), 100, 0, 100, "#phi", "Centrality", "dN^{Tracks}/(d#phi)");
  AddHistogram2D<TH2D>("hTrackEta", "Track angular distribution in #eta", "LEGO2", 100, -2.5, 2.5, 100, 0, 100, "#eta", "Centrality", "dN^{Tracks}/(d#eta)");
  AddHistogram2D<TH2D>("hTrackPhiEta", "Track angular distribution #phi/#eta", "COLZ", 180, 0., 2*TMath::Pi(), 100, -2.5, 2.5, "#phi", "#eta", "dN^{Tracks}/d#phi d#eta");

  AddHistogram2D<TH2D>("hTrackEtaPt", "Track angular distribution in #eta vs. p_{T}", "LEGO2", 100, -2.5, 2.5, 300, 0., 300., "#eta", "p_{T} (GeV/c)", "dN^{Tracks}/(d#eta dp_{T})");
  AddHistogram2D<TH2D>("hTrackPhiPt", "Track angular distribution in #phi vs. p_{T}", "LEGO2", 180, 0, 2*TMath::Pi(), 300, 0., 300., "#phi", "p_{T} (GeV/c)", "dN^{Tracks}/(d#phi dp_{T})");
}


//________________________________________________________________________
void AliAnalysisTaskJetExtractor::ExecOnce()
{
  AliAnalysisTaskEmcalJet::ExecOnce();
  if (fTruthParticleArrayName != "")
    fTruthParticleArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTruthParticleArrayName.Data()));

  // ### Prepare vertexer
  if(fJetTree->GetSaveSecondaryVertices())
  {
    if(!fVertexerCuts)
      AliFatal("VertexerCuts not given but secondary vertex calculation turned on.");
    fVtxTagger = new AliHFJetsTaggingVertex();
    fVtxTagger->SetCuts(fVertexerCuts);
  }

  // ### Save tree extraction percentage to histogram
  std::vector<Float_t> extractionPtBins      = fJetTree->GetExtractionPercentagePtBins();
  std::vector<Float_t> extractionPercentages = fJetTree->GetExtractionPercentages();

  for(Int_t i=0; i<static_cast<Int_t>(extractionPercentages.size()); i++)
  {
    Double_t percentage = extractionPercentages[i];
    for(Int_t pt=static_cast<Int_t>(extractionPtBins[i*2]); pt<static_cast<Int_t>(extractionPtBins[i*2+1]); pt++)
    {
      FillHistogram("hExtractionPercentage", pt, percentage);
    }
  }

  // ### Add HIJING Ncoll histogram (in case we need it)
  if(MCEvent())
  {
    if(dynamic_cast<AliGenHijingEventHeader*>(MCEvent()->GenEventHeader()))
      AddHistogram1D<TH1D>("hHijingNcoll", "N_{coll} from HIJING", "e1p", 1000, 0., 5000, "N_{coll}", "dN^{Events}/dN^{N_{coll}}");
  }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetExtractor::Run()
{
  // ################################### EVENT SELECTION

  if(!IsTriggerTrackInEvent())
    return kFALSE;
  if(gRandom->Rndm() >= fEventPercentage)
    return kFALSE;

  // ################################### EVENT PROPERTIES

  FillEventControlHistograms();

  // Load vertex if possible
  Double_t vtxX = 0;
  Double_t vtxY = 0;
  Double_t vtxZ = 0;
  Long64_t eventID = 0;
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
  if(eventIDHeader)
    eventID = eventIDHeader->GetEventIdAsLong();

  // If available, get Ncoll from HIJING
  if(MCEvent())
  {
    if(AliGenHijingEventHeader* mcHeader = dynamic_cast<AliGenHijingEventHeader*>(MCEvent()->GenEventHeader()))
    {
      Double_t ncoll = mcHeader->NN() + mcHeader->NNw() + mcHeader->NwN() + mcHeader->NwNw();
      FillHistogram("hHijingNcoll", ncoll);
    }
  }

  // ################################### MAIN JET LOOP
  fJetsCont->ResetCurrentID();
  while(AliEmcalJet *jet = fJetsCont->GetNextAcceptJet())
  {
    FillJetControlHistograms(jet);

    Double_t matchedJetPt = 0;
    Double_t truePtFraction = 0;
    Int_t currentJetType_HM = 0;
    Int_t currentJetType_PM = 0;
    Int_t currentJetType_IC = 0;
    if(fJetTree->GetSaveMCInformation())
    {
      // Get jet type from MC (hadron matching, parton matching definition - for HF jets)
      GetJetType(jet, currentJetType_HM, currentJetType_PM, currentJetType_IC);
      // Get true pT estimators
      GetJetTruePt(jet, matchedJetPt, truePtFraction);
    }

    // ### CONSTITUENT LOOP: Retrieve PID values + impact parameters
    std::vector<Float_t> vecSigITS; std::vector<Float_t> vecSigTPC; std::vector<Float_t> vecSigTRD; std::vector<Float_t> vecSigTOF; std::vector<Short_t> vecRecoPID; std::vector<Short_t> vecTruePID;
    std::vector<Float_t> vec_d0; std::vector<Float_t> vec_d0cov; std::vector<Float_t> vec_z0; std::vector<Float_t> vec_z0cov;

    if(fJetTree->GetSaveConstituentPID() || fJetTree->GetSaveConstituentsIP())
      for(Int_t i = 0; i < jet->GetNumberOfTracks(); i++)
      {
        AliVParticle* particle = static_cast<AliVParticle*>(jet->TrackAt(i, fTracksCont->GetArray()));
        if(!particle) continue;
        if(fJetTree->GetSaveConstituentPID())
        {
          Float_t sigITS = 0; Float_t sigTPC = 0; Float_t sigTOF = 0; Float_t sigTRD = 0; Short_t recoPID = 0; Short_t truePID = 0;
          AddPIDInformation(particle, sigITS, sigTPC, sigTOF, sigTRD, recoPID, truePID);
          vecSigITS.push_back(sigITS); vecSigTPC.push_back(sigTPC); vecSigTOF.push_back(sigTOF); vecSigTRD.push_back(sigTRD); vecRecoPID.push_back(recoPID); vecTruePID.push_back(truePID);
        }
        if(fJetTree->GetSaveConstituentsIP())
        {
          Float_t d0 = 0; Float_t d0cov = 0; Float_t z0 = 0; Float_t z0cov = 0;
          GetTrackImpactParameters(myVertex, dynamic_cast<AliAODTrack*>(particle), d0, d0cov, z0, z0cov);
          vec_d0.push_back(d0); vec_d0cov.push_back(d0cov); vec_z0.push_back(z0); vec_z0cov.push_back(z0cov);
        }
      }

    // Reconstruct secondary vertices
    std::vector<Float_t> secVtx_X; std::vector<Float_t> secVtx_Y; std::vector<Float_t> secVtx_Z; std::vector<Float_t> secVtx_Mass; std::vector<Float_t> secVtx_Lxy; std::vector<Float_t> secVtx_SigmaLxy; std::vector<Float_t> secVtx_Chi2; std::vector<Float_t> secVtx_Dispersion;
    if(fJetTree->GetSaveSecondaryVertices())
      ReconstructSecondaryVertices(myVertex, jet, secVtx_X, secVtx_Y, secVtx_Z, secVtx_Mass, secVtx_Lxy, secVtx_SigmaLxy, secVtx_Chi2, secVtx_Dispersion);

    // Now change trigger tracks eta/phi to dEta/dPhi relative to jet
    std::vector<Float_t> triggerTracks_dEta(fTriggerTracks_Eta);
    std::vector<Float_t> triggerTracks_dPhi(fTriggerTracks_Phi);
    for(UInt_t i=0; i<triggerTracks_dEta.size(); i++)
    {
      triggerTracks_dEta[i] = jet->Eta() - fTriggerTracks_Eta[i];
      triggerTracks_dPhi[i] = TMath::Min(TMath::Abs(jet->Phi() - fTriggerTracks_Phi[i]), TMath::TwoPi() - TMath::Abs(jet->Phi() - fTriggerTracks_Phi[i]));
      if(  ((TMath::Abs(jet->Phi() - fTriggerTracks_Phi[i]) <= TMath::Pi()) && (jet->Phi() - fTriggerTracks_Phi[i] > 0)) 
        || ((TMath::Abs(jet->Phi() - fTriggerTracks_Phi[i]) > TMath::Pi())  && (jet->Phi() - fTriggerTracks_Phi[i] < 0)) )
        triggerTracks_dPhi[i] = -triggerTracks_dPhi[i];
    }

    // Fill jet to tree
    Bool_t accepted = fJetTree->AddJetToTree(jet, fJetsCont->GetRhoVal(), vtxX, vtxY, vtxZ, fCent, eventID, InputEvent()->GetMagneticField(), fTracksCont,
              currentJetType_PM,currentJetType_HM,currentJetType_IC,matchedJetPt,truePtFraction,fPtHard,
              vecSigITS, vecSigTPC, vecSigTOF, vecSigTRD, vecRecoPID, vecTruePID,
              fTriggerTracks_Pt, triggerTracks_dEta, triggerTracks_dPhi,
              vec_d0, vec_z0, vec_d0cov, vec_z0cov,
              secVtx_X, secVtx_Y, secVtx_Z, secVtx_Mass, secVtx_Lxy, secVtx_SigmaLxy, secVtx_Chi2, secVtx_Dispersion);

    if(accepted)
      FillHistogram("hJetPtExtracted", jet->Pt() - fJetsCont->GetRhoVal()*jet->Area(), fCent);
  }

  // ################################### PARTICLE PROPERTIES
  fTracksCont->ResetCurrentID();
  while(AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle()))
    FillTrackControlHistograms(track);

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetExtractor::IsTriggerTrackInEvent()
{
  // Cut for trigger track requirement
  if(fEventCut_TriggerTrackMinPt || fEventCut_TriggerTrackMaxPt)
  {
    // Clear vector of trigger tracks
    fTriggerTracks_Pt.clear();
    fTriggerTracks_Eta.clear();
    fTriggerTracks_Phi.clear();

    // Go through all tracks and check whether trigger tracks can be found
    fTracksCont->ResetCurrentID();
    while(AliVTrack *track = static_cast<AliVTrack*>(fTracksCont->GetNextAcceptParticle()))
    {
      if( (track->GetLabel() >= fEventCut_TriggerTrackMinLabel) && (track->GetLabel() < fEventCut_TriggerTrackMaxLabel) )
        if( (track->Pt() >= fEventCut_TriggerTrackMinPt) && (track->Pt() < fEventCut_TriggerTrackMaxPt) )
        {
          fTriggerTracks_Pt.push_back(track->Pt());
          fTriggerTracks_Eta.push_back(track->Eta());
          fTriggerTracks_Phi.push_back(track->Phi());
        }
    }
    // No particle has been found that fulfills requirement -> Do not accept event
    if(fTriggerTracks_Pt.size() == 0)
      return kFALSE;
  }
  return kTRUE;
}


//________________________________________________________________________
void AliAnalysisTaskJetExtractor::GetJetTruePt(AliEmcalJet* jet, Double_t& matchedJetPt, Double_t& truePtFraction)
{
  // #################################################################################
  // ##### METHOD 1: If fTruthJetsArrayName is set, a matching jet is searched for
  Double_t     bestMatchDeltaR = 999.;
  if(fTruthJetsArrayName != "")
  {
    // "True" background
    AliRhoParameter* rho = static_cast<AliRhoParameter*>(InputEvent()->FindListObject(fTruthJetsRhoName.Data()));
    Double_t trueRho = 0;
    if(rho)
     trueRho = rho->GetVal();

    TClonesArray* truthArray = static_cast<TClonesArray*>(InputEvent()->FindListObject(Form("%s", fTruthJetsArrayName.Data())));

    // Loop over all true jets to find the best match
    matchedJetPt = 0;
    if(truthArray)
      for(Int_t i=0; i<truthArray->GetEntries(); i++)
      {
        AliEmcalJet* truthJet = static_cast<AliEmcalJet*>(truthArray->At(i));
        if(truthJet->Pt() < 0.15)
          continue;

        Double_t deltaEta = (truthJet->Eta()-jet->Eta());
        Double_t deltaPhi = TMath::Min(TMath::Abs(truthJet->Phi()-jet->Phi()),TMath::TwoPi() - TMath::Abs(truthJet->Phi()-jet->Phi()));
        Double_t deltaR = TMath::Sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);

        // Cut jets too far away
        if (deltaR > fTrueJetMatchingRadius)
          continue;

        // Search for the best match
        if(deltaR < bestMatchDeltaR)
        {
          bestMatchDeltaR = deltaR;
          matchedJetPt = truthJet->Pt() - truthJet->Area()* trueRho;
        }
      }
  }

  // #################################################################################
  // ##### METHOD 2: Calculate fraction of "true" pT -- pT which is not from a toy
  Double_t pt_nonMC = 0.;
  Double_t pt_all   = 0.;
  truePtFraction = 0;

  for(Int_t iConst = 0; iConst < jet->GetNumberOfTracks(); iConst++)
  {
    // Loop over all valid jet constituents
    AliVParticle* particle = static_cast<AliVParticle*>(jet->TrackAt(iConst, fTracksCont->GetArray()));
    if(!particle) continue;
    if(particle->Pt() < 1e-6) continue;

    // Particles marked w/ labels within label range are considered from toy
    if( (particle->GetLabel() >= fTruthMinLabel) && (particle->GetLabel() < fTruthMaxLabel))
      pt_nonMC += particle->Pt();
    pt_all += particle->Pt();
  }
  if(pt_all)
    truePtFraction = (pt_nonMC/pt_all);

}


//________________________________________________________________________
void AliAnalysisTaskJetExtractor::GetJetType(AliEmcalJet* jet, Int_t& typeHM, Int_t& typePM, Int_t& typeIC)
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
    // No HF jet (according to hadron matching) -- now also separate jets in udg (1) and s-jets (3)
    if(IsStrangeJet(jet))
      typeHM = 3;
    else
      typeHM = 1;
  }
  else {
    Int_t pdg = TMath::Abs(parton[1]->PdgCode());
    if(fVtxTagger->IsDMeson(pdg)) typeHM = 4;
    else if (fVtxTagger->IsBMeson(pdg)) typeHM = 5;
  }

  // Set flavour of AliEmcalJet object (set ith bit while i corresponds to type)
  if(fSetEmcalJetFlavour)
    jet->AddFlavourTag(static_cast<Int_t>(TMath::Power(2, typeHM)));


  const AliEmcalPythiaInfo* partonsInfo = GetPythiaInfo();
  typeIC = 0;
  if (partonsInfo)
  {
    // Get primary partons directions 
    Double_t parton1phi = partonsInfo->GetPartonPhi6();
    Double_t parton1eta = partonsInfo->GetPartonEta6();
    Double_t parton2phi = partonsInfo->GetPartonPhi7();
    Double_t parton2eta = partonsInfo->GetPartonEta7();


    Double_t delta1Eta = (parton1eta-jet->Eta());
    Double_t delta1Phi = TMath::Min(TMath::Abs(parton1phi-jet->Phi()),TMath::TwoPi() - TMath::Abs(parton1phi-jet->Phi()));
    Double_t delta1R   = TMath::Sqrt(delta1Eta*delta1Eta + delta1Phi*delta1Phi);
    Double_t delta2Eta = (parton2eta-jet->Eta());
    Double_t delta2Phi = TMath::Min(TMath::Abs(parton2phi-jet->Phi()),TMath::TwoPi() - TMath::Abs(parton2phi-jet->Phi()));
    Double_t delta2R   = TMath::Sqrt(delta2Eta*delta2Eta + delta2Phi*delta2Phi);

    // Check if one of the partons if closer than matching criterion
    Bool_t matched = (delta1R < fJetsCont->GetJetRadius()/2.) || (delta2R < fJetsCont->GetJetRadius()/2.);

    // Matching criterion fulfilled -> Set flag to closest
    if(matched)
    {
      if(delta1R < delta2R)
        typeIC = partonsInfo->GetPartonFlag6();
      else
        typeIC = partonsInfo->GetPartonFlag7();
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskJetExtractor::GetTrackImpactParameters(const AliVVertex* vtx, AliAODTrack* track, Float_t& d0, Float_t& d0cov, Float_t& z0, Float_t& z0cov)
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
void AliAnalysisTaskJetExtractor::ReconstructSecondaryVertices(const AliVVertex* primVtx, const AliEmcalJet* jet, std::vector<Float_t>& secVtx_X, std::vector<Float_t>& secVtx_Y, std::vector<Float_t>& secVtx_Z, std::vector<Float_t>& secVtx_Mass, std::vector<Float_t>& secVtx_Lxy, std::vector<Float_t>& secVtx_SigmaLxy, std::vector<Float_t>& secVtx_Chi2, std::vector<Float_t>& secVtx_Dispersion)
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

    secVtx_X.push_back(secVtx->GetX()); secVtx_Y.push_back(secVtx->GetY()); secVtx_Z.push_back(secVtx->GetZ()); secVtx_Chi2.push_back(secVtx->GetChi2perNDF());
    secVtx_Dispersion.push_back(dispersion); secVtx_Mass.push_back(mass); secVtx_Lxy.push_back(Lxy); secVtx_SigmaLxy.push_back(sigmaLxy); 
  }

  if(fCalculateSecondaryVertices)
  {
    secVertexArr->Clear();
    delete secVertexArr;
  }
  delete esdVtx;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetExtractor::IsStrangeJet(AliEmcalJet* jet)
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
void AliAnalysisTaskJetExtractor::FillEventControlHistograms()
{
  // ### Event control plots
  FillHistogram("hTrackCount", fTracksCont->GetNAcceptedParticles(), fCent);
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
    Double_t tmpRandConeEta = fJetsCont->GetJetEtaMin() + gRandom->Rndm()*TMath::Abs(fJetsCont->GetJetEtaMax()-fJetsCont->GetJetEtaMin());
    Double_t tmpRandConePhi = gRandom->Rndm()*TMath::TwoPi();
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
void AliAnalysisTaskJetExtractor::FillTrackControlHistograms(AliVTrack* track)
{
  FillHistogram("hTrackPt", track->Pt(), fCent);
  FillHistogram("hTrackPhi", track->Phi(), fCent);
  FillHistogram("hTrackEta", track->Eta(), fCent);
  FillHistogram("hTrackEtaPt", track->Eta(), track->Pt());
  FillHistogram("hTrackPhiPt", track->Phi(), track->Pt());
  FillHistogram("hTrackPhiEta", track->Phi(), track->Eta());
}

//________________________________________________________________________
void AliAnalysisTaskJetExtractor::AddPIDInformation(AliVParticle* particle, Float_t& sigITS, Float_t& sigTPC, Float_t& sigTOF, Float_t& sigTRD, Short_t& recoPID, Short_t& truePID)
{
  truePID = 9;
  if(!particle) return;

  // If we have AODs, retrieve particle PID signals
  AliAODTrack* aodtrack = dynamic_cast<AliAODTrack*>(particle);

  if(aodtrack)
  {
    // Get AOD value from reco
    recoPID  = aodtrack->GetMostProbablePID();
    AliAODPid* pidObj = aodtrack->GetDetPid();
    if(!pidObj)
      return;

    sigITS = pidObj->GetITSsignal();
    sigTPC = pidObj->GetTPCsignal();
    sigTOF = pidObj->GetTOFsignal();
    sigTRD = pidObj->GetTRDsignal();
  }

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
          truePID = 4;
        else if (TMath::Abs(mcParticle->PdgCode()) == 211) // pion
          truePID = 2;
        else if (TMath::Abs(mcParticle->PdgCode()) == 321) // kaon
          truePID = 3;
        else if (TMath::Abs(mcParticle->PdgCode()) == 11) // electron
          truePID = 0;
        else if (TMath::Abs(mcParticle->PdgCode()) == 13) // muon
          truePID = 1;
        else if (TMath::Abs(mcParticle->PdgCode()) == 700201) // deuteron
          truePID = 5;
        else if (TMath::Abs(mcParticle->PdgCode()) == 700301) // triton
          truePID = 6;
        else if (TMath::Abs(mcParticle->PdgCode()) == 700302) // He3
          truePID = 7;
        else if (TMath::Abs(mcParticle->PdgCode()) == 700202) // alpha
          truePID = 8;
        else
          truePID = 9;

        break;
      }
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskJetExtractor::PrintConfig()
{
  // Print properties for extraction
  // to be implemented later
}


//########################################################################
// HELPERS
//########################################################################

//________________________________________________________________________
inline Bool_t AliAnalysisTaskJetExtractor::IsTrackInCone(AliVParticle* track, Double_t eta, Double_t phi, Double_t radius)
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
