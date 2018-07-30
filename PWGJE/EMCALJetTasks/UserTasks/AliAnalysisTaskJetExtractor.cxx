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
#include "AliAODVertex.h"
#include "AliAODPid.h"

#include "AliVTrack.h"
#include "AliVHeader.h"
#include "AliEmcalJet.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliTrackContainer.h"
#include "AliAODTrack.h"
#include "AliVParticle.h"
#include "TRandom3.h"
#include "AliAnalysisTaskEmcalJet.h"

#include "AliAnalysisTaskJetExtractor.h"

/// \cond CLASSIMP
ClassImp(AliEmcalJetTree)
/// \endcond

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskJetExtractor)
/// \endcond

//________________________________________________________________________
AliEmcalJetTree::AliEmcalJetTree() : TNamed("CustomTree", "CustomTree"), fJetTree(0), fInitialized(0), fSaveEventProperties(0), fSaveConstituents(0), fSaveConstituentsIP(0), fSaveConstituentPID(0), fSaveJetShapes(0), fSaveMCInformation(0), fSaveSecondaryVertices(0), fExtractionPercentages(), fExtractionPercentagePtBins()
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
AliEmcalJetTree::AliEmcalJetTree(const char* name) : TNamed(name, name), fJetTree(0), fInitialized(0), fSaveEventProperties(0), fSaveConstituents(0), fSaveConstituentsIP(0), fSaveConstituentPID(0), fSaveJetShapes(0), fSaveMCInformation(0), fSaveSecondaryVertices(0), fExtractionPercentages(), fExtractionPercentagePtBins()
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
  AliTrackContainer* trackCont, Int_t motherParton, Int_t motherHadron, Float_t truePt, Float_t ptHard,
  Float_t* trackPID_ITS, Float_t* trackPID_TPC, Float_t* trackPID_TOF, Float_t* trackPID_TRD, Short_t* trackPID_Reco, Short_t* trackPID_Truth,
  Float_t* trackIP_d0, Float_t* trackIP_z0, Float_t* trackIP_d0cov, Float_t* trackIP_z0cov,
  Int_t numSecVertices, Float_t* secVtx_X, Float_t* secVtx_Y, Float_t* secVtx_Z, Float_t* secVtx_Mass, Float_t* secVtx_Lxy, Float_t* secVtx_SigmaLxy, Float_t* secVtx_Chi2, Float_t* secVtx_Dispersion)
{
  // ############ 
  // Approach: Fill buffer and then add to tree
  //
  if(!fInitialized)
    AliFatal("Tree is not initialized.");

  fBuffer_JetPt                                   = jet->Pt() - bgrdDensity*jet->Area();

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
    fBuffer_Jet_MC_TruePt = truePt;
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
    fJetTree->Branch("Jet_Const_CoVIPz",&dummy,"Jet_Const_CovIPz[Jet_NumConstituents]/F");

    fJetTree->Branch("Jet_Const_ProdVtx_X",&dummy,"Jet_Const_ProdVtx_X[Jet_NumConstituents]/F");
    fJetTree->Branch("Jet_Const_ProdVtx_Y",&dummy,"Jet_Const_ProdVtx_Y[Jet_NumConstituents]/F");
    fJetTree->Branch("Jet_Const_ProdVtx_Z",&dummy,"Jet_Const_ProdVtx_Z[Jet_NumConstituents]/F");
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
    fJetTree->Branch("Jet_MC_TruePt",&fBuffer_Jet_MC_TruePt,"Jet_MC_TruePt/F");
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

  fInitialized = kTRUE;
}




//________________________________________________________________________
AliAnalysisTaskJetExtractor::AliAnalysisTaskJetExtractor() :
  AliAnalysisTaskEmcalJet("AliAnalysisTaskJetExtractor", kTRUE),
  fJetTree(0),
  fJetsCont(0),
  fTracksCont(0)
{
  SetMakeGeneralHistograms(kTRUE);
  fJetTree = new AliEmcalJetTree(GetName());
  fJetTree->SetSaveEventProperties(kTRUE);
  fJetTree->SetSaveConstituents(kTRUE);
  fJetTree->SetSaveJetShapes(kTRUE);
}

//________________________________________________________________________
AliAnalysisTaskJetExtractor::AliAnalysisTaskJetExtractor(const char *name) :
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fJetTree(0),
  fJetsCont(0),
  fTracksCont(0)
{
  SetMakeGeneralHistograms(kTRUE);
  fJetTree = new AliEmcalJetTree(GetName());
  fJetTree->SetSaveEventProperties(kTRUE);
  fJetTree->SetSaveConstituents(kTRUE);
  fJetTree->SetSaveJetShapes(kTRUE);
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
void AliAnalysisTaskJetExtractor::ExecOnce()
{
  AliAnalysisTaskEmcalJet::ExecOnce();
  // ### Initialize the jet tree (settings must all be given at this stage)
  fJetTree->InitializeTree();
  fOutput->Add(fJetTree->GetTreePointer());
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetExtractor::Run()
{
  // ################################### EVENT PROPERTIES
  FillEventControlHistograms();

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


  // ################################### MAIN JET LOOP
  fJetsCont->ResetCurrentID();
  while(AliEmcalJet *jet = fJetsCont->GetNextAcceptJet())
  {
    FillJetControlHistograms(jet);

    // ### CONSTITUENT LOOP: Retrieve PID values
    std::vector<Float_t> vecSigITS; std::vector<Float_t> vecSigTPC; std::vector<Float_t> vecSigTRD; std::vector<Float_t> vecSigTOF; std::vector<Short_t> vecRecoPID;
    for(Int_t i = 0; i < jet->GetNumberOfTracks(); i++)
    {
      AliVParticle* particle = static_cast<AliVParticle*>(jet->TrackAt(i, fTracksCont->GetArray()));
      if(!particle) continue;
      Float_t sigITS = 0; Float_t sigTPC = 0; Float_t sigTOF = 0; Float_t sigTRD = 0; Short_t recoPID = 0;
      AddPIDInformation(particle, sigITS, sigTPC, sigTOF, sigTRD, recoPID);
      vecSigITS.push_back(sigITS); vecSigTPC.push_back(sigTPC); vecSigTOF.push_back(sigTOF); vecSigTRD.push_back(sigTRD); vecRecoPID.push_back(recoPID);
    }

    // Fill jet to tree
    Bool_t accepted = fJetTree->AddJetToTree(jet, fJetsCont->GetRhoVal(), vtxX, vtxY, vtxZ, fCent, eventID, 0, fTracksCont, 0,0,0,0, vecSigITS, vecSigTPC, vecSigTOF, vecSigTRD, vecRecoPID);
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
void AliAnalysisTaskJetExtractor::AddPIDInformation(AliVParticle* particle, Float_t& sigITS, Float_t& sigTPC, Float_t& sigTOF, Float_t& sigTRD, Short_t& recoPID)
{
  // If we have AODs, retrieve particle PID signals
  AliAODTrack* aodtrack = dynamic_cast<AliAODTrack*>(particle);

  if(!aodtrack)
    return;

  // Get AOD value from reco
  recoPID  = aodtrack->GetMostProbablePID();
  AliAODPid* pidObj = aodtrack->GetDetPid();
  sigITS = pidObj->GetITSsignal();
  sigTPC = pidObj->GetTPCsignal();
  sigTOF = pidObj->GetTOFsignal();
  sigTRD = pidObj->GetTRDsignal();
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
