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

#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TTree.h>
#include <TList.h>
#include <TLorentzVector.h>
#include <TNamed.h>
#include <TGrid.h>
#include <TFile.h>
#include <TSystem.h>
//#include <TTimeStamp.h>

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
  #include <TPython.h>
#endif

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
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliAODTrack.h"
#include "AliVParticle.h"
#include "TRandom3.h"
#include "AliEmcalPythiaInfo.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliAnalysisManager.h"
#include "AliEmcalContainerUtils.h"
#include "AliFJWrapper.h"


#include "AliGenHepMCEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliHFJetsTaggingVertex.h"
#include "AliRDHFJetsCutsVertex.h"

#include "AliAnalysisTaskJetExtractor.h"

/// \cond CLASSIMP
ClassImp(AliEmcalJetTree)
/// \endcond

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskJetExtractor)
/// \endcond

//________________________________________________________________________
AliEmcalJetTree::AliEmcalJetTree() : TNamed("CustomTree", "CustomTree"), fJetTree(0), fInitialized(0), fExtractionPercentages(), fExtractionPercentagePtBins(), fExtractionJetTypes_HM(), fExtractionJetTypes_PM()
{
  // For these arrays, we need to reserve memory
  fBuffer_Track_Pt         = new Float_t[kMaxNumConstituents];
  fBuffer_Track_Eta        = new Float_t[kMaxNumConstituents];
  fBuffer_Track_Phi        = new Float_t[kMaxNumConstituents];
  fBuffer_Track_Charge     = new Float_t[kMaxNumConstituents];
  fBuffer_Track_Label      = new Int_t  [kMaxNumConstituents];
  fBuffer_Track_ProdVtx_X  = new Float_t[kMaxNumConstituents];
  fBuffer_Track_ProdVtx_Y  = new Float_t[kMaxNumConstituents];
  fBuffer_Track_ProdVtx_Z  = new Float_t[kMaxNumConstituents];

  fBuffer_Cluster_Pt       = new Float_t[kMaxNumConstituents];
  fBuffer_Cluster_E        = new Float_t[kMaxNumConstituents];
  fBuffer_Cluster_Eta      = new Float_t[kMaxNumConstituents];
  fBuffer_Cluster_Phi      = new Float_t[kMaxNumConstituents];
  fBuffer_Cluster_M02      = new Float_t[kMaxNumConstituents];
  fBuffer_Cluster_Time     = new Float_t[kMaxNumConstituents];
  fBuffer_Cluster_Label    = new Int_t[kMaxNumConstituents];
}

//________________________________________________________________________
AliEmcalJetTree::AliEmcalJetTree(const char* name) : TNamed(name, name), fJetTree(0), fInitialized(0), fExtractionPercentages(), fExtractionPercentagePtBins(), fExtractionJetTypes_HM(), fExtractionJetTypes_PM()
{
  // For these arrays, we need to reserve memory
  fBuffer_Track_Pt         = new Float_t[kMaxNumConstituents];
  fBuffer_Track_Eta        = new Float_t[kMaxNumConstituents];
  fBuffer_Track_Phi        = new Float_t[kMaxNumConstituents];
  fBuffer_Track_Charge     = new Float_t[kMaxNumConstituents];
  fBuffer_Track_Label      = new Int_t  [kMaxNumConstituents];
  fBuffer_Track_ProdVtx_X  = new Float_t[kMaxNumConstituents];
  fBuffer_Track_ProdVtx_Y  = new Float_t[kMaxNumConstituents];
  fBuffer_Track_ProdVtx_Z  = new Float_t[kMaxNumConstituents];

  fBuffer_Cluster_Pt       = new Float_t[kMaxNumConstituents];
  fBuffer_Cluster_E        = new Float_t[kMaxNumConstituents];
  fBuffer_Cluster_Eta      = new Float_t[kMaxNumConstituents];
  fBuffer_Cluster_Phi      = new Float_t[kMaxNumConstituents];
  fBuffer_Cluster_M02      = new Float_t[kMaxNumConstituents];
  fBuffer_Cluster_Time     = new Float_t[kMaxNumConstituents];
  fBuffer_Cluster_Label    = new Int_t[kMaxNumConstituents];
}

//________________________________________________________________________
Bool_t AliEmcalJetTree::AddJetToTree(AliEmcalJet* jet, Bool_t saveConstituents, Bool_t saveConstituentsIP, Bool_t saveCaloClusters, Double_t* vertex, Float_t rho, Float_t rhoMass, Float_t centrality, Int_t multiplicity, Long64_t eventID, Float_t magField)
{
  if(!fInitialized)
    AliFatal("Tree is not initialized.");

  fBuffer_JetPt                                   = jet->Pt() - rho*jet->Area();

  // Check if jet type is contained in extraction list
  if( (fExtractionJetTypes_PM.size() || fExtractionJetTypes_HM.size()) &&
      (std::find(fExtractionJetTypes_PM.begin(), fExtractionJetTypes_PM.end(), fBuffer_Jet_MC_MotherParton) == fExtractionJetTypes_PM.end()) &&
      (std::find(fExtractionJetTypes_HM.begin(), fExtractionJetTypes_HM.end(), fBuffer_Jet_MC_MotherHadron) == fExtractionJetTypes_HM.end()) )
    return kFALSE;

  // Check acceptance percentage for the given jet and discard statistically on demand
  Bool_t inPtRange = kFALSE;
  for(size_t i = 0; i<fExtractionPercentages.size(); i++)
  {
    if(fBuffer_JetPt>=fExtractionPercentagePtBins[i*2] && fBuffer_JetPt<fExtractionPercentagePtBins[i*2+1])
    {
      inPtRange = kTRUE;
      if(fRandomGenerator->Rndm() >= fExtractionPercentages[i])
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
  fBuffer_Event_BackgroundDensity               = rho;
  fBuffer_Event_BackgroundDensityMass           = rhoMass;
  fBuffer_Event_Vertex_X                        = vertex ? vertex[0] : 0;
  fBuffer_Event_Vertex_Y                        = vertex ? vertex[1] : 0;
  fBuffer_Event_Vertex_Z                        = vertex ? vertex[2] : 0;
  fBuffer_Event_Centrality                      = centrality;
  fBuffer_Event_Multiplicity                    = multiplicity;
  fBuffer_Event_ID                              = eventID;
  fBuffer_Event_MagneticField                   = magField;

  // Extract basic constituent track properties directly from AliEmcalJet object
  fBuffer_NumTracks = 0;
  if(saveConstituents || saveConstituentsIP)
    for(Int_t i = 0; i < jet->GetNumberOfParticleConstituents(); i++)
    {
      const AliVParticle* particle = jet->GetParticleConstituents()[i].GetParticle();
      if(!particle) continue;

      if(saveConstituents)
      {
        fBuffer_Track_Pt[fBuffer_NumTracks] = particle->Pt();
        fBuffer_Track_Eta[fBuffer_NumTracks] = particle->Eta();
        fBuffer_Track_Phi[fBuffer_NumTracks] = particle->Phi();
        fBuffer_Track_Charge[fBuffer_NumTracks] = particle->Charge();
        fBuffer_Track_Label[fBuffer_NumTracks] = particle->GetLabel();
      }
      if(saveConstituentsIP) // track production vertices
      {
        fBuffer_Track_ProdVtx_X[fBuffer_NumTracks] = particle->Xv();
        fBuffer_Track_ProdVtx_Y[fBuffer_NumTracks] = particle->Yv();
        fBuffer_Track_ProdVtx_Z[fBuffer_NumTracks] = particle->Zv();
      }
      fBuffer_NumTracks++;
    }

  // Extract basic constituent cluster properties directly from AliEmcalJet object
  fBuffer_NumClusters = 0;
  if(saveCaloClusters)
    for(Int_t i = 0; i < jet->GetNumberOfClusterConstituents(); i++)
    {
      const AliVCluster* cluster = jet->GetClusterConstituents()[i].GetCluster();
      if(!cluster) continue;

      // #### Retrieve cluster pT
      TLorentzVector clusterMomentum;
      cluster->GetMomentum(clusterMomentum, vertex);
      // ####

      fBuffer_Cluster_Pt[fBuffer_NumClusters] = clusterMomentum.Perp();
      fBuffer_Cluster_E[fBuffer_NumClusters] = cluster->E();
      fBuffer_Cluster_Eta[fBuffer_NumClusters] = clusterMomentum.Eta();
      fBuffer_Cluster_Phi[fBuffer_NumClusters] = clusterMomentum.Phi();
      fBuffer_Cluster_M02[fBuffer_NumClusters] = cluster->GetM02();
      fBuffer_Cluster_Time[fBuffer_NumClusters] = cluster->GetTOF();
      fBuffer_Cluster_Label[fBuffer_NumClusters] = cluster->GetLabel();
      fBuffer_NumClusters++;
    }


  // Add all buffers to tree
  fJetTree->Fill();

  return kTRUE;
}

//________________________________________________________________________
void AliEmcalJetTree::FillBuffer_TriggerTracks(std::vector<Float_t>& triggerTrackPt, std::vector<Float_t>& triggerTrackDeltaEta, std::vector<Float_t>& triggerTrackDeltaPhi)
{
  fBuffer_NumTriggerTracks = triggerTrackPt.size();
  fJetTree->SetBranchAddress("Jet_TriggerTrack_Pt", triggerTrackPt.data());
  fJetTree->SetBranchAddress("Jet_TriggerTrack_dEta", triggerTrackDeltaEta.data());
  fJetTree->SetBranchAddress("Jet_TriggerTrack_dPhi", triggerTrackDeltaPhi.data());
}

//________________________________________________________________________
void AliEmcalJetTree::FillBuffer_ImpactParameters(std::vector<Float_t>& trackIP_d0, std::vector<Float_t>& trackIP_z0, std::vector<Float_t>& trackIP_d0cov, std::vector<Float_t>& trackIP_z0cov)
{
  fJetTree->SetBranchAddress("Jet_Track_CovIPd", trackIP_d0cov.data());
  fJetTree->SetBranchAddress("Jet_Track_CovIPz", trackIP_z0cov.data());
  fJetTree->SetBranchAddress("Jet_Track_IPd", trackIP_d0.data());
  fJetTree->SetBranchAddress("Jet_Track_IPz", trackIP_z0.data());
}

//________________________________________________________________________
void AliEmcalJetTree::FillBuffer_MonteCarlo(Int_t motherParton, Int_t motherHadron, Int_t partonInitialCollision, Float_t matchedJetDistance_Det, Float_t matchedJetPt_Det, Float_t matchedJetMass_Det, Float_t matchedJetAngularity_Det, Float_t matchedJetpTD_Det, Float_t matchedJetDistance_Part, Float_t matchedJetPt_Part, Float_t matchedJetMass_Part, Float_t matchedJetAngularity_Part, Float_t matchedJetpTD_Part, Float_t truePtFraction, Float_t truePtFraction_PartLevel, Float_t ptHard, Float_t eventWeight, Float_t impactParameter)
{
  fBuffer_Jet_MC_MotherParton = motherParton;
  fBuffer_Jet_MC_MotherHadron = motherHadron;
  fBuffer_Jet_MC_MotherIC = partonInitialCollision;
  fBuffer_Jet_MC_MatchedDetLevelJet_Distance = matchedJetDistance_Det;
  fBuffer_Jet_MC_MatchedDetLevelJet_Pt = matchedJetPt_Det;
  fBuffer_Jet_MC_MatchedDetLevelJet_Mass = matchedJetMass_Det;
  fBuffer_Jet_MC_MatchedDetLevelJet_Angularity = matchedJetAngularity_Det;
  fBuffer_Jet_MC_MatchedDetLevelJet_pTD = matchedJetpTD_Det;

  fBuffer_Jet_MC_MatchedPartLevelJet_Distance = matchedJetDistance_Part;
  fBuffer_Jet_MC_MatchedPartLevelJet_Pt = matchedJetPt_Part;
  fBuffer_Jet_MC_MatchedPartLevelJet_Mass = matchedJetMass_Part;
  fBuffer_Jet_MC_MatchedPartLevelJet_Angularity = matchedJetAngularity_Part;
  fBuffer_Jet_MC_MatchedPartLevelJet_pTD = matchedJetpTD_Part;

  fBuffer_Jet_MC_TruePtFraction = truePtFraction;
  fBuffer_Jet_MC_TruePtFraction_PartLevel = truePtFraction_PartLevel;


  fBuffer_Event_PtHard = ptHard;
  fBuffer_Event_Weight = eventWeight;
  fBuffer_Event_ImpactParameter = impactParameter;
}
//________________________________________________________________________
void AliEmcalJetTree::FillBuffer_PID(std::vector<Float_t>& trackPID_ITS, std::vector<Float_t>& trackPID_TPC, std::vector<Float_t>& trackPID_TOF, std::vector<Float_t>& trackPID_TRD, std::vector<Short_t>& trackPID_Reco, std::vector<Int_t>& trackPID_Truth)
{
  fJetTree->SetBranchAddress("Jet_Track_PID_ITS", trackPID_ITS.data());
  fJetTree->SetBranchAddress("Jet_Track_PID_TPC", trackPID_TPC.data());
  fJetTree->SetBranchAddress("Jet_Track_PID_TOF", trackPID_TOF.data());
  fJetTree->SetBranchAddress("Jet_Track_PID_TRD", trackPID_TRD.data());
  fJetTree->SetBranchAddress("Jet_Track_PID_Reconstructed", trackPID_Reco.data());
  if(trackPID_Truth.data())
    fJetTree->SetBranchAddress("Jet_Track_PID_Truth", trackPID_Truth.data());
}

//________________________________________________________________________
void AliEmcalJetTree::FillBuffer_JetShapes(AliEmcalJet* jet, Double_t leSub_noCorr, Double_t angularity, Double_t momentumDispersion, Double_t trackPtMean, Double_t trackPtMedian)
{
  fBuffer_Shape_Mass_NoCorr = jet->M();
  fBuffer_Shape_Mass_DerivCorr_1 = jet->GetShapeProperties()->GetFirstOrderSubtracted();
  fBuffer_Shape_Mass_DerivCorr_2 = jet->GetShapeProperties()->GetSecondOrderSubtracted();
  fBuffer_Shape_pTD_DerivCorr_1  = jet->GetShapeProperties()->GetFirstOrderSubtractedpTD();
  fBuffer_Shape_pTD_DerivCorr_2  = jet->GetShapeProperties()->GetSecondOrderSubtractedpTD();
  fBuffer_Shape_Angularity_NoCorr = angularity;
  fBuffer_Shape_Angularity_DerivCorr_1  = jet->GetShapeProperties()->GetFirstOrderSubtractedAngularity();
  fBuffer_Shape_Angularity_DerivCorr_2  = jet->GetShapeProperties()->GetSecondOrderSubtractedAngularity();
  fBuffer_Shape_Circularity_DerivCorr_1  = jet->GetShapeProperties()->GetFirstOrderSubtractedCircularity();
  fBuffer_Shape_Circularity_DerivCorr_2  = jet->GetShapeProperties()->GetSecondOrderSubtractedCircularity();
  fBuffer_Shape_LeSub_DerivCorr = jet->GetShapeProperties()->GetSecondOrderSubtractedLeSub();
  fBuffer_Shape_Sigma2_DerivCorr_1 = jet->GetShapeProperties()->GetFirstOrderSubtractedSigma2();
  fBuffer_Shape_Sigma2_DerivCorr_2 = jet->GetShapeProperties()->GetSecondOrderSubtractedSigma2();
  fBuffer_Shape_NumTracks_DerivCorr = jet->GetShapeProperties()->GetSecondOrderSubtractedConstituent();
  fBuffer_Shape_LeSub_NoCorr = leSub_noCorr;
  fBuffer_Shape_MomentumDispersion = momentumDispersion;
  fBuffer_Shape_TrackPtMean = trackPtMean;
  fBuffer_Shape_TrackPtMedian = trackPtMedian;
}

//________________________________________________________________________
void AliEmcalJetTree::FillBuffer_Splittings(std::vector<Float_t>& splittings_radiatorE, std::vector<Float_t>& splittings_kT, std::vector<Float_t>& splittings_theta, Bool_t saveSecondaryVertices, std::vector<Int_t>& splittings_secVtx_rank, std::vector<Int_t>& splittings_secVtx_index)
{
  fBuffer_NumSplittings = splittings_radiatorE.size();
  fJetTree->SetBranchAddress("Jet_Splitting_RadiatorE", splittings_radiatorE.data());
  fJetTree->SetBranchAddress("Jet_Splitting_kT", splittings_kT.data());
  fJetTree->SetBranchAddress("Jet_Splitting_Theta", splittings_theta.data());
  if(saveSecondaryVertices)
  {
    fJetTree->SetBranchAddress("Jet_Splitting_SecVtx_Rank", splittings_secVtx_rank.data());
    fJetTree->SetBranchAddress("Jet_Splitting_SecVtx_Index", splittings_secVtx_index.data());
  }
}

//________________________________________________________________________
void AliEmcalJetTree::FillBuffer_SecVertices(std::vector<Float_t>& secVtx_X, std::vector<Float_t>& secVtx_Y, std::vector<Float_t>& secVtx_Z, std::vector<Float_t>& secVtx_Mass, std::vector<Float_t>& secVtx_Lxy, std::vector<Float_t>& secVtx_SigmaLxy, std::vector<Float_t>& secVtx_Chi2, std::vector<Float_t>& secVtx_Dispersion)
{
  fBuffer_NumSecVertices = secVtx_X.size();
  fJetTree->SetBranchAddress("Jet_SecVtx_X", secVtx_X.data());
  fJetTree->SetBranchAddress("Jet_SecVtx_Y", secVtx_Y.data());
  fJetTree->SetBranchAddress("Jet_SecVtx_Z", secVtx_Z.data());
  fJetTree->SetBranchAddress("Jet_SecVtx_Mass", secVtx_Mass.data());
  fJetTree->SetBranchAddress("Jet_SecVtx_Lxy", secVtx_Lxy.data());
  fJetTree->SetBranchAddress("Jet_SecVtx_SigmaLxy", secVtx_SigmaLxy.data());
  fJetTree->SetBranchAddress("Jet_SecVtx_Chi2", secVtx_Chi2.data());
  fJetTree->SetBranchAddress("Jet_SecVtx_Dispersion", secVtx_Dispersion.data());
}

//________________________________________________________________________
void AliEmcalJetTree::InitializeTree(Bool_t saveCaloClusters, Bool_t saveMCInformation, Bool_t saveMatchedJets_Det, Bool_t saveMatchedJets_Part, Bool_t saveConstituents, Bool_t saveConstituentsIP, Bool_t saveConstituentPID, Bool_t saveJetShapes, Bool_t saveSplittings, Bool_t saveSecondaryVertices, Bool_t saveTriggerTracks)
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
  fJetTree->Branch("Jet_NumTracks",&fBuffer_NumTracks,"Jet_NumTracks/I");
  if(saveCaloClusters)
    fJetTree->Branch("Jet_NumClusters",&fBuffer_NumClusters,"Jet_NumClusters/I");

  fJetTree->Branch("Event_BackgroundDensity",&fBuffer_Event_BackgroundDensity,"Event_BackgroundDensity/F");
  fJetTree->Branch("Event_BackgroundDensityMass",&fBuffer_Event_BackgroundDensityMass,"Event_BackgroundDensityMass/F");
  fJetTree->Branch("Event_Vertex_X",&fBuffer_Event_Vertex_X,"Event_Vertex_X/F");
  fJetTree->Branch("Event_Vertex_Y",&fBuffer_Event_Vertex_Y,"Event_Vertex_Y/F");
  fJetTree->Branch("Event_Vertex_Z",&fBuffer_Event_Vertex_Z,"Event_Vertex_Z/F");
  fJetTree->Branch("Event_Centrality",&fBuffer_Event_Centrality,"Event_Centrality/F");
  fJetTree->Branch("Event_Multiplicity",&fBuffer_Event_Multiplicity,"Event_Multiplicity/I");
  fJetTree->Branch("Event_ID",&fBuffer_Event_ID,"Event_ID/L");
  fJetTree->Branch("Event_MagneticField",&fBuffer_Event_MagneticField,"Event_MagneticField/F");

  if(saveMCInformation)
  {
    fJetTree->Branch("Event_PtHard",&fBuffer_Event_PtHard,"Event_PtHard/F");
    fJetTree->Branch("Event_Weight",&fBuffer_Event_Weight,"Event_Weight/F");
    fJetTree->Branch("Event_ImpactParameter",&fBuffer_Event_ImpactParameter,"Event_ImpactParameter/F");
  }

  if(saveConstituents)
  {
    fJetTree->Branch("Jet_Track_Pt",fBuffer_Track_Pt,"Jet_Track_Pt[Jet_NumTracks]/F");
    fJetTree->Branch("Jet_Track_Phi",fBuffer_Track_Phi,"Jet_Track_Phi[Jet_NumTracks]/F");
    fJetTree->Branch("Jet_Track_Eta",fBuffer_Track_Eta,"Jet_Track_Eta[Jet_NumTracks]/F");
    fJetTree->Branch("Jet_Track_Charge",fBuffer_Track_Charge,"Jet_Track_Charge[Jet_NumTracks]/F");
    if(saveMCInformation)
      fJetTree->Branch("Jet_Track_Label",fBuffer_Track_Label,"Jet_Track_Label[Jet_NumTracks]/I");
  }

  if(saveCaloClusters)
  {
    fJetTree->Branch("Jet_Cluster_Pt",fBuffer_Cluster_Pt,"Jet_Cluster_Pt[Jet_NumClusters]/F");
    fJetTree->Branch("Jet_Cluster_E",fBuffer_Cluster_E,"Jet_Cluster_Pt[Jet_NumClusters]/F");
    fJetTree->Branch("Jet_Cluster_Phi",fBuffer_Cluster_Phi,"Jet_Cluster_Phi[Jet_NumClusters]/F");
    fJetTree->Branch("Jet_Cluster_Eta",fBuffer_Cluster_Eta,"Jet_Cluster_Eta[Jet_NumClusters]/F");
    fJetTree->Branch("Jet_Cluster_M02",fBuffer_Cluster_M02,"Jet_Cluster_M02[Jet_NumClusters]/F");
    fJetTree->Branch("Jet_Cluster_Time",fBuffer_Cluster_Time,"Jet_Cluster_Time[Jet_NumClusters]/F");
    if(saveMCInformation)
      fJetTree->Branch("Jet_Cluster_Label",fBuffer_Cluster_Label,"Jet_Cluster_Label[Jet_NumClusters]/I");
  }

  if(saveConstituentsIP)
  {
    fJetTree->Branch("Jet_Track_IPd",&dummy,"Jet_Track_IPd[Jet_NumTracks]/F");
    fJetTree->Branch("Jet_Track_IPz",&dummy,"Jet_Track_IPz[Jet_NumTracks]/F");
    fJetTree->Branch("Jet_Track_CovIPd",&dummy,"Jet_Track_CovIPd[Jet_NumTracks]/F");
    fJetTree->Branch("Jet_Track_CovIPz",&dummy,"Jet_Track_CovIPz[Jet_NumTracks]/F");

    fJetTree->Branch("Jet_Track_ProdVtx_X",fBuffer_Track_ProdVtx_X,"Jet_Track_ProdVtx_X[Jet_NumTracks]/F");
    fJetTree->Branch("Jet_Track_ProdVtx_Y",fBuffer_Track_ProdVtx_Y,"Jet_Track_ProdVtx_Y[Jet_NumTracks]/F");
    fJetTree->Branch("Jet_Track_ProdVtx_Z",fBuffer_Track_ProdVtx_Z,"Jet_Track_ProdVtx_Z[Jet_NumTracks]/F");
  }

  if(saveConstituentPID)
  {
    fJetTree->Branch("Jet_Track_PID_ITS",&dummy,"Jet_Track_PID_ITS[Jet_NumTracks]/F");
    fJetTree->Branch("Jet_Track_PID_TPC",&dummy,"Jet_Track_PID_TPC[Jet_NumTracks]/F");
    fJetTree->Branch("Jet_Track_PID_TOF",&dummy,"Jet_Track_PID_TOF[Jet_NumTracks]/F");
    fJetTree->Branch("Jet_Track_PID_TRD",&dummy,"Jet_Track_PID_TRD[Jet_NumTracks]/F");

    fJetTree->Branch("Jet_Track_PID_Reconstructed",&dummy,"Jet_Track_PID_Reconstructed[Jet_NumTracks]/S");
    if(saveMCInformation)
      fJetTree->Branch("Jet_Track_PID_Truth",&dummy,"Jet_Track_PID_Truth[Jet_NumTracks]/I");
  }

  if(saveJetShapes)
  {
    fJetTree->Branch("Jet_Shape_Mass_NoCorr",&fBuffer_Shape_Mass_NoCorr,"Jet_Shape_Mass_NoCorr/F");
    fJetTree->Branch("Jet_Shape_Mass_DerivCorr_1",&fBuffer_Shape_Mass_DerivCorr_1,"Jet_Shape_Mass_DerivCorr_1/F");
    fJetTree->Branch("Jet_Shape_Mass_DerivCorr_2",&fBuffer_Shape_Mass_DerivCorr_2,"Jet_Shape_Mass_DerivCorr_2/F");
    fJetTree->Branch("Jet_Shape_pTD_DerivCorr_1",&fBuffer_Shape_pTD_DerivCorr_1,"Jet_Shape_pTD_DerivCorr_1/F");
    fJetTree->Branch("Jet_Shape_pTD_DerivCorr_2",&fBuffer_Shape_pTD_DerivCorr_2,"Jet_Shape_pTD_DerivCorr_2/F");
    fJetTree->Branch("Jet_Shape_LeSub_NoCorr",&fBuffer_Shape_LeSub_NoCorr,"Jet_Shape_LeSub_NoCorr/F");
    fJetTree->Branch("Jet_Shape_LeSub_DerivCorr",&fBuffer_Shape_LeSub_DerivCorr,"Jet_Shape_LeSub_DerivCorr/F");
    fJetTree->Branch("Jet_Shape_Angularity",&fBuffer_Shape_Angularity_NoCorr,"Jet_Shape_Angularity/F");
    fJetTree->Branch("Jet_Shape_Angularity_DerivCorr_1",&fBuffer_Shape_Angularity_DerivCorr_1,"Jet_Shape_Angularity_DerivCorr_1/F");
    fJetTree->Branch("Jet_Shape_Angularity_DerivCorr_2",&fBuffer_Shape_Angularity_DerivCorr_2,"Jet_Shape_Angularity_DerivCorr_2/F");
    fJetTree->Branch("Jet_Shape_Circularity_DerivCorr_1",&fBuffer_Shape_Circularity_DerivCorr_1,"Jet_Shape_Circularity_DerivCorr_1/F");
    fJetTree->Branch("Jet_Shape_Circularity_DerivCorr_2",&fBuffer_Shape_Circularity_DerivCorr_2,"Jet_Shape_Circularity_DerivCorr_2/F");
    fJetTree->Branch("Jet_Shape_Sigma2_DerivCorr_1",&fBuffer_Shape_Sigma2_DerivCorr_1,"Jet_Shape_Sigma2_DerivCorr_1/F");
    fJetTree->Branch("Jet_Shape_Sigma2_DerivCorr_2",&fBuffer_Shape_Sigma2_DerivCorr_2,"Jet_Shape_Sigma2_DerivCorr_2/F");
    fJetTree->Branch("Jet_Shape_NumTracks_DerivCorr",&fBuffer_Shape_NumTracks_DerivCorr,"Jet_Shape_NumTracks_DerivCorr/F");
    fJetTree->Branch("Jet_Shape_MomentumDispersion",&fBuffer_Shape_MomentumDispersion,"Jet_Shape_MomentumDispersion/F");
    fJetTree->Branch("Jet_Shape_TrackPtMean",&fBuffer_Shape_TrackPtMean,"Jet_Shape_TrackPtMean/F");
    fJetTree->Branch("Jet_Shape_TrackPtMedian",&fBuffer_Shape_TrackPtMedian,"Jet_Shape_TrackPtMedian/F");
  }

  if(saveSplittings)
  {
    fJetTree->Branch("Jet_NumSplittings",&fBuffer_NumSplittings,"Jet_NumSplittings/I");
    fJetTree->Branch("Jet_Splitting_Theta",&dummy,"Jet_Splitting_Theta[Jet_NumSplittings]/F");
    fJetTree->Branch("Jet_Splitting_RadiatorE",&dummy,"Jet_Splitting_RadiatorE[Jet_NumSplittings]/F");
    fJetTree->Branch("Jet_Splitting_kT",&dummy,"Jet_Splitting_kT[Jet_NumSplittings]/F");
    if(saveSecondaryVertices)
    {
      fJetTree->Branch("Jet_Splitting_SecVtx_Rank",&dummy,"Jet_Splitting_SecVtx_Rank[Jet_NumSplittings]/I");
      fJetTree->Branch("Jet_Splitting_SecVtx_Index",&dummy,"Jet_Splitting_SecVtx_Index[Jet_NumSplittings]/I");
    }
  }

  if(saveMCInformation)
  {
    fJetTree->Branch("Jet_MC_MotherParton",&fBuffer_Jet_MC_MotherParton,"Jet_MC_MotherParton/I");
    fJetTree->Branch("Jet_MC_MotherHadron",&fBuffer_Jet_MC_MotherHadron,"Jet_MC_MotherHadron/I");
    fJetTree->Branch("Jet_MC_MotherIC",&fBuffer_Jet_MC_MotherIC,"Jet_MC_MotherIC/I");

    if(saveMatchedJets_Det)
    {
      fJetTree->Branch("Jet_MC_MatchedDetLevelJet_Distance",&fBuffer_Jet_MC_MatchedDetLevelJet_Distance,"Jet_MC_MatchedDetLevelJet_Distance/F");
      fJetTree->Branch("Jet_MC_MatchedDetLevelJet_Pt",&fBuffer_Jet_MC_MatchedDetLevelJet_Pt,"Jet_MC_MatchedDetLevelJet_Pt/F");
      fJetTree->Branch("Jet_MC_MatchedDetLevelJet_Mass",&fBuffer_Jet_MC_MatchedDetLevelJet_Mass,"Jet_MC_MatchedDetLevelJet_Mass/F");
      fJetTree->Branch("Jet_MC_MatchedDetLevelJet_Angularity",&fBuffer_Jet_MC_MatchedDetLevelJet_Angularity,"Jet_MC_MatchedDetLevelJet_Angularity/F");
      fJetTree->Branch("Jet_MC_MatchedDetLevelJet_pTD",&fBuffer_Jet_MC_MatchedDetLevelJet_pTD,"Jet_MC_MatchedDetLevelJet_pTD/F");
    }
    if(saveMatchedJets_Part)
    {
      fJetTree->Branch("Jet_MC_MatchedPartLevelJet_Distance",&fBuffer_Jet_MC_MatchedPartLevelJet_Distance,"Jet_MC_MatchedPartLevelJet_Distance/F");
      fJetTree->Branch("Jet_MC_MatchedPartLevelJet_Pt",&fBuffer_Jet_MC_MatchedPartLevelJet_Pt,"Jet_MC_MatchedPartLevelJet_Pt/F");
      fJetTree->Branch("Jet_MC_MatchedPartLevelJet_Mass",&fBuffer_Jet_MC_MatchedPartLevelJet_Mass,"Jet_MC_MatchedPartLevelJet_Mass/F");
      fJetTree->Branch("Jet_MC_MatchedPartLevelJet_Angularity",&fBuffer_Jet_MC_MatchedPartLevelJet_Angularity,"Jet_MC_MatchedPartLevelJet_Angularity/F");
      fJetTree->Branch("Jet_MC_MatchedPartLevelJet_pTD",&fBuffer_Jet_MC_MatchedPartLevelJet_pTD,"Jet_MC_MatchedPartLevelJet_pTD/F");
    }

    fJetTree->Branch("Jet_MC_TruePtFraction",&fBuffer_Jet_MC_TruePtFraction,"Jet_MC_TruePtFraction/F");
    fJetTree->Branch("Jet_MC_TruePtFraction_PartLevel",&fBuffer_Jet_MC_TruePtFraction_PartLevel,"Jet_MC_TruePtFraction_PartLevel/F");
  }


  if(saveSecondaryVertices)
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
  if(saveTriggerTracks)
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
  fSaveConstituents(0), fSaveConstituentsIP(0), fSaveConstituentPID(0), fSaveJetShapes(0), fSaveJetSplittings(0), fSaveMCInformation(0), fSaveSecondaryVertices(0), fSaveTriggerTracks(0), fSaveCaloClusters(0),
  fEventPercentage(1.0),
  fEventCut_TriggerTrackMinPt(0),
  fEventCut_TriggerTrackMaxPt(0),
  fEventCut_TriggerTrackMinLabel(-9999999),
  fEventCut_TriggerTrackMaxLabel(+9999999),
  fEventCut_TriggerTrackOrigin(-1),
  fTruthMinLabel(0),
  fTruthMaxLabel(100000),
  fHadronMatchingRadius(0.4),
  fJetMatchingRadius(0.3),
  fJetMatchingSharedPtFraction(0.5),
  fMCParticleArrayName("mcparticles"),
  fNeedEmbedClusterContainer(0),
  fRandomSeed(0),
  fRandomSeedCones(0),
  fVertexerCuts(0),
  fSecondaryVertexMaxChi2(1e10),
  fSecondaryVertexMaxDispersion(0.05),
  fCustomStartupScript(),
  fSetEmcalJetFlavour(0),
  fSaveTrackPDGCode(kTRUE),
  fJetTree(0),
  fEventWeight(0.),
  fImpactParameter(0.),
  fMultiplicity(0),
  fTriggerTracks_Pt(),
  fTriggerTracks_Eta(),
  fTriggerTracks_Phi(),
  fMCParticleArray(0),
  fRandomGenerator(0),
  fRandomGeneratorCones(0),
  fVtxTagger(0),
  fIsEmbeddedEvent(kFALSE),
  fDoDetLevelMatching(kFALSE),
  fDoPartLevelMatching(kFALSE),
  fSimpleSecVertices()
{
  fRandomGenerator = new TRandom3();
  fRandomGeneratorCones = new TRandom3();

  SetMakeGeneralHistograms(kTRUE);
  fJetTree = new AliEmcalJetTree(GetName());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskJetExtractor::AliAnalysisTaskJetExtractor(const char *name) :
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fSaveConstituents(0), fSaveConstituentsIP(0), fSaveConstituentPID(0), fSaveJetShapes(0), fSaveJetSplittings(0), fSaveMCInformation(0), fSaveSecondaryVertices(0), fSaveTriggerTracks(0), fSaveCaloClusters(0),
  fEventPercentage(1.0),
  fEventCut_TriggerTrackMinPt(0),
  fEventCut_TriggerTrackMaxPt(0),
  fEventCut_TriggerTrackMinLabel(-9999999),
  fEventCut_TriggerTrackMaxLabel(+9999999),
  fEventCut_TriggerTrackOrigin(-1),
  fTruthMinLabel(0),
  fTruthMaxLabel(100000),
  fHadronMatchingRadius(0.4),
  fJetMatchingRadius(0.3),
  fJetMatchingSharedPtFraction(0.5),
  fMCParticleArrayName("mcparticles"),
  fNeedEmbedClusterContainer(0),
  fRandomSeed(0),
  fRandomSeedCones(0),
  fVertexerCuts(0),
  fSecondaryVertexMaxChi2(1e10),
  fSecondaryVertexMaxDispersion(0.05),
  fCustomStartupScript(),
  fSetEmcalJetFlavour(0),
  fSaveTrackPDGCode(kTRUE),
  fJetTree(0),
  fEventWeight(0.),
  fImpactParameter(0.),
  fMultiplicity(0),
  fTriggerTracks_Pt(),
  fTriggerTracks_Eta(),
  fTriggerTracks_Phi(),
  fMCParticleArray(0),
  fRandomGenerator(0),
  fRandomGeneratorCones(0),
  fVtxTagger(0),
  fIsEmbeddedEvent(kFALSE),
  fDoDetLevelMatching(kFALSE),
  fDoPartLevelMatching(kFALSE),
  fSimpleSecVertices()
{
  fRandomGenerator = new TRandom3();
  fRandomGeneratorCones = new TRandom3();

  SetMakeGeneralHistograms(kTRUE);
  fJetTree = new AliEmcalJetTree(GetName());
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
  if(!GetJetContainer(0))
    AliFatal("Jet input container not found!");
  GetJetContainer(0)->PrintCuts();
  if(!GetParticleContainer(0))
    AliFatal("At least one particle input container needs to be added");
  if(!GetClusterContainer(0) && fSaveCaloClusters)
    AliFatal("Cluster input container not found although cluster extraction demanded!");

  fRandomGenerator->SetSeed(fRandomSeed);
  fRandomGeneratorCones->SetSeed(fRandomSeedCones);

  
  // Activate saving of trigger tracks if this is demanded
  if(fEventCut_TriggerTrackMinPt || fEventCut_TriggerTrackMaxPt)
    fSaveTriggerTracks = kTRUE;

  // Use the jet containers from the task to activate matching
  if (GetJetContainer(1)) {
    fDoDetLevelMatching = kTRUE;
  }
  if (GetJetContainer(2)) {
    fDoPartLevelMatching = kTRUE;
  }
  // ### Initialize the jet tree (settings must all be given at this stage)
  fJetTree->SetRandomGenerator(fRandomGenerator);
  fJetTree->InitializeTree(fSaveCaloClusters, fSaveMCInformation, fDoDetLevelMatching, fDoPartLevelMatching, fSaveConstituents, fSaveConstituentsIP, fSaveConstituentPID, fSaveJetShapes, fSaveJetSplittings, fSaveSecondaryVertices, fSaveTriggerTracks);
  OpenFile(2);
  PostData(2, fJetTree->GetTreePointer());

  // ### Add control histograms (already some created in base task)
  AddHistogram2D<TH2D>("hTrackCount", "Number of tracks in acceptance vs. centrality", "COLZ", 500, 0., 5000., 100, 0, 100, "N tracks","Centrality", "dN^{Events}/dN^{Tracks}");
  AddHistogram2D<TH2D>("hBackgroundPt", "Background p_{T} distribution", "", 150, 0., 150., 100, 0, 100, "Background p_{T} (GeV/c)", "Centrality", "dN^{Events}/dp_{T}");

  AddHistogram2D<TH2D>("hJetPtRaw", "Jets p_{T} distribution (raw)", "COLZ", 300, 0., 300., 100, 0, 100, "p_{T, jet} (GeV/c)", "Centrality", "dN^{Jets}/dp_{T}");
  AddHistogram2D<TH2D>("hJetPt", "Jets p_{T} distribution (background subtracted)", "COLZ", 400, -100., 300., 100, 0, 100, "p_{T, jet} (GeV/c)", "Centrality", "dN^{Jets}/dp_{T}");
  AddHistogram2D<TH2D>("hJetPtExtracted", "Extracted jets p_{T} distribution (background subtracted)", "COLZ", 400, -100., 300., 100, 0, 100, "p_{T, jet} (GeV/c)", "Centrality", "dN^{Jets}/dp_{T}");
  AddHistogram2D<TH2D>("hJetPhiEta", "Jet angular distribution #phi/#eta", "COLZ", 180, 0., 2*TMath::Pi(), 100, -2.5, 2.5, "#phi", "#eta", "dN^{Jets}/d#phi d#eta");
  AddHistogram2D<TH2D>("hJetArea", "Jet area", "COLZ", 200, 0., 2., 100, 0, 100, "Jet A", "Centrality", "dN^{Jets}/dA");
  AddHistogram2D<TH2D>("hConstituentPt", "Jet constituent p_{T} distribution", "COLZ", 400, 0., 300., 100, 0, 100, "p_{T, const} (GeV/c)", "Centrality", "dN^{Const}/dp_{T}");
  AddHistogram2D<TH2D>("hConstituentPhiEta", "Jet constituent relative #phi/#eta distribution", "COLZ", 120, -0.6, 0.6, 120, -0.6, 0.6, "#Delta#phi", "#Delta#eta", "dN^{Const}/d#phi d#eta");
  AddHistogram1D<TH1D>("hExtractionPercentage", "Extracted jets p_{T} distribution (background subtracted)", "e1p", 400, -100., 600., "p_{T, jet} (GeV/c)", "Extracted percentage");

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


  // ### If there is an embedding track container, set flag that an embedded event is used
  for(Int_t iCont=0; iCont<fParticleCollArray.GetEntriesFast(); iCont++)
    if (GetParticleContainer(iCont)->GetIsEmbedding())
      fIsEmbeddedEvent = kTRUE;

  // ### Need to explicitly tell jet container to load rho mass object
  GetJetContainer(0)->LoadRhoMass(InputEvent());

  if (fMCParticleArrayName != "")
  {
    if(!fIsEmbeddedEvent)
      fMCParticleArray = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fMCParticleArrayName.Data()));
    else
    {
      // In case of embedding, the MC particle array needs to be fetched differently
      AliVEvent* event = AliEmcalContainerUtils::GetEvent(InputEvent(), kTRUE);
      fMCParticleArray = dynamic_cast<TClonesArray*>(event->FindListObject(fMCParticleArrayName.Data()));
    }
  }

  // ### Prepare vertexer
  if(fSaveSecondaryVertices)
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

  // ### Add HIJING/JEWEL histograms (in case we need it)
  if(MCEvent())
  {
    if(dynamic_cast<AliGenHijingEventHeader*>(MCEvent()->GenEventHeader()))
      AddHistogram1D<TH1D>("hHIJING_Ncoll", "N_{coll} from HIJING", "e1p", 1000, 0., 5000, "N_{coll}", "dN^{Events}/dN^{N_{coll}}");

    if(dynamic_cast<AliGenHepMCEventHeader*>(MCEvent()->GenEventHeader()))
    {
      AddHistogram1D<TH1D>("hJEWEL_NProduced", "NProduced from JEWEL", "e1p", 1000, 0., 1000, "Nproduced", "dN^{Events}/dN^{Nproduced}");
      AddHistogram1D<TH1D>("hJEWEL_ImpactParameter", "Impact parameter from JEWEL", "e1p", 1000, 0., 1, "IP", "dN^{Events}/dN^{IP}");
      TProfile* evweight = new TProfile("hJEWEL_EventWeight", "Event weight from JEWEL", 1, 0,1);
      fOutput->Add(evweight);
    }
  }

  // ### Execute shell script at startup
  if(!fCustomStartupScript.IsNull())
  {
    TGrid::Connect("alien://");
    TFile::Cp(fCustomStartupScript.Data(), "./myScript.sh");
    gSystem->Exec("bash ./myScript.sh");
  }

  PrintConfig();

}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetExtractor::Run()
{
  // ################################### EVENT SELECTION
  // For debugging
  //auto startTime = std::chrono::high_resolution_clock::now();
  //auto elapsedTime = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - startTime).count();
  //std::cout << Form("%s, Time (start): %E s", GetName(), elapsedTime/100000000.) << std::endl;

  if(!IsTriggerTrackInEvent())
    return kFALSE;
  if(fRandomGenerator->Rndm() >= fEventPercentage)
    return kFALSE;

  // ################################### EVENT PROPERTIES
  FillEventControlHistograms();

  // Load vertex if possible
  Long64_t eventID = 0;
  const AliVVertex* myVertex = InputEvent()->GetPrimaryVertex();
  if(!myVertex && MCEvent())
    myVertex = MCEvent()->GetPrimaryVertex();
  Double_t vtx[3] = {0, 0, 0};
  if(myVertex)
  {
    vtx[0] = myVertex->GetX(); vtx[1] = myVertex->GetY(); vtx[2] = myVertex->GetZ();
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
      FillHistogram("hHIJING_Ncoll", ncoll);
    }
    if(AliGenHepMCEventHeader* mcHeader = dynamic_cast<AliGenHepMCEventHeader*>(MCEvent()->GenEventHeader()))
    {
      fEventWeight = mcHeader->EventWeight();
      fImpactParameter = mcHeader->impact_parameter();

      TProfile* tmpHist = static_cast<TProfile*>(fOutput->FindObject("hJEWEL_EventWeight"));
      tmpHist->Fill(0., fEventWeight);
      FillHistogram("hJEWEL_NProduced", mcHeader->NProduced());
      FillHistogram("hJEWEL_ImpactParameter", fImpactParameter);
    }
  }

  // Perform the matching before the main jet loop 
  if(fDoDetLevelMatching) DoJetMatching();
  
  // ################################### MAIN JET LOOP
  GetJetContainer(0)->ResetCurrentID();
  Int_t jetCount = 0;
  while(AliEmcalJet *jet = GetJetContainer(0)->GetNextAcceptJet())
  {
    // LOCAL BUFFER (Will be set for each jet, only added to tree if jet is accepted )
    std::vector<Float_t> vecSigITS; std::vector<Float_t> vecSigTPC; std::vector<Float_t> vecSigTRD; std::vector<Float_t> vecSigTOF; std::vector<Short_t> vecRecoPID; std::vector<Int_t> vecTruePID;
    std::vector<Float_t> vec_d0; std::vector<Float_t> vec_d0cov; std::vector<Float_t> vec_z0; std::vector<Float_t> vec_z0cov;
    std::vector<Float_t> secVtx_X; std::vector<Float_t> secVtx_Y; std::vector<Float_t> secVtx_Z; std::vector<Float_t> secVtx_Mass; std::vector<Float_t> secVtx_Lxy; std::vector<Float_t> secVtx_SigmaLxy; std::vector<Float_t> secVtx_Chi2; std::vector<Float_t> secVtx_Dispersion;
    std::vector<Float_t> triggerTracks_dEta(fTriggerTracks_Eta);
    std::vector<Float_t> triggerTracks_dPhi(fTriggerTracks_Phi);
    std::vector<Float_t> splittings_radiatorE; std::vector<Float_t> splittings_kT; std::vector<Float_t> splittings_theta; std::vector<Int_t> splittings_secVtx_rank; std::vector<Int_t> splittings_secVtx_index;

    FillJetControlHistograms(jet);
    if(fSaveMCInformation)
    {
      Double_t matchedJetDistance_Det = 0;
      Double_t matchedJetPt_Det = 0;
      Double_t matchedJetMass_Det = 0;
      Double_t matchedJetAngularity_Det = 0;
      Double_t matchedJetpTD_Det = 0;
      Double_t matchedJetDistance_Part = 0;
      Double_t matchedJetPt_Part = 0;
      Double_t matchedJetMass_Part = 0;
      Double_t matchedJetAngularity_Part = 0;
      Double_t matchedJetpTD_Part = 0;
      Double_t truePtFraction = 0;
      Double_t truePtFraction_PartLevel = 0;
      Int_t currentJetType_HM = 0;
      Int_t currentJetType_PM = 0;
      Int_t currentJetType_IC = 0;

      // Get jet type from MC (hadron matching, parton matching definition - for HF jets)
      GetJetType(jet, currentJetType_HM, currentJetType_PM, currentJetType_IC);
      // Get true estimators: for pt, jet mass, ...
      GetTrueJetPtFraction(jet, truePtFraction, truePtFraction_PartLevel);
      
      GetMatchedJetObservables(jet, matchedJetPt_Det, matchedJetPt_Part, matchedJetDistance_Det, matchedJetDistance_Part, matchedJetMass_Det, matchedJetMass_Part, matchedJetAngularity_Det, matchedJetAngularity_Part, matchedJetpTD_Det,matchedJetpTD_Part);

      
      fJetTree->FillBuffer_MonteCarlo(currentJetType_PM,currentJetType_HM,currentJetType_IC,
                                      matchedJetDistance_Det,matchedJetPt_Det,matchedJetMass_Det,matchedJetAngularity_Det,matchedJetpTD_Det,
                                      matchedJetDistance_Part,matchedJetPt_Part,matchedJetMass_Part,matchedJetAngularity_Part,matchedJetpTD_Part,
                                      truePtFraction,truePtFraction_PartLevel,fPtHard,fEventWeight,fImpactParameter);
    }

    // ### CONSTITUENT LOOP: Retrieve PID values + impact parameters
    if(fSaveConstituentPID || fSaveConstituentsIP)
    {
      for(Int_t i = 0; i < jet->GetNumberOfParticleConstituents(); i++)
      {
        const AliVParticle* particle = jet->GetParticleConstituents()[i].GetParticle();
        if(!particle) continue;
        if(fSaveConstituentPID)
        {
          Float_t sigITS = 0; Float_t sigTPC = 0; Float_t sigTOF = 0; Float_t sigTRD = 0; Short_t recoPID = 0; Int_t truePID = 0;
          AddPIDInformation(const_cast<AliVParticle*>(particle), sigITS, sigTPC, sigTOF, sigTRD, recoPID, truePID);
          vecSigITS.push_back(sigITS); vecSigTPC.push_back(sigTPC); vecSigTOF.push_back(sigTOF); vecSigTRD.push_back(sigTRD); vecRecoPID.push_back(recoPID); vecTruePID.push_back(truePID);
          fJetTree->FillBuffer_PID(vecSigITS, vecSigTPC, vecSigTOF, vecSigTRD, vecRecoPID, vecTruePID);
        }
        if(fSaveConstituentsIP)
        {
          Float_t d0 = 0; Float_t d0cov = 0; Float_t z0 = 0; Float_t z0cov = 0;
          GetTrackImpactParameters(myVertex, dynamic_cast<const AliAODTrack*>(particle), d0, d0cov, z0, z0cov);
          vec_d0.push_back(d0); vec_d0cov.push_back(d0cov); vec_z0.push_back(z0); vec_z0cov.push_back(z0cov);
          fJetTree->FillBuffer_ImpactParameters(vec_d0, vec_z0, vec_d0cov, vec_z0cov);
        }
      }
    }

    // Reconstruct secondary vertices
    if(fSaveSecondaryVertices)
    {
      ReconstructSecondaryVertices(myVertex, jet, secVtx_X, secVtx_Y, secVtx_Z, secVtx_Mass, secVtx_Lxy, secVtx_SigmaLxy, secVtx_Chi2, secVtx_Dispersion);
      fJetTree->FillBuffer_SecVertices(secVtx_X, secVtx_Y, secVtx_Z, secVtx_Mass, secVtx_Lxy, secVtx_SigmaLxy, secVtx_Chi2, secVtx_Dispersion);
    }

    // Now change trigger tracks eta/phi to dEta/dPhi relative to jet and save them
    if(fSaveTriggerTracks)
    {
      for(UInt_t i=0; i<triggerTracks_dEta.size(); i++)
      {
        triggerTracks_dEta[i] = jet->Eta() - fTriggerTracks_Eta[i];
        triggerTracks_dPhi[i] = TMath::Min(TMath::Abs(jet->Phi() - fTriggerTracks_Phi[i]), TMath::TwoPi() - TMath::Abs(jet->Phi() - fTriggerTracks_Phi[i]));
        if(  ((TMath::Abs(jet->Phi() - fTriggerTracks_Phi[i]) <= TMath::Pi()) && (jet->Phi() - fTriggerTracks_Phi[i] > 0)) 
          || ((TMath::Abs(jet->Phi() - fTriggerTracks_Phi[i]) > TMath::Pi())  && (jet->Phi() - fTriggerTracks_Phi[i] < 0)) )
          triggerTracks_dPhi[i] = -triggerTracks_dPhi[i];
      }
      fJetTree->FillBuffer_TriggerTracks(fTriggerTracks_Pt, triggerTracks_dEta, triggerTracks_dPhi);
    }

    if(fSaveJetShapes)
    {
      // Calculate jet shapes and set them in the tree (some are retrieved in the tree itself)
      Double_t leSub_noCorr = 0;
      Double_t angularity = 0;
      Double_t momentumDispersion = 0;
      Double_t trackPtMean = 0;
      Double_t trackPtMedian = 0;
      CalculateJetShapes(jet, leSub_noCorr, angularity, momentumDispersion, trackPtMean, trackPtMedian);
      fJetTree->FillBuffer_JetShapes(jet, leSub_noCorr, angularity, momentumDispersion, trackPtMean, trackPtMedian);
    }

    if(fSaveJetSplittings)
    {
      GetJetSplittings(jet, splittings_radiatorE, splittings_kT, splittings_theta, splittings_secVtx_rank, splittings_secVtx_index);
      fJetTree->FillBuffer_Splittings(splittings_radiatorE, splittings_kT, splittings_theta, fSaveSecondaryVertices, splittings_secVtx_rank, splittings_secVtx_index);
    }

    // Fill jet to tree (here adding the minimum properties)
    Bool_t accepted = fJetTree->AddJetToTree(jet, fSaveConstituents, fSaveConstituentsIP, fSaveCaloClusters, vtx, GetJetContainer(0)->GetRhoVal(), GetJetContainer(0)->GetRhoMassVal(), fCent, fMultiplicity, eventID, InputEvent()->GetMagneticField());
    if(accepted)
      FillHistogram("hJetPtExtracted", jet->Pt() - GetJetContainer(0)->GetRhoVal()*jet->Area(), fCent);
    jetCount++;
  }

  // ################################### PARTICLE PROPERTIES
  for(Int_t iCont=0; iCont<fParticleCollArray.GetEntriesFast(); iCont++)
  {
    GetParticleContainer(iCont)->ResetCurrentID();
    while(AliVTrack *track = static_cast<AliVTrack*>(GetParticleContainer(iCont)->GetNextAcceptParticle()))
      FillTrackControlHistograms(track);
  }

  // For debugging
  //elapsedTime = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::high_resolution_clock::now() - startTime).count();
  //std::cout << Form("Time (end): %E s", elapsedTime/100000000.) << std::endl;
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

    // Go through all tracks, in all attached track containers, and check whether trigger tracks can be found
    for(Int_t iCont=0; iCont<fParticleCollArray.GetEntriesFast(); iCont++)
    {
      AliParticleContainer* partCont = GetParticleContainer(iCont);
      partCont->ResetCurrentID();
      if((fEventCut_TriggerTrackOrigin == 1) && !partCont->GetIsEmbedding())
        continue;
      else if((fEventCut_TriggerTrackOrigin == 2) && partCont->GetIsEmbedding())
        continue;

      while(AliVTrack *track = static_cast<AliVTrack*>(partCont->GetNextAcceptParticle()))
      {
        if( (track->GetLabel() >= fEventCut_TriggerTrackMinLabel) && (track->GetLabel() < fEventCut_TriggerTrackMaxLabel) )
          if( (track->Pt() >= fEventCut_TriggerTrackMinPt) && (track->Pt() < fEventCut_TriggerTrackMaxPt) )
          {
            fTriggerTracks_Pt.push_back(track->Pt());
            fTriggerTracks_Eta.push_back(track->Eta());
            fTriggerTracks_Phi.push_back(track->Phi());
          }
      }
    }
    // No particle has been found that fulfills requirement -> Do not accept event
    if(fTriggerTracks_Pt.size() == 0)
      return kFALSE;
  }
  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskJetExtractor::CalculateJetShapes(AliEmcalJet* jet, Double_t& leSub_noCorr, Double_t& angularity, Double_t& momentumDispersion, Double_t& trackPtMean, Double_t& trackPtMedian)
{
  // #### Calculate mean, median of constituents, radial moment (angularity), momentum dispersion, leSub (no correction)
  Double_t jetLeadingHadronPt = -999.;
  Double_t jetSubleadingHadronPt = -999.;
  Double_t jetSummedPt = 0;
  Double_t jetSummedPt2 = 0;
  trackPtMean = 0;
  trackPtMedian = 0;
  angularity = 0;
  momentumDispersion = 0;
  std::vector<PWG::JETFW::AliEmcalParticleJetConstituent> tracks_sorted = jet->GetParticleConstituents();
  std::sort(tracks_sorted.rbegin(), tracks_sorted.rend());
  Int_t     numTracks = tracks_sorted.size();
  if(!numTracks) return;
  Double_t* constPts = new Double_t[numTracks];

  // Loop over all constituents and do jet shape calculations
  for (Int_t i=0;i<numTracks;i++)
  {
    const AliVParticle* particle = tracks_sorted[i].GetParticle();
    trackPtMean += particle->Pt();
    constPts[i] = particle->Pt();
    if(particle->Pt() > jetLeadingHadronPt)
    {
      jetSubleadingHadronPt = jetLeadingHadronPt;
      jetLeadingHadronPt = particle->Pt();
    }
    else if(particle->Pt() > jetSubleadingHadronPt)
      jetSubleadingHadronPt = particle->Pt();

    Double_t deltaR = GetDistance(particle->Eta(), jet->Eta(), particle->Phi(), jet->Phi());
    jetSummedPt += particle->Pt();
    jetSummedPt2 += particle->Pt()*particle->Pt();
    angularity += particle->Pt() * deltaR;
  }

  if(numTracks)
  {
    trackPtMean   /= numTracks;
    trackPtMedian = TMath::Median(numTracks, constPts);
  }

  if(numTracks > 1)
    leSub_noCorr = jetLeadingHadronPt - jetSubleadingHadronPt;
  else
    leSub_noCorr = jetLeadingHadronPt;

  if(jetSummedPt)
  {
    momentumDispersion = TMath::Sqrt(jetSummedPt2)/jetSummedPt;
    angularity /= jetSummedPt;
  }
}

//________________________________________________________________________
void AliAnalysisTaskJetExtractor::GetTrueJetPtFraction(AliEmcalJet* jet, Double_t& truePtFraction, Double_t& truePtFraction_mcparticles)
{
  // #################################################################################
  // ##### FRACTION OF TRUE PT IN JET: Defined as "not from toy"
  Double_t pt_truth = 0.;
  Double_t pt_truth_mcparticles = 0.;
  Double_t pt_all   = 0.;
  truePtFraction = 0;
  truePtFraction_mcparticles = 0;

  // ### Loop over all tracks constituents
  for(Int_t iConst = 0; iConst < jet->GetNumberOfParticleConstituents(); iConst++)
  {
    const AliVParticle* particle = jet->GetParticleConstituents()[iConst].GetParticle();
    if(!particle) continue;
    if(particle->Pt() < 1e-6) continue;

    // Particles marked w/ labels within label range OR explicitly set as embedded tracks are considered to be from truth
    if (  (fIsEmbeddedEvent && jet->GetParticleConstituents()[iConst].IsFromEmbeddedEvent()) ||
          (!fIsEmbeddedEvent && ((particle->GetLabel() >= fTruthMinLabel) && (particle->GetLabel() < fTruthMaxLabel)))  ){
            pt_truth += particle->Pt();
    }
    pt_all += particle->Pt();
  }
  // Get the Primary Vertex
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

  Double_t primVtx[3] = {vtxX, vtxY, vtxZ};

  // ### Loop over all cluster constituents
  for(Int_t iConst = 0; iConst < jet->GetNumberOfClusterConstituents(); iConst++)
  {
    const AliVCluster* cluster = jet->GetClusterConstituents()[iConst].GetCluster();
    if(!cluster) continue;
    if(cluster->E() < 1e-6) continue;

    // #### Retrieve cluster pT
    TLorentzVector clusterMomentum;
    cluster->GetMomentum(clusterMomentum, primVtx);
    // Get the total pT from clusters in the jet
    Double_t ClusterPt = clusterMomentum.Perp();
    pt_all += ClusterPt;
  } // end loop over clusters

  // Calculate the true pt fraction for the case where we have clusters
  if (GetClusterContainer(0)){
    // get the cluster container from the input event corresponding
    AliClusterContainer* clusterCont = 0;
    if(fNeedEmbedClusterContainer){
      // if hybrid event, take external cluster container
      clusterCont = GetClusterContainer(1);
    }
    else{
      clusterCont = GetClusterContainer(0);
    }
    // loop over clusters in the input event
    // check to make sure that we are getting a cluster container
    if(clusterCont){
      for(int k=0; k< clusterCont->GetNClusters(); k++){
        const AliVCluster* cluster = clusterCont->GetAcceptCluster(k);
        if (!cluster)continue;
        TLorentzVector clusterMomentum;
        cluster->GetMomentum(clusterMomentum, primVtx);
        Double_t ClusterPt = clusterMomentum.Perp();
        if(IsClusterInCone(clusterMomentum, jet->Eta(), jet->Phi(), GetJetContainer(0)->GetJetRadius()) ){
          pt_truth += ClusterPt;
        }
      } // end loop over clusters
    }// end if clusterCont
  }

  // ### Loop over all primary (charged) MC particles and check if they have a corresponding track/cluster
  //     Correspondence is checked geometrically, sum of matched particles pT is truth
  Bool_t fulljets = (GetJetContainer(0)->GetJetType() == AliJetContainer::kFullJet);
  Double_t jetRadius = GetJetContainer(0)->GetJetRadius();
  if(fMCParticleArray)
    for(Int_t iPart=0; iPart<fMCParticleArray->GetEntriesFast();iPart++)
    {
      AliAODMCParticle* part = (AliAODMCParticle*)fMCParticleArray->At(iPart);
      if(!part) continue;
      if(!part->IsPhysicalPrimary()) continue;
      if(!fulljets && !part->Charge()) continue;
      if(part->Pt() < 1e-6) continue;

      if(IsTrackInCone(part, jet->Eta(), jet->Phi(), jetRadius))
        pt_truth_mcparticles += part->Pt();
    }

  if(pt_all)
  {
    truePtFraction = (pt_truth/pt_all);
    truePtFraction_mcparticles = (pt_truth_mcparticles/pt_all);
  }
}
//________________________________________________________________________
void AliAnalysisTaskJetExtractor::DoJetMatching(){
   // Perform the matching before the main jet loop (Only if we use the tagger for matching)                                              
  AliJetContainer * jetsHybrid = GetJetContainer(0);
  AliJetContainer * jetsDetLevel = GetJetContainer(1);
  AliJetContainer * jetsPartLevel = GetJetContainer(2);

  
  // Now, begin the actual matching.
  // Hybrid <-> det first
  AliDebugStream(2) << "Matching hybrid to detector level jets.\n";
  // First, we reset the tagging
  for(auto j : jetsHybrid->all()){
    j->ResetMatching();
  }
  for(auto j : jetsDetLevel->all()){
    j->ResetMatching();
  }
  // Next, we perform the matching
  PerformGeometricalJetMatching(*jetsHybrid, *jetsDetLevel, fJetMatchingRadius);
  // Now, begin the next matching stage
  // det <-> particle
  AliDebugStream(2) << "Matching detector level to particle level jets.\n";
  // First, we reset the tagging. We need to reset the det matching again to ensure
  // that it doesn't accidentally keep some latent matches to the hybrid jets.
  for(auto j : jetsDetLevel->all()){
    j->ResetMatching();
  }
  // if we do not need to do particle level matching return
  if(!fDoPartLevelMatching) return;
  
  for(auto j : jetsPartLevel->all()){
    j->ResetMatching();
  }
  // Next, we perform the matching
  PerformGeometricalJetMatching(*jetsHybrid, *jetsDetLevel, fJetMatchingRadius);
  // Now, begin the next matching stage 
  // det <-> particle
  AliDebugStream(2) << "Matching detector level to particle level jets.\n";
  // First, we reset the tagging. We need to reset the det matching again to ensure
  // that it doesn't accidentally keep some latent matches to the hybrid jets.
  for(auto j : jetsDetLevel->all()){
    j->ResetMatching();
  }
  for(auto j : jetsPartLevel->all()){
    j->ResetMatching();
  }
  // Next, we perform the matching
  PerformGeometricalJetMatching(*jetsDetLevel, *jetsPartLevel, fJetMatchingRadius);
}

//________________________________________________________________________
bool AliAnalysisTaskJetExtractor::PerformGeometricalJetMatching(AliJetContainer& contBase,
                                    AliJetContainer& contTag, double maxDist) 
{
  // Note that this function is also utilized in /PWGJE/EMCALJetTasks/UserTasks/AliAnalysisTaskEmcalJetHPerformance.cxx. For more details, see this file.
  // Setup
  const Int_t kNacceptedBase = contBase.GetNAcceptedJets(), kNacceptedTag = contTag.GetNAcceptedJets();
  if (!(kNacceptedBase && kNacceptedTag)) {
    return false;
  }

  // Build up vectors of jet pointers to use when assigning the closest jets.
  // The storages are needed later for applying the tagging, in order to avoid multiple occurrence of jet selection
  std::vector<AliEmcalJet*> jetsBase(kNacceptedBase), jetsTag(kNacceptedTag);

  int countBase(0), countTag(0);
  for (auto jb : contBase.accepted()) {
    jetsBase[countBase] = jb;
    countBase++;
  }
  for (auto jt : contTag.accepted()) {
    jetsTag[countTag] = jt;
    countTag++;
  }

  TArrayI faMatchIndexTag(kNacceptedBase), faMatchIndexBase(kNacceptedTag);
  faMatchIndexBase.Reset(-1);
  faMatchIndexTag.Reset(-1);

  // find the closest distance to the base jet
  countBase = 0;
  for (auto jet1 : contBase.accepted()) {
    double distance = maxDist;

    // Loop over all accepted jets and brute force search for the closest jet.
    // NOTE: current_index() returns the jet index in the underlying array, not
    //       the index within the accepted jets that are returned.
    int contTagAcceptedIndex = 0;
    for (auto jet2 : contTag.accepted()) {
      double dR = jet1->DeltaR(jet2);
      if (dR < distance && dR < maxDist) {
        faMatchIndexTag[countBase] = contTagAcceptedIndex;
        distance = dR;
      }
      contTagAcceptedIndex++;
    }


    countBase++;
  }

  // other way around
  countTag = 0;
  for (auto jet1 : contTag.accepted()) {
    double distance = maxDist;

    // Loop over all accepted jets and brute force search for the closest jet.
    // NOTE: current_index() returns the jet index in the underlying array, not
    //       the index within the accepted jets that are returned.
    int contBaseAcceptedIndex = 0;
    for (auto jet2 : contBase.accepted()) {
      double dR = jet1->DeltaR(jet2);
      if (dR < distance && dR < maxDist) {
        faMatchIndexBase[countTag] = contBaseAcceptedIndex;
        distance = dR;
      }
      contBaseAcceptedIndex++;
    }
    countTag++;
  }

  // check for "true" correlations
  // these are pairs where the base jet is the closest to the tag jet and vice versa
  // As the lists are linear a loop over the outer base jet is sufficient.
  AliDebugStream(1) << "Starting true jet loop: nbase(" << kNacceptedBase << "), ntag(" << kNacceptedTag << ")\n";
  for (int ibase = 0; ibase < kNacceptedBase; ibase++) {
    AliDebugStream(2) << "base jet " << ibase << ": match index in tag jet container " << faMatchIndexTag[ibase]
             << "\n";
    if (faMatchIndexTag[ibase] > -1) {
      AliDebugStream(2) << "tag jet " << faMatchIndexTag[ibase] << ": matched base jet " << faMatchIndexBase[faMatchIndexTag[ibase]] << "\n";
    }
    // We have a true correlation where each jet points to the other.
    if (faMatchIndexTag[ibase] > -1 && faMatchIndexBase[faMatchIndexTag[ibase]] == ibase) {
      AliDebugStream(2) << "found a true match \n";
      AliEmcalJet *jetBase = jetsBase[ibase], *jetTag = jetsTag[faMatchIndexTag[ibase]];
      // We have a valid pair of matched jets, so set the closest jet properties.
      if (jetBase && jetTag) {
        Double_t dR = jetBase->DeltaR(jetTag);
        jetBase->SetClosestJet(jetTag, dR);
        jetTag->SetClosestJet(jetBase, dR);
      }
    }
  }
  return true;
}
//________________________________________________________________________
void AliAnalysisTaskJetExtractor::GetMatchedJetObservables(AliEmcalJet* jet, Double_t& detJetPt, Double_t& partJetPt, Double_t& detJetDistance, Double_t& partJetDistance, Double_t& detJetMass, Double_t& partJetMass, Double_t& detJetAngularity, Double_t& partJetAngularity, Double_t& detJetpTD, Double_t& partJetpTD)
{

  
  // Get the Matched Observables                                                                                                                   
  AliJetContainer * hybridJetCont = GetJetContainer(0);
  AliJetContainer * detJetCont    = GetJetContainer(1);
  AliJetContainer * partJetCont   = GetJetContainer(2);
  
  // hybrid to detector level matching 
  AliEmcalJet * jet2 = jet->ClosestJet();
  if (!jet2) { return; }// if there is no match return.
  Double_t ptJet2 = jet2->Pt() - detJetCont->GetRhoVal() * jet2->Area();
  // This will retrieve the fraction of jet2's momentum in jet1.
  Double_t fraction = hybridJetCont->GetFractionSharedPt(jet);
  if (fraction < fJetMatchingSharedPtFraction) {
    return;
  }

  // Set temp jet shape parameters
  Double_t leSub_noCorr  = 0;
  Double_t trackPtMean   = 0;
  Double_t trackPtMedian = 0;

  // if we are not doing the particle level match, fill all observables here and return
  if(!fDoPartLevelMatching){
    detJetPt          = ptJet2;
    detJetDistance    = jet->DeltaR(jet2);
    detJetMass        = jet2->M();

    // Get Jet Shape parameters
    detJetAngularity       = 0;
    detJetpTD              = 0;
    CalculateJetShapes(jet2, leSub_noCorr, detJetAngularity, detJetpTD, trackPtMean, trackPtMedian);
    return; 
  }
  
  // detector to particle matching
  AliEmcalJet * jet3 = jet2->ClosestJet();
  if(!jet3){return;}
  Double_t ptJet3 = jet3->Pt() - partJetCont->GetRhoVal() * jet3->Area(); 
  // make no fraction cut on the detector to particle

  // if we are doing particle level matching, require that there is both a particle and detector level match
  detJetPt  = ptJet2;
  partJetPt = ptJet3;
  detJetDistance  = jet->DeltaR(jet2);
  partJetDistance = jet2->DeltaR(jet3);
  detJetMass  = jet2->M();
  partJetMass = jet3->M();

  // Get Jet Shape Parameters
  detJetAngularity          = 0;
  detJetpTD                 = 0;
  partJetAngularity         = 0;
  partJetpTD                = 0; 

  CalculateJetShapes(jet2, leSub_noCorr, detJetAngularity, detJetpTD, trackPtMean, trackPtMedian);

  // reset the temp variables
  leSub_noCorr     = 0;
  trackPtMean      = 0;
  trackPtMedian    = 0;

  CalculateJetShapes(jet3, leSub_noCorr, partJetAngularity, partJetpTD, trackPtMean, trackPtMedian);
}


  

//________________________________________________________________________
void AliAnalysisTaskJetExtractor::GetJetType(AliEmcalJet* jet, Int_t& typeHM, Int_t& typePM, Int_t& typeIC)
{
  if(!fMCParticleArray)
    return;

  typeHM = 0;
  typePM = 0;

  AliAODMCParticle* parton[2];
  parton[0] = (AliAODMCParticle*) fVtxTagger->IsMCJetParton(fMCParticleArray, jet, fHadronMatchingRadius);  // method 1 (parton matching)
  parton[1] = (AliAODMCParticle*) fVtxTagger->IsMCJetMeson(fMCParticleArray, jet, fHadronMatchingRadius);   // method 2 (hadron matching)

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

    Double_t delta1R   = GetDistance(parton1eta, jet->Eta(), parton1phi, jet->Phi());
    Double_t delta2R   = GetDistance(parton2eta, jet->Eta(), parton2phi, jet->Phi());

    // Check if one of the partons if closer than matching criterion
    Bool_t matched = (delta1R < GetJetContainer(0)->GetJetRadius()/2.) || (delta2R < GetJetContainer(0)->GetJetRadius()/2.);

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
void AliAnalysisTaskJetExtractor::GetTrackImpactParameters(const AliVVertex* vtx, const AliAODTrack* track, Float_t& d0, Float_t& d0cov, Float_t& z0, Float_t& z0cov)
{
  if (track)
  {
    AliAODTrack localTrack(*track);
    Double_t d0rphiz[2],covd0[3];
    Bool_t isDCA=localTrack.PropagateToDCA(vtx,InputEvent()->GetMagneticField(),3.0,d0rphiz,covd0);
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
  //###########################################################################
  // ********* Calculate secondary vertices
  // Derived from AliAnalysisTaskEmcalJetBtagSV
  secVertexArr = new TClonesArray("AliAODVertex");
  Int_t nDauRejCount = 0;
  Int_t nVtx = fVtxTagger->FindVertices(jet,
                                        GetParticleContainer(0)->GetArray(),
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

  // Loop over all potential secondary vertices
  fSimpleSecVertices.clear();
  for(Int_t i=0; i<secVertexArr->GetEntriesFast(); i++)
  {
    AliAODVertex* secVtx = (AliAODVertex*)(secVertexArr->UncheckedAt(i));

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

    // Internally, save sec. vertices to a list
    // Each secondary vertex is reconstructed from 3 prongs
    SimpleSecondaryVertex vtx;
    vtx.fIndex = secVtx_X.size();
    vtx.fLxy = TMath::Abs(Lxy);
    vtx.fDaughter1 = static_cast<AliVParticle*>(secVtx->GetDaughter(0));
    vtx.fDaughter2 = static_cast<AliVParticle*>(secVtx->GetDaughter(1));
    vtx.fDaughter3 = static_cast<AliVParticle*>(secVtx->GetDaughter(2));
    fSimpleSecVertices.push_back(vtx);

    secVtx_X.push_back(secVtx->GetX()); secVtx_Y.push_back(secVtx->GetY()); secVtx_Z.push_back(secVtx->GetZ()); secVtx_Chi2.push_back(secVtx->GetChi2perNDF());
    secVtx_Dispersion.push_back(dispersion); secVtx_Mass.push_back(mass); secVtx_Lxy.push_back(Lxy); secVtx_SigmaLxy.push_back(sigmaLxy); 
  }

  // Sort simple sec. vertices w/ descending Lxy
  std::sort(fSimpleSecVertices.rbegin(), fSimpleSecVertices.rend(), [](const SimpleSecondaryVertex a, const SimpleSecondaryVertex b) {return a.fLxy<b.fLxy;});
  if(fSimpleSecVertices.size() > 10) fSimpleSecVertices.resize(10);

  secVertexArr->Clear();
  delete secVertexArr;
  delete esdVtx;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskJetExtractor::IsStrangeJet(AliEmcalJet* jet)
{
  // Do hadron matching jet type tagging using mcparticles
  // ... if not explicitly deactivated
  if (fMCParticleArray)
  {
    for(Int_t i=0; i<fMCParticleArray->GetEntriesFast();i++)
    {
      AliAODMCParticle* part = (AliAODMCParticle*)fMCParticleArray->At(i);
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
  fMultiplicity = 0;
  for(Int_t iCont=0; iCont<fParticleCollArray.GetEntriesFast(); iCont++)
    fMultiplicity += GetParticleContainer(iCont)->GetNAcceptedParticles();

  // ### Event control plots
  FillHistogram("hTrackCount", fMultiplicity, fCent);
  FillHistogram("hBackgroundPt", GetJetContainer(0)->GetRhoVal(), fCent);
}

//________________________________________________________________________
void AliAnalysisTaskJetExtractor::FillJetControlHistograms(AliEmcalJet* jet)
{
  // ### Jet control plots
  AliJetContainer* jetContainer = GetJetContainer(0);
  FillHistogram("hJetPtRaw", jet->Pt(), fCent); 
  FillHistogram("hJetPt", jet->Pt() - jetContainer->GetRhoVal()*jet->Area(), fCent);
  FillHistogram("hJetPhiEta", jet->Phi(), jet->Eta());
  FillHistogram("hJetArea", jet->Area(), fCent);

  // ### Jet constituent plots
  for(Int_t i = 0; i < jet->GetNumberOfParticleConstituents(); i++)
  {
    const AliVParticle* particle = jet->GetParticleConstituents()[i].GetParticle();
    if(!particle) continue;

    // Constituent eta/phi (relative to jet)
    Double_t deltaEta = jet->Eta() - particle->Eta();
    Double_t deltaPhi = TMath::Min(TMath::Abs(jet->Phi() - particle->Phi()), TMath::TwoPi() - TMath::Abs(jet->Phi() - particle->Phi()));
    if(!((jet->Phi() - particle->Phi() < 0) && (jet->Phi() - particle->Phi() <= TMath::Pi())))
    deltaPhi = -deltaPhi;

    FillHistogram("hConstituentPt", particle->Pt(), fCent);
    FillHistogram("hConstituentPhiEta", deltaPhi, deltaEta);
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
void AliAnalysisTaskJetExtractor::AddPIDInformation(const AliVParticle* particle, Float_t& sigITS, Float_t& sigTPC, Float_t& sigTOF, Float_t& sigTRD, Short_t& recoPID, Int_t& truePID)
{
  truePID = 9;
  if(!particle) return;

  // If we have AODs, retrieve particle PID signals
  const AliAODTrack* aodtrack = dynamic_cast<const AliAODTrack*>(particle);

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
  if(fMCParticleArray)
  {
    for(Int_t i=0; i<fMCParticleArray->GetEntriesFast();i++)
    {
      AliAODMCParticle* mcParticle = dynamic_cast<AliAODMCParticle*>(fMCParticleArray->At(i));
      if(!mcParticle) continue;

      if (mcParticle->GetLabel() == particle->GetLabel())
      {
        if(fSaveTrackPDGCode)
        {
          truePID = mcParticle->PdgCode();
        }
        else // Use same convention as for PID in AODs
        {
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
        }

        break;
      }
    }
  }
}


//________________________________________________________________________
void AliAnalysisTaskJetExtractor::GetJetSplittings(AliEmcalJet* jet, std::vector<Float_t>& splittings_radiatorE, std::vector<Float_t>& splittings_kT, std::vector<Float_t>& splittings_theta, std::vector<Int_t>& splittings_secVtx_rank, std::vector<Int_t>& splittings_secVtx_index)
{
  // ### Adapted from code in AliAnalysisTaskDmesonJetsSub ###
  // Define jet reclusterizer
  fastjet::JetAlgorithm   jetAlgo(fastjet::cambridge_algorithm);
  Double_t                jetRadius_CA = 1.0;
  fastjet::JetDefinition  jetDefinition(jetAlgo, jetRadius_CA,static_cast<fastjet::RecombinationScheme>(0), fastjet::Best);

  try{
    // Convert jet constituents to vector of fastjet::PseudoJet
    std::vector<fastjet::PseudoJet> particles;
    for(Int_t iConst=0; iConst<jet->GetNumberOfParticleConstituents(); iConst++)
    {
      const AliVParticle* constituent = jet->GetParticleConstituents()[iConst].GetParticle();
      Double_t p[3];
      constituent->PxPyPz(p);
      fastjet::PseudoJet pseudoJet = fastjet::PseudoJet(p[0], p[1], p[2], constituent->E());

      // If constituent is part of the N most significant sec. vertices, mark it
      for(UInt_t iVtx=0; iVtx<fSimpleSecVertices.size(); iVtx++)
      {
        if((constituent == fSimpleSecVertices[iVtx].fDaughter1) || (constituent == fSimpleSecVertices[iVtx].fDaughter2) || (constituent == fSimpleSecVertices[iVtx].fDaughter3))
        {
          pseudoJet.set_user_index(iVtx); // user_index is now index in temp sec. vtx vector
          break;
        }
      }
      particles.push_back(pseudoJet);
    }

    // Perform jet reclustering
    fastjet::ClusterSequence        clusterSeq_CA(particles, jetDefinition);
    std::vector<fastjet::PseudoJet> jets_CA = clusterSeq_CA.inclusive_jets(0);
    jets_CA = sorted_by_pt(jets_CA);
    fastjet::PseudoJet radiator = jets_CA[0];
    fastjet::PseudoJet leadingSubJet;
    fastjet::PseudoJet subleadingSubJet;

    // Iterate through the splitting history of the CA clusterization
    while(radiator.has_parents(leadingSubJet,subleadingSubJet))
    {
      if(leadingSubJet.perp() < subleadingSubJet.perp())
        std::swap(leadingSubJet,subleadingSubJet);

      // Angle theta
      Float_t theta = leadingSubJet.delta_R(subleadingSubJet);
      // Radiator energy
      Float_t radiatorEnergy = leadingSubJet.e()+subleadingSubJet.e();
      // kT
      Float_t kT = subleadingSubJet.perp()*theta;

      // Go through leading subjet constituents and check if there are tracks belonging to one of the ten most significant sec. vertices
      Int_t secVtx_rank  = -1; // rank  = nth most displaced
      Int_t secVtx_index = -1; // index = index in sec. vertex array
      if(fSimpleSecVertices.size())
      {
        std::vector<fastjet::PseudoJet> leadingConsts = leadingSubJet.constituents();
        for(UInt_t iLeadingConst=0; iLeadingConst<leadingConsts.size(); iLeadingConst++)
        {
          if(leadingConsts[iLeadingConst].user_index()>=0)
          {
            secVtx_rank  = leadingConsts[iLeadingConst].user_index();
            secVtx_index = fSimpleSecVertices[leadingConsts[iLeadingConst].user_index()].fIndex;
          }
        }
      }

      // Now add splitting properties to result vectors
      splittings_radiatorE.push_back(radiatorEnergy);
      splittings_theta.push_back(theta);
      splittings_kT.push_back(kT);
      splittings_secVtx_rank.push_back(secVtx_rank);
      splittings_secVtx_index.push_back(secVtx_index);

      // Continue with leadingSubJet as new radiator
      radiator=leadingSubJet;
    }
  } catch (fastjet::Error) { /*return -1;*/ }
}


//________________________________________________________________________
void AliAnalysisTaskJetExtractor::PrintConfig()
{
  std::cout << "#########################################\n";
  std::cout << "Settings for extractor task " << GetName() << std::endl;
  std::cout << std::endl;

  std::cout << "### Event cuts ###" << std::endl;
  std::cout << std::endl;
  if(fEventCut_TriggerTrackMinPt || fEventCut_TriggerTrackMaxPt)
   std::cout << Form("* Trigger track requirement: pT=%3.1f-%3.1f GeV/c, labels=%d-%d", fEventCut_TriggerTrackMinPt, fEventCut_TriggerTrackMaxPt, fEventCut_TriggerTrackMinLabel, fEventCut_TriggerTrackMaxLabel) << std::endl;
  if(fEventCut_TriggerTrackOrigin>0)
   std::cout << Form("* Trigger track event cut: %s", (fEventCut_TriggerTrackOrigin==1)? "From embedded event" : "Not from embedded event") << std::endl;
  if(fEventPercentage < 1)
   std::cout << Form("* Artificially lowered event efficiency: %f", fEventPercentage) << std::endl;
  if(fRandomSeed || fRandomSeedCones)
   std::cout << Form("* Random seeds explicitly set: %lu (for event/jet eff.), %lu (RCs)", fRandomSeed, fRandomSeedCones) << std::endl;
  std::cout << std::endl;

  std::cout << "### Jet properties ###" << std::endl;
  std::cout << "* Jet array name = " << GetJetContainer(0)->GetName() << std::endl;
  std::cout << "* Rho name = " << GetJetContainer(0)->GetRhoName() << std::endl;
  std::cout << "* Rho mass name = " << GetJetContainer(0)->GetRhoMassName() << std::endl;
  std::cout << std::endl;

  std::cout << "### Tree extraction properties ###" << std::endl;
  std::vector<Float_t> extractionPtBins      = fJetTree->GetExtractionPercentagePtBins();
  std::vector<Float_t> extractionPercentages = fJetTree->GetExtractionPercentages();
  std::vector<Int_t>   extractionHM          = fJetTree->GetExtractionJetTypes_HM();
  std::vector<Int_t>   extractionPM          = fJetTree->GetExtractionJetTypes_PM();
  if(extractionPercentages.size())
  {
    std::cout << "* Extraction bins: (";
    for(Int_t i=0; i<static_cast<Int_t>(extractionPercentages.size()-1); i++)
      std::cout << extractionPtBins[i*2] << "-" << extractionPtBins[i*2+1] << " = " << extractionPercentages[i] << ", ";
    std::cout << extractionPtBins[(extractionPercentages.size()-1)*2] << "-" << extractionPtBins[(extractionPercentages.size()-1)*2+1] << " = " << extractionPercentages[(extractionPercentages.size()-1)];

    std::cout << ")" << std::endl;
  }
  if(extractionHM.size())
  {
    std::cout << "* Extraction of hadron-matched jets with types: (";
    for(Int_t i=0; i<static_cast<Int_t>(extractionHM.size()-1); i++)
      std::cout << extractionHM[i] << ", ";
    std::cout << extractionHM[extractionHM.size()-1];
    std::cout << ")" << std::endl;
  }
  if(extractionPM.size())
  {
    std::cout << "* Extraction of parton-matched jets with types: (";
    for(Int_t i=0; i<static_cast<Int_t>(extractionPM.size()-1); i++)
      std::cout << extractionPM[i] << ", ";
    std::cout << extractionPM[extractionPM.size()-1];
    std::cout << ")" << std::endl;
  }
  std::cout << std::endl;

  std::cout << "### Extracted data groups ###" << std::endl;
  std::cout << "* Basic event properties (ID, vertex, centrality, bgrd. densities, ...)" << std::endl;
  if(fSaveConstituents)
    std::cout << "* Jet constituents, basic properties (pt, eta, phi, charge, ...)" << std::endl;
  if(fSaveConstituentsIP)
    std::cout << "* Jet constituents, IPs" << std::endl;
  if(fSaveConstituentPID)
    std::cout << "* Jet constituents, PID" << std::endl;
  if(fSaveCaloClusters)
    std::cout << "* Jet calorimeter clusters" << std::endl;
  if(fSaveMCInformation)
    std::cout << "* MC information (origin, matched jets, ...)" << std::endl;
  if(fSaveSecondaryVertices)
    std::cout << "* Secondary vertices" << std::endl;
  if(fSaveJetShapes)
    std::cout << "* Jet shapes (jet mass, LeSub, pTD, ...)" << std::endl;
  if(fSaveJetSplittings)
    std::cout << "* Jet splittings (kT, theta, E from iterative CA reclustering)" << std::endl;
  if(fSaveTriggerTracks)
    std::cout << "* Trigger tracks" << std::endl;
  std::cout << std::endl;
  std::cout << "### Further settings ###" << std::endl;
  if(fIsEmbeddedEvent)
    std::cout << Form("* EMCal embedding framework will be used (at least on container has IsEmbedded() true)") << std::endl;
  if(fDoDetLevelMatching)
    std::cout << Form("* Detector level jet matching active with containger: %s  *", GetJetContainer(1)->GetName()) << std::endl;
  if(fDoPartLevelMatching)
    std::cout << Form("* Particle level jet matching active with container: %s *", GetJetContainer(2)->GetName()) << std::endl;
  if(fMCParticleArray)
    std::cout << Form("* Particle level information available (for jet origin calculation, particle code): %s", fMCParticleArrayName.Data()) << std::endl;
  if(extractionHM.size())
    std::cout << Form("* Hadronic b/c-jet matching radius=%3.3f", fHadronMatchingRadius) << std::endl;
  if(fSaveMCInformation)
    std::cout << Form("* True particle label range for pt fraction=%d-%d", fTruthMinLabel, fTruthMaxLabel) << std::endl;
  if(fSaveMCInformation && fSaveTrackPDGCode && fSaveConstituentPID)
    std::cout << Form("* Particle PID codes will be PDG codes") << std::endl;
  else if(fSaveMCInformation && !fSaveTrackPDGCode && fSaveConstituentPID)
    std::cout << Form("* Particle PID codes will follow AliAODPid convention") << std::endl;
  if(fSaveMCInformation && fSetEmcalJetFlavour)
    std::cout << Form("* AliEmcalJet flavour tag is set from HF-jet tagging") << std::endl;

  std::cout << "#########################################\n\n";
}


//########################################################################
// HELPERS
//########################################################################

//________________________________________________________________________
inline Bool_t AliAnalysisTaskJetExtractor::IsTrackInCone(const AliVParticle* track, Double_t eta, Double_t phi, Double_t radius)
{
  Double_t deltaR = GetDistance(track->Eta(), eta, track->Phi(), phi);
  if(deltaR <= radius)
    return kTRUE;

  return kFALSE;
}

//________________________________________________________________________
inline Bool_t AliAnalysisTaskJetExtractor::IsClusterInCone( TLorentzVector clusterMomentum, Double_t eta, Double_t phi,Double_t radius)
{
  Double_t deltaR = GetDistance(clusterMomentum.Eta(), eta, clusterMomentum.Phi(), phi);
  if(deltaR <= radius)
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

// ### ADDTASK MACRO
//________________________________________________________________________
AliAnalysisTaskJetExtractor* AliAnalysisTaskJetExtractor::AddTaskJetExtractor(TString trackArray, TString clusterArray, TString jetArray, TString rhoObject, Double_t jetRadius, AliRDHFJetsCutsVertex* vertexerCuts, const char* taskNameSuffix)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  Double_t    minJetEta          = 0.5;
  Double_t    minJetPt           = 0.15;
  Double_t    minTrackPt         = 0.15;
  Double_t    minJetAreaPerc     = 0.557;
  TString     suffix             = "";
  if(taskNameSuffix != 0)
    suffix = taskNameSuffix;

  // ###### Task name
  TString name("AliAnalysisTaskJetExtractor");
  if (jetArray != "") {
    name += "_";
    name += jetArray;
  }
  if (rhoObject != "") {
    name += "_";
    name += rhoObject;
  }
  if (suffix != "") {
    name += "_";
    name += suffix;
  }

  // ###### Setup task with default settings
  AliAnalysisTaskJetExtractor* myTask = new AliAnalysisTaskJetExtractor(name);
  myTask->SetNeedEmcalGeom(kFALSE);
  myTask->SetVzRange(-10.,10.);

  // Particle container and track pt cut
  AliParticleContainer* trackCont = 0;
  if(trackArray == "mcparticles")
    trackCont = myTask->AddMCParticleContainer(trackArray);
  else if(trackArray =="mctracks")
    trackCont = myTask->AddParticleContainer(trackArray);
  else
    trackCont = myTask->AddTrackContainer(trackArray);
  trackCont->SetParticlePtCut(minTrackPt);

  // Particle container and track pt cut
  AliClusterContainer* clusterCont = 0;
  if(clusterArray != ""){
    clusterCont = myTask->AddClusterContainer(clusterArray);
    clusterCont->SetClusECut(0.);
    clusterCont->SetClusPtCut(0.);
    clusterCont->SetClusNonLinCorrEnergyCut(0.);
    clusterCont->SetClusHadCorrEnergyCut(0.30);
    clusterCont->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }


  // Secondary vertex cuts (default settings from PWGHF)
  // (can be overwritten by using myTask->SetVertexerCuts(cuts) from outside macro)
  AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts", "default");
  esdTrackCuts->SetRequireSigmaToVertex(kFALSE);
  esdTrackCuts->SetMinNClustersTPC(90);
  esdTrackCuts->SetMaxChi2PerClusterTPC(4);
  esdTrackCuts->SetRequireTPCRefit(kTRUE);
  esdTrackCuts->SetRequireITSRefit(kTRUE);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
  esdTrackCuts->SetMinDCAToVertexXY(0.);
  esdTrackCuts->SetEtaRange(-0.8, 0.8);
  esdTrackCuts->SetPtRange(1.0, 1.e10);
  vertexerCuts->AddTrackCuts(esdTrackCuts);
  vertexerCuts->SetNprongs(3);
  vertexerCuts->SetMinPtHardestTrack(1.0);//default 0.3
  vertexerCuts->SetSecVtxWithKF(kFALSE);//default with StrLinMinDist
  vertexerCuts->SetImpParCut(0.);//default 0
  vertexerCuts->SetDistPrimSec(0.);//default 0
  vertexerCuts->SetCospCut(-1);//default -1
  myTask->SetVertexerCuts(vertexerCuts);

  // Jet container
  AliJetContainer *jetCont = myTask->AddJetContainer(jetArray,6,jetRadius);
  if (jetCont) {
    jetCont->SetRhoName(rhoObject);
    jetCont->SetPercAreaCut(minJetAreaPerc);
    jetCont->SetJetPtCut(minJetPt);
    jetCont->SetLeadingHadronType(0);
    jetCont->SetPtBiasJetTrack(minTrackPt);
    jetCont->SetJetEtaLimits(-minJetEta, +minJetEta);
    jetCont->ConnectParticleContainer(trackCont);
    if(clusterCont)
      jetCont->ConnectClusterContainer(clusterCont);
    jetCont->SetMaxTrackPt(1000);
  }

  mgr->AddTask(myTask);

  // ###### Connect inputs/outputs
  mgr->ConnectInput  (myTask, 0,  mgr->GetCommonInputContainer() );
  mgr->ConnectOutput (myTask, 1, mgr->CreateContainer(Form("%s_histos", name.Data()), AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, Form("%s:ChargedJetsHadronCF", mgr->GetCommonFileName())) );
  mgr->ConnectOutput (myTask, 2, mgr->CreateContainer(Form("%s_tree", name.Data()), TTree::Class(), AliAnalysisManager::kOutputContainer, mgr->GetCommonFileName()) );

  return myTask;
}
