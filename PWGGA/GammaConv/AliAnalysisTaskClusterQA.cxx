/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Nicolas Schmidt                                               *
 * Version 1.0                                                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// Class used to create a tree for QA of clusters
//---------------------------------------------
////////////////////////////////////////////////

#include "AliAnalysisTaskClusterQA.h"
#include "TChain.h"
#include "TRandom.h"
#include "AliAnalysisManager.h"
#include "TParticle.h"
#include "TVectorF.h"
#include "AliPIDResponse.h"
#include "TFile.h"
#include "AliESDtrackCuts.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliAODEvent.h"

class iostream;

using namespace std;

ClassImp(AliAnalysisTaskClusterQA)

//________________________________________________________________________
AliAnalysisTaskClusterQA::AliAnalysisTaskClusterQA() : AliAnalysisTaskSE(),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fCorrTaskSetting(""),
  fConversionCuts(NULL),
  fEventCuts(NULL),
  fClusterCutsEMC(NULL),
  fClusterCutsDMC(NULL),
  fMesonCuts(NULL),
  fMinNLMCut(1),
  fMaxNLMCut(1),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fWeightJetJetMC(1),
  fGeomEMCAL(NULL),
  fClusterTree(NULL),
  fIsHeavyIon(kFALSE),
  ffillTree(-100),
  ffillHistograms(kFALSE),
  fOutputList(NULL),
  fIsMC(kFALSE),
  fCorrectForNonlinearity(kFALSE),
  fSaveEventProperties(0),
  fSaveCells(0),
  fSaveSurroundingCells(0),
  fSaveTracks(0),
  fConeRadius(0),
  fMinTrackPt(0),
  fMinClusterEnergy(0),
  fSaveMCInformation(0),
  fExtractionPercentages(),
  fExtractionPercentagePtBins(),
  fBuffer_ClusterE(0),
  fBuffer_ClusterPhi(0),
  fBuffer_ClusterEta(0),
  fBuffer_ClusterIsEMCAL(kFALSE),
  fBuffer_MC_Cluster_Flag(0),
  fBuffer_ClusterNumCells(0),
  fBuffer_LeadingCell_ID(0),
  fBuffer_LeadingCell_E(0),
  fBuffer_LeadingCell_Eta(0),
  fBuffer_LeadingCell_Phi(0),
  fBuffer_ClusterM02(0),
  fBuffer_ClusterM20(0),
  fBuffer_Event_Vertex_X(0),
  fBuffer_Event_Vertex_Y(0),
  fBuffer_Event_Vertex_Z(0),
  fBuffer_Event_Multiplicity(0),
  fBuffer_Event_NumActiveCells(0),
  fBuffer_Cells_ID(0),
  fBuffer_Cells_E(0),
  fBuffer_Cells_RelativeEta(0),
  fBuffer_Cells_RelativePhi(0),
  fBuffer_Surrounding_NCells(0),
  fBuffer_Surrounding_Cells_ID(0),
  fBuffer_Surrounding_Cells_R(0),
  fBuffer_Surrounding_Cells_E(0),
  fBuffer_Surrounding_Cells_RelativeEta(0),
  fBuffer_Surrounding_Cells_RelativePhi(0),
  fBuffer_Surrounding_NTracks(10),
  fBuffer_Surrounding_Tracks_R(0),
  fBuffer_Surrounding_Tracks_Pt(0),
  fBuffer_Surrounding_Tracks_RelativeEta(0),
  fBuffer_Surrounding_Tracks_RelativePhi(0),
  fBuffer_Cluster_MC_Label(-10)
{
  fBuffer_Cells_ID                = new Int_t[kMaxActiveCells];
  fBuffer_Cells_E                 = new Float_t[kMaxActiveCells];
  fBuffer_Cells_RelativeEta               = new Float_t[kMaxActiveCells];
  fBuffer_Cells_RelativePhi               = new Float_t[kMaxActiveCells];
  fBuffer_Surrounding_Cells_ID              = new Int_t[kMaxActiveCells];
  fBuffer_Surrounding_Cells_R               = new Float_t[kMaxActiveCells];
  fBuffer_Surrounding_Cells_E               = new Float_t[kMaxActiveCells];
  fBuffer_Surrounding_Cells_RelativeEta           = new Float_t[kMaxActiveCells];
  fBuffer_Surrounding_Cells_RelativePhi             = new Float_t[kMaxActiveCells];
  fBuffer_Surrounding_Tracks_R           = new Float_t[kMaxNTracks];
  fBuffer_Surrounding_Tracks_Pt           = new Float_t[kMaxNTracks];
  fBuffer_Surrounding_Tracks_RelativeEta           = new Float_t[kMaxNTracks];
  fBuffer_Surrounding_Tracks_RelativePhi             = new Float_t[kMaxNTracks];
}

AliAnalysisTaskClusterQA::AliAnalysisTaskClusterQA(const char *name) : AliAnalysisTaskSE(name),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fCorrTaskSetting(""),
  fConversionCuts(NULL),
  fEventCuts(NULL),
  fClusterCutsEMC(NULL),
  fClusterCutsDMC(NULL),
  fMesonCuts(NULL),
  fMinNLMCut(1),
  fMaxNLMCut(1),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fWeightJetJetMC(1),
  fGeomEMCAL(NULL),
  fClusterTree(NULL),
  fIsHeavyIon(kFALSE),
  ffillTree(-100),
  ffillHistograms(kFALSE),
  fOutputList(NULL),
  fIsMC(kFALSE),
  fCorrectForNonlinearity(kFALSE),
  fSaveEventProperties(0),
  fSaveCells(0),
  fSaveSurroundingCells(0),
  fSaveTracks(0),
  fConeRadius(0),
  fMinTrackPt(0),
  fMinClusterEnergy(0),
  fSaveMCInformation(0),
  fExtractionPercentages(),
  fExtractionPercentagePtBins(),
  fBuffer_ClusterE(0),
  fBuffer_ClusterPhi(0),
  fBuffer_ClusterEta(0),
  fBuffer_ClusterIsEMCAL(kFALSE),
  fBuffer_MC_Cluster_Flag(0),
  fBuffer_ClusterNumCells(0),
  fBuffer_LeadingCell_ID(0),
  fBuffer_LeadingCell_E(0),
  fBuffer_LeadingCell_Eta(0),
  fBuffer_LeadingCell_Phi(0),
  fBuffer_ClusterM02(0),
  fBuffer_ClusterM20(0),
  fBuffer_Event_Vertex_X(0),
  fBuffer_Event_Vertex_Y(0),
  fBuffer_Event_Vertex_Z(0),
  fBuffer_Event_Multiplicity(0),
  fBuffer_Event_NumActiveCells(0),
  fBuffer_Cells_ID(0),
  fBuffer_Cells_E(0),
  fBuffer_Cells_RelativeEta(0),
  fBuffer_Cells_RelativePhi(0),
  fBuffer_Surrounding_NCells(0),
  fBuffer_Surrounding_Cells_ID(0),
  fBuffer_Surrounding_Cells_R(0),
  fBuffer_Surrounding_Cells_E(0),
  fBuffer_Surrounding_Cells_RelativeEta(0),
  fBuffer_Surrounding_Cells_RelativePhi(0),
  fBuffer_Surrounding_NTracks(10),
  fBuffer_Surrounding_Tracks_R(0),
  fBuffer_Surrounding_Tracks_Pt(0),
  fBuffer_Surrounding_Tracks_RelativeEta(0),
  fBuffer_Surrounding_Tracks_RelativePhi(0),
  fBuffer_Cluster_MC_Label(-10)
{
  fBuffer_Cells_ID                = new Int_t[kMaxActiveCells];
  fBuffer_Cells_E                 = new Float_t[kMaxActiveCells];
  fBuffer_Cells_RelativeEta               = new Float_t[kMaxActiveCells];
  fBuffer_Cells_RelativePhi               = new Float_t[kMaxActiveCells];
  fBuffer_Surrounding_Cells_ID              = new Int_t[kMaxActiveCells];
  fBuffer_Surrounding_Cells_R               = new Float_t[kMaxActiveCells];
  fBuffer_Surrounding_Cells_E               = new Float_t[kMaxActiveCells];
  fBuffer_Surrounding_Cells_RelativeEta           = new Float_t[kMaxActiveCells];
  fBuffer_Surrounding_Cells_RelativePhi             = new Float_t[kMaxActiveCells];
  fBuffer_Surrounding_Tracks_R           = new Float_t[kMaxNTracks];
  fBuffer_Surrounding_Tracks_Pt           = new Float_t[kMaxNTracks];
  fBuffer_Surrounding_Tracks_RelativeEta           = new Float_t[kMaxNTracks];
  fBuffer_Surrounding_Tracks_RelativePhi             = new Float_t[kMaxNTracks];
  // Default constructor

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskClusterQA::~AliAnalysisTaskClusterQA()
{
  // default deconstructor
  
}
//________________________________________________________________________
void AliAnalysisTaskClusterQA::UserCreateOutputObjects()
{
  // Create User Output Objects

  if(fOutputList != NULL){
    delete fOutputList;
    fOutputList                         = NULL;
  }
  if(fOutputList == NULL){
    fOutputList                         = new TList();
    fOutputList->SetOwner(kTRUE);
  }
  
  if(ffillHistograms){
    
    if(((AliConvEventCuts*)fEventCuts)->GetCutHistograms()){
      fOutputList->Add(((AliConvEventCuts*)fEventCuts)->GetCutHistograms());
    }
    
    if(((AliCaloPhotonCuts*)fClusterCutsEMC)->GetCutHistograms()){
      fOutputList->Add(((AliCaloPhotonCuts*)fClusterCutsEMC)->GetCutHistograms());
    }
    if(((AliCaloPhotonCuts*)fClusterCutsDMC)->GetCutHistograms()){
      fOutputList->Add(((AliCaloPhotonCuts*)fClusterCutsDMC)->GetCutHistograms());
    }
    if(((AliConversionMesonCuts*)fMesonCuts)->GetCutHistograms()){
      fOutputList->Add(((AliConversionMesonCuts*)fMesonCuts)->GetCutHistograms());
    }
    PostData(1, fOutputList);
  }
  
  fClusterTree = new TTree("ClusterQA","ClusterQA");

  fClusterTree->Branch("Cluster_E",                         &fBuffer_ClusterE,                        "Cluster_E/F");
  fClusterTree->Branch("Cluster_Eta",                       &fBuffer_ClusterEta,                      "Cluster_Eta/F");
  fClusterTree->Branch("Cluster_Phi",                       &fBuffer_ClusterPhi,                      "Cluster_Phi/F");
  fClusterTree->Branch("Cluster_IsEMCAL",                   &fBuffer_ClusterIsEMCAL,                  "Cluster_IsEMCAL/O");
  fClusterTree->Branch("Cluster_NumCells",                  &fBuffer_ClusterNumCells,                 "Cluster_NumCells/I");
  fClusterTree->Branch("Cluster_LeadingCell_ID",            &fBuffer_LeadingCell_ID,                  "Cluster_LeadingCell_ID/I");
  fClusterTree->Branch("Cluster_LeadingCell_E",             &fBuffer_LeadingCell_E,                   "Cluster_LeadingCell_E/F");
  fClusterTree->Branch("Cluster_LeadingCell_Eta",           &fBuffer_LeadingCell_Eta,                 "Cluster_LeadingCell_Eta/F");
  fClusterTree->Branch("Cluster_LeadingCell_Phi",           &fBuffer_LeadingCell_Phi,                 "Cluster_LeadingCell_Phi/F");
  fClusterTree->Branch("Cluster_M02",                       &fBuffer_ClusterM02,                      "Cluster_M02/F");
  fClusterTree->Branch("Cluster_M20",                       &fBuffer_ClusterM20,                      "Cluster_M20/F");
  if(fSaveEventProperties)
  {
    fClusterTree->Branch("Event_Vertex_X",                  &fBuffer_Event_Vertex_X,                  "Event_Vertex_X/F");
    fClusterTree->Branch("Event_Vertex_Y",                  &fBuffer_Event_Vertex_Y,                  "Event_Vertex_Y/F");
    fClusterTree->Branch("Event_Vertex_Z",                  &fBuffer_Event_Vertex_Z,                  "Event_Vertex_Z/F");
    fClusterTree->Branch("Event_Multiplicity",              &fBuffer_Event_Multiplicity,              "Event_Multiplicity/F");
    fClusterTree->Branch("Event_NumActiveCells",            &fBuffer_Event_NumActiveCells,            "Event_NumActiveCells/I");
  }
  
  if(fSaveCells)
  {
    fClusterTree->Branch("Cluster_Cells_ID",                fBuffer_Cells_ID,                         "Cluster_Cells_ID[Cluster_NumCells]/I");
    fClusterTree->Branch("Cluster_Cells_E",                 fBuffer_Cells_E,                          "Cluster_Cells_E[Cluster_NumCells]/F");
    fClusterTree->Branch("Cluster_Cells_RelativeEta",       fBuffer_Cells_RelativeEta,                "Cluster_Cells_RelativeEta[Cluster_NumCells]/F");
    fClusterTree->Branch("Cluster_Cells_RelativePhi",       fBuffer_Cells_RelativePhi,                "Cluster_Cells_RelativePhi[Cluster_NumCells]/F");
  }
  
  if(fSaveSurroundingCells)
  {
    fClusterTree->Branch("Surrounding_NCells",             &fBuffer_Surrounding_NCells,               "Surrounding_NCells/I");
    fClusterTree->Branch("Surrounding_Cells_ID",            fBuffer_Surrounding_Cells_ID,             "Surrounding_Cells_ID[Surrounding_NCells]/I");
    fClusterTree->Branch("Surrounding_Cells_R",             fBuffer_Surrounding_Cells_R,              "Surrounding_Cells_R[Surrounding_NCells]/F");
    fClusterTree->Branch("Surrounding_Cells_E",             fBuffer_Surrounding_Cells_E,              "Surrounding_Cells_E[Surrounding_NCells]/F");
    fClusterTree->Branch("Surrounding_Cells_RelativeEta",   fBuffer_Surrounding_Cells_RelativeEta,    "Surrounding_Cells_RelativeEta[Surrounding_NCells]/F");
    fClusterTree->Branch("Surrounding_Cells_RelativePhi",   fBuffer_Surrounding_Cells_RelativePhi,    "Surrounding_Cells_RelativePhi[Surrounding_NCells]/F");
  }
  if(fSaveTracks)
  {
    fClusterTree->Branch("Surrounding_NTracks",             &fBuffer_Surrounding_NTracks,             "Surrounding_NTracks/I");
    fClusterTree->Branch("Surrounding_Tracks_R",            fBuffer_Surrounding_Tracks_R,             "Surrounding_Tracks_R[Surrounding_NTracks]/F");
    fClusterTree->Branch("Surrounding_Tracks_Pt",           fBuffer_Surrounding_Tracks_Pt,            "Surrounding_Tracks_Pt[Surrounding_NTracks]/F");
    fClusterTree->Branch("Surrounding_Tracks_RelativeEta",  fBuffer_Surrounding_Tracks_RelativeEta,   "Surrounding_Tracks_RelativeEta[Surrounding_NTracks]/F");
    fClusterTree->Branch("Surrounding_Tracks_RelativePhi",  fBuffer_Surrounding_Tracks_RelativePhi,   "Surrounding_Tracks_RelativePhi[Surrounding_NTracks]/F");
  }

  if(fSaveMCInformation)
  {
    fClusterTree->Branch("Cluster_MC_Label",                &fBuffer_Cluster_MC_Label,                         "Cluster_MC_Label/I");
  }

  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  OpenFile(2);
  PostData(2, fClusterTree);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskClusterQA::Notify()
{
  if (fEventCuts->GetPeriodEnum() == AliConvEventCuts::kNoPeriod && ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum() != AliConvEventCuts::kNoPeriod){
    fEventCuts->SetPeriodEnumExplicit(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum());
  } else if (fEventCuts->GetPeriodEnum() == AliConvEventCuts::kNoPeriod ){
    fEventCuts->SetPeriodEnum(fV0Reader->GetPeriodName());
  }
  
  return kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskClusterQA::UserExec(Option_t *){

  Int_t eventQuality                  = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
  if(eventQuality != 0){// Event Not Accepted
    return;
  }
  fInputEvent                         = InputEvent();
  if(fIsMC) fMCEvent                  = MCEvent();

  Int_t eventNotAccepted              = fEventCuts->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon,kFALSE);
  if(eventNotAccepted) return; // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1

  fGeomEMCAL                          = AliEMCALGeometry::GetInstance();
  if(!fGeomEMCAL){ AliFatal("EMCal geometry not initialized!");}
  
  Int_t nclus                         = 0;
  TClonesArray * arrClustersProcess   = NULL;
  // fNCurrentClusterBasic             = 0;
  if(!fCorrTaskSetting.CompareTo("")){
    nclus = fInputEvent->GetNumberOfCaloClusters();
  } else {
    arrClustersProcess                = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(Form("%sClustersBranch",fCorrTaskSetting.Data())));
    if(!arrClustersProcess)
      AliFatal(Form("%sClustersBranch was not found in AliAnalysisTaskGammaCalo! Check the correction framework settings!",fCorrTaskSetting.Data()));
    nclus                             = arrClustersProcess->GetEntries();
  }

  if(nclus == 0)  return;

  ((AliCaloPhotonCuts*)fClusterCutsEMC)->FillHistogramsExtendedQA(fInputEvent,fIsMC);
  // ((AliCaloPhotonCuts*)fClusterCutsDMC)->FillHistogramsExtendedQA(fInputEvent,fIsMC);

  // match tracks to clusters
  ((AliCaloPhotonCuts*)fClusterCutsEMC)->MatchTracksToClusters(fInputEvent,fWeightJetJetMC,kTRUE, fMCEvent);
  // ((AliCaloPhotonCuts*)fClusterCutsDMC)->MatchTracksToClusters(fInputEvent,fWeightJetJetMC,kTRUE, fMCEvent);

  // vertex
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);

  map<Long_t,Int_t> mapIsClusterAccepted;
  map<Long_t,Int_t> mapIsClusterAcceptedWithoutTrackMatch;
  // Loop over EMCal clusters
  for(Long_t i = 0; i < nclus; i++){
    Double_t tempClusterWeight              =  fWeightJetJetMC;
    AliVCluster* clus                       = NULL;
    if(fInputEvent->IsA()==AliESDEvent::Class()){
      if(arrClustersProcess)
        clus                                = new AliESDCaloCluster(*(AliESDCaloCluster*)arrClustersProcess->At(i));
      else
        clus                                = new AliESDCaloCluster(*(AliESDCaloCluster*)fInputEvent->GetCaloCluster(i));
    } else if(fInputEvent->IsA()==AliAODEvent::Class()){
      if(arrClustersProcess)
        clus                                = new AliAODCaloCluster(*(AliAODCaloCluster*)arrClustersProcess->At(i));
      else
        clus                                = new AliAODCaloCluster(*(AliAODCaloCluster*)fInputEvent->GetCaloCluster(i));
    }
    if(!clus) continue;


    if ( !clus->IsEMCAL()){ // for PHOS: cluster->GetType() == AliVCluster::kPHOSNeutral
      delete clus;
      continue;
    }
    if ( clus->E() < fMinClusterEnergy){
      delete clus;
      continue;
    }

    if(!((AliCaloPhotonCuts*)fClusterCutsEMC)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i) && !((AliCaloPhotonCuts*)fClusterCutsDMC)->ClusterIsSelected(clus,fInputEvent,fMCEvent,fIsMC, tempClusterWeight,i)){
      delete clus;
      continue;
    }
    ProcessQATreeCluster(fInputEvent,clus,i);
  }
  PostData(1, fOutputList);
}

///________________________________________________________________________
void AliAnalysisTaskClusterQA::ProcessQATreeCluster(AliVEvent *event, AliVCluster* cluster, Long_t indexCluster){

  Float_t clusPos[3]                      = { 0,0,0 };
  cluster->GetPosition(clusPos);
  TVector3 clusterVector(clusPos[0],clusPos[1],clusPos[2]);
  Double_t etaCluster                     = clusterVector.Eta();
  Double_t phiCluster                     = clusterVector.Phi();
  if (phiCluster < 0) phiCluster          += 2*TMath::Pi();

  // Vertex position x, y, z
  fBuffer_Event_Vertex_X                  = fInputEvent->GetPrimaryVertex()->GetX();
  fBuffer_Event_Vertex_Y                  = fInputEvent->GetPrimaryVertex()->GetY();
  fBuffer_Event_Vertex_Z                  = fInputEvent->GetPrimaryVertex()->GetZ();

  // V0-based multiplicity of the event
  fBuffer_Event_Multiplicity              = fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C();

  // Properties of the current cluster
  fBuffer_ClusterE                        = cluster->E();

  if(etaCluster < 0.67 && etaCluster > -0.67){
    if(phiCluster < 3.2 && phiCluster > 1.4){
      fBuffer_ClusterIsEMCAL = kTRUE;
    }
  }
  if(etaCluster < 0.67 && etaCluster > -0.67){
    if(etaCluster > 0.23 || etaCluster < -0.23){
      if(phiCluster < 5.6 && phiCluster > 4.6){
        fBuffer_ClusterIsEMCAL = kFALSE;
      }
    }
  }

  fBuffer_ClusterPhi                      = phiCluster;
  fBuffer_ClusterEta                      = etaCluster;
  fBuffer_ClusterM02                      = cluster->GetM02();
  fBuffer_ClusterM20                      = cluster->GetM20();

  // Get all cells from the event
  AliVCaloCells* cells                    = NULL;
  if(cluster->IsEMCAL())
    cells                                 = event->GetEMCALCells();
  else
    return;

  // Determine all cells in the event
  fBuffer_Event_NumActiveCells            = cells->GetNumberOfCells();

  // Get the number of cells from the current cluster
  Int_t nCellCluster = cluster->GetNCells();
  fBuffer_ClusterNumCells                 = nCellCluster;

  // Find the leading cell in the cluster and its position
  fBuffer_LeadingCell_ID                  = FindLargestCellInCluster(cluster,fInputEvent);
  fBuffer_LeadingCell_E                  = cells->GetCellAmplitude(fBuffer_LeadingCell_ID);
  Float_t leadcelleta = 0.;
  Float_t leadcellphi = 0.;
  fGeomEMCAL->EtaPhiFromIndex(fBuffer_LeadingCell_ID, leadcelleta, leadcellphi);
  if ( leadcellphi < 0 ) leadcellphi+=TMath::TwoPi();
  fBuffer_LeadingCell_Eta                = leadcelleta;
  fBuffer_LeadingCell_Phi              = leadcellphi;
  // Treat the remaining cells of the cluster and get their relative position compared to the leading cell
  for(Int_t iCell=0;iCell<nCellCluster;iCell++){
    if(iCell<100){ // maximum number of cells for a cluster is set to 100
      fBuffer_Cells_ID[iCell]               = cluster->GetCellAbsId(iCell);
      fBuffer_Cells_E[iCell]                = cells->GetCellAmplitude(cluster->GetCellAbsId(iCell));
      Float_t celleta = 0.;
      Float_t cellphi = 0.;
      fGeomEMCAL->EtaPhiFromIndex(fBuffer_Cells_ID[iCell], celleta, cellphi);
      if ( cellphi < 0 ) cellphi+=TMath::TwoPi();
      fBuffer_Cells_RelativeEta[iCell]      = leadcelleta-celleta;
      fBuffer_Cells_RelativePhi[iCell]   = leadcellphi-cellphi;
    }
  }
  Int_t nActiveCellsSurroundingInR = 0;
  // fill surrounding cell buffer
  for(Int_t aCell=0;aCell<cells->GetNumberOfCells();aCell++){
    // Define necessary variables
    Short_t cellNumber                    = 0;
    Double_t cellAmplitude = 0,  cellTime = 0, cellEFrac = 0;
    Int_t cellMCLabel = 0;
    Float_t surrcelleta = 0.;
    Float_t surrcellphi = 0.;
    // Get Cell ID
    cells->GetCell(aCell,cellNumber,cellAmplitude,cellTime,cellMCLabel,cellEFrac);

    // Get eta and phi for the surounding cells
    fGeomEMCAL->EtaPhiFromIndex(cellNumber, surrcelleta, surrcellphi);
    if ( surrcellphi < 0 ) surrcellphi+=TMath::TwoPi();
    Double_t dR2 = pow(leadcelleta-surrcelleta,2) + pow(leadcellphi-surrcellphi,2);
    // Select those cells that are within fConeRadius of the leading cluster cell
    if( dR2 < fConeRadius){
      fBuffer_Surrounding_Cells_E[aCell]              = cells->GetCellAmplitude(cellNumber);
      fBuffer_Surrounding_Cells_ID[aCell]             = cellNumber;
      fBuffer_Surrounding_Cells_R[aCell]             = dR2;
      fBuffer_Surrounding_Cells_RelativeEta[aCell]      = leadcelleta-surrcelleta;
      fBuffer_Surrounding_Cells_RelativePhi[aCell]   = leadcellphi-surrcellphi;
      nActiveCellsSurroundingInR+=1;
    }
  }
  fBuffer_Surrounding_NCells = nActiveCellsSurroundingInR;
  if(fIsMC) fBuffer_Cluster_MC_Label = MakePhotonCandidates(cluster, cells,indexCluster);

  ProcessTracksAndMatching(cluster,indexCluster);
  // Add everything to the tree
  if (fClusterTree) fClusterTree->Fill();
}

//________________________________________________________________________
void AliAnalysisTaskClusterQA::Terminate(Option_t *){

}

//________________________________________________________________________
Int_t AliAnalysisTaskClusterQA::FindLargestCellInCluster(AliVCluster* cluster, AliVEvent* event){

  const Int_t nCells      = cluster->GetNCells();
  AliVCaloCells* cells    = NULL;

  if(cluster->IsEMCAL())
    cells = event->GetEMCALCells();
  else if(cluster->IsPHOS())
    cells = event->GetPHOSCells();

  Float_t eMax            = 0.;
  Int_t idMax             = -1;

  if (nCells < 1) return idMax;
  for (Int_t iCell = 0;iCell < nCells;iCell++){
    Int_t cellAbsID       = cluster->GetCellsAbsId()[iCell];
    if (cells->GetCellAmplitude(cellAbsID)> eMax){
      eMax                = cells->GetCellAmplitude(cellAbsID);
      idMax               = cellAbsID;
    }
  }
  return idMax;
}

//________________________________________________________________________
void AliAnalysisTaskClusterQA::GetRowAndColumnFromAbsCellID(Int_t cellIndex, Int_t& row, Int_t& column){
  Int_t nSupMod=0, nModule=0, nIphi=0, nIeta=0;
  row=0;
  column=0;
  // Get SM number and relative row/column for SM
  fGeomEMCAL->GetCellIndex(cellIndex, nSupMod,nModule,nIphi,nIeta);
  fGeomEMCAL->GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, row,column);
  row    += nSupMod/2 * 24;
  column += nSupMod%2 * 48;
}

//________________________________________________________________________
Int_t AliAnalysisTaskClusterQA::MakePhotonCandidates(AliVCluster* clus, AliVCaloCells* cells, Long_t indexCluster){
  
  Double_t vertex[3] = {0};
  InputEvent()->GetPrimaryVertex()->GetXYZ(vertex);
  
  TLorentzVector clusterVector;
  clus->GetMomentum(clusterVector,vertex);

  TLorentzVector* tmpvec = new TLorentzVector();
  tmpvec->SetPxPyPzE(clusterVector.Px(),clusterVector.Py(),clusterVector.Pz(),clusterVector.E());

  // convert to AODConversionPhoton
  AliAODConversionPhoton *PhotonCandidate=new AliAODConversionPhoton(tmpvec);
  if(!PhotonCandidate){
    delete clus;
    delete tmpvec;
    if (PhotonCandidate) delete PhotonCandidate;
    return -9;
  }
  // Flag Photon as CaloPhoton
  PhotonCandidate->SetIsCaloPhoton();
  PhotonCandidate->SetCaloClusterRef(indexCluster);

  // get MC label
  if(fIsMC> 0){
    Int_t* mclabelsCluster  = clus->GetLabels();
    Int_t nValidClusters    = 0;
    if (clus->GetNLabels()>0){
      for (Int_t k =0; k< (Int_t)clus->GetNLabels(); k++){
        if (mclabelsCluster[k]>0){
          if (nValidClusters< 50)PhotonCandidate->SetCaloPhotonMCLabel(nValidClusters,mclabelsCluster[k]);
            nValidClusters++;
        }
      }
    }
    PhotonCandidate->SetNCaloPhotonMCLabels(nValidClusters);
  }
    
  AliAODCaloCluster* clusSub1 = new AliAODCaloCluster();
  AliAODCaloCluster* clusSub2 = new AliAODCaloCluster();


  // split clusters according to their shares in the cluster (NLM == 1) needs to be treated differently
  if (fMinNLMCut == 1 && fMaxNLMCut == 1 ){
    Int_t absCellIdFirst    = ((AliCaloPhotonCuts*)fClusterCutsEMC)->FindLargestCellInCluster(clus, fInputEvent);
    Int_t absCellIdSecond   = ((AliCaloPhotonCuts*)fClusterCutsEMC)->FindSecondLargestCellInCluster(clus, fInputEvent);

    ((AliCaloPhotonCuts*)fClusterCutsEMC)->SplitEnergy(absCellIdFirst, absCellIdSecond, clus, fInputEvent, fIsMC, clusSub1, clusSub2);
  } else if (fMinNLMCut > 1 ){
    const Int_t   nc = clus->GetNCells();
    Int_t   absCellIdList[nc];
    
    ((AliCaloPhotonCuts*)fClusterCutsEMC)->SplitEnergy(absCellIdList[0], absCellIdList[1], clus, fInputEvent, fIsMC, clusSub1, clusSub2);
  }
  
  // TLorentzvector with sub cluster 1
  TLorentzVector clusterVector1;
  clusSub1->GetMomentum(clusterVector1,vertex);
  TLorentzVector* tmpvec1 = new TLorentzVector();
  tmpvec1->SetPxPyPzE(clusterVector1.Px(),clusterVector1.Py(),clusterVector1.Pz(),clusterVector1.E());
  // convert to AODConversionPhoton
  AliAODConversionPhoton *PhotonCandidate1=new AliAODConversionPhoton(tmpvec1);
  if(!PhotonCandidate1){
    delete clusSub1;
    delete tmpvec1;
    return -9;
  }
  // Flag Photon as CaloPhoton
  PhotonCandidate1->SetIsCaloPhoton();
  // TLorentzvector with sub cluster 2
  TLorentzVector clusterVector2;
  clusSub2->GetMomentum(clusterVector2,vertex);
  TLorentzVector* tmpvec2 = new TLorentzVector();
  tmpvec2->SetPxPyPzE(clusterVector2.Px(),clusterVector2.Py(),clusterVector2.Pz(),clusterVector2.E());
  // convert to AODConversionPhoton
  AliAODConversionPhoton *PhotonCandidate2=new AliAODConversionPhoton(tmpvec2);
  if(!PhotonCandidate2){
    delete clusSub2;
    delete tmpvec2;
    return -9;
  }
  // Flag Photon as CaloPhoton
  PhotonCandidate2->SetIsCaloPhoton();
  Int_t mclabel = -3;
  // create pi0 candidate
  AliAODConversionMother *pi0cand = new AliAODConversionMother(PhotonCandidate1,PhotonCandidate2);
  if((((AliConversionMesonCuts*)fMesonCuts)->MesonIsSelected(pi0cand,kTRUE,fEventCuts->GetEtaShift()))){
    if(fIsMC> 0 && PhotonCandidate && PhotonCandidate1 && PhotonCandidate2 && fSaveMCInformation){
      // if(fInputEvent->IsA()==AliESDEvent::Class())
        mclabel = ProcessTrueClusterCandidates(PhotonCandidate, clus->GetM02(), PhotonCandidate1, PhotonCandidate2);
      // if(fInputEvent->IsA()==AliAODEvent::Class())
        // ProcessTrueClusterCandidatesAOD(PhotonCandidate, clus->GetM02(), PhotonCandidate1, PhotonCandidate2);
        return mclabel;
      }
  } else {
    return -7;
  }
  return -1;
}

//________________________________________________________________________
void AliAnalysisTaskClusterQA::ProcessTracksAndMatching(AliVCluster* clus, Long_t indexCluster){


  Int_t nTracksInR = 0;
  Int_t nModules = fGeomEMCAL->GetNumberOfSuperModules();

  AliESDEvent *esdev = dynamic_cast<AliESDEvent*>(fInputEvent);
  AliESDtrackCuts *EsdTrackCuts = 0x0;
  // Using standard function for setting Cuts
  static int prevRun = -1;
  // Using standard function for setting Cuts
  if (esdev){
    Int_t runNumber = fInputEvent->GetRunNumber();
    if (prevRun!=runNumber) {
      delete EsdTrackCuts;
      EsdTrackCuts = 0;
      prevRun = runNumber;
    }
    // if LHC11a or earlier or if LHC13g or if LHC12a-i -> use 2010 cuts
    if( (runNumber<=146860) || (runNumber>=197470 && runNumber<=197692) || (runNumber>=172440 && runNumber<=193766) ){
      EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();
    // else if run2 data use 2015 PbPb cuts
    }else if (runNumber>=209122){
      // hard coded track cuts for the moment, because AliESDtrackCuts::GetStandardITSTPCTrackCuts2015PbPb() gives spams warnings
      EsdTrackCuts = new AliESDtrackCuts();
      EsdTrackCuts->AliESDtrackCuts::SetMinNCrossedRowsTPC(70);
      EsdTrackCuts->AliESDtrackCuts::SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
      EsdTrackCuts->AliESDtrackCuts::SetMaxChi2PerClusterTPC(4);
      EsdTrackCuts->AliESDtrackCuts::SetAcceptKinkDaughters(kFALSE);
      EsdTrackCuts->AliESDtrackCuts::SetRequireTPCRefit(kTRUE);
      // ITS
      EsdTrackCuts->AliESDtrackCuts::SetRequireITSRefit(kTRUE);
      EsdTrackCuts->AliESDtrackCuts::SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                                              AliESDtrackCuts::kAny);
      EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
      EsdTrackCuts->AliESDtrackCuts::SetMaxChi2TPCConstrainedGlobal(36);
      EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexZ(2);
      EsdTrackCuts->AliESDtrackCuts::SetDCAToVertex2D(kFALSE);
      EsdTrackCuts->AliESDtrackCuts::SetRequireSigmaToVertex(kFALSE);
      EsdTrackCuts->AliESDtrackCuts::SetMaxChi2PerClusterITS(36);
    // else use 2011 version of track cuts
    }else{
      EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
    }
    EsdTrackCuts->SetMaxDCAToVertexZ(2);
    EsdTrackCuts->SetEtaRange(-0.8, 0.8);
    EsdTrackCuts->SetPtRange(0.15);
  }

  for (Int_t itr=0;itr<fInputEvent->GetNumberOfTracks();itr++){
    AliVTrack *inTrack = esdev->GetTrack(itr);
    if(!inTrack) continue;
    if(inTrack->Pt()<fMinTrackPt) continue;
    AliESDtrack *esdt = dynamic_cast<AliESDtrack*>(inTrack);
    // check track cuts
    if(!EsdTrackCuts->AcceptTrack(esdt)) continue;
    const AliExternalTrackParam *in = esdt->GetInnerParam();
    if (!in){AliError("Could not get InnerParam of Track, continue");continue;}
    AliExternalTrackParam *trackParam =new AliExternalTrackParam(*in);
    AliExternalTrackParam emcParam(*trackParam);
    Float_t eta, phi, pt;

    //propagate tracks to emc surfaces
    if (!AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(&emcParam, 440., 0.139, 20., eta, phi, pt)) {
      delete trackParam;
      continue;
    }
    if( TMath::Abs(eta) > 0.75 ) {
      delete trackParam;
      continue;
    }
    // Save some time and memory in case of no DCal present
    if( nModules < 13 && ( phi < 70*TMath::DegToRad() || phi > 190*TMath::DegToRad())){
      delete trackParam;
      continue;
    }
    // Save some time and memory in case of run2
    if( nModules > 12 ){
      if (( phi < 70*TMath::DegToRad() || phi > 190*TMath::DegToRad()) && ( phi < 250*TMath::DegToRad() || phi > 340*TMath::DegToRad())){
        delete trackParam;
        continue;
      }
    }
    Float_t dEta=-999, dPhi=-999;
    Double_t trkPos[3] = {0.,0.,0.};
    if (!emcParam.GetXYZ(trkPos)){
      delete trackParam;
      continue;
    }

    AliExternalTrackParam trackParamTmp(emcParam);//Retrieve the starting point every time before the extrapolation
    if(!AliEMCALRecoUtils::ExtrapolateTrackToCluster(&trackParamTmp, clus, 0.139, 5., dEta, dPhi)) continue;

    Float_t dR2 = dPhi*dPhi + dEta*dEta;
    if(dR2 < fConeRadius){
      fBuffer_Surrounding_Tracks_R[nTracksInR]=dR2;
      fBuffer_Surrounding_Tracks_Pt[nTracksInR]=inTrack->Pt();
      fBuffer_Surrounding_Tracks_RelativeEta[nTracksInR]=dEta;
      fBuffer_Surrounding_Tracks_RelativePhi[nTracksInR]=dPhi;
      nTracksInR+=1;
    }
  }
  fBuffer_Surrounding_NTracks = nTracksInR;
  if(nTracksInR==0){
    fBuffer_Surrounding_Tracks_R[nTracksInR]=-1;
    fBuffer_Surrounding_Tracks_Pt[nTracksInR]=-1;
    fBuffer_Surrounding_Tracks_RelativeEta[nTracksInR]=-1;
    fBuffer_Surrounding_Tracks_RelativePhi[nTracksInR]=-1;
  }

  if(EsdTrackCuts){
    delete EsdTrackCuts;
    EsdTrackCuts=0x0;
  }
}

//________________________________________________________________________
Int_t AliAnalysisTaskClusterQA::ProcessTrueClusterCandidates(AliAODConversionPhoton *TrueClusterCandidate, Float_t m02,
                                    AliAODConversionPhoton *TrueSubClusterCandidate1,
                                    AliAODConversionPhoton *TrueSubClusterCandidate2)
{
  Int_t mcLabelCluster = -5;
  const AliVVertex* primVtxMC   = fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX   = primVtxMC->GetX();
  Double_t mcProdVtxY   = primVtxMC->GetY();
  Double_t mcProdVtxZ   = primVtxMC->GetZ();

  TParticle *Photon = NULL;
  if (!TrueClusterCandidate->GetIsCaloPhoton()) AliFatal("CaloPhotonFlag has not been set task will abort");
  if (TrueClusterCandidate->GetCaloPhotonMCLabel(0) < 0) return mcLabelCluster;
  if (TrueClusterCandidate->GetNCaloPhotonMCLabels()>0) Photon = fMCEvent->Particle(TrueClusterCandidate->GetCaloPhotonMCLabel(0));
    else return mcLabelCluster;
    
  if(Photon == NULL){
    return mcLabelCluster;
  }
  AliAODConversionMother *mesoncand = new AliAODConversionMother(TrueSubClusterCandidate1,TrueSubClusterCandidate2);
  Bool_t mesonIsSelected            = (((AliConversionMesonCuts*)fMesonCuts)->MesonIsSelected(mesoncand,kTRUE,fEventCuts->GetEtaShift()));
  if (!mesonIsSelected){
    delete mesoncand;
    return mcLabelCluster;
  }

  TrueClusterCandidate->SetCaloPhotonMCFlags(fMCEvent, kFALSE);

  Int_t clusterClass    = 0;
  Bool_t isPrimary      = fEventCuts->IsConversionPrimaryESD( fMCEvent, TrueClusterCandidate->GetCaloPhotonMCLabel(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);

  // cluster classification:
  // 1    - nice merged cluster (2 gamma | contributions from 2 gamma) from pi0/eta
  // 2    - contribution from only 1 partner (1 gamma, 1 fully coverted gamma) from pi0/eta
  // 3    - contribution from part of 1 partner (1 electron) from pi0/eta
  Long_t motherLab = -1;
  if (TrueClusterCandidate->IsMerged() || TrueClusterCandidate->IsMergedPartConv()){
      clusterClass    = 1;
      motherLab       = TrueClusterCandidate->GetCaloPhotonMotherMCLabel(0);
  } else if (TrueClusterCandidate->GetNCaloPhotonMotherMCLabels()> 0){
    if (TrueClusterCandidate->IsLargestComponentElectron() || TrueClusterCandidate->IsLargestComponentPhoton()){
      if (TrueClusterCandidate->GetCaloPhotonMotherMCLabel(0) > -1 && (fMCEvent->Particle(TrueClusterCandidate->GetCaloPhotonMotherMCLabel(0))->GetPdgCode() == 111 || fMCEvent->Particle(TrueClusterCandidate->GetCaloPhotonMotherMCLabel(0))->GetPdgCode() == 221) ){
        if ( TrueClusterCandidate->IsConversion() && !TrueClusterCandidate->IsConversionFullyContained() ){
          clusterClass  = 3;
          motherLab       = TrueClusterCandidate->GetCaloPhotonMotherMCLabel(0);
        } else {
          clusterClass  = 2;
          motherLab       = TrueClusterCandidate->GetCaloPhotonMotherMCLabel(0);
        }
      }
    } else if (TrueClusterCandidate->IsSubLeadingEM()){
      if (TrueClusterCandidate->GetNCaloPhotonMotherMCLabels()> 1){
        if ( TrueClusterCandidate->GetCaloPhotonMotherMCLabel(1) > -1){
          if (fMCEvent->Particle(TrueClusterCandidate->GetCaloPhotonMotherMCLabel(1))->GetPdgCode() == 111 || fMCEvent->Particle(TrueClusterCandidate->GetCaloPhotonMotherMCLabel(1))->GetPdgCode() == 221 ){
            clusterClass  = 2;
            motherLab       = TrueClusterCandidate->GetCaloPhotonMotherMCLabel(1);
          }
        }
      }
    } else {
      motherLab       = TrueClusterCandidate->GetCaloPhotonMotherMCLabel(0);
    }
  }

  // Get Mother particle
  TParticle *mother = NULL;
  Int_t motherPDG   = -1;
  if (motherLab > -1)
     mother           = fMCEvent->Particle(motherLab);
  if (mother)
      motherPDG = TMath::Abs(mother->GetPdgCode());

  //
  if (clusterClass == 1 || clusterClass == 2 || clusterClass == 3 ){
    // separate different components
    if (clusterClass == 1 && TrueClusterCandidate->IsMerged()){
      if (motherPDG == 111){
        mcLabelCluster = 10; // NOTE merged pi0
        if (!isPrimary && m02 >= 0 && m02 <= 4.8 )
          mcLabelCluster = 11;  // NOTE secondary merged pi0
      }
      if (motherPDG == 221)
        mcLabelCluster = 12;  // NOTE merged eta
    } else if (clusterClass == 1 && TrueClusterCandidate->IsMergedPartConv()){
      if (motherPDG == 111){
        mcLabelCluster = 13;  // NOTE merged pi0 with one converted gamma
        if (!isPrimary && m02 >= 0 && m02 <= 4.8 )
          mcLabelCluster = 14;  // NOTE merged secondary pi0 with one converted gamma
      }
      if (motherPDG == 221)
        mcLabelCluster = 15;  // NOTE merged eta with one converted gamma
    } else if (clusterClass == 2){
      if (motherPDG == 111){
        mcLabelCluster = 20;  // NOTE decay photon from pi0
        if (!isPrimary && m02 >= 0 && m02 <= 4.8 )
          mcLabelCluster = 21;  // NOTE decay photon from secondary pi0
      }
      if (motherPDG == 221)
        mcLabelCluster = 22;  // NOTE decay photon from eta
    } else if (clusterClass == 3){
      if (motherPDG == 111) {
        mcLabelCluster = 30;  // NOTE electron from decayed pi0
        if (!isPrimary && m02 >= 0 && m02 <= 4.8 )
          mcLabelCluster = 31;  // NOTE electron from decayed secondary pi0
      }
      if (motherPDG == 221)
        mcLabelCluster = 32;  // NOTE electron from decayed eta
    }
  // leading particle is a photon or the conversion is fully contained and its not from pi0 || eta
  } else if (TrueClusterCandidate->IsLargestComponentPhoton() || TrueClusterCandidate->IsConversionFullyContained()
               || TrueClusterCandidate->IsElectronFromFragPhoton()){
    if (motherLab == -1)
      mcLabelCluster = 40; // NOTE  direct photon
    else
      mcLabelCluster = 41;
  // leading particle is an electron and its not from pi0 || eta and no electron from fragmentation photon conversion
  } else if (TrueClusterCandidate->IsLargestComponentElectron() &&  !TrueClusterCandidate->IsElectronFromFragPhoton()){
    mcLabelCluster = 50; // NOTE cluster from hadron
  } else {
    mcLabelCluster = 60; // NOTE BG cluster
  }
  
  delete mesoncand;
  return mcLabelCluster;
}