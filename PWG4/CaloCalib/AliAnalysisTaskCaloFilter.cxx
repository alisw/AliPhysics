/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
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

/* $Id: AliAnalysisTaskCaloFilter.cxx $ */

//////////////////////////////////////////////////////////
// Filter the ESDCaloClusters and ESDCaloCells of EMCAL,
// PHOS or both, creating the corresponing AODCaloClusters
// and AODCaloCells.
// Keep also the AODHeader information and the vertex.
// Needed for calorimeter calibration.
// Copy of AliAnalysisTaskESDfilter.
// Author: Gustavo Conesa Balbastre (INFN - Frascati)
//////////////////////////////////////////////////////////

#include "AliAnalysisTaskCaloFilter.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliLog.h"

ClassImp(AliAnalysisTaskCaloFilter)
  
////////////////////////////////////////////////////////////////////////

AliAnalysisTaskCaloFilter::AliAnalysisTaskCaloFilter():
  AliAnalysisTaskSE(),
  fCalorimeter("EMCAL PHOS")
{
  // Default constructor
}

//__________________________________________________
AliAnalysisTaskCaloFilter::AliAnalysisTaskCaloFilter(const char* name):
  AliAnalysisTaskSE(name),
  fCalorimeter("EMCAL PHOS")
{
  // Constructor
}

//__________________________________________________
void AliAnalysisTaskCaloFilter::CreateAODFromAOD()
{

  // Copy AOD header, vertex, CaloClusters and CaloCells to output AOD
  
  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(InputEvent());
  
  // set arrays and pointers
  Float_t posF[3];
  Double_t pos[3];
  
  Double_t covVtx[6];
  
  for (Int_t i = 0; i < 6; i++)  covVtx[i] = 0.;
  
  // Update the header
   AliAODHeader* headerin = aod->GetHeader();
  AliAODHeader* header = AODEvent()->GetHeader();
  header->SetRunNumber(headerin->GetRunNumber());
  header->SetBunchCrossNumber(headerin->GetBunchCrossNumber());
  header->SetOrbitNumber(headerin->GetOrbitNumber());
  header->SetPeriodNumber(headerin->GetPeriodNumber());
  header->SetEventType(headerin->GetEventType());
  header->SetMuonMagFieldScale(headerin->GetMuonMagFieldScale()); // FIXME
  header->SetCentrality(headerin->GetCentrality());        // FIXME
  
  
  header->SetTriggerMask(headerin->GetTriggerMask()); 
  header->SetTriggerCluster(headerin->GetTriggerCluster());
  header->SetMagneticField(headerin->GetMagneticField());
  header->SetZDCN1Energy(headerin->GetZDCN1Energy());
  header->SetZDCP1Energy(headerin->GetZDCP1Energy());
  header->SetZDCN2Energy(headerin->GetZDCN2Energy());
  header->SetZDCP2Energy(headerin->GetZDCP2Energy());
  header->SetZDCEMEnergy(headerin->GetZDCEMEnergy(0),headerin->GetZDCEMEnergy(1));
  Float_t diamxy[2]={aod->GetDiamondX(),aod->GetDiamondY()};
  Float_t diamcov[3]; aod->GetDiamondCovXY(diamcov);
  header->SetDiamond(diamxy,diamcov);
  //
  //
  Int_t nVertices = 1 ;/* = prim. vtx*/;
  Int_t nCaloClus = aod->GetNCaloClusters();
  
  AODEvent()->ResetStd(0, nVertices, 0, 0, 0, nCaloClus, 0, 0);
  
  // Access to the AOD container of vertices
  TClonesArray &vertices = *(AODEvent()->GetVertices());
  Int_t jVertices=0;
  
  // Add primary vertex. The primary tracks will be defined
  // after the loops on the composite objects (V0, cascades, kinks)
  const AliAODVertex *vtx = aod->GetPrimaryVertex();
  
  vtx->GetXYZ(pos); // position
  vtx->GetCovMatrix(covVtx); //covariance matrix
  
  AliAODVertex * primary = new(vertices[jVertices++])
    AliAODVertex(pos, covVtx, vtx->GetChi2perNDF(), NULL, -1, AliAODVertex::kPrimary);
  primary->SetName(vtx->GetName());
  primary->SetTitle(vtx->GetTitle());
  
  // Access to the AOD container of clusters
  TClonesArray &caloClusters = *(AODEvent()->GetCaloClusters());
  Int_t jClusters=0;
  
  for (Int_t iClust=0; iClust<nCaloClus; ++iClust) {
    
    AliAODCaloCluster * cluster = aod->GetCaloCluster(iClust);
    
    //Check which calorimeter information we want to keep.
    if     (fCalorimeter.Contains("PHOS")  && !fCalorimeter.Contains("EMCAL") && cluster->IsEMCALCluster()) continue ;
    else if(fCalorimeter.Contains("EMCAL") && !fCalorimeter.Contains("PHOS")  && cluster->IsPHOSCluster())  continue ;
    
    Int_t id       = cluster->GetID();
    Float_t energy = cluster->E();
    cluster->GetPosition(posF);
    Char_t ttype   = cluster->GetType(); 
    
    
    AliAODCaloCluster *caloCluster = new(caloClusters[jClusters++]) 
      AliAODCaloCluster(id,
			0,
			0x0,
			energy,
			posF,
			NULL,
			ttype);
    
    caloCluster->SetCaloCluster(cluster->GetDistToBadChannel(),
				cluster->GetDispersion(),
				cluster->GetM20(), cluster->GetM02(),
				cluster->GetEmcCpvDistance(),  
				cluster->GetNExMax(),cluster->GetTOF()) ;
    
    caloCluster->SetPIDFromESD(cluster->PID());
    caloCluster->SetNCells(cluster->GetNCells());
    caloCluster->SetCellsAbsId(cluster->GetCellsAbsId());
    caloCluster->SetCellsAmplitudeFraction(cluster->GetCellsAmplitudeFraction());
    
  } 
  caloClusters.Expand(jClusters); // resize TObjArray to 'remove' slots for pseudo clusters	 
  // end of loop on calo clusters
  
  // fill EMCAL cell info
  if (fCalorimeter.Contains("EMCAL") && aod->GetEMCALCells()) { // protection against missing ESD information
    AliAODCaloCells &aodinEMcells = *(aod->GetEMCALCells());
    Int_t nEMcell = aodinEMcells.GetNumberOfCells() ;
    
    AliAODCaloCells &aodEMcells = *(AODEvent()->GetEMCALCells());
    aodEMcells.CreateContainer(nEMcell);
    aodEMcells.SetType(AliAODCaloCells::kEMCAL);
    for (Int_t iCell = 0; iCell < nEMcell; iCell++) {      
      aodEMcells.SetCell(iCell,aodinEMcells.GetCellNumber(iCell),aodinEMcells.GetAmplitude(iCell));
    }
    aodEMcells.Sort();
  }
  
  // fill PHOS cell info
  if (fCalorimeter.Contains("PHOS") && aod->GetPHOSCells()) { // protection against missing ESD information
    AliAODCaloCells &aodinPHcells = *(aod->GetPHOSCells());
    Int_t nPHcell = aodinPHcells.GetNumberOfCells() ;
    
    AliAODCaloCells &aodPHcells = *(AODEvent()->GetPHOSCells());
    aodPHcells.CreateContainer(nPHcell);
    aodPHcells.SetType(AliAODCaloCells::kPHOS);
    for (Int_t iCell = 0; iCell < nPHcell; iCell++) {      
      aodPHcells.SetCell(iCell,aodinPHcells.GetCellNumber(iCell),aodinPHcells.GetAmplitude(iCell));
    }
    aodPHcells.Sort();
  }
  
  
  return;
}

//__________________________________________________
void AliAnalysisTaskCaloFilter::CreateAODFromESD()
{

  // Copy ESD header, vertex, CaloClusters and CaloCells to output AOD
  
  AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent());
  
  // set arrays and pointers
  Float_t posF[3];
  Double_t pos[3];
  
  Double_t covVtx[6];
  
  for (Int_t i = 0; i < 6; i++)  covVtx[i] = 0.;
  
  // Update the header
  
  AliAODHeader* header = AODEvent()->GetHeader();
  header->SetRunNumber(esd->GetRunNumber());
  header->SetBunchCrossNumber(esd->GetBunchCrossNumber());
  header->SetOrbitNumber(esd->GetOrbitNumber());
  header->SetPeriodNumber(esd->GetPeriodNumber());
  header->SetEventType(esd->GetEventType());
  header->SetMuonMagFieldScale(-999.); // FIXME
  header->SetCentrality(-999.);        // FIXME
  
  
  header->SetTriggerMask(esd->GetTriggerMask()); 
  header->SetTriggerCluster(esd->GetTriggerCluster());
  header->SetMagneticField(esd->GetMagneticField());
  header->SetZDCN1Energy(esd->GetZDCN1Energy());
  header->SetZDCP1Energy(esd->GetZDCP1Energy());
  header->SetZDCN2Energy(esd->GetZDCN2Energy());
  header->SetZDCP2Energy(esd->GetZDCP2Energy());
  header->SetZDCEMEnergy(esd->GetZDCEMEnergy(0),esd->GetZDCEMEnergy(1));
  Float_t diamxy[2]={esd->GetDiamondX(),esd->GetDiamondY()};
  Float_t diamcov[3]; esd->GetDiamondCovXY(diamcov);
  header->SetDiamond(diamxy,diamcov);
  //
  //
  Int_t nVertices = 1 ;/* = prim. vtx*/;
  Int_t nCaloClus = esd->GetNumberOfCaloClusters();
  
  AODEvent()->ResetStd(0, nVertices, 0, 0, 0, nCaloClus, 0, 0);
  
  // Access to the AOD container of vertices
  TClonesArray &vertices = *(AODEvent()->GetVertices());
  Int_t jVertices=0;
  
  // Add primary vertex. The primary tracks will be defined
  // after the loops on the composite objects (V0, cascades, kinks)
  const AliESDVertex *vtx = esd->GetPrimaryVertex();
  
  vtx->GetXYZ(pos); // position
  vtx->GetCovMatrix(covVtx); //covariance matrix
  
  AliAODVertex * primary = new(vertices[jVertices++])
    AliAODVertex(pos, covVtx, vtx->GetChi2toNDF(), NULL, -1, AliAODVertex::kPrimary);
  primary->SetName(vtx->GetName());
  primary->SetTitle(vtx->GetTitle());
  
  // Access to the AOD container of clusters
  TClonesArray &caloClusters = *(AODEvent()->GetCaloClusters());
  Int_t jClusters=0;
  
  for (Int_t iClust=0; iClust<nCaloClus; ++iClust) {
    
    AliESDCaloCluster * cluster = esd->GetCaloCluster(iClust);
    
    //Check which calorimeter information we want to keep.
    if     (fCalorimeter.Contains("PHOS")  && !fCalorimeter.Contains("EMCAL") && cluster->IsEMCAL()) continue ;
    else if(fCalorimeter.Contains("EMCAL") && !fCalorimeter.Contains("PHOS")  && cluster->IsPHOS())  continue ;
    
    Int_t id       = cluster->GetID();
    Float_t energy = cluster->E();
    cluster->GetPosition(posF);
    Char_t ttype   = AliAODCluster::kUndef; 
    
    if (cluster->GetClusterType() == AliESDCaloCluster::kPHOSCluster) {
      ttype=AliAODCluster::kPHOSNeutral;
    } 
    else if (cluster->GetClusterType() == AliESDCaloCluster::kEMCALClusterv1) {
      ttype = AliAODCluster::kEMCALClusterv1;
    }
    
    
    AliAODCaloCluster *caloCluster = new(caloClusters[jClusters++]) 
      AliAODCaloCluster(id,
			0,
			0x0,
			energy,
			posF,
			NULL,
			ttype);
    
    caloCluster->SetCaloCluster(cluster->GetDistanceToBadChannel(),
				cluster->GetClusterDisp(),
				cluster->GetM20(), cluster->GetM02(),
				cluster->GetEmcCpvDistance(),  
				cluster->GetNExMax(),cluster->GetTOF()) ;
    
    caloCluster->SetPIDFromESD(cluster->GetPid());
    caloCluster->SetNCells(cluster->GetNCells());
    caloCluster->SetCellsAbsId(cluster->GetCellsAbsId());
    caloCluster->SetCellsAmplitudeFraction(cluster->GetCellsAmplitudeFraction());
    
  } 
  caloClusters.Expand(jClusters); // resize TObjArray to 'remove' slots for pseudo clusters	 
  // end of loop on calo clusters
  
  // fill EMCAL cell info
  if (fCalorimeter.Contains("EMCAL") && esd->GetEMCALCells()) { // protection against missing ESD information
    AliESDCaloCells &esdEMcells = *(esd->GetEMCALCells());
    Int_t nEMcell = esdEMcells.GetNumberOfCells() ;
    
    AliAODCaloCells &aodEMcells = *(AODEvent()->GetEMCALCells());
    aodEMcells.CreateContainer(nEMcell);
    aodEMcells.SetType(AliAODCaloCells::kEMCAL);
    for (Int_t iCell = 0; iCell < nEMcell; iCell++) {      
      aodEMcells.SetCell(iCell,esdEMcells.GetCellNumber(iCell),esdEMcells.GetAmplitude(iCell));
    }
    aodEMcells.Sort();
  }
  
  // fill PHOS cell info
  if (fCalorimeter.Contains("PHOS") && esd->GetPHOSCells()) { // protection against missing ESD information
    AliESDCaloCells &esdPHcells = *(esd->GetPHOSCells());
    Int_t nPHcell = esdPHcells.GetNumberOfCells() ;
    
    AliAODCaloCells &aodPHcells = *(AODEvent()->GetPHOSCells());
    aodPHcells.CreateContainer(nPHcell);
    aodPHcells.SetType(AliAODCaloCells::kPHOS);
    for (Int_t iCell = 0; iCell < nPHcell; iCell++) {      
      aodPHcells.SetCell(iCell,esdPHcells.GetCellNumber(iCell),esdPHcells.GetAmplitude(iCell));
    }
    aodPHcells.Sort();
  }
  
  
  return;
}

//__________________________________________________
void AliAnalysisTaskCaloFilter::Init()
{
  // Initialization
  if (fDebug > 1) AliInfo("Init() \n");
  
}

//__________________________________________________
void AliAnalysisTaskCaloFilter::UserCreateOutputObjects()
{
  // Create the output container
}

//__________________________________________________
void AliAnalysisTaskCaloFilter::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event
  //
  
  Long64_t ientry = Entry();
  //  if (fDebug > 0) 
  printf("CaloFilter: Analysing event # %5d\n", (Int_t) ientry);

  const char * dataevent =InputEvent()->GetName();
  if(!strcmp(dataevent,"AliAODEvent"))   CreateAODFromAOD();  
  
  else if(!strcmp(dataevent,"AliESDEvent"))  CreateAODFromESD();
  else {
    AliFatal(Form("Unknown event type %s, ABORT",dataevent));
    
  }

}

//__________________________________________________
void AliAnalysisTaskCaloFilter::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
    if (fDebug > 1) printf("AnalysisCaloFilter: Terminate() \n");
}

