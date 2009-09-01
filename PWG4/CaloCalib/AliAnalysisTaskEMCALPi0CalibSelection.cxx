/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Boris Polishchuk                                               *
 * Adapted to AOD reading by Gustavo Conesa  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//---------------------------------------------------------------------------// 
//                                                                           //
// Fill histograms (one per cell) with two-cluster invariant mass            //
// using calibration coefficients of the previous iteration.                 //
// Histogram for a given cell is filled if the most energy of one cluster    //
// is deposited in this cell and the other cluster could be anywherein EMCAL.//
//                                                                           //
//---------------------------------------------------------------------------//

//#include <cstdlib>
//#include <Riostream.h>
// Root 
#include "TLorentzVector.h"
//#include "TVector3.h"
#include "TRefArray.h"
#include "TList.h"
#include "TH1F.h"

// AliRoot
#include "AliAnalysisTaskEMCALPi0CalibSelection.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliEMCALGeometry.h"
#include "AliAODCaloCluster.h"
#include "AliAODCaloCells.h"
#include "AliEMCALAodCluster.h"
#include "AliEMCALCalibData.h"

ClassImp(AliAnalysisTaskEMCALPi0CalibSelection)

//__________________________________________________
AliAnalysisTaskEMCALPi0CalibSelection::AliAnalysisTaskEMCALPi0CalibSelection() :
AliAnalysisTaskSE(),fEMCALGeo(0x0),fCalibData(0x0), fEmin(0.), fLogWeight(4.5), fCopyAOD(kFALSE), 
fEMCALGeoName("EMCAL_COMPLETE"), fOldData(kFALSE),fOutputContainer(0x0),fHmgg(0x0)
{
  //Default constructor.

  for(Int_t iMod=0; iMod < 12; iMod++) {
    for(Int_t iX=0; iX<24; iX++) {
      for(Int_t iZ=0; iZ<48; iZ++) {
	fHmpi0[iMod][iX][iZ]=0;
      }
    } 
  }
		
}

//__________________________________________________
AliAnalysisTaskEMCALPi0CalibSelection::AliAnalysisTaskEMCALPi0CalibSelection(const char* name) :
  AliAnalysisTaskSE(name),fEMCALGeo(0x0),fCalibData(0x0), fEmin(0.), fLogWeight(4.5), fCopyAOD(kFALSE), 
  fEMCALGeoName("EMCAL_COMPLETE"), fOldData(kFALSE), fOutputContainer(0x0),fHmgg(0x0)
{
  //Named constructor which should be used.
  
  DefineOutput(1,TList::Class());
  
  for(Int_t iMod=0; iMod < 12; iMod++) {
    for(Int_t iX=0; iX<24; iX++) {
      for(Int_t iZ=0; iZ<48; iZ++) {
	fHmpi0[iMod][iX][iZ]=0;
      }
    } 
  }

}

//__________________________________________________
AliAnalysisTaskEMCALPi0CalibSelection::~AliAnalysisTaskEMCALPi0CalibSelection()
{
  //Destructor.
  
  if(fOutputContainer){
    fOutputContainer->Delete() ; 
    delete fOutputContainer ;
  }
	
  if(fCalibData)  delete fCalibData;
  if(fEMCALGeo)   delete fEMCALGeo;
	
}


//__________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelection::CreateAODFromAOD()
{
  // Copy AOD header, vertex, CaloClusters and CaloCells to output AOD
  
  AliAODEvent* aod = dynamic_cast<AliAODEvent*>(InputEvent());
  
  // set arrays and pointers
  Float_t posF[3];
  Double_t pos[3];
  
  Double_t covVtx[6];
  
  for (Int_t i = 0; i < 6; i++)  covVtx[i] = 0.;
  
  // Update the header
  AliAODHeader*headerin = aod->GetHeader();
  AliAODHeader* header = AODEvent()->GetHeader();
  header->SetRunNumber(headerin->GetRunNumber());
  header->SetBunchCrossNumber(headerin->GetBunchCrossNumber());
  header->SetOrbitNumber(headerin->GetOrbitNumber());
  header->SetPeriodNumber(headerin->GetPeriodNumber());
  header->SetEventType(headerin->GetEventType());
  header->SetMuonMagFieldScale(headerin->GetMuonMagFieldScale());
  header->SetCentrality(headerin->GetCentrality()); 
  
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
    
    //Check if it is a EMCAL cluster
    if(!cluster->IsEMCALCluster())  continue ;
    
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
	if(!fOldData){
		caloCluster->SetCellsAmplitudeFraction(cluster->GetCellsAmplitudeFraction());
	}
	else{
		Double_t *fraction = cluster->GetCellsAmplitudeFraction();
		for(Int_t i = 0; i < cluster->GetNCells() ; i++) fraction[i] = 1;
		caloCluster->SetCellsAmplitudeFraction(fraction);
	}
    
  } 

  caloClusters.Expand(jClusters); // resize TObjArray to 'remove' slots for pseudo clusters	 
  // end of loop on calo clusters
  
  // fill EMCAL cell info
  if (aod->GetEMCALCells()) { // protection against missing AOD information
	    AliAODCaloCells &aodinEMcells = *(aod->GetEMCALCells());
    Int_t nEMcell = aodinEMcells.GetNumberOfCells() ;
    
    AliAODCaloCells &aodEMcells = *(AODEvent()->GetEMCALCells());
    aodEMcells.CreateContainer(nEMcell);
    aodEMcells.SetType(AliAODCaloCells::kEMCAL);
	
	Double_t calibFactor = 1;
	if(fOldData) calibFactor = 0.0153;

    for (Int_t iCell = 0; iCell < nEMcell; iCell++) {      
      aodEMcells.SetCell(iCell,aodinEMcells.GetCellNumber(iCell),aodinEMcells.GetAmplitude(iCell)*calibFactor);
    }
    aodEMcells.Sort();
	  
  }
  
}

//__________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelection::CreateAODFromESD()
{
  
  // Copy Header, Vertex, CaloClusters and CaloCells from ESDs to AODs
  
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
    if(!cluster->IsEMCAL())  continue ;
    
    Int_t id       = cluster->GetID();
    Float_t energy = cluster->E();
    cluster->GetPosition(posF);
    
    AliAODCaloCluster *caloCluster = new(caloClusters[jClusters++]) 
      AliAODCaloCluster(id,
			0,
			0x0,
			energy,
			posF,
			NULL,
			AliAODCluster::kEMCALClusterv1);
    
    caloCluster->SetCaloCluster(cluster->GetDistanceToBadChannel(),
				cluster->GetClusterDisp(),
				cluster->GetM20(), cluster->GetM02(),
				cluster->GetEmcCpvDistance(),  
				cluster->GetNExMax(),cluster->GetTOF()) ;
    
    caloCluster->SetPIDFromESD(cluster->GetPid());
    caloCluster->SetNCells(cluster->GetNCells());
    caloCluster->SetCellsAbsId(cluster->GetCellsAbsId());
	if(!fOldData){
		  caloCluster->SetCellsAmplitudeFraction(cluster->GetCellsAmplitudeFraction());
	}
	else{
		  Double_t *fraction = cluster->GetCellsAmplitudeFraction();
		  for(Int_t i = 0; i < cluster->GetNCells() ; i++) fraction[i] = 1;
		  caloCluster->SetCellsAmplitudeFraction(fraction);
	  }    
  } 
  caloClusters.Expand(jClusters); // resize TObjArray to 'remove' slots for pseudo clusters	 
  // end of loop on calo clusters
  
  // fill EMCAL cell info

  if( esd->GetEMCALCells()) { // protection against missing ESD information
    AliESDCaloCells &esdEMcells = *(esd->GetEMCALCells());
	Int_t nEMcell = esdEMcells.GetNumberOfCells() ;
    
    AliAODCaloCells &aodEMcells = *(AODEvent()->GetEMCALCells());
    aodEMcells.CreateContainer(nEMcell);
    aodEMcells.SetType(AliAODCaloCells::kEMCAL);  
	  
	Double_t calibFactor = 1;
	if(fOldData) calibFactor = 0.0153;
   
	for (Int_t iCell = 0; iCell < nEMcell; iCell++) {      
      aodEMcells.SetCell(iCell,esdEMcells.GetCellNumber(iCell),esdEMcells.GetAmplitude(iCell)*calibFactor);
    }
    aodEMcells.Sort();
	  
  }

}

//__________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelection::UserCreateOutputObjects()
{
  //Create output container
  fOutputContainer = new TList();
  
  char hname[128], htitl[128];
  
	for(Int_t iMod=0; iMod < AliEMCALGeoParams::fgkEMCALModules; iMod++) {
    for(Int_t iX=0; iX<AliEMCALGeoParams::fgkEMCALCols; iX++) {
      for(Int_t iZ=0; iZ<AliEMCALGeoParams::fgkEMCALRows; iZ++) {
	sprintf(hname,"%d_%d_%d",iMod,iX,iZ);
	sprintf(htitl,"Two-gamma inv. mass for mod %d, cell (%d,%d)",iMod,iX,iZ);
	fHmpi0[iMod][iX][iZ] = new TH1F(hname,htitl,100,0.,300.);
	fOutputContainer->Add(fHmpi0[iMod][iX][iZ]);
      }
    }
  }

  fHmgg = new TH1F("hmgg","2-cluster invariant mass",100,0.,300.);
  fOutputContainer->Add(fHmgg);
	
  fCalibData = new AliEMCALCalibData();
	
  if(!fCalibData)
    AliFatal("Calibration parameters not found in CDB!");
  
  printf("Get Geometry : %s\n", fEMCALGeoName.Data());
  fEMCALGeo =  AliEMCALGeometry::GetInstance(fEMCALGeoName) ;	

}

//__________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelection::UserExec(Option_t* /* option */)
{
  //Analysis per event.
  if(DebugLevel() > 1) printf("AliAnalysisTaskEMCALPi0CalibSelection <<< Event %d >>>\n",(Int_t)Entry());
  
  AliAODEvent* aod = 0x0;
  if(!strcmp(InputEvent()->GetName(),"AliAODEvent")) {
    //Input are ESDs
    aod = dynamic_cast<AliAODEvent*>(InputEvent());
    // Create new AOD with only CaloClusters and CaloCells
    if(fCopyAOD) CreateAODFromAOD();
  }
  else  if(!strcmp(InputEvent()->GetName(),"AliESDEvent")) {
    //Input are ESDs
    aod = AODEvent();
    // Create AOD with CaloClusters and use it as input.
    // If filtering task is already executed, this is not needed.
    if(fCopyAOD) CreateAODFromESD();
  }
  else {
    printf("AliAnalysisTaskEMCALPi0CalibSelection: Unknown event type, STOP!\n");
    abort();
  }	
  
  Double_t v[] = {aod->GetVertex(0)->GetX(),aod->GetVertex(0)->GetY(),aod->GetVertex(0)->GetZ()}; //to check!!
  //aod->GetVertex()->GetXYZ(v) ;
  //TVector3 vtx(v); 
  
  //if(DebugLevel() > 1) printf("AliAnalysisTaskEMCALPi0CalibSelection Vertex: (%.3f,%.3f,%.3f)\n",vtx.X(),vtx.Y(),vtx.Z());
  if(DebugLevel() > 1) printf("AliAnalysisTaskEMCALPi0CalibSelection Vertex: (%.3f,%.3f,%.3f)\n",v[0],v[1],v[2]);
  
  Int_t runNum = aod->GetRunNumber();
  if(DebugLevel() > 1) printf("Run number: %d\n",runNum);
  
  //Get the matrix with geometry information
  //Still not implemented in AOD, just a workaround to be able to work at least with ESDs	
  if(!strcmp(InputEvent()->GetName(),"AliAODEvent")) {
    if(DebugLevel() > 1) 
      printf("AliAnalysisTaskEMCALPi0CalibSelection Use ideal geometry, values geometry matrix not kept in AODs.\n");
  }
  else{	
    if(DebugLevel() > 1) printf("AliAnalysisTaskEMCALPi0CalibSelection Load Misaligned matrices. \n");
    AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent()) ;
    for(Int_t mod=0; mod < 12; mod++){ 
      if(esd->GetEMCALMatrix(mod)) fEMCALGeo->SetMisalMatrix(esd->GetEMCALMatrix(mod),mod) ;
    }
  }
  
  if(DebugLevel() > 1) printf("AliAnalysisTaskEMCALPi0CalibSelection Will use fLogWeight %.3f .\n",fLogWeight);
  
  //AliEMCALPID pid;
  
  Int_t maxId;
  Int_t iSupMod = -1;
  Int_t iTower  = -1;
  Int_t iIphi   = -1;
  Int_t iIeta   = -1;
  Int_t iphi   = -1;
  Int_t ieta   = -1;
  
  TLorentzVector p1;
  TLorentzVector p2;
  TLorentzVector p12;
  
  
  TRefArray * caloClustersArr  = new TRefArray();
  aod->GetEMCALClusters(caloClustersArr);
  
  const Int_t kNumberOfEMCALClusters   = caloClustersArr->GetEntries() ;
  if(DebugLevel() > 1) printf("AliAnalysisTaskEMCALPi0CalibSelection CaloClusters: %d\n", kNumberOfEMCALClusters);
  
  // EMCAL cells
  AliAODCaloCells *emCells = aod->GetEMCALCells();
	
  // Check if is old data not modified for calibration and change the fraction factor and calibrate amplitudes
  // Only for aliroot older than release 17 tag 3
  if(fOldData && !fCopyAOD){
	  //Calibrate amplitudes
	  AliAODCaloCells &aodEMcells = *(AODEvent()->GetEMCALCells());
	  for (Int_t iCell = 0; iCell < aodEMcells.GetNumberOfCells(); iCell++) {      
			aodEMcells.SetCell(iCell,aodEMcells.GetCellNumber(iCell),aodEMcells.GetAmplitude(iCell)*0.0153);
	  }
	  
	  //Change fraction
	  for(Int_t iClu=0; iClu<kNumberOfEMCALClusters; iClu++) {
		  AliAODCaloCluster *cl = (AliAODCaloCluster *) caloClustersArr->At(iClu);
		  if(!cl->IsEMCALCluster()) continue; // EMCAL cluster!
		  Double_t *fraction = cl->GetCellsAmplitudeFraction();
		  for(Int_t i = 0; i < cl->GetNCells() ; i++) fraction[i] = 1;
		  cl->SetCellsAmplitudeFraction(fraction);	  
	  }
  }
	
  // loop over EMCAL clusters
  for(Int_t iClu=0; iClu<kNumberOfEMCALClusters; iClu++) {
    
    AliAODCaloCluster *c1 = (AliAODCaloCluster *) caloClustersArr->At(iClu);
    if(!c1->IsEMCALCluster()) continue; // EMCAL cluster!
		  
	Float_t e1i = c1->E();   // cluster energy before correction   
	if(e1i<fEmin) continue;

	if(DebugLevel() > 2){ 
		printf("Std  : i %d, E %f, dispersion %f, m02 %f, m20 %f\n",iClu,e1i, c1->GetDispersion(),c1->GetM02(),c1->GetM20());
		Double_t pos[]={0,0,0};
		c1->GetPosition(pos);
		printf("Std  : i %d, x %f, y %f, z %f\n",iClu, pos[0], pos[1], pos[2]);
	}
	  
    AliEMCALAodCluster clu1(*c1);
    clu1.Recalibrate(fCalibData, emCells, fEMCALGeoName);
	  printf("recalibrated, now calculate the rest\n");
    clu1.EvalAll(fLogWeight, fEMCALGeoName);
    //clu1.EnergyCorrection(&pid) ;
	  
	if(DebugLevel() > 2){ 
		printf("Recal: i %d, E %f, dispersion %f, m02 %f, m20 %f\n",iClu,clu1.E(), clu1.GetDispersion(),clu1.GetM02(),clu1.GetM20());    
		Double_t pos2[]={0,0,0};
		clu1.GetPosition(pos2);
		printf("Recal: i %d, x %f, y %f, z %f\n",iClu, pos2[0], pos2[1], pos2[2]);
	}
	  
	clu1.GetMomentum(p1,v);

	MaxEnergyCellPos(emCells,&clu1,maxId);
    
    //Get from the absid the supermodule, tower and eta/phi numbers
    fEMCALGeo->GetCellIndex(maxId,iSupMod,iTower,iIphi,iIeta); 
    //Gives SuperModule and Tower numbers
    fEMCALGeo->GetCellPhiEtaIndexInSModule(iSupMod,iTower,
					   iIphi, iIeta,iphi,ieta);    
    
    Float_t e1ii = clu1.E(); // cluster energy after correction
    
    for (Int_t jClu=iClu; jClu<kNumberOfEMCALClusters; jClu++) {
      AliAODCaloCluster *c2 = (AliAODCaloCluster *) caloClustersArr->At(jClu);
      if(!c2->IsEMCALCluster())   continue; // EMCAL cluster!
      if(c2->IsEqual(c1)) continue;
      
      Float_t e2i = c2->E();
      if(e2i<fEmin) continue;
      
      AliEMCALAodCluster clu2(*c2);
      //printf("i2 %d, E %f\n",iClu,e2i);
      clu2.Recalibrate(fCalibData, emCells,fEMCALGeoName);
      clu2.EvalAll(fLogWeight,fEMCALGeoName);
      //clu2.EnergyCorrection(&pid) ;
     // printf("i2 %d, Erc %f\n",iClu,clu2.E());

      clu2.GetMomentum(p2,v);
      Float_t e2ii = clu2.E();
      
      p12 = p1+p2;
      Float_t invmass = p12.M()*1000; 

      fHmgg->Fill(invmass); 

      //printf("iSM %d, ieta %d, iphi %d,  mass %f  \n",iSupMod, ieta, iphi, invmass);
      if(invmass < 300)
	  fHmpi0[iSupMod][ieta][iphi]->Fill(invmass);
      
      
      if(DebugLevel() > 1) printf("AliAnalysisTaskEMCALPi0CalibSelection Mass in (SM%d,%d,%d): %.3f GeV  E1_i=%f E1_ii=%f  E2_i=%f E2_ii=%f\n",
				  iSupMod,iphi,ieta,p12.M(),e1i,e1ii,e2i,e2ii);
    }
    
  } // end of loop over EMCAL clusters
  
  delete caloClustersArr;
  PostData(1,fOutputContainer);
}

//__________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelection::MaxEnergyCellPos(AliAODCaloCells* const cells, AliAODCaloCluster* const clu, Int_t& maxId)
{
  //For a given CaloCluster calculates the absId of the cell 
  //with maximum energy deposit.
  
  Double_t eMax = -111;
  
  for (Int_t iDig=0; iDig< clu->GetNCells(); iDig++) {
    Int_t cellAbsId = clu->GetCellAbsId(iDig);
    Double_t eCell  = cells->GetCellAmplitude(cellAbsId)*clu->GetCellAmplitudeFraction(iDig);
    if(eCell>eMax)  { 
      eMax = eCell; 
      maxId = cellAbsId;
    }
  }

}

//__________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelection::SetCalibCorrections(AliEMCALCalibData* const cdata)
{
  //Set new correction factors (~1) to calibration coefficients, delete previous.

   if(fCalibData) delete fCalibData;
   fCalibData = cdata;
	
}

