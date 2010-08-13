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
#include "AliEMCALGeometry.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
//#include "AliEMCALAodCluster.h"
//#include "AliEMCALCalibData.h"

ClassImp(AliAnalysisTaskEMCALPi0CalibSelection)


//__________________________________________________
AliAnalysisTaskEMCALPi0CalibSelection::AliAnalysisTaskEMCALPi0CalibSelection(const char* name) :
  AliAnalysisTaskSE(name),fEMCALGeo(0x0),//fCalibData(0x0), 
  fEmin(0.5), fEmax(15.), fMinNCells(2), fGroupNCells(0),
  fLogWeight(4.5), fCopyAOD(kFALSE), fEMCALGeoName("EMCAL_FIRSTYEAR"),     
  fRemoveBadChannels(kFALSE),fEMCALBadChannelMap(0x0),
  fRecalibration(kFALSE),fEMCALRecalibrationFactors(), 
  fNbins(300), fMinBin(0.), fMaxBin(300.),
  fOutputContainer(0x0),fHmgg(0x0), fhNEvents(0x0),fCuts(0x0)
{
  //Named constructor which should be used.
  
  for(Int_t iMod=0; iMod < 12; iMod++) {
    for(Int_t iX=0; iX<24; iX++) {
      for(Int_t iZ=0; iZ<48; iZ++) {
	fHmpi0[iMod][iX][iZ]=0;
      }
    } 
  }
  
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());  // will contain cuts or local params

}

//__________________________________________________
AliAnalysisTaskEMCALPi0CalibSelection::~AliAnalysisTaskEMCALPi0CalibSelection()
{
  //Destructor.
  
  if(fOutputContainer){
    fOutputContainer->Delete() ; 
    delete fOutputContainer ;
  }
	
  //if(fCalibData)  delete fCalibData;
  if(fEMCALGeo)   delete fEMCALGeo;
	
	if(fEMCALBadChannelMap) { 
		fEMCALBadChannelMap->Clear();
		delete  fEMCALBadChannelMap;
	}
	
	
	if(fEMCALRecalibrationFactors) { 
		fEMCALRecalibrationFactors->Clear();
		delete  fEMCALRecalibrationFactors;
	}	
}

//_____________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelection::LocalInit()
{
	// Local Initialization
	
	// Create cuts/param objects and publish to slot
	
	char onePar[255] ;
	fCuts = new TList();

	sprintf(onePar,"Custer cuts: %2.2f < E < %2.2f GeV; min number of cells %d;", fEmin,fEmax, fMinNCells) ;
	fCuts->Add(new TObjString(onePar));
	sprintf(onePar,"Group %d cells;", fGroupNCells) ;
	fCuts->Add(new TObjString(onePar));
	sprintf(onePar,"Histograms: bins %d; energy range: %2.2f < E < %2.2f GeV;",fNbins,fMinBin,fMaxBin) ;
	fCuts->Add(new TObjString(onePar));
	sprintf(onePar,"Switchs: Remove Bad Channels? %d; Copy AODs? %d;  Recalibrate %d?",fRemoveBadChannels,fCopyAOD,fRecalibration) ;
	fCuts->Add(new TObjString(onePar));
	sprintf(onePar,"EMCAL Geometry name: < %s >",fEMCALGeoName.Data()) ;
	fCuts->Add(new TObjString(onePar));

	// Post Data
	PostData(2, fCuts);
	
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
  Int_t nCaloClus = aod->GetNumberOfCaloClusters();
  
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
  printf("nCaloClus %d\n",nCaloClus);
  
  for (Int_t iClust=0; iClust<nCaloClus; ++iClust) {
    
    AliAODCaloCluster * cluster = aod->GetCaloCluster(iClust);
    
    //Check if it is a EMCAL cluster
    if(!cluster->IsEMCAL())  continue ;
    printf("EMCAL cluster %d, ncells %d\n",iClust, cluster->GetNCells());
    if(ClusterContainsBadChannel(cluster->GetCellsAbsId(), cluster->GetNCells())) continue;	
    printf("copy\n");
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
    
    caloCluster->SetCaloCluster(cluster->GetDistanceToBadChannel(),
                                cluster->GetDispersion(),
                                cluster->GetM20(), cluster->GetM02(),
                                cluster->GetEmcCpvDistance(),  
                                cluster->GetNExMax(),cluster->GetTOF()) ;
    
    caloCluster->SetPIDFromESD(cluster->GetPID());
    caloCluster->SetNCells(cluster->GetNCells());
    caloCluster->SetCellsAbsId(cluster->GetCellsAbsId());
    
    caloCluster->SetCellsAmplitudeFraction(cluster->GetCellsAmplitudeFraction());
    
  } 
  
  caloClusters.Expand(jClusters); // resize TObjArray	 
  // end of loop on calo clusters
  
  // fill EMCAL cell info
  if (aod->GetEMCALCells()) { // protection against missing AOD information
    AliAODCaloCells &aodinEMcells = *(aod->GetEMCALCells());
    Int_t nEMcell = aodinEMcells.GetNumberOfCells() ;
    
    AliAODCaloCells &aodEMcells = *(AODEvent()->GetEMCALCells());
    aodEMcells.CreateContainer(nEMcell);
    aodEMcells.SetType(AliAODCaloCells::kEMCALCell);
    
    Double_t calibFactor = 1;
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
  printf("nCaloClus %d\n",nCaloClus);
  
  for (Int_t iClust=0; iClust<nCaloClus; ++iClust) {
    
    AliESDCaloCluster * cluster = esd->GetCaloCluster(iClust);
    
    //Check which calorimeter information we want to keep.
    if(!cluster->IsEMCAL())  continue ;
    printf("EMCAL cluster %d, ncells %d\n",iClust, cluster->GetNCells());
    
    if(ClusterContainsBadChannel(cluster->GetCellsAbsId(), cluster->GetNCells())) continue;	
    printf("copy\n");
    
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
                      AliVCluster::kEMCALClusterv1);
    
    caloCluster->SetCaloCluster(cluster->GetDistanceToBadChannel(),
                                cluster->GetDispersion(),
                                cluster->GetM20(), cluster->GetM02(),
                                cluster->GetEmcCpvDistance(),  
                                cluster->GetNExMax(),cluster->GetTOF()) ;
    
    caloCluster->SetPIDFromESD(cluster->GetPID());
    caloCluster->SetNCells(cluster->GetNCells());
    caloCluster->SetCellsAbsId(cluster->GetCellsAbsId());
    caloCluster->SetCellsAmplitudeFraction(cluster->GetCellsAmplitudeFraction());
    
  } 
  
  caloClusters.Expand(jClusters); // resize TObjArray
  // end of loop on calo clusters
  
  // fill EMCAL cell info
  
  if( esd->GetEMCALCells()) { // protection against missing ESD information
    AliESDCaloCells &esdEMcells = *(esd->GetEMCALCells());
    Int_t nEMcell = esdEMcells.GetNumberOfCells() ;
    
    AliAODCaloCells &aodEMcells = *(AODEvent()->GetEMCALCells());
    aodEMcells.CreateContainer(nEMcell);
    aodEMcells.SetType(AliAODCaloCells::kEMCALCell);  
	  
    Double_t calibFactor = 1;   
    for (Int_t iCell = 0; iCell < nEMcell; iCell++) {      
      aodEMcells.SetCell(iCell,esdEMcells.GetCellNumber(iCell),esdEMcells.GetAmplitude(iCell)*calibFactor);
    }
    aodEMcells.Sort();
	  
  }

}

//__________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelection::UserCreateOutputObjects()
{
  //Create output container, init geometry and calibration
	
  fEMCALGeo =  AliEMCALGeometry::GetInstance(fEMCALGeoName) ;	
	
  fOutputContainer = new TList();
  
  char hname[128], htitl[128];
  
  for(Int_t iMod=0; iMod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); iMod++) {
    for(Int_t iRow=0; iRow<AliEMCALGeoParams::fgkEMCALRows; iRow++) {
      for(Int_t iCol=0; iCol<AliEMCALGeoParams::fgkEMCALCols; iCol++) {
	sprintf(hname,"%d_%d_%d",iMod,iCol,iRow);
	sprintf(htitl,"Two-gamma inv. mass for super mod %d, cell(col,row)=(%d,%d)",iMod,iCol,iRow);
	fHmpi0[iMod][iCol][iRow] = new TH1F(hname,htitl,fNbins,fMinBin,fMaxBin);
	fOutputContainer->Add(fHmpi0[iMod][iCol][iRow]);
      }
    }
  }

  fHmgg = new TH1F("hmgg","2-cluster invariant mass",fNbins,fMinBin,fMaxBin);
  fOutputContainer->Add(fHmgg);
  
  fhNEvents        = new TH1I("hNEvents", "Number of analyzed events"   , 1 , 0 , 1  ) ;
  fOutputContainer->Add(fhNEvents);

//  fCalibData = new AliEMCALCalibData();
		
  PostData(1,fOutputContainer);

}

//__________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelection::UserExec(Option_t* /* option */)
{
  //Analysis per event.
  if(DebugLevel() > 1) printf("AliAnalysisTaskEMCALPi0CalibSelection <<< Event %d >>>\n",(Int_t)Entry());
  
  fhNEvents->Fill(0); //Event analyzed
  
  AliAODEvent* aod = 0x0;
  Bool_t kAOD = kFALSE;
  if(!strcmp(InputEvent()->GetName(),"AliAODEvent")) kAOD=kTRUE;
  Bool_t kESD = kFALSE;
  if(!strcmp(InputEvent()->GetName(),"AliESDEvent")) kESD=kTRUE;

  if(kAOD){
    //Input are ESDs
    aod = dynamic_cast<AliAODEvent*>(InputEvent());
    // Create new AOD with only CaloClusters and CaloCells
    if(fCopyAOD) CreateAODFromAOD();
  }
  else  if(kESD) {
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
  
  Double_t v[3];// = {aod->GetVertex(0)->GetX(),aod->GetVertex(0)->GetY(),aod->GetVertex(0)->GetZ()}; //to check!!
  aod->GetPrimaryVertex()->GetXYZ(v) ;
  //TVector3 vtx(v); 
  
  //if(DebugLevel() > 1) printf("AliAnalysisTaskEMCALPi0CalibSelection Vertex: (%.3f,%.3f,%.3f)\n",vtx.X(),vtx.Y(),vtx.Z());
  if(DebugLevel() > 1) printf("AliAnalysisTaskEMCALPi0CalibSelection Vertex: (%.3f,%.3f,%.3f)\n",v[0],v[1],v[2]);
  
  Int_t runNum = aod->GetRunNumber();
  if(DebugLevel() > 1) printf("Run number: %d\n",runNum);
  
  //FIXME Not neede the matrices for the moment MEFIX
  //Get the matrix with geometry information
  //Still not implemented in AOD, just a workaround to be able to work at least with ESDs	
//  if(!strcmp(InputEvent()->GetName(),"AliAODEvent")) {
//    if(DebugLevel() > 1) 
//      printf("AliAnalysisTaskEMCALPi0CalibSelection Use ideal geometry, values geometry matrix not kept in AODs.\n");
//  }
//  else{	
//    if(DebugLevel() > 1) printf("AliAnalysisTaskEMCALPi0CalibSelection Load Misaligned matrices. \n");
//    AliESDEvent* esd = dynamic_cast<AliESDEvent*>(InputEvent()) ;
//    for(Int_t mod=0; mod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++){ 
//      if(esd->GetEMCALMatrix(mod)) fEMCALGeo->SetMisalMatrix(esd->GetEMCALMatrix(mod),mod) ;
//    }
//  }
  
  if(DebugLevel() > 1) printf("AliAnalysisTaskEMCALPi0CalibSelection Will use fLogWeight %.3f .\n",fLogWeight);
    
  Int_t iSupMod1 = -1;
  Int_t iphi1    = -1;
  Int_t ieta1    = -1;
  Int_t iSupMod2 = -1;
  Int_t iphi2    = -1;
  Int_t ieta2    = -1;
	
  TLorentzVector p1;
  TLorentzVector p2;
  TLorentzVector p12;
  
  TRefArray * caloClustersArr  = new TRefArray();
  aod->GetEMCALClusters(caloClustersArr);
  
  const Int_t kNumberOfEMCALClusters   = caloClustersArr->GetEntries() ;
  if(DebugLevel() > 1) printf("AliAnalysisTaskEMCALPi0CalibSelection - N CaloClusters: %d \n", kNumberOfEMCALClusters);
  
  // EMCAL cells
  AliAODCaloCells *emCells = aod->GetEMCALCells();
  
  // loop over EMCAL clusters
  for(Int_t iClu=0; iClu<kNumberOfEMCALClusters; iClu++) {
    
    AliAODCaloCluster *c1 = (AliAODCaloCluster *) caloClustersArr->At(iClu);
    if(!fCopyAOD && ClusterContainsBadChannel(c1->GetCellsAbsId(), c1->GetNCells())) continue;	
    
    Float_t e1i = c1->E();   // cluster energy before correction   
    if(e1i < fEmin) continue;
    else if(e1i > fEmax) continue;
    else if (c1->GetNCells() < fMinNCells) continue; 
    
    if(DebugLevel() > 2)
      { 
	printf("Std  : i %d, E %f, dispersion %f, m02 %f, m20 %f\n",iClu,e1i, c1->GetDispersion(),c1->GetM02(),c1->GetM20());
	Float_t pos[]={0,0,0};
	c1->GetPosition(pos);
	printf("Std  : i %d, x %f, y %f, z %f\n",iClu, pos[0], pos[1], pos[2]);
      }
    
    //AliEMCALAodCluster clu1(*c1);
    //clu1.Recalibrate(fCalibData, emCells, fEMCALGeoName);
    //clu1.EvalEnergy();
    //clu1.EvalAll(fLogWeight, fEMCALGeoName);
    if(IsRecalibrationOn())	{
      Float_t energy = RecalibrateClusterEnergy(c1, emCells);
      //clu1.SetE(energy);
      c1->SetE(energy);
    }
    
    //Float_t e1ii = clu1.E(); // cluster energy after correction
    Float_t e1ii = c1->E(); // cluster energy after correction
    
    if(DebugLevel() > 2)
      { 
	//printf("Recal: i %d, E %f, dispersion %f, m02 %f, m20 %f\n",iClu,e1ii, clu1.GetDispersion(),clu1.GetM02(),clu1.GetM20()); 
	printf("Recal: i %d, E %f, dispersion %f, m02 %f, m20 %f\n",iClu,e1ii, c1->GetDispersion(),c1->GetM02(),c1->GetM20());    
	Float_t pos2[]={0,0,0};
	//clu1.GetPosition(pos2);
	c1->GetPosition(pos2);
	printf("Recal: i %d, x %f, y %f, z %f\n",iClu, pos2[0], pos2[1], pos2[2]);
      }
    
    //clu1.GetMomentum(p1,v);
    c1->GetMomentum(p1,v);
    
    //Get tower with maximum energy and fill in the end the pi0 histogram for this cell.
    //MaxEnergyCellPos(emCells,&clu1,iSupMod1,ieta1,iphi1);
    MaxEnergyCellPos(emCells,c1,iSupMod1,ieta1,iphi1);
    
    for (Int_t jClu=iClu; jClu<kNumberOfEMCALClusters; jClu++) {
      AliAODCaloCluster *c2 = (AliAODCaloCluster *) caloClustersArr->At(jClu);
      if(c2->IsEqual(c1)) continue;
      if(!fCopyAOD && ClusterContainsBadChannel(c2->GetCellsAbsId(), c2->GetNCells())) continue;	
      
      Float_t e2i = c2->E();
      if(e2i < fEmin) continue;
      else if (e2i > fEmax) continue;
      else if (c2->GetNCells() < fMinNCells) continue; 
      
      //AliEMCALAodCluster clu2(*c2);
      //clu2.Recalibrate(fCalibData, emCells,fEMCALGeoName);
      //clu2.EvalEnergy();
      //clu2.EvalAll(fLogWeight,fEMCALGeoName);
      if(IsRecalibrationOn())	{
	Float_t energy = RecalibrateClusterEnergy(c2, emCells);
	//clu2.SetE(energy);
	c2->SetE(energy);
      }
      
      Float_t e2ii = c2->E();
      
      //clu2.GetMomentum(p2,v);
      c2->GetMomentum(p2,v);
      
      //Get tower with maximum energy and fill in the end the pi0 histogram for this cell.
      //MaxEnergyCellPos(emCells,&clu2,iSupMod2,ieta2,iphi2);
      MaxEnergyCellPos(emCells,c2,iSupMod2,ieta2,iphi2);
      
      p12 = p1+p2;
      Float_t invmass = p12.M()*1000; 
      
      if(invmass < fMaxBin && invmass > fMinBin){
	fHmgg->Fill(invmass); 
	
	if (fGroupNCells == 0){
	  fHmpi0[iSupMod1][ieta1][iphi1]->Fill(invmass);
	  fHmpi0[iSupMod2][ieta2][iphi2]->Fill(invmass);
	}	
	else  {
	  //printf("Regroup N %d, eta1 %d, phi1 %d, eta2 %d, phi2 %d \n",fGroupNCells, ieta1, iphi1, ieta2, iphi2);
	  for (Int_t i = -fGroupNCells; i < fGroupNCells+1; i++) {
	    for (Int_t j = -fGroupNCells; j < fGroupNCells+1; j++) {
	      //printf("\t i %d, j %d\n",i,j);
	      if((ieta1+i >= 0) && (iphi1+j >= 0) && (ieta1+i < 48) && (iphi1+j < 24)){
		//printf("\t \t eta1+i %d, phi1+j %d\n", ieta1+i, iphi1+j);
		fHmpi0[iSupMod1][ieta1+i][iphi1+j]->Fill(invmass);
	      }
	      if((ieta2+i >= 0) && (iphi2+j >= 0) && (ieta2+i < 48) && (iphi2+j < 24)){
		//printf("\t \t eta2+i %d, phi2+j %d\n", ieta2+i, iphi2+j);
		fHmpi0[iSupMod2][ieta2+i][iphi2+j]->Fill(invmass);
	      }
	    }// j loop
	  }//i loop
	}//group cells
	
	if(DebugLevel() > 1) printf("AliAnalysisTaskEMCALPi0CalibSelection Mass in (SM%d,%d,%d) and  (SM%d,%d,%d): %.3f GeV  E1_i=%f E1_ii=%f  E2_i=%f E2_ii=%f\n",
				    iSupMod1,iphi1,ieta1,iSupMod2,iphi2,ieta2,p12.M(),e1i,e1ii,e2i,e2ii);
      }
      
    }
    
  } // end of loop over EMCAL clusters
  
  delete caloClustersArr;
  
  PostData(1,fOutputContainer);
  
}

//__________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelection::MaxEnergyCellPos(AliAODCaloCells* const cells, AliAODCaloCluster* const clu, Int_t& iSupMod, Int_t& ieta, Int_t& iphi)
{
  //For a given CaloCluster calculates the absId of the cell 
  //with maximum energy deposit.
  
  Double_t eMax       = -1.;
  Double_t eCell      = -1.;
  Float_t  fraction   = 1.;
  Int_t    cellAbsId  = -1;
  Float_t recalFactor = 1.;
	
  Int_t maxId   = -1;
  Int_t iTower  = -1;
  Int_t iIphi   = -1;
  Int_t iIeta   = -1;
	
  for (Int_t iDig=0; iDig< clu->GetNCells(); iDig++) {
    cellAbsId = clu->GetCellAbsId(iDig);
    fraction  = clu->GetCellAmplitudeFraction(iDig);
    if(fraction < 1e-4) fraction = 1.; // in case unfolding is off
    if(IsRecalibrationOn()) {
      Int_t imodrc   = -1, iphirc  = -1, ietarc  =-1;
      Int_t iTowerrc = -1, iIphirc = -1, iIetarc =-1;
      fEMCALGeo->GetCellIndex(cellAbsId,imodrc,iTowerrc,iIphirc,iIetarc); 
      fEMCALGeo->GetCellPhiEtaIndexInSModule(imodrc,iTowerrc,iIphirc, iIetarc,iphirc,ietarc);			
      recalFactor = GetEMCALChannelRecalibrationFactor(imodrc,ietarc,iphirc);
    }
    eCell  = cells->GetCellAmplitude(cellAbsId)*fraction*recalFactor;
    //printf("Max cell? cell %d, amplitude org %f, fraction %f, recalibration %f, amplitude new %f \n",cellAbsId, cells->GetCellAmplitude(cellAbsId), fraction, recalFactor, eCell) ;
    
    if(eCell > eMax)  { 
      eMax  = eCell; 
      maxId = cellAbsId;
      //printf("\t new max: cell %d, e %f, ecell %f\n",maxId, eMax,eCell);
    }
  }
  
  //Get from the absid the supermodule, tower and eta/phi numbers
  fEMCALGeo->GetCellIndex(maxId,iSupMod,iTower,iIphi,iIeta); 
  //Gives SuperModule and Tower numbers
  fEMCALGeo->GetCellPhiEtaIndexInSModule(iSupMod,iTower,
					 iIphi, iIeta,iphi,ieta);    
  
  //printf("\t Max : cell %d, iSupMod %d, ieta %d, iphi %d \n",maxId,iSupMod, ieta,iphi);
  
}

//__________________________________________________
//void AliAnalysisTaskEMCALPi0CalibSelection::SetCalibCorrections(AliEMCALCalibData* const cdata)
//{
//  //Set new correction factors (~1) to calibration coefficients, delete previous.
//
//   if(fCalibData) delete fCalibData;
//   fCalibData = cdata;
//	
//}

//________________________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelection::InitEMCALBadChannelStatusMap(){
	//Init EMCAL bad channels map
	if(DebugLevel() > 0 )printf("AliAnalysisTaskEMCALPi0CalibSelection::InitEMCALBadChannelStatusMap()\n");
	//In order to avoid rewriting the same histograms
	Bool_t oldStatus = TH1::AddDirectoryStatus();
	TH1::AddDirectory(kFALSE);
	
	fEMCALBadChannelMap = new TObjArray(12);
	//TH2F * hTemp = new  TH2I("EMCALBadChannelMap","EMCAL SuperModule bad channel map", 48, 0, 48, 24, 0, 24);
	for (int i = 0; i < 12; i++) {
		fEMCALBadChannelMap->Add(new TH2I(Form("EMCALBadChannelMap_Mod%d",i),Form("EMCALBadChannelMap_Mod%d",i), 48, 0, 48, 24, 0, 24));
		//fEMCALBadChannelMap->Add((TH2I*)hTemp->Clone(Form("EMCALBadChannelMap_Mod%d",i)));
	}
	
	//delete hTemp;
	
	fEMCALBadChannelMap->SetOwner(kTRUE);
	fEMCALBadChannelMap->Compress();
	
	//In order to avoid rewriting the same histograms
	TH1::AddDirectory(oldStatus);		
}


//_________________________________________________________________________________________________________
Bool_t AliAnalysisTaskEMCALPi0CalibSelection::ClusterContainsBadChannel(UShort_t* cellList, Int_t nCells){
	// Check that in the cluster cells, there is no bad channel of those stored 
	// in fEMCALBadChannelMap or fPHOSBadChannelMap
	
	if(!fRemoveBadChannels)  return kFALSE;
	if(!fEMCALBadChannelMap) return kFALSE;
	
	Int_t icol = -1;
	Int_t irow = -1;
	Int_t imod = -1;
	for(Int_t iCell = 0; iCell<nCells; iCell++){
		
		//Get the column and row
			Int_t iTower = -1, iIphi = -1, iIeta = -1; 
			fEMCALGeo->GetCellIndex(cellList[iCell],imod,iTower,iIphi,iIeta); 
			if(fEMCALBadChannelMap->GetEntries() <= imod) continue;
			fEMCALGeo->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,irow,icol);			
			if(GetEMCALChannelStatus(imod, icol, irow))return kTRUE;
		
	}// cell cluster loop
	
	return kFALSE;
	
}


//________________________________________________________________
void AliAnalysisTaskEMCALPi0CalibSelection::InitEMCALRecalibrationFactors(){
	//Init EMCAL recalibration factors
	if(DebugLevel() > 0 )printf("AliAnalysisTaskEMCALPi0CalibSelection::InitEMCALRecalibrationFactors()\n");
	//In order to avoid rewriting the same histograms
	Bool_t oldStatus = TH1::AddDirectoryStatus();
	TH1::AddDirectory(kFALSE);
	
	fEMCALRecalibrationFactors = new TObjArray(12);
	for (int i = 0; i < 12; i++) fEMCALRecalibrationFactors->Add(new TH2F(Form("EMCALRecalFactors_SM%d",i),Form("EMCALRecalFactors_SM%d",i),  48, 0, 48, 24, 0, 24));
	//Init the histograms with 1
	for (Int_t sm = 0; sm < 12; sm++) {
		for (Int_t i = 0; i < 48; i++) {
			for (Int_t j = 0; j < 24; j++) {
				SetEMCALChannelRecalibrationFactor(sm,i,j,1.);
			}
		}
	}
	fEMCALRecalibrationFactors->SetOwner(kTRUE);
	fEMCALRecalibrationFactors->Compress();
	
	//In order to avoid rewriting the same histograms
	TH1::AddDirectory(oldStatus);		
}

//________________________________________________________________
Float_t AliAnalysisTaskEMCALPi0CalibSelection::RecalibrateClusterEnergy(AliAODCaloCluster * cluster, AliAODCaloCells * cells){
	// Recalibrate the cluster energy, considering the recalibration map and the energy of the cells that compose the cluster.
	// AOD case
	
	if(!cells) {
		printf("AliAnalysisTaskEMCALPi0CalibSelection::RecalibrateClusterEnergy(AOD) - Cells pointer does not exist, stop!");
		abort();
	}
	
	//Get the cluster number of cells and list of absId, check what kind of cluster do we have.
	UShort_t * index    = cluster->GetCellsAbsId() ;
	Double_t * fraction = cluster->GetCellsAmplitudeFraction() ;
	Int_t ncells = cluster->GetNCells();
	
	//Initialize some used variables
	Float_t energy = 0;
	Int_t absId    = -1;
    Int_t icol = -1, irow = -1, imod=1;
	Float_t factor = 1, frac = 0;
	
	//Loop on the cells, get the cell amplitude and recalibration factor, multiply and and to the new energy
	for(Int_t icell = 0; icell < ncells; icell++){
		absId = index[icell];
		frac =  fraction[icell];
		if(frac < 1e-5) frac = 1; //in case of EMCAL, this is set as 0 since unfolding is off
		Int_t iTower = -1, iIphi = -1, iIeta = -1; 
		fEMCALGeo->GetCellIndex(absId,imod,iTower,iIphi,iIeta); 
		if(fEMCALRecalibrationFactors->GetEntries() <= imod) continue;
		fEMCALGeo->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,irow,icol);			
		factor = GetEMCALChannelRecalibrationFactor(imod,icol,irow);
		if(DebugLevel()>2)
			printf("AliAnalysisTaskEMCALPi0CalibSelection::RecalibrateClusterEnergy - recalibrate cell: module %d, col %d, row %d, cell fraction %f,recalibration factor %f, cell energy %f\n",
				   imod,icol,irow,frac,factor,cells->GetCellAmplitude(absId));
		
		energy += cells->GetCellAmplitude(absId)*factor*frac;
	}
	
	if(DebugLevel()>1)
		printf("AliAnalysisTaskEMCALPi0CalibSelection::RecalibrateClusterEnergy - Energy before %f, after %f\n",cluster->E(),energy);
	
	return energy;
	
}

