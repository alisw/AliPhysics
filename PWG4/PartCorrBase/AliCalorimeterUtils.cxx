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
/* $Id:  $ */

//_________________________________________________________________________
// Class utility for Calorimeter specific selection methods                ///
//
//
//
//-- Author: Gustavo Conesa (LPSC-Grenoble) 
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---
#include "TGeoManager.h"

//---- ANALYSIS system ----
#include "AliCalorimeterUtils.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliAODPWG4Particle.h"
#include "AliAODCaloCluster.h"
#include "AliESDCaloCluster.h"
#include "AliAODCaloCells.h"
#include "AliESDCaloCells.h"


ClassImp(AliCalorimeterUtils)
  
  
//____________________________________________________________________________
  AliCalorimeterUtils::AliCalorimeterUtils() : 
    TObject(), fDebug(0), 
    fEMCALGeoName("EMCAL_COMPLETE"),fPHOSGeoName("PHOSgeo"), 
    fEMCALGeo(0x0), fPHOSGeo(0x0), 
    fEMCALGeoMatrixSet(kFALSE), fPHOSGeoMatrixSet(kFALSE), 
    fRemoveBadChannels(kFALSE),
    fEMCALBadChannelMap(new TObjArray()),fPHOSBadChannelMap(new TObjArray()), 
    fNCellsFromEMCALBorder(0), fNCellsFromPHOSBorder(0), 
    fNoEMCALBorderAtEta0(kFALSE)
{
  //Ctor
  
  //Initialize parameters
  InitParameters();
}

//____________________________________________________________________________
AliCalorimeterUtils::AliCalorimeterUtils(const AliCalorimeterUtils & calo) :   
  TObject(calo), fDebug(calo.fDebug),
  fEMCALGeoName(calo.fEMCALGeoName), fPHOSGeoName(calo.fPHOSGeoName),
  fEMCALGeo(new AliEMCALGeoUtils(*calo.fEMCALGeo)), 
  fPHOSGeo(new AliPHOSGeoUtils(*calo.fPHOSGeo)),
  fEMCALGeoMatrixSet(calo.fEMCALGeoMatrixSet), 
  fPHOSGeoMatrixSet(calo.fPHOSGeoMatrixSet),
  fRemoveBadChannels(calo.fRemoveBadChannels),
  fEMCALBadChannelMap(new TObjArray(*calo.fEMCALBadChannelMap)),
  fPHOSBadChannelMap(new TObjArray(*calo.fPHOSBadChannelMap)),
  fNCellsFromEMCALBorder(calo.fNCellsFromEMCALBorder),
  fNCellsFromPHOSBorder(calo.fNCellsFromPHOSBorder), 
  fNoEMCALBorderAtEta0(calo.fNoEMCALBorderAtEta0)
{
  // cpy ctor  
}


//_________________________________
AliCalorimeterUtils::~AliCalorimeterUtils() {
  //Dtor
  
  //if(fPHOSGeo)  delete fPHOSGeo  ;
  if(fEMCALGeo) delete fEMCALGeo ;
	
  if(fEMCALBadChannelMap) { 
    fEMCALBadChannelMap->Clear();
    delete  fEMCALBadChannelMap;
  }
  if(fPHOSBadChannelMap) { 
    fPHOSBadChannelMap->Clear();
    delete  fPHOSBadChannelMap;
  }
	
}

//_______________________________________________________________
Bool_t AliCalorimeterUtils::CheckCellFiducialRegion(AliAODCaloCluster* cluster, AliAODCaloCells* cells) const {
	// Given the list of AbsId of the cluster, get the maximum cell and 
	// check if there are fNCellsFromBorder from the calorimeter border
	
	//If the distance to the border is 0 or negative just exit accept all clusters
	if(cells->GetType()==AliAODCaloCells::kEMCAL && fNCellsFromEMCALBorder <= 0 ) return kTRUE;
	if(cells->GetType()==AliAODCaloCells::kPHOS  && fNCellsFromPHOSBorder  <= 0 ) return kTRUE;

	//Find cells with maximum amplitude
	Int_t absIdMax	= -1;
	Float_t ampMax  = -1;
	for(Int_t i = 0; i < cluster->GetNCells() ; i++){
		Int_t absId = cluster->GetCellAbsId(i) ;
		Float_t amp	= cells->GetCellAmplitude(absId);
		if(amp > ampMax){
			ampMax   = amp;
			absIdMax = absId;
		}
	}
	
	if(fDebug > 1)
		printf("AliCalorimeterUtils::CheckCellFiducialRegion(AOD) - Cluster Max AbsId %d, Cell Energy %2.2f, Cluster Energy %2.2f\n", 
		   absIdMax, ampMax, cluster->E());
	
	if(absIdMax==-1) return kFALSE;
	
	//Check if the cell is close to the borders:
	Bool_t okrow = kFALSE;
	Bool_t okcol = kFALSE;

	if(cells->GetType()==AliAODCaloCells::kEMCAL){
	
		Int_t iTower = -1, iIphi = -1, iIeta = -1, iphi = -1, ieta = -1, iSM = -1; 
		fEMCALGeo->GetCellIndex(absIdMax,iSM,iTower,iIphi,iIeta); 
		fEMCALGeo->GetCellPhiEtaIndexInSModule(iSM,iTower,iIphi, iIeta,iphi,ieta);
		
		//Check rows/phi
		if(iSM < 10){
			if(iphi >= fNCellsFromEMCALBorder && iphi < 24-fNCellsFromEMCALBorder) okrow =kTRUE; 
	    }
		else{
			if(iphi >= fNCellsFromEMCALBorder && iphi < 12-fNCellsFromEMCALBorder) okrow =kTRUE; 
		}
		
		//Check collumns/eta
		if(!fNoEMCALBorderAtEta0){
			if(ieta  > fNCellsFromEMCALBorder && ieta < 48-fNCellsFromEMCALBorder) okcol =kTRUE; 
		}
		else{
			if(iSM%2==0){
				if(ieta >= fNCellsFromEMCALBorder)     okcol = kTRUE;	
			}
			else {
				if(ieta <  48-fNCellsFromEMCALBorder)  okcol = kTRUE;	
			}
		}//eta 0 not checked
		if(fDebug > 1)
		{
			printf("AliCalorimeterUtils::CheckCellFiducialRegion(AOD) - EMCAL Cluster in %d cells fiducial volume: ieta %d, iphi %d, SM %d ?",
				   fNCellsFromEMCALBorder, ieta, iphi, iSM);
			if (okcol && okrow ) printf(" YES \n");
			else  printf(" NO: column ok? %d, row ok? %d \n",okcol,okrow);
		}
	}//EMCAL
	else if(cells->GetType()==AliAODCaloCells::kPHOS){
		Int_t relId[4];
		Int_t irow = -1, icol = -1;
		fPHOSGeo->AbsToRelNumbering(absIdMax,relId);
		irow = relId[2];
		icol = relId[3];
		//imod = relId[0]-1;
		if(irow >= fNCellsFromPHOSBorder && irow < 64-fNCellsFromPHOSBorder) okrow =kTRUE; 
		if(icol >= fNCellsFromPHOSBorder && icol < 56-fNCellsFromPHOSBorder) okcol =kTRUE; 
		if(fDebug > 1)
		{
			printf("AliCalorimeterUtils::CheckCellFiducialRegion(AOD) - PHOS Cluster in %d cells fiducial volume: icol %d, irow %d, Module %d?",
				   fNCellsFromPHOSBorder, icol, irow, relId[0]-1);
			if (okcol && okrow ) printf(" YES \n");
			else  printf(" NO: column ok? %d, row ok? %d \n",okcol,okrow);
		}
	}//PHOS
	
	if (okcol && okrow) return kTRUE; 
	else                return kFALSE;
	
}	

//_______________________________________________________________
Bool_t AliCalorimeterUtils::CheckCellFiducialRegion(AliESDCaloCluster* cluster, AliESDCaloCells* cells) const {
	// Given the list of AbsId of the cluster, get the maximum cell and 
	// check if there are fNCellsFromBorder from the calorimeter border
	
	//If the distance to the border is 0 or negative just exit accept all clusters
	if(cells->GetType()==AliESDCaloCells::kEMCALCell && fNCellsFromEMCALBorder <= 0 ) return kTRUE;
	if(cells->GetType()==AliESDCaloCells::kPHOSCell  && fNCellsFromPHOSBorder  <= 0 ) return kTRUE;
	
	//Find cell with maximum amplitude
	Int_t absIdMax	= -1;
	Float_t ampMax  = -1;
	for(Int_t i = 0; i < cluster->GetNCells() ; i++){
		Int_t absId = cluster->GetCellAbsId(i) ;
		Float_t amp	= cells->GetCellAmplitude(absId);
		if(amp > ampMax){
			ampMax   = amp;
			absIdMax = absId;
		}
	}
	
	if(fDebug > 1)
		printf("AliCalorimeterUtils::CheckCellFiducialRegion(ESD) - Cluster Max AbsId %d, Cell Energy %2.2f, Cluster Energy %2.2f\n", 
						 absIdMax, ampMax, cluster->E());
	if(absIdMax==-1) return kFALSE;
	
	//Check if the cell is close to the borders:
	Bool_t okrow = kFALSE;
	Bool_t okcol = kFALSE;
	
	if(cells->GetType()==AliESDCaloCells::kEMCALCell){
		
		Int_t iTower = -1, iIphi = -1, iIeta = -1, iphi = -1, ieta = -1, iSM = -1; 
		fEMCALGeo->GetCellIndex(absIdMax,iSM,iTower,iIphi,iIeta); 
		fEMCALGeo->GetCellPhiEtaIndexInSModule(iSM,iTower,iIphi, iIeta,iphi,ieta);
		
		//Check rows/phi
		if(iSM < 10){
			if(iphi >= fNCellsFromEMCALBorder && iphi < 24-fNCellsFromEMCALBorder) okrow =kTRUE; 
	    }
		else{
			if(iphi >= fNCellsFromEMCALBorder && iphi < 12-fNCellsFromEMCALBorder) okrow =kTRUE; 
		}
		
		//Check collumns/eta
		if(!fNoEMCALBorderAtEta0){
			if(ieta  > fNCellsFromEMCALBorder && ieta < 48-fNCellsFromEMCALBorder) okcol =kTRUE; 
		}
		else{
			if(iSM%2==0){
				if(ieta >= fNCellsFromEMCALBorder)     okcol = kTRUE;	
			}
			else {
				if(ieta <  48-fNCellsFromEMCALBorder)  okcol = kTRUE;	
			}
		}//eta 0 not checked
		if(fDebug > 1)
		{
			printf("AliCalorimeterUtils::CheckCellFiducialRegion(ESD) - EMCAL Cluster in %d cells fiducial volume: ieta %d, iphi %d, SM %d ?",
				   fNCellsFromEMCALBorder, ieta, iphi, iSM);
			if (okcol && okrow ) printf(" YES \n");
			else  printf(" NO: column ok? %d, row ok? %d \n",okcol,okrow);
		}
	}//EMCAL
	else if(cells->GetType()==AliESDCaloCells::kPHOSCell){
		Int_t relId[4];
		Int_t irow = -1, icol = -1;
		fPHOSGeo->AbsToRelNumbering(absIdMax,relId);
		irow = relId[2];
		icol = relId[3];
		//imod = relId[0]-1;
		if(irow >= fNCellsFromPHOSBorder && irow < 64-fNCellsFromPHOSBorder) okrow =kTRUE; 
		if(icol >= fNCellsFromPHOSBorder && icol < 56-fNCellsFromPHOSBorder) okcol =kTRUE; 
		if(fDebug > 1)
		{
			printf("AliCalorimeterUtils::CheckCellFiducialRegion(ESD) - PHOS Cluster in %d cells fiducial volume: icol %d, irow %d, Module %d ?",
				   fNCellsFromPHOSBorder, icol, irow,relId[0]-1);
			if (okcol && okrow ) printf(" YES \n");
			else  printf(" NO: column ok? %d, row ok? %d \n",okcol,okrow);
		}
	}//PHOS
	
	if (okcol && okrow) return kTRUE; 
	else                return kFALSE;
	
}	


//_________________________________________________________________________________________________________
Bool_t AliCalorimeterUtils::ClusterContainsBadChannel(TString calorimeter,UShort_t* cellList, Int_t nCells){
	// Check that in the cluster cells, there is no bad channel of those stored 
	// in fEMCALBadChannelMap or fPHOSBadChannelMap
	
	if (!fRemoveBadChannels) return kFALSE;
	
	if(calorimeter == "EMCAL" && !fEMCALBadChannelMap->GetEntries()) return kFALSE;
	if(calorimeter == "PHOS"  && !fPHOSBadChannelMap ->GetEntries()) return kFALSE;

	Int_t icol = -1;
	Int_t irow = -1;
	Int_t imod = -1;
	for(Int_t iCell = 0; iCell<nCells; iCell++){
	
		//Get the column and row
		if(calorimeter == "EMCAL"){
			Int_t iTower = -1, iIphi = -1, iIeta = -1; 
			fEMCALGeo->GetCellIndex(cellList[iCell],imod,iTower,iIphi,iIeta); 
			if(fEMCALBadChannelMap->GetEntries() <= imod) continue;
			fEMCALGeo->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,irow,icol);			
			if(GetEMCALChannelStatus(imod, icol, irow))return kTRUE;
		}
		else if(calorimeter=="PHOS"){
			Int_t    relId[4];
			fPHOSGeo->AbsToRelNumbering(cellList[iCell],relId);
			irow = relId[2];
			icol = relId[3];
			imod = relId[0]-1;
			if(fPHOSBadChannelMap->GetEntries() <= imod)continue;
			if(GetPHOSChannelStatus(imod, icol, irow)) return kTRUE;
		}
		else return kFALSE;
		
	}// cell cluster loop

	return kFALSE;

}

//____________________________________________________________________________________________________________________________________________________
Int_t AliCalorimeterUtils::GetModuleNumber(AliAODPWG4Particle * particle, AliVEvent * inputEvent) const
{
	//Get the EMCAL/PHOS module number that corresponds to this particle
	
	Int_t absId = -1;
	if(particle->GetDetector()=="EMCAL"){
		fEMCALGeo->GetAbsCellIdFromEtaPhi(particle->Eta(),particle->Phi(), absId);
		if(fDebug > 2) 
		  printf("AliCalorimeterUtils::GetModuleNumber(PWG4AOD) - EMCAL: cluster eta %f, phi %f, absid %d, SuperModule %d\n",
			 particle->Eta(), particle->Phi()*TMath::RadToDeg(),absId, fEMCALGeo->GetSuperModuleNumber(absId));
		return fEMCALGeo->GetSuperModuleNumber(absId) ;
	}//EMCAL
	else if(particle->GetDetector()=="PHOS"){
		if(!strcmp(inputEvent->GetName(),"AliESDEvent"))   {
			AliESDCaloCluster *cluster = ((AliESDEvent*)inputEvent)->GetCaloCluster(particle->GetCaloLabel(0));
			return GetModuleNumber(cluster);
		}//ESDs
		else{
			AliAODCaloCluster *cluster = ((AliAODEvent*)inputEvent)->GetCaloCluster(particle->GetCaloLabel(0));
			return GetModuleNumber(cluster);
		}//AODs
	}//PHOS
	
	return -1;
}

//____________________________________________________________________________________________________________________________________________________
Int_t AliCalorimeterUtils::GetModuleNumber(AliAODCaloCluster * cluster) const
{
	//Get the EMCAL/PHOS module number that corresponds to this cluster, input are AODs
	TLorentzVector lv;
	Double_t v[]={0.,0.,0.}; //not necessary to pass the real vertex.
	cluster->GetMomentum(lv,v);
	Float_t phi = lv.Phi();
	if(phi < 0) phi+=TMath::TwoPi();	
	Int_t absId = -1;
	if(cluster->IsEMCALCluster()){
		fEMCALGeo->GetAbsCellIdFromEtaPhi(lv.Eta(),phi, absId);
		if(fDebug > 2) 
		  printf("AliCalorimeterUtils::GetModuleNumber(ESD) - EMCAL: cluster eta %f, phi %f, absid %d, SuperModule %d\n",
			 lv.Eta(), phi*TMath::RadToDeg(),absId, fEMCALGeo->GetSuperModuleNumber(absId));
		return fEMCALGeo->GetSuperModuleNumber(absId) ;
	}//EMCAL
	else if(cluster->IsPHOSCluster()) {
		Int_t    relId[4];
		if ( cluster->GetNCells() > 0) {
			absId = cluster->GetCellAbsId(0);
			if(fDebug > 2) 
				printf("AliCalorimeterUtils::GetModuleNumber(AOD) - PHOS: cluster eta %f, phi %f, e %f, absId %d\n",
					   lv.Eta(), phi*TMath::RadToDeg(), lv.E(), absId);
		}
		else return -1;
		
		if ( absId >= 0) {
			fPHOSGeo->AbsToRelNumbering(absId,relId);
			if(fDebug > 2) 
			  printf("AliCalorimeterUtils::GetModuleNumber(AOD) - PHOS: Module %d\n",relId[0]-1);
			return relId[0]-1;
		}
		else return -1;
	}//PHOS
	
	return -1;
}

//____________________________________________________________________________________________________________________________________________________
Int_t AliCalorimeterUtils::GetModuleNumber(AliESDCaloCluster * cluster) const
{
	//Get the EMCAL/PHOS module number that corresponds to this cluster, input are ESDs
	TLorentzVector lv;
	Double_t v[]={0.,0.,0.}; //not necessary to pass the real vertex.
	cluster->GetMomentum(lv,v);
	Float_t phi = lv.Phi();
	if(phi < 0) phi+=TMath::TwoPi();	
	Int_t absId = -1;
	if(cluster->IsEMCAL()){
		fEMCALGeo->GetAbsCellIdFromEtaPhi(lv.Eta(),phi, absId);
		if(fDebug > 2) 
		  printf("AliCalorimeterUtils::GetModuleNumber(ESD) - EMCAL: cluster eta %f, phi %f, absid %d, SuperModule %d\n",
			 lv.Eta(), phi*TMath::RadToDeg(),absId, fEMCALGeo->GetSuperModuleNumber(absId));
		return fEMCALGeo->GetSuperModuleNumber(absId) ;
	}//EMCAL
	else if(cluster->IsPHOS()){
		Int_t    relId[4];
		if ( cluster->GetNCells() > 0) {
			absId = cluster->GetCellAbsId(0);
			if(fDebug > 2) 
			  printf("AliCalorimeterUtils::GetModuleNumber(ESD) - PHOS: cluster eta %f, phi %f, e %f, absId %d\n",
				 lv.Eta(), phi*TMath::RadToDeg(), lv.E(), absId);
		}
		else return -1;
		
		if ( absId >= 0) {
			fPHOSGeo->AbsToRelNumbering(absId,relId);
			if(fDebug > 2) 
			  printf("AliCalorimeterUtils::GetModuleNumber(ESD) - PHOS: Module %d\n",relId[0]-1);
			return relId[0]-1;
		}
		else return -1;
	}//PHOS
	
	return -1;
}


//_____________________________________________________________________________________________________________
Int_t AliCalorimeterUtils::GetModuleNumberCellIndexes(const Int_t absId, const TString calo, Int_t & icol, Int_t & irow, Int_t & iRCU) const
{
	//Get the EMCAL/PHOS module, columns, row and RCU number that corresponds to this absId
	Int_t imod = -1;
	if ( absId >= 0) {
		if(calo=="EMCAL"){
			Int_t iTower = -1, iIphi = -1, iIeta = -1; 
			fEMCALGeo->GetCellIndex(absId,imod,iTower,iIphi,iIeta); 
			fEMCALGeo->GetCellPhiEtaIndexInSModule(imod,iTower,iIphi, iIeta,irow,icol);
			
			//RCU0
			if (0<=irow&&irow<8) iRCU=0; // first cable row
			else if (8<=irow&&irow<16 && 0<=icol&&icol<24) iRCU=0; // first half; 
			//second cable row
			//RCU1
			else if(8<=irow&&irow<16 && 24<=icol&&icol<48) iRCU=1; // second half; 
			//second cable row
			else if(16<=irow&&irow<24) iRCU=1; // third cable row
			
			if (imod%2==1) iRCU = 1 - iRCU; // swap for odd=C side, to allow us to cable both sides the same
			if (iRCU<0) {
				printf("AliCalorimeterUtils::GetModuleNumberCellIndexes() - Wrong EMCAL RCU number = %d\n", iRCU);
				abort();
			}			
			
			return imod ;
		}//EMCAL
		else{//PHOS
			Int_t    relId[4];
			fPHOSGeo->AbsToRelNumbering(absId,relId);
			irow = relId[2];
			icol = relId[3];
			imod = relId[0]-1;
			iRCU= (Int_t)(relId[2]-1)/16 ;
			//Int_t iBranch= (Int_t)(relid[3]-1)/28 ; //0 to 1
			if (iRCU >= 4) {
				printf("AliCalorimeterUtils::GetModuleNumberCellIndexes() - Wrong PHOS RCU number = %d\n", iRCU);
				abort();
			}			
			return imod;
		}//PHOS	
	}
	
	return -1;
}

//_______________________________________________________________
//void AliCalorimeterUtils::Init()
//{
//	//Init reader. Method to be called in AliAnaPartCorrMaker
//	
//	fEMCALBadChannelMap->SetName(Form("EMCALBadMap_%s",fTaskName.Data()));
//	fPHOSBadChannelMap->SetName(Form("PHOSBadMap_%s",fTaskName.Data()));
//}

//_______________________________________________________________
void AliCalorimeterUtils::InitParameters()
{
  //Initialize the parameters of the analysis.
  fEMCALGeoName = "EMCAL_COMPLETE";
  fPHOSGeoName  = "PHOSgeo";
	
  if(gGeoManager) {// geoManager was set
	if(fDebug > 2)printf("AliCalorimeterUtils::InitParameters() - Geometry manager available\n");
	fEMCALGeoMatrixSet = kTRUE;	 
	fPHOSGeoMatrixSet  = kTRUE;	 
  }
  else{
	fEMCALGeoMatrixSet = kFALSE;
	fPHOSGeoMatrixSet  = kFALSE;
  }
		
  fRemoveBadChannels   = kFALSE;
	
  fNCellsFromEMCALBorder   = 0;
  fNCellsFromPHOSBorder    = 0;
  fNoEMCALBorderAtEta0     = kFALSE;
}

//________________________________________________________________
void AliCalorimeterUtils::InitEMCALBadChannelStatusMap(){
  //Init EMCAL bad channels map
   if(fDebug > 0 )printf("AliCalorimeterUtils::InitEMCALBadChannelStatusMap()\n");
  //In order to avoid rewriting the same histograms
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  //for (int i = 0; i < 12; i++) fEMCALBadChannelMap->Add(new TH2I(Form("EMCALBadChannelMap_SM%d_%s",i,fTaskName.Data()),Form("EMCALBadChannelMap_SM%d",i),  48, 0, 48, 24, 0, 24));
  for (int i = 0; i < 12; i++) fEMCALBadChannelMap->Add(new TH2I(Form("EMCALBadChannelMap_SM%d",i),Form("EMCALBadChannelMap_SM%d",i),  48, 0, 48, 24, 0, 24));

  fEMCALBadChannelMap->SetOwner(kTRUE);
  fEMCALBadChannelMap->Compress();
  
  //In order to avoid rewriting the same histograms
  TH1::AddDirectory(oldStatus);		
}

//________________________________________________________________
void AliCalorimeterUtils::InitPHOSBadChannelStatusMap(){
  //Init PHOS bad channels map
  if(fDebug > 0 )printf("AliCalorimeterUtils::InitPHOSBadChannelStatusMap()\n");
  //In order to avoid rewriting the same histograms
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  //for (int i = 0; i < 5; i++)fPHOSBadChannelMap->Add(new TH2I(Form("PHOSBadChannelMap_Mod%d_%s",i,fTaskName.Data()),Form("PHOSBadChannelMap_Mod%d",i), 56, 0, 56, 64, 0, 64));
  for (int i = 0; i < 5; i++)fPHOSBadChannelMap->Add(new TH2I(Form("PHOSBadChannelMap_Mod%d",i),Form("PHOSBadChannelMap_Mod%d",i), 56, 0, 56, 64, 0, 64));

  fPHOSBadChannelMap->SetOwner(kTRUE);
  fPHOSBadChannelMap->Compress();
  
  //In order to avoid rewriting the same histograms
  TH1::AddDirectory(oldStatus);		
}

//________________________________________________________________
void AliCalorimeterUtils::InitEMCALGeometry()
{
	//Initialize EMCAL geometry if it did not exist previously
	if (!fEMCALGeo){
		fEMCALGeo = new AliEMCALGeoUtils(fEMCALGeoName);
		if(fDebug > 0){
			printf("AliCalorimeterUtils::InitEMCALGeometry()");
			if (!gGeoManager) printf(" - Careful!, gGeoManager not loaded, load misalign matrices");
			printf("\n");
		}
	}
}

//________________________________________________________________
void AliCalorimeterUtils::InitPHOSGeometry()
{
	//Initialize PHOS geometry if it did not exist previously
	if (!fPHOSGeo){
		fPHOSGeo = new AliPHOSGeoUtils(fPHOSGeoName); 
		if(fDebug > 0){
			printf("AliCalorimeterUtils::InitPHOSGeometry()");
			if (!gGeoManager) printf(" - Careful!, gGeoManager not loaded, load misalign matrices");
			printf("\n");
		}	
	}	
}

//________________________________________________________________
void AliCalorimeterUtils::Print(const Option_t * opt) const
{

  //Print some relevant parameters set for the analysis
  if(! opt)
    return;

  printf("***** Print: %s %s ******\n", GetName(), GetTitle() ) ;
  printf("Remove Clusters with bad channels? %d\n",fRemoveBadChannels);
  printf("Remove Clusters with max cell at less than %d cells from EMCAL border and %d cells from PHOS border\n",
	 fNCellsFromEMCALBorder, fNCellsFromPHOSBorder);
  if(fNoEMCALBorderAtEta0) printf("Do not remove EMCAL clusters at Eta = 0\n");

  printf("    \n") ;
} 

//________________________________________________________________
void AliCalorimeterUtils::SetGeometryTransformationMatrices(AliVEvent* inputEvent) 
{
  //Set the calorimeters transformation matrices
	
  //Get the EMCAL transformation geometry matrices from ESD 
  if (!gGeoManager && fEMCALGeo) {//&& !fEMCALGeoMatrixSet) { FIXME
	  if(fDebug > 1) 
		  printf(" AliCalorimeterUtils::SetGeometryTransformationMatrices() - Load EMCAL misalignment matrices. \n");
	  if(!strcmp(inputEvent->GetName(),"AliESDEvent"))  {
			for(Int_t mod=0; mod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++){ 
				if(((AliESDEvent*)inputEvent)->GetEMCALMatrix(mod)) {
					//printf("EMCAL: mod %d, matrix %p\n",mod, ((AliESDEvent*)inputEvent)->GetEMCALMatrix(mod));
					fEMCALGeo->SetMisalMatrix(((AliESDEvent*)inputEvent)->GetEMCALMatrix(mod),mod) ;
					fEMCALGeoMatrixSet = kTRUE;//At least one, so good
				}
			}// loop over super modules	
		}//ESD as input
		else {
			if(fDebug > 1)
				printf("AliCalorimeterUtils::SetGeometryTransformationMatrices() - Setting of EMCAL transformation matrixes for AODs not implemented yet. \n Import geometry.root file\n");
				}//AOD as input
  }//EMCAL geo && no geoManager
	
	//Get the PHOS transformation geometry matrices from ESD 
	if (!gGeoManager && fPHOSGeo && !fPHOSGeoMatrixSet) {
		if(fDebug > 1) 
			printf(" AliCalorimeterUtils::SetGeometryTransformationMatrices() - Load PHOS misalignment matrices. \n");
			if(!strcmp(inputEvent->GetName(),"AliESDEvent"))  {
				for(Int_t mod=0; mod < 5; mod++){ 
					if(((AliESDEvent*)inputEvent)->GetPHOSMatrix(mod)) {
						//printf("PHOS: mod %d, matrix %p\n",mod, ((AliESDEvent*)inputEvent)->GetPHOSMatrix(mod));
						fPHOSGeo->SetMisalMatrix(((AliESDEvent*)inputEvent)->GetPHOSMatrix(mod),mod) ;
						fPHOSGeoMatrixSet  = kTRUE; //At least one so good
					}
				}// loop over modules	
			}//ESD as input
			else {
				if(fDebug > 1) 
					printf("AliCalorimeterUtils::SetGeometryTransformationMatrices() - Setting of EMCAL transformation matrixes for AODs not implemented yet. \n Import geometry.root file\n");
					}//AOD as input
	}//PHOS geo	and  geoManager was not set

}

