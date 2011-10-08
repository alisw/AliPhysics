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
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliAODPWG4Particle.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "AliMixedEvent.h"

ClassImp(AliCalorimeterUtils)
  
  
//____________________________________________________________________________
  AliCalorimeterUtils::AliCalorimeterUtils() : 
    TObject(), fDebug(0), 
    fEMCALGeoName("EMCAL_COMPLETEV1"),fPHOSGeoName("PHOSgeo"), 
    fEMCALGeo(0x0),                   fPHOSGeo(0x0), 
    fEMCALGeoMatrixSet(kFALSE),       fPHOSGeoMatrixSet(kFALSE), 
    fLoadEMCALMatrices(kFALSE),       fLoadPHOSMatrices(kFALSE),
    fRemoveBadChannels(kFALSE),       fPHOSBadChannelMap(0x0), 
    fNCellsFromPHOSBorder(0),
    fNMaskCellColumns(0),             fMaskCellColumns(0x0),
    fRecalibration(kFALSE),           fPHOSRecalibrationFactors(),
    fEMCALRecoUtils(new AliEMCALRecoUtils),
    fRecalculatePosition(kFALSE),     fCorrectELinearity(kFALSE),
    fRecalculateMatching(kFALSE),
    fCutR(20),                        fCutZ(20),
    fCutEta(20),                      fCutPhi(20)
{
  //Ctor
  
  //Initialize parameters
  InitParameters();
  for(Int_t i = 0; i < 10; i++) fEMCALMatrix[i] = 0 ;
  for(Int_t i = 0; i < 5 ; i++) fPHOSMatrix[i]  = 0 ;
  
}

//_________________________________
AliCalorimeterUtils::~AliCalorimeterUtils() {
  //Dtor
  
  //if(fPHOSGeo)  delete fPHOSGeo  ;
  if(fEMCALGeo) delete fEMCALGeo ;
  
  if(fPHOSBadChannelMap) { 
    fPHOSBadChannelMap->Clear();
    delete  fPHOSBadChannelMap;
  }
	
  if(fPHOSRecalibrationFactors) { 
    fPHOSRecalibrationFactors->Clear();
    delete  fPHOSRecalibrationFactors;
  }
	
  if(fEMCALRecoUtils)   delete fEMCALRecoUtils ;
  if(fNMaskCellColumns) delete [] fMaskCellColumns;
  
}

//_______________________________________________________________
Bool_t AliCalorimeterUtils::CheckCellFiducialRegion(AliVCluster* cluster, AliVCaloCells* cells, AliVEvent * event, Int_t iev) const {
  
	// Given the list of AbsId of the cluster, get the maximum cell and 
	// check if there are fNCellsFromBorder from the calorimeter border
	
  //If the distance to the border is 0 or negative just exit accept all clusters
	if(cells->GetType()==AliVCaloCells::kEMCALCell && fEMCALRecoUtils->GetNumberOfCellsFromEMCALBorder() <= 0 ) return kTRUE;
	if(cells->GetType()==AliVCaloCells::kPHOSCell  && fNCellsFromPHOSBorder  <= 0 ) return kTRUE;
  
  Int_t absIdMax	= -1;
	Float_t ampMax  = -1;
	
  AliMixedEvent * mixEvent = dynamic_cast<AliMixedEvent*> (event);
  Int_t nMixedEvents = 0 ; 
  Int_t * cellsCumul = NULL ;
  Int_t numberOfCells = 0 ;  
  if (mixEvent){
    nMixedEvents = mixEvent->GetNumberOfEvents() ; 
    if (cells->GetType()==AliVCaloCells::kEMCALCell) {
      cellsCumul =  mixEvent->GetEMCALCellsCumul() ; 
      numberOfCells = mixEvent->GetNumberOfEMCALCells() ;
    } 
    
    else if (cells->GetType()==AliVCaloCells::kPHOSCell) {
      cellsCumul =  mixEvent->GetPHOSCellsCumul() ; 
      numberOfCells = mixEvent->GetNumberOfPHOSCells() ;
    } 
    
    if(cellsCumul){
      
      Int_t startCell = cellsCumul[iev] ; 
      Int_t endCell   = (iev+1 < nMixedEvents)?cellsCumul[iev+1]:numberOfCells;
      //Find cells with maximum amplitude
      for(Int_t i = 0; i < cluster->GetNCells() ; i++){
        Int_t absId = cluster->GetCellAbsId(i) ;
        for (Int_t j = startCell; j < endCell ;  j++) {
          Short_t cellNumber; 
          Double_t amp ; 
          Double_t time; 
          cells->GetCell(j, cellNumber, amp, time) ; 
          if (absId == cellNumber) {
            if(amp > ampMax){
              ampMax   = amp;
              absIdMax = absId;
            }        
          }
        }
      }//loop on cluster cells
    }// cells cumul available
    else {
      printf("AliCalorimeterUtils::CheckCellFiducialRegion() - CellsCumul is NULL!!!\n");
      abort();
    }
  } else {//Normal SE Events
    for(Int_t i = 0; i < cluster->GetNCells() ; i++){
      Int_t absId = cluster->GetCellAbsId(i) ;
      Float_t amp	= cells->GetCellAmplitude(absId);
      if(amp > ampMax){
        ampMax   = amp;
        absIdMax = absId;
      }
    }
  }
	
	if(fDebug > 1)
		printf("AliCalorimeterUtils::CheckCellFiducialRegion() - Cluster Max AbsId %d, Cell Energy %2.2f, Cluster Energy %2.2f\n", 
           absIdMax, ampMax, cluster->E());
	
	if(absIdMax==-1) return kFALSE;
	
	//Check if the cell is close to the borders:
	Bool_t okrow = kFALSE;
	Bool_t okcol = kFALSE;
  
	if(cells->GetType()==AliVCaloCells::kEMCALCell){
    
		Int_t iTower = -1, iIphi = -1, iIeta = -1, iphi = -1, ieta = -1, iSM = -1; 
		fEMCALGeo->GetCellIndex(absIdMax,iSM,iTower,iIphi,iIeta); 
		fEMCALGeo->GetCellPhiEtaIndexInSModule(iSM,iTower,iIphi, iIeta,iphi,ieta);
		if(iSM < 0 || iphi < 0 || ieta < 0 ) {
      Fatal("CheckCellFidutialRegion","Negative value for super module: %d, or cell ieta: %d, or cell iphi: %d, check EMCAL geometry name\n",iSM,ieta,iphi);
    }
    
		//Check rows/phi
    Int_t nborder = fEMCALRecoUtils->GetNumberOfCellsFromEMCALBorder();
		if(iSM < 10){
			if(iphi >= nborder && iphi < 24-nborder) okrow =kTRUE; 
    }
		else{
			if(iphi >= nborder && iphi < 12-nborder) okrow =kTRUE; 
		}
		
		//Check columns/eta
		if(!fEMCALRecoUtils->IsEMCALNoBorderAtEta0()){
			if(ieta  > nborder && ieta < 48-nborder) okcol =kTRUE; 
		}
		else{
			if(iSM%2==0){
				if(ieta >= nborder)     okcol = kTRUE;	
			}
			else {
				if(ieta <  48-nborder)  okcol = kTRUE;	
			}
		}//eta 0 not checked
		if(fDebug > 1)
		{
			printf("AliCalorimeterUtils::CheckCellFiducialRegion() - EMCAL Cluster in %d cells fiducial volume: ieta %d, iphi %d, SM %d ?",
             nborder, ieta, iphi, iSM);
			if (okcol && okrow ) printf(" YES \n");
			else  printf(" NO: column ok? %d, row ok? %d \n",okcol,okrow);
		}
	}//EMCAL
	else if(cells->GetType()==AliVCaloCells::kPHOSCell){
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
			printf("AliCalorimeterUtils::CheckCellFiducialRegion() - PHOS Cluster in %d cells fiducial volume: icol %d, irow %d, Module %d?",
             fNCellsFromPHOSBorder, icol, irow, relId[0]-1);
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
	//printf("fEMCALBadChannelMap %p, fPHOSBadChannelMap %p \n",fEMCALBadChannelMap,fPHOSBadChannelMap);
	if(calorimeter == "EMCAL" && !fEMCALRecoUtils->GetEMCALChannelStatusMap(0)) return kFALSE;
	if(calorimeter == "PHOS"  && !fPHOSBadChannelMap)  return kFALSE;
  
	Int_t icol = -1;
	Int_t irow = -1;
	Int_t imod = -1;
	for(Int_t iCell = 0; iCell<nCells; iCell++){
    
		//Get the column and row
		if(calorimeter == "EMCAL"){
      return fEMCALRecoUtils->ClusterContainsBadChannel((AliEMCALGeometry*)fEMCALGeo,cellList,nCells);
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
void AliCalorimeterUtils::CorrectClusterEnergy(AliVCluster *clus){
  // Correct cluster energy non linearity
  
  clus->SetE(fEMCALRecoUtils->CorrectClusterEnergyLinearity(clus));

}

//____________________________________________________________________________________________________________________________________________________
Int_t  AliCalorimeterUtils::GetMaxEnergyCell(AliVCaloCells* cells, AliVCluster* clu, Float_t & clusterFraction) const {
  
  //For a given CaloCluster gets the absId of the cell 
  //with maximum energy deposit.
  
  if( !clu || !cells ){
    AliInfo("Cluster or cells pointer is null!");
    return -1;
  }
  
  Double_t eMax        =-1.;
  Double_t eTot        = 0.;
  Double_t eCell       =-1.;
  Float_t  fraction    = 1.;
  Float_t  recalFactor = 1.;
  Int_t    cellAbsId   =-1 , absId =-1 ;
  Int_t    iSupMod     =-1 , ieta  =-1 , iphi = -1, iRCU = -1;
  
  TString           calo = "EMCAL";
  if(clu->IsPHOS()) calo = "PHOS";
  
  for (Int_t iDig=0; iDig< clu->GetNCells(); iDig++) {
    
    cellAbsId = clu->GetCellAbsId(iDig);
    
    fraction  = clu->GetCellAmplitudeFraction(iDig);
    if(fraction < 1e-4) fraction = 1.; // in case unfolding is off
    
    iSupMod = GetModuleNumberCellIndexes(cellAbsId, calo, ieta, iphi, iRCU);
    
    if(IsRecalibrationOn()) {
      if(calo=="EMCAL") recalFactor = GetEMCALChannelRecalibrationFactor(iSupMod,ieta,iphi);
      else              recalFactor = GetPHOSChannelRecalibrationFactor(iSupMod,ieta,iphi);
    }
    
    eCell  = cells->GetCellAmplitude(cellAbsId)*fraction*recalFactor;
    
    if(eCell > eMax)  { 
      eMax  = eCell; 
      absId = cellAbsId;
    }
    
    eTot+=eCell;
    
  }// cell loop
  
  clusterFraction = (eTot-eMax)/eTot; //Do not use cluster energy in case it was corrected for non linearity.
  
  return absId;
  
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
    // In case we use the MC reader, the input are TParticles, 
    // in this case use the corresponing method in PHOS Geometry to get the particle.
    if(strcmp(inputEvent->ClassName(), "AliMCEvent") == 0 )
    {
      Int_t mod =-1;
      Double_t z = 0., x=0.;
      TParticle* primary = 0x0;
      AliStack * stack = ((AliMCEvent*)inputEvent)->Stack();
      if(stack) {
        primary = stack->Particle(particle->GetCaloLabel(0));
      }
      else {
        Fatal("GetModuleNumber(PWG4AOD)", "Stack not available, stop!");
      }
      
      if(primary){
        fPHOSGeo->ImpactOnEmc(primary,mod,z,x) ;
      }
      else{
        Fatal("GetModuleNumber(PWG4AOD)", "Primary not available, stop!");
      }
      return mod;
    }
    // Input are ESDs or AODs, get the PHOS module number like this.
    else{
      //FIXME
      //AliVCluster *cluster = inputEvent->GetCaloCluster(particle->GetCaloLabel(0));
      //return GetModuleNumber(cluster);
      //MEFIX
      return -1;
    }
	}//PHOS
	
	return -1;
}

//____________________________________________________________________________________________________________________________________________________
Int_t AliCalorimeterUtils::GetModuleNumber(AliVCluster * cluster) const
{
	//Get the EMCAL/PHOS module number that corresponds to this cluster
	TLorentzVector lv;
	Double_t v[]={0.,0.,0.}; //not necessary to pass the real vertex.
  if(!cluster){
    if(fDebug > 1) printf("AliCalorimeterUtils::GetModuleNumber() - NUL Cluster, please check!!!");
    return -1;
  }
	cluster->GetMomentum(lv,v);
	Float_t phi = lv.Phi();
	if(phi < 0) phi+=TMath::TwoPi();	
	Int_t absId = -1;
	if(cluster->IsEMCAL()){
		fEMCALGeo->GetAbsCellIdFromEtaPhi(lv.Eta(),phi, absId);
		if(fDebug > 2) 
		  printf("AliCalorimeterUtils::GetModuleNumber() - EMCAL: cluster eta %f, phi %f, absid %d, SuperModule %d\n",
             lv.Eta(), phi*TMath::RadToDeg(),absId, fEMCALGeo->GetSuperModuleNumber(absId));
		return fEMCALGeo->GetSuperModuleNumber(absId) ;
	}//EMCAL
	else if(cluster->IsPHOS()) {
		Int_t    relId[4];
		if ( cluster->GetNCells() > 0) {
			absId = cluster->GetCellAbsId(0);
			if(fDebug > 2) 
				printf("AliCalorimeterUtils::GetModuleNumber() - PHOS: cluster eta %f, phi %f, e %f, absId %d\n",
               lv.Eta(), phi*TMath::RadToDeg(), lv.E(), absId);
		}
		else return -1;
		
		if ( absId >= 0) {
			fPHOSGeo->AbsToRelNumbering(absId,relId);
			if(fDebug > 2) 
			  printf("AliCalorimeterUtils::GetModuleNumber() - PHOS: Module %d\n",relId[0]-1);
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
      if(imod < 0 || irow < 0 || icol < 0 ) {
        Fatal("GetModuleNumberCellIndexes()","Negative value for super module: %d, or cell icol: %d, or cell irow: %d, check EMCAL geometry name\n",imod,icol,irow);
      }
      
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
				Fatal("GetModuleNumberCellIndexes()","Wrong EMCAL RCU number = %d\n", iRCU);
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
				Fatal("GetModuleNumberCellIndexes()","Wrong PHOS RCU number = %d\n", iRCU);
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
  fEMCALGeoName         = "EMCAL_COMPLETEV1";
  fPHOSGeoName          = "PHOSgeo";
  fEMCALGeoMatrixSet    = kFALSE;
  fPHOSGeoMatrixSet     = kFALSE;
  fRemoveBadChannels    = kFALSE;
  fNCellsFromPHOSBorder = 0;
  
  //  fMaskCellColumns = new Int_t[fNMaskCellColumns];
  //  fMaskCellColumns[0] = 6 ;  fMaskCellColumns[1] = 7 ;  fMaskCellColumns[2] = 8 ; 
  //  fMaskCellColumns[3] = 35;  fMaskCellColumns[4] = 36;  fMaskCellColumns[5] = 37; 
  //  fMaskCellColumns[6] = 12+AliEMCALGeoParams::fgkEMCALCols; fMaskCellColumns[7] = 13+AliEMCALGeoParams::fgkEMCALCols;
  //  fMaskCellColumns[8] = 40+AliEMCALGeoParams::fgkEMCALCols; fMaskCellColumns[9] = 41+AliEMCALGeoParams::fgkEMCALCols; 
  //  fMaskCellColumns[10]= 42+AliEMCALGeoParams::fgkEMCALCols; 
  
}


//________________________________________________________________
void AliCalorimeterUtils::InitPHOSBadChannelStatusMap(){
  //Init PHOS bad channels map
  if(fDebug > 0 )printf("AliCalorimeterUtils::InitPHOSBadChannelStatusMap()\n");
  //In order to avoid rewriting the same histograms
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  
  fPHOSBadChannelMap = new TObjArray(5);	
  for (int i = 0; i < 5; i++)fPHOSBadChannelMap->Add(new TH2I(Form("PHOSBadChannelMap_Mod%d",i),Form("PHOSBadChannelMap_Mod%d",i), 56, 0, 56, 64, 0, 64));
  
  fPHOSBadChannelMap->SetOwner(kTRUE);
  fPHOSBadChannelMap->Compress();
  
  //In order to avoid rewriting the same histograms
  TH1::AddDirectory(oldStatus);		
}

//________________________________________________________________
void AliCalorimeterUtils::InitPHOSRecalibrationFactors(){
	//Init EMCAL recalibration factors
	if(fDebug > 0 )printf("AliCalorimeterUtils::InitPHOSRecalibrationFactors()\n");
	//In order to avoid rewriting the same histograms
	Bool_t oldStatus = TH1::AddDirectoryStatus();
	TH1::AddDirectory(kFALSE);
  
	fPHOSRecalibrationFactors = new TObjArray(5);
	for (int i = 0; i < 5; i++)fPHOSRecalibrationFactors->Add(new TH2F(Form("PHOSRecalFactors_Mod%d",i),Form("PHOSRecalFactors_Mod%d",i), 56, 0, 56, 64, 0, 64));
	//Init the histograms with 1
	for (Int_t m = 0; m < 5; m++) {
		for (Int_t i = 0; i < 56; i++) {
			for (Int_t j = 0; j < 64; j++) {
				SetPHOSChannelRecalibrationFactor(m,i,j,1.);
			}
		}
	}
	fPHOSRecalibrationFactors->SetOwner(kTRUE);
	fPHOSRecalibrationFactors->Compress();
	
	//In order to avoid rewriting the same histograms
	TH1::AddDirectory(oldStatus);		
}


//________________________________________________________________
void AliCalorimeterUtils::InitEMCALGeometry()
{
	//Initialize EMCAL geometry if it did not exist previously
	if (!fEMCALGeo){
		fEMCALGeo = AliEMCALGeometry::GetInstance(fEMCALGeoName);
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
         fEMCALRecoUtils->GetNumberOfCellsFromEMCALBorder(), fNCellsFromPHOSBorder);
  if(fEMCALRecoUtils->IsEMCALNoBorderAtEta0()) printf("Do not remove EMCAL clusters at Eta = 0\n");
  printf("Recalibrate Clusters? %d\n",fRecalibration);
  printf("Recalculate Clusters Position? %d\n",fRecalculatePosition);
  printf("Recalculate Clusters Energy? %d\n",fCorrectELinearity);
  printf("Matching criteria: dR < %2.2f[cm], dZ < %2.2f[cm]\n",fCutR,fCutZ);
  
  printf("    \n") ;
} 

//__________________________________________________
Bool_t AliCalorimeterUtils::MaskFrameCluster(const Int_t iSM,  const Int_t ieta) const {
  //Check if cell is in one of the regions where we have significant amount 
  //of material in front. Only EMCAL
  
  Int_t icol = ieta;
  if(iSM%2) icol+=48; // Impair SM, shift index [0-47] to [48-96]
  
  if (fNMaskCellColumns && fMaskCellColumns) {
    for (Int_t imask = 0; imask < fNMaskCellColumns; imask++) {
      if(icol==fMaskCellColumns[imask]) return kTRUE;
    }
  }
  
  return kFALSE;
  
}


//________________________________________________________________
Float_t AliCalorimeterUtils::RecalibrateClusterEnergy(AliVCluster * cluster, AliVCaloCells * cells){
	// Recalibrate the cluster energy, considering the recalibration map and the energy of the cells that compose the cluster.
  
  //Initialize some used variables
	Float_t energy   = 0;
	Int_t absId      = -1;
	Int_t icol = -1, irow = -1, iRCU = -1, module=1;
	Float_t factor = 1, frac = 0;  
  
	if(cells) {
    
    //Get the cluster number of cells and list of absId, check what kind of cluster do we have.
    UShort_t * index    = cluster->GetCellsAbsId() ;
    Double_t * fraction = cluster->GetCellsAmplitudeFraction() ;
    Int_t ncells     = cluster->GetNCells();	
    TString calo     = "EMCAL";
    if(cluster->IsPHOS()) calo = "PHOS";
    
    //Loop on the cells, get the cell amplitude and recalibration factor, multiply and and to the new energy
    for(Int_t icell = 0; icell < ncells; icell++){
      absId = index[icell];
      frac =  fraction[icell];
      if(frac < 1e-3) frac = 1; //in case of EMCAL, this is set as 0, not used.
      module = GetModuleNumberCellIndexes(absId,calo,icol,irow,iRCU);
      if(cluster->IsPHOS()) factor = GetPHOSChannelRecalibrationFactor (module,icol,irow);
      else                  factor = GetEMCALChannelRecalibrationFactor(module,icol,irow);
      if(fDebug>2)
        printf("AliCalorimeterUtils::RecalibrateClusterEnergy() - recalibrate cell: %s, module %d, col %d, row %d, cell fraction %f, recalibration factor %f, cell energy %f\n", 
               calo.Data(),module,icol,irow,frac,factor,cells->GetCellAmplitude(absId));
      
      energy += cells->GetCellAmplitude(absId)*factor*frac;
    }
    
    if(fDebug>1)
      printf("AliCalorimeterUtils::RecalibrateClusterEnergy() - Energy before %f, after %f\n",cluster->E(),energy);
    
	}// cells available
  else{
    Fatal("RecalibrateClusterEnergy()","Cells pointer does not exist!");
  }
  
	return energy;
}

//________________________________________________________________
void AliCalorimeterUtils::SetGeometryTransformationMatrices(AliVEvent* inputEvent) 
{
  //Set the calorimeters transformation matrices
  
  //Get the EMCAL transformation geometry matrices from ESD 
  if(!fEMCALGeoMatrixSet && fEMCALGeo){
    if(fLoadEMCALMatrices){
      printf("AliCalorimeterUtils::SetGeometryTransformationMatrices() - Load user defined geometry matrices\n");
      for(Int_t mod=0; mod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++){
        if(fEMCALMatrix[mod]){
          if(fDebug > 1) 
            fEMCALMatrix[mod]->Print();
          fEMCALGeo->SetMisalMatrix(fEMCALMatrix[mod],mod) ;  
        }
      }//SM loop
      fEMCALGeoMatrixSet = kTRUE;//At least one, so good
      
    }//Load matrices
    else if (!gGeoManager) { 
      
      if(fDebug > 1) 
        printf(" AliCalorimeterUtils::SetGeometryTransformationMatrices() - Load EMCAL misalignment matrices. \n");
      if(!strcmp(inputEvent->GetName(),"AliESDEvent"))  {
        for(Int_t mod=0; mod < (fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++){ 
          //printf("Load ESD matrix %d, %p\n",mod,((AliESDEvent*)inputEvent)->GetEMCALMatrix(mod));
          if(((AliESDEvent*)inputEvent)->GetEMCALMatrix(mod)) {
            fEMCALGeo->SetMisalMatrix(((AliESDEvent*)inputEvent)->GetEMCALMatrix(mod),mod) ;
          }
        }// loop over super modules	
        fEMCALGeoMatrixSet = kTRUE;//At least one, so good
        
      }//ESD as input
      else {
        if(fDebug > 1)
          printf("AliCalorimeterUtils::SetGeometryTransformationMatrices() - Setting of EMCAL transformation matrixes for AODs not implemented yet. \n Import geometry.root file\n");
      }//AOD as input
    }//Get matrix from data
    else if(gGeoManager){
      fEMCALGeoMatrixSet = kTRUE;
    }
  }//EMCAL geo && no geoManager
	
	//Get the PHOS transformation geometry matrices from ESD 
  if(!fPHOSGeoMatrixSet && fPHOSGeo){
    if(fLoadPHOSMatrices){
      printf("AliCalorimeterUtils::SetGeometryTransformationMatrices() - Load user defined geometry matrices\n");
      for(Int_t mod=0; mod < 5; mod++){
        if(fPHOSMatrix[mod]){
          if(fDebug > 1) 
            fPHOSMatrix[mod]->Print();
          fPHOSGeo->SetMisalMatrix(fPHOSMatrix[mod],mod) ;  
        }
      }//SM loop
      fPHOSGeoMatrixSet = kTRUE;//At least one, so good
    }//Load matrices
    else if (!gGeoManager) { 
      if(fDebug > 1) 
        printf(" AliCalorimeterUtils::SetGeometryTransformationMatrices() - Load PHOS misalignment matrices. \n");
			if(!strcmp(inputEvent->GetName(),"AliESDEvent"))  {
				for(Int_t mod=0; mod < 5; mod++){ 
					if( ((AliESDEvent*)inputEvent)->GetPHOSMatrix(mod)) {
						//printf("PHOS: mod %d, matrix %p\n",mod, ((AliESDEvent*)inputEvent)->GetPHOSMatrix(mod));
						fPHOSGeo->SetMisalMatrix( ((AliESDEvent*)inputEvent)->GetPHOSMatrix(mod),mod) ;
					}
				}// loop over modules
        fPHOSGeoMatrixSet  = kTRUE; //At least one so good
			}//ESD as input
			else {
				if(fDebug > 1) 
					printf("AliCalorimeterUtils::SetGeometryTransformationMatrices() - Setting of EMCAL transformation matrixes for AODs not implemented yet. \n Import geometry.root file\n");
      }//AOD as input
    }// get matrix from data
    else if(gGeoManager){
      fPHOSGeoMatrixSet = kTRUE;
    }
	}//PHOS geo	and  geoManager was not set
}

//________________________________________________________________
void AliCalorimeterUtils::RecalculateClusterPosition(AliVCaloCells* cells, AliVCluster* clu){
  
  //Recalculate EMCAL cluster position
  
  fEMCALRecoUtils->RecalculateClusterPosition((AliEMCALGeometry*)fEMCALGeo, cells,clu);
  
}
