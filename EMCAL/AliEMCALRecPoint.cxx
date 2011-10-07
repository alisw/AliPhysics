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
/* $Id$ */
//_________________________________________________________________________
//  Reconstructed Points for the EMCAL
//  A RecPoint is a cluster of digits
//  
//  
//*-- Author: Yves Schutz (SUBATECH)
//*-- Author: Dmitri Peressounko (RRC KI & SUBATECH)
//*-- Author: Heather Gray (LBL) merged AliEMCALRecPoint and AliEMCALTowerRecPoint 02/04

// --- ROOT system ---
#include "TPad.h"
#include "TGraph.h"
#include "TPaveText.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TGeoMatrix.h"
#include "TGeoManager.h"
#include "TGeoPhysicalNode.h"
#include "TRandom.h"

// --- Standard library ---
#include <Riostream.h>

// --- AliRoot header files ---
//#include "AliGenerator.h"
class AliGenerator;
class AliEMCAL;
#include "AliLog.h"
#include "AliGeomManager.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALHit.h"
#include "AliEMCALDigit.h"
#include "AliEMCALRecPoint.h"
#include "AliCaloCalibPedestal.h"
#include "AliEMCALGeoParams.h"

ClassImp(AliEMCALRecPoint)

//____________________________________________________________________________
AliEMCALRecPoint::AliEMCALRecPoint()
  : AliCluster(), fGeomPtr(0),
    fAmp(0), fIndexInList(-1), //to be set when the point is already stored
    fGlobPos(0,0,0),fLocPos(0,0,0), 
    fMaxDigit(100), fMulDigit(0), fMaxTrack(200),
    fMulTrack(0), fDigitsList(0), fTracksList(0),
    fClusterType(-1), fCoreEnergy(0), fDispersion(0),
    fEnergyList(0), fAbsIdList(0),
    fTime(0.), fNExMax(0), fCoreRadius(10),  //HG check this 
    fDETracksList(0), fMulParent(0), fMaxParent(0),
    fParentsList(0), fDEParentsList(0), fSuperModuleNumber(0),
    fDigitIndMax(-1), fDistToBadTower(-1), fSharedCluster(kFALSE)
{
  // ctor
  fGeomPtr = AliEMCALGeometry::GetInstance();
  
  fLambda[0] = 0;
  fLambda[1] = 0;

}

//____________________________________________________________________________
AliEMCALRecPoint::AliEMCALRecPoint(const char *) 
  : AliCluster(), fGeomPtr(0),
    fAmp(0), fIndexInList(-1), //to be set when the point is already stored
    fGlobPos(0,0,0), fLocPos(0,0,0),
    fMaxDigit(100), fMulDigit(0), fMaxTrack(1000), fMulTrack(0),
    fDigitsList(new Int_t[fMaxDigit]), fTracksList(new Int_t[fMaxTrack]),
    fClusterType(-1), fCoreEnergy(0), fDispersion(0),
    fEnergyList(new Float_t[fMaxDigit]), 
    fAbsIdList(new Int_t[fMaxDigit]), fTime(-1.), fNExMax(0), fCoreRadius(10),
    fDETracksList(new Float_t[fMaxTrack]), fMulParent(0), fMaxParent(1000),
    fParentsList(new Int_t[fMaxParent]), fDEParentsList(new Float_t[fMaxParent]),
    fSuperModuleNumber(0), fDigitIndMax(-1), fDistToBadTower(-1),fSharedCluster(kFALSE)
{
  // ctor
  for (Int_t i = 0; i < fMaxTrack; i++)
    fDETracksList[i] = 0;
  for (Int_t i = 0; i < fMaxParent; i++) {
    fParentsList[i] = -1;
    fDEParentsList[i] = 0;
  }

  fGeomPtr = AliEMCALGeometry::GetInstance();
  fLambda[0] = 0;
  fLambda[1] = 0;
}

//____________________________________________________________________________
AliEMCALRecPoint::AliEMCALRecPoint(const AliEMCALRecPoint & rp) 
  : AliCluster(rp), fGeomPtr(rp.fGeomPtr),
    fAmp(rp.fAmp), fIndexInList(rp.fIndexInList),
    fGlobPos(rp.fGlobPos),fLocPos(rp.fLocPos),
    fMaxDigit(rp.fMaxDigit), fMulDigit(rp.fMulDigit),
    fMaxTrack(rp.fMaxTrack), fMulTrack(rp.fMaxTrack),
    fDigitsList(new Int_t[rp.fMaxDigit]), fTracksList(new Int_t[rp.fMaxTrack]),
    fClusterType(rp.fClusterType), fCoreEnergy(rp.fCoreEnergy), 
    fDispersion(rp.fDispersion),
    fEnergyList(new Float_t[rp.fMaxDigit]),  
    fAbsIdList(new Int_t[rp.fMaxDigit]), fTime(rp.fTime), fNExMax(rp.fNExMax),fCoreRadius(rp.fCoreRadius),
    fDETracksList(new Float_t[rp.fMaxTrack]), fMulParent(rp.fMulParent), 
    fMaxParent(rp.fMaxParent), fParentsList(new Int_t[rp.fMaxParent]), 
    fDEParentsList(new Float_t[rp.fMaxParent]),
    fSuperModuleNumber(rp.fSuperModuleNumber), fDigitIndMax(rp.fDigitIndMax), 
    fDistToBadTower(rp.fDistToBadTower), fSharedCluster(rp.fSharedCluster)
{
  //copy ctor
  fLambda[0] = rp.fLambda[0];
  fLambda[1] = rp.fLambda[1];

  for(Int_t i = 0; i < rp.fMulDigit; i++) {
    fEnergyList[i] = rp.fEnergyList[i];
    fAbsIdList[i]  = rp.fAbsIdList[i];
  }

  for(Int_t i = 0; i < rp.fMulTrack; i++) fDETracksList[i] = rp.fDETracksList[i];

  for(Int_t i = 0; i < rp.fMulParent; i++) {
    fParentsList[i] = rp.fParentsList[i];
    fDEParentsList[i] = rp.fDEParentsList[i];
  }

}
//____________________________________________________________________________
AliEMCALRecPoint::~AliEMCALRecPoint()
{
  // dtor
  if ( fEnergyList )
    delete[] fEnergyList ; 
  if ( fAbsIdList )
    delete[] fAbsIdList ; 
   if ( fDETracksList)
    delete[] fDETracksList;
   if ( fParentsList)
    delete[] fParentsList;
   if ( fDEParentsList)
    delete[] fDEParentsList;
	
   delete [] fDigitsList ;
   delete [] fTracksList ;
}

//____________________________________________________________________________
AliEMCALRecPoint& AliEMCALRecPoint::operator= (const AliEMCALRecPoint &rp)
{
  // assignment operator

  if(&rp == this) return *this;

  fGeomPtr = rp.fGeomPtr;
  fAmp = rp.fAmp;
  fIndexInList = rp.fIndexInList;
  fGlobPos = rp.fGlobPos;
  fLocPos  = rp.fLocPos;
  fMaxDigit = rp.fMaxDigit;
  fMulDigit = rp.fMulDigit;
  fMaxTrack = rp.fMaxTrack;
  fMulTrack = rp.fMaxTrack;
  for(Int_t i = 0; i<fMaxDigit; i++) fDigitsList[i] = rp.fDigitsList[i];
  for(Int_t i = 0; i<fMaxTrack; i++) fTracksList[i] = rp.fTracksList[i];
  fClusterType = rp.fClusterType;
  fCoreEnergy  = rp.fCoreEnergy; 
  fDispersion  = rp.fDispersion;
  for(Int_t i = 0; i<fMaxDigit; i++) {
    fEnergyList[i] = rp.fEnergyList[i];
    fAbsIdList[i]  = rp.fAbsIdList[i];
  }
  fTime = rp.fTime;
  fNExMax = rp.fNExMax;
  fCoreRadius = rp.fCoreRadius;
  for(Int_t i = 0; i < fMaxTrack; i++) fDETracksList[i] = rp.fDETracksList[i];
  fMulParent = rp.fMulParent;
  fMaxParent = rp.fMaxParent;
  for(Int_t i = 0; i < fMaxParent; i++) {
    fParentsList[i] = rp.fParentsList[i]; 
    fDEParentsList[i] = rp.fDEParentsList[i];
  }
  fSuperModuleNumber = rp.fSuperModuleNumber;
  fDigitIndMax = rp.fDigitIndMax;

  fLambda[0] = rp.fLambda[0];
  fLambda[1] = rp.fLambda[1];
	
  fDistToBadTower = rp.fDistToBadTower;
  fSharedCluster  = rp.fSharedCluster;
	
  return *this;

}

//____________________________________________________________________________
void AliEMCALRecPoint::AddDigit(AliEMCALDigit & digit, const Float_t energy, const Bool_t shared)
{
  // Adds a digit to the RecPoint
  // and accumulates the total amplitude and the multiplicity 
  
  if(fEnergyList == 0)
    fEnergyList =  new Float_t[fMaxDigit]; 
 
  if(fAbsIdList == 0) {
    fAbsIdList  =  new Int_t  [fMaxDigit];
  }

  if ( fMulDigit >= fMaxDigit ) { // increase the size of the lists 
    fMaxDigit*=2 ; 
    Int_t   * tempo   = new Int_t  [fMaxDigit]; 
    Float_t * tempoE  = new Float_t[fMaxDigit];
    Int_t   * tempoId = new Int_t  [fMaxDigit]; 

    Int_t index ;     
    for ( index = 0 ; index < fMulDigit ; index++ ){
      tempo  [index] = fDigitsList[index] ;
      tempoE [index] = fEnergyList[index] ; 
      tempoId[index] = fAbsIdList [index] ; 
    }
    
    delete [] fDigitsList ;
    delete [] fEnergyList ;
    delete [] fAbsIdList ;

    fDigitsList = tempo;
    fEnergyList = tempoE; 
    fAbsIdList  = tempoId;
  } // if
  
  fDigitsList[fMulDigit]   = digit.GetIndexInList()  ; 
  fEnergyList[fMulDigit]   = energy ;
  fAbsIdList [fMulDigit]   = digit.GetId();
  fMulDigit++ ; 
  fAmp += energy ; 
	
  if(shared) fSharedCluster = kTRUE;
}
//____________________________________________________________________________
Bool_t AliEMCALRecPoint::AreNeighbours(AliEMCALDigit * digit1, AliEMCALDigit * digit2 ) const
{
  // Tells if (true) or not (false) two digits are neighbours
  // A neighbour is defined as being two digits which share a corner
  // ONLY USED IN CASE OF UNFOLDING 
	
  Bool_t areNeighbours = kFALSE ;
  Int_t nSupMod=0, nModule=0, nIphi=0, nIeta=0;
  Int_t nSupMod1=0, nModule1=0, nIphi1=0, nIeta1=0;
  Int_t relid1[2] , relid2[2] ; // ieta, iphi
  Int_t rowdiff=0, coldiff=0;

  areNeighbours = kFALSE ;

  fGeomPtr->GetCellIndex(digit1->GetId(), nSupMod,nModule,nIphi,nIeta);
  fGeomPtr->GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, relid1[0],relid1[1]);

  fGeomPtr->GetCellIndex(digit2->GetId(), nSupMod1,nModule1,nIphi1,nIeta1);
  fGeomPtr->GetCellPhiEtaIndexInSModule(nSupMod1,nModule1,nIphi1,nIeta1, relid2[0],relid2[1]);
  
  // In case of a shared cluster, index of SM in C side, columns start at 48 and ends at 48*2-1
  // C Side impair SM, nSupMod%2=1; A side pair SM nSupMod%2=0
  if(fSharedCluster){
    if(nSupMod1%2) relid1[1]+=AliEMCALGeoParams::fgkEMCALCols;
    else           relid2[1]+=AliEMCALGeoParams::fgkEMCALCols;
  }
	
  rowdiff = TMath::Abs( relid1[0] - relid2[0] ) ;  
  coldiff = TMath::Abs( relid1[1] - relid2[1] ) ;  

  if (( coldiff <= 1 )  && ( rowdiff <= 1 ) && (coldiff + rowdiff > 0)) 
  areNeighbours = kTRUE ;
  
  return areNeighbours;
}

//____________________________________________________________________________
Int_t AliEMCALRecPoint::Compare(const TObject * obj) const
{
  // Compares two RecPoints according to their position in the EMCAL modules

  Float_t delta = 1 ; //Width of "Sorting row". 
	
  Int_t rv = 2 ; 

  AliEMCALRecPoint * clu = (AliEMCALRecPoint *)obj ; 

  TVector3 locpos1; 
  GetLocalPosition(locpos1);
  TVector3 locpos2;  
  clu->GetLocalPosition(locpos2);  

  Int_t rowdif = (Int_t)(TMath::Ceil(locpos1.X()/delta)-TMath::Ceil(locpos2.X()/delta)) ;
  if (rowdif> 0) 
    rv = 1 ;
  else if(rowdif < 0) 
    rv = -1 ;
  else if(locpos1.Y()>locpos2.Y()) 
    rv = -1 ;
  else 
    rv = 1 ; 

  return rv ; 
}

//___________________________________________________________________________
 void AliEMCALRecPoint::Draw(Option_t *option)
 {
   // Draw this AliEMCALRecPoint with its current attributes
   
   AppendPad(option);
 }

//____________________________________________________________________________
void AliEMCALRecPoint::EvalAll(Float_t logWeight,TClonesArray * digits, const Bool_t justClusters) 
{
  // Evaluates cluster parameters
	
  // First calculate the index of digit with maximum amplitude and get 
  // the supermodule number where it sits.
    
  fDigitIndMax       = GetMaximalEnergyIndex();
  fSuperModuleNumber = fGeomPtr->GetSuperModuleNumber(GetAbsIdMaxDigit());
  
  //Evaluate global and local position
  EvalGlobalPosition(logWeight, digits) ;
  EvalLocalPosition(logWeight, digits) ;
	
  //Evaluate shower parameters
  EvalElipsAxis(logWeight, digits) ;
  EvalDispersion(logWeight, digits) ;

  //EvalCoreEnergy(logWeight, digits);
  EvalTime(digits) ;
  EvalPrimaries(digits) ;
  EvalParents(digits);
	
  //Called last because it sets the global position of the cluster?
  //Do not call it when recalculating clusters out of standard reconstruction
  if(!justClusters){ 
    EvalLocal2TrackingCSTransform();
  }

}

//____________________________________________________________________________
void  AliEMCALRecPoint::EvalDispersion(Float_t logWeight, TClonesArray * digits)
{
  // Calculates the dispersion of the shower at the origin of the RecPoint
  // in cell units - Nov 16,2006

  Double_t d = 0., wtot = 0., w = 0.;
  Int_t iDigit=0, nstat=0;
  AliEMCALDigit * digit=0;
	
  // Calculates the dispersion in cell units 
  Double_t etai, phii, etaMean=0.0, phiMean=0.0; 
  int nSupMod=0, nModule=0, nIphi=0, nIeta=0;
  int iphi=0, ieta=0;
  // Calculate mean values
  for(iDigit=0; iDigit < fMulDigit; iDigit++) {
    digit = (AliEMCALDigit *) digits->At(fDigitsList[iDigit])  ;

    if (fAmp>0 && fEnergyList[iDigit]>0) {
      fGeomPtr->GetCellIndex(digit->GetId(), nSupMod,nModule,nIphi,nIeta);
      fGeomPtr->GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, iphi,ieta);
	
      // In case of a shared cluster, index of SM in C side, columns start at 48 and ends at 48*2
      // C Side impair SM, nSupMod%2=1; A side pair SM nSupMod%2=0
      if(fSharedCluster && nSupMod%2) ieta+=AliEMCALGeoParams::fgkEMCALCols;
		
      etai=(Double_t)ieta;
      phii=(Double_t)iphi;
      w = TMath::Max(0.,logWeight+TMath::Log(fEnergyList[iDigit]/fAmp ) ) ;

      if(w>0.0) {
        phiMean += phii*w;
        etaMean += etai*w;
        wtot    += w;
      }
    }
  }
  if (wtot>0) {
    phiMean /= wtot ;
    etaMean /= wtot ;
  } else AliError(Form("Wrong weight %f\n", wtot));

  // Calculate dispersion
  for(iDigit=0; iDigit < fMulDigit; iDigit++) {
    digit = (AliEMCALDigit *) digits->At(fDigitsList[iDigit])  ;

    if (fAmp>0 && fEnergyList[iDigit]>0) {
      fGeomPtr->GetCellIndex(digit->GetId(), nSupMod,nModule,nIphi,nIeta);
      fGeomPtr->GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, iphi,ieta);
		
      // In case of a shared cluster, index of SM in C side, columns start at 48 and ends at 48*2
      // C Side impair SM, nSupMod%2=1; A side pair SM, nSupMod%2=0
      if(fSharedCluster && nSupMod%2) ieta+=AliEMCALGeoParams::fgkEMCALCols;
      
      etai=(Double_t)ieta;
      phii=(Double_t)iphi;
      w = TMath::Max(0.,logWeight+TMath::Log(fEnergyList[iDigit]/fAmp ) ) ;

      if(w>0.0) {
        nstat++;
        d += w*((etai-etaMean)*(etai-etaMean)+(phii-phiMean)*(phii-phiMean));
      }
    }
  }
  
  if ( wtot > 0 && nstat>1) d /= wtot ;
  else                      d = 0. ; 

  fDispersion = TMath::Sqrt(d) ;
  //printf("AliEMCALRecPoint::EvalDispersion() : Dispersion %f \n",fDispersion);
}

//____________________________________________________________________________
void AliEMCALRecPoint::EvalDistanceToBadChannels(AliCaloCalibPedestal* caloped)
{
  //For each EMC rec. point set the distance to the nearest bad channel.
  //AliInfo(Form("%d bad channel(s) found.\n", caloped->GetDeadTowerCount()));
  //It is done in cell units and not in global or local position as before (Sept 2010)
  
  if(!caloped->GetDeadTowerCount()) return;
  
  //Get channels map of the supermodule where the cluster is.
  TH2D* hMap  = caloped->GetDeadMap(fSuperModuleNumber);
  
  Int_t dRrow, dReta;	
  Float_t  minDist = 10000.;
  Float_t  dist    = 0.;
  Int_t nSupMod, nModule;
  Int_t nIphi, nIeta;
  Int_t iphi, ieta;
  fDigitIndMax  = GetMaximalEnergyIndex();
  fGeomPtr->GetCellIndex(fAbsIdList[fDigitIndMax], nSupMod,nModule,nIphi,nIeta);
  fGeomPtr->GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, iphi,ieta);
  
  //Loop on tower status map 
  for(Int_t irow = 0; irow < AliEMCALGeoParams::fgkEMCALRows; irow++){
    for(Int_t icol = 0; icol < AliEMCALGeoParams::fgkEMCALCols; icol++){
      //Check if tower is bad.
      if(hMap->GetBinContent(icol,irow)==AliCaloCalibPedestal::kAlive) continue;
      //printf("AliEMCALRecPoint::EvalDistanceToBadChannels() - Bad channel in SM %d, col %d, row %d\n",iSM,icol, irow);
      
      dRrow=TMath::Abs(irow-iphi);
      dReta=TMath::Abs(icol-ieta);
      dist=TMath::Sqrt(dRrow*dRrow+dReta*dReta);
      if(dist < minDist) minDist = dist;
      
    }
  }
  
  //In case the cluster is shared by 2 SuperModules, need to check the map of the second Super Module
  if (fSharedCluster) {
    TH2D* hMap2 = 0;
    Int_t nSupMod2 = -1;
    
    //The only possible combinations are (0,1), (2,3) ... (10,11)
    if(fSuperModuleNumber%2) nSupMod2 = fSuperModuleNumber-1;
    else                     nSupMod2 = fSuperModuleNumber+1;
    hMap2  = caloped->GetDeadMap(nSupMod2);
    
    //Loop on tower status map of second super module
    for(Int_t irow = 0; irow < AliEMCALGeoParams::fgkEMCALRows; irow++){
      for(Int_t icol = 0; icol < AliEMCALGeoParams::fgkEMCALCols; icol++){
	//Check if tower is bad.
	if(hMap2->GetBinContent(icol,irow)==AliCaloCalibPedestal::kAlive) continue;
	//printf("AliEMCALRecPoint::EvalDistanceToBadChannels() - Bad channel in SM %d, col %d, row %d\n",iSM,icol, irow);
        dRrow=TMath::Abs(irow-iphi);
	
        if(fSuperModuleNumber%2) {
	  dReta=TMath::Abs(icol-(AliEMCALGeoParams::fgkEMCALCols+ieta));
       }
       else {
	 dReta=TMath::Abs(AliEMCALGeoParams::fgkEMCALCols+icol-ieta);
       }                    
        
       dist=TMath::Sqrt(dRrow*dRrow+dReta*dReta);
        if(dist < minDist) minDist = dist;        

      }
    }
    
  }// shared cluster in 2 SuperModules
  
  fDistToBadTower = minDist;
  //printf("AliEMCALRecPoint::EvalDistanceToBadChannel() - Distance to Bad is %f cm, shared cluster? %d \n",fDistToBadTower,fSharedCluster);
}


//____________________________________________________________________________
void AliEMCALRecPoint::EvalLocalPosition(Float_t logWeight, TClonesArray * digits)
{
  // Calculates the center of gravity in the local EMCAL-module coordinates 
  //  Info("Print", " logWeight %f : cluster energy %f ", logWeight, fAmp); // for testing
  
  AliEMCALDigit * digit=0;
  Int_t i=0, nstat=0;
  
  Double_t dist  = TmaxInCm(Double_t(fAmp));
  //Int_t	idMax = GetAbsIdMaxDigit();// idMax is not used at all in RelPosCellInSModule, why use it?
  
  Double_t clXYZ[3]={0.,0.,0.}, clRmsXYZ[3]={0.,0.,0.}, xyzi[3], wtot=0., w=0.;
  
  //printf(" dist : %f  e : %f \n", dist, fAmp);
  for(Int_t iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = dynamic_cast<AliEMCALDigit *>(digits->At(fDigitsList[iDigit])) ;
    
    if(!digit) {
      AliError("No Digit!!");
      continue;
    }
    
    fGeomPtr->RelPosCellInSModule(digit->GetId(), dist, xyzi[0], xyzi[1], xyzi[2]);
    
    //Temporal patch, due to mapping problem, need to swap "y" in one of the 2 SM, although no effect in position calculation. GCB 05/2010
    if(fSharedCluster && fSuperModuleNumber != fGeomPtr->GetSuperModuleNumber(digit->GetId())) xyzi[1]*=-1;
    
    //printf("EvalLocalPosition Cell:  Id %i, SM %i : dist %f Local x,y,z %f %f %f \n", 
    //		digit->GetId(), fGeomPtr->GetSuperModuleNumber(digit->GetId()), dist, xyzi[0], xyzi[1], xyzi[2]);
    
    if(logWeight > 0.0)  w = TMath::Max( 0., logWeight + TMath::Log( fEnergyList[iDigit] / fAmp ));
    else  w = fEnergyList[iDigit]; // just energy
    
    if(w>0.0) {
      wtot += w ;
      nstat++;
      for(i=0; i<3; i++ ) {
	clXYZ[i]    += (w*xyzi[i]);
	clRmsXYZ[i] += (w*xyzi[i]*xyzi[i]);
      }
    }
  }
  //  cout << " wtot " << wtot << endl;
  if ( wtot > 0 ) { 
    //    xRMS   = TMath::Sqrt(x2m - xMean*xMean);
    for(i=0; i<3; i++ ) {
      clXYZ[i] /= wtot;
      if(nstat>1) {
	clRmsXYZ[i] /= (wtot*wtot);  
	clRmsXYZ[i]  =  clRmsXYZ[i] - clXYZ[i]*clXYZ[i];
	if(clRmsXYZ[i] > 0.0) {
	  clRmsXYZ[i] =  TMath::Sqrt(clRmsXYZ[i]);
	} else      clRmsXYZ[i] = 0;
      } else        clRmsXYZ[i] = 0;
    }    
  } else {
    for(i=0; i<3; i++ ) {
      clXYZ[i] = clRmsXYZ[i] = -1.;
    }
  }
  
  //	// Cluster of one single digit, smear the position to avoid discrete position
  //	// smear x and z with +- 3 cm to uniform (avoid discrete effects). Tower size is approx 6 cm.
  //	// Rndm generates a number in ]0,1]
  //	if (fMulDigit==1) { 
  //	  clXYZ[0] += fGeomPtr->GetPhiTileSize()*(0.5 - gRandom->Rndm()); 
  //	  clXYZ[2] += fGeomPtr->GetEtaTileSize()*(0.5 - gRandom->Rndm()); 
  //	}
  
  //Set position in local vector
  fLocPos.SetX(clXYZ[0]);
  fLocPos.SetY(clXYZ[1]);
  fLocPos.SetZ(clXYZ[2]);
  
  if (gDebug==2)
    printf("EvalLocalPosition Cluster: Local (x,y,z) = (%f,%f,%f) \n", fLocPos.X(), fLocPos.Y(), fLocPos.Z()) ; 
  
}


//____________________________________________________________________________
void AliEMCALRecPoint::EvalGlobalPosition(Float_t logWeight, TClonesArray * digits)
{
  // Calculates the center of gravity in the global ALICE coordinates 
  //  Info("Print", " logWeight %f : cluster energy %f ", logWeight, fAmp); // for testing
  
  AliEMCALDigit * digit=0;
  Int_t i=0, nstat=0;
	
  Double_t dist  = TmaxInCm(Double_t(fAmp));
  //Int_t	idMax = GetAbsIdMaxDigit();// idMax is not used at all in RelPosCellInSModule, why use it?
	
  Double_t clXYZ[3]={0.,0.,0.}, clRmsXYZ[3]={0.,0.,0.}, lxyzi[3], xyzi[3], wtot=0., w=0.;

  //printf(" dist : %f  e : %f \n", dist, fAmp);
  for(Int_t iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = dynamic_cast<AliEMCALDigit *>(digits->At(fDigitsList[iDigit])) ;

    if(!digit) {
      AliError("No Digit!!");
      continue;
    }    
    
    //Get the local coordinates of the cell
    fGeomPtr->RelPosCellInSModule(digit->GetId(), dist, lxyzi[0], lxyzi[1], lxyzi[2]);
    
    //Now get the global coordinate
    fGeomPtr->GetGlobal(lxyzi,xyzi, fGeomPtr->GetSuperModuleNumber(digit->GetId()));
    //TVector3 pos(xyzi[0], xyzi[1], xyzi[2]);
    //printf("EvalGlobalPosition Cell:  Id %i, SM %i : dist %f Local (x,y,z) = (%f %f %f), eta %f, phi%f \n", 
    //	   digit->GetId(), fGeomPtr->GetSuperModuleNumber(digit->GetId()),dist, xyzi[0], xyzi[1], xyzi[2],pos.Eta(),pos.Phi()*TMath::RadToDeg());
	  
    if(logWeight > 0.0)  w = TMath::Max( 0., logWeight + TMath::Log( fEnergyList[iDigit] / fAmp ));
    else  w = fEnergyList[iDigit]; // just energy

    if(w>0.0) {
      wtot += w ;
      nstat++;
      for(i=0; i<3; i++ ) {
        clXYZ[i]    += (w*xyzi[i]);
        clRmsXYZ[i] += (w*xyzi[i]*xyzi[i]);
      }
    }
  }
  //  cout << " wtot " << wtot << endl;
  if ( wtot > 0 ) { 
    //    xRMS   = TMath::Sqrt(x2m - xMean*xMean);
    for(i=0; i<3; i++ ) {
      clXYZ[i] /= wtot;
      if(nstat>1) {
        clRmsXYZ[i] /= (wtot*wtot);  
        clRmsXYZ[i]  =  clRmsXYZ[i] - clXYZ[i]*clXYZ[i];
	if(clRmsXYZ[i] > 0.0) {
          clRmsXYZ[i] =  TMath::Sqrt(clRmsXYZ[i]);
        } else      clRmsXYZ[i] = 0;
      } else        clRmsXYZ[i] = 0;
    }    
  } else {
    for(i=0; i<3; i++ ) {
      clXYZ[i] = clRmsXYZ[i] = -1.;
    }
  }

//  // Cluster of one single digit, smear the position to avoid discrete position
//  // smear x and z with +- 3 cm to uniform (avoid discrete effects). Tower size is approx 6 cm.
//  // Rndm generates a number in ]0,1]
//  if (fMulDigit==1) { 
//    clXYZ[0] += fGeomPtr->GetPhiTileSize()*(0.5 - gRandom->Rndm()); 
//    clXYZ[2] += fGeomPtr->GetEtaTileSize()*(0.5 - gRandom->Rndm()); 	
//  }
	
  //Set position in global vector
  fGlobPos.SetX(clXYZ[0]);
  fGlobPos.SetY(clXYZ[1]);
  fGlobPos.SetZ(clXYZ[2]);
		
  if (gDebug==2)
	printf("EvalGlobalPosition Cluster: (x ,y ,z) = (%f,%f,%f), eta %f,phi %f\n", 
		   fGlobPos.X(), fGlobPos.Y(), fGlobPos.Z(),fGlobPos.Eta(),fGlobPos.Phi()*TMath::RadToDeg()) ; 
}

//____________________________________________________________________________
void AliEMCALRecPoint::EvalLocalPositionFit(Double_t deff, Double_t logWeight, 
Double_t phiSlope, TClonesArray * digits)
{
  // Evaluates local position of clusters in SM
  
  Double_t ycorr=0;
  AliEMCALDigit *digit=0;
  Int_t i=0, nstat=0;
  Double_t clXYZ[3]={0.,0.,0.}, clRmsXYZ[3]={0.,0.,0.}, xyzi[3], wtot=0., w=0.; 

  Double_t dist  = TmaxInCm(Double_t(fAmp));
  //Int_t	idMax = GetAbsIdMaxDigit();// idMax is not used at all in RelPosCellInSModule, why use it?
	
  for(Int_t iDigit=0; iDigit<digits->GetEntries(); iDigit++) {
    digit = dynamic_cast<AliEMCALDigit *>(digits->At(fDigitsList[iDigit])) ;
    if(digit){
      dist = deff;
      //fGeomPtr->RelPosCellInSModule(digit->GetId(), idMax, dist, xyzi[0], xyzi[1], xyzi[2]);
      fGeomPtr->RelPosCellInSModule(digit->GetId(), dist, xyzi[0], xyzi[1], xyzi[2]);
      
      if(logWeight > 0.0)  w = TMath::Max( 0., logWeight + TMath::Log( fEnergyList[iDigit] / fAmp ));
      else                 w = fEnergyList[iDigit]; // just energy
      
      if(w>0.0) {
        wtot += w ;
        nstat++;
        for(i=0; i<3; i++ ) {
          clXYZ[i]    += (w*xyzi[i]);
          clRmsXYZ[i] += (w*xyzi[i]*xyzi[i]);
        }
      }
    }else AliError("Digit null");
  }//loop
  //  cout << " wtot " << wtot << endl;
  if ( wtot > 0 ) { 
    //    xRMS   = TMath::Sqrt(x2m - xMean*xMean);
    for(i=0; i<3; i++ ) {
      clXYZ[i] /= wtot;
      if(nstat>1) {
        clRmsXYZ[i] /= (wtot*wtot);  
        clRmsXYZ[i]  =  clRmsXYZ[i] - clXYZ[i]*clXYZ[i];
	if(clRmsXYZ[i] > 0.0) {
          clRmsXYZ[i] =  TMath::Sqrt(clRmsXYZ[i]);
        } else      clRmsXYZ[i] = 0;
      } else        clRmsXYZ[i] = 0;
    }    
  } else {
    for(i=0; i<3; i++ ) {
      clXYZ[i] = clRmsXYZ[i] = -1.;
    }
  }
  // clRmsXYZ[i] ??
  if(phiSlope != 0.0 && logWeight > 0.0 && wtot) { 
    // Correction in phi direction (y - coords here); Aug 16;
    // May be put to global level or seperate method
    ycorr = clXYZ[1] * (1. + phiSlope);
    //printf(" y %f : ycorr %f : slope %f \n", clXYZ[1], ycorr, phiSlope); 
    clXYZ[1] = ycorr;
  }
	
  fLocPos.SetX(clXYZ[0]);
  fLocPos.SetY(clXYZ[1]);
  fLocPos.SetZ(clXYZ[2]);
    
//  if (gDebug==2)
//    printf("EvalLocalPosition: eta,phi,r = %f,%f,%f", fLocPos.X(), fLocPos.Y(), fLocPos.Z()) ; 
}

//_____________________________________________________________________________
Bool_t AliEMCALRecPoint::EvalLocalPosition2(TClonesArray * digits, TArrayD &ed)
{
  // Evaluated local position of rec.point using digits 
  // and parametrisation of w0 and deff
  //printf(" <I> AliEMCALRecPoint::EvalLocalPosition2() \n"); 
  return AliEMCALRecPoint::EvalLocalPositionFromDigits(digits, ed, fLocPos);
}

//_____________________________________________________________________________
Bool_t AliEMCALRecPoint::EvalLocalPositionFromDigits(TClonesArray *digits, TArrayD &ed, TVector3 &locPos)
{
  // Used when digits should be recalibrated
  Double_t deff=0, w0=0, esum=0;
  Int_t iDigit=0;
  //  AliEMCALDigit *digit;

  if(ed.GetSize() && (digits->GetEntries()!=ed.GetSize())) return kFALSE;

  // Calculate sum energy of digits
  esum = 0.0;
  for(iDigit=0; iDigit<ed.GetSize(); iDigit++) esum += ed[iDigit];

  GetDeffW0(esum, deff, w0);
  
  return EvalLocalPositionFromDigits(esum, deff, w0, digits, ed, locPos); 
}

//_____________________________________________________________________________
Bool_t AliEMCALRecPoint::EvalLocalPositionFromDigits(const Double_t esum, const Double_t deff, const Double_t w0, TClonesArray *digits, TArrayD &ed, TVector3 &locPos)
{
  //Evaluate position of digits in supermodule.
  AliEMCALDigit *digit=0;

  Int_t i=0, nstat=0;
  Double_t clXYZ[3]={0.,0.,0.}, xyzi[3], wtot=0., w=0.; 
  //Int_t	idMax = GetAbsIdMaxDigit();// idMax is not used at all in RelPosCellInSModule, why use it?
	
  // Get pointer to EMCAL geometry
  // (can't use fGeomPtr in static method)
  AliEMCALGeometry* geo = AliEMCALGeometry::GetInstance(); 

  for(Int_t iDigit=0; iDigit<digits->GetEntries(); iDigit++) {
    digit = dynamic_cast<AliEMCALDigit *>(digits->At(iDigit));
    if(digit){
      //geo->RelPosCellInSModule(digit->GetId(), idMax, deff, xyzi[0], xyzi[1], xyzi[2]);
      geo->RelPosCellInSModule(digit->GetId(), deff, xyzi[0], xyzi[1], xyzi[2]);
      
      if(w0 > 0.0)  w = TMath::Max( 0., w0 + TMath::Log(ed[iDigit] / esum));
      else          w = ed[iDigit]; // just energy
      
      if(w>0.0) {
        wtot += w ;
        nstat++;
        for(i=0; i<3; i++ ) {
          clXYZ[i] += (w*xyzi[i]);
        }
      }
    }else AliError("Digit null");
  }//loop
  //  cout << " wtot " << wtot << endl;
  if (wtot > 0) { 
    for(i=0; i<3; i++ ) {
      clXYZ[i] /= wtot;
    }
    locPos.SetX(clXYZ[0]);
    locPos.SetY(clXYZ[1]);
    locPos.SetZ(clXYZ[2]);
    return kTRUE;
  } else {
    return kFALSE;
  }

}

//_____________________________________________________________________________
void AliEMCALRecPoint::GetDeffW0(const Double_t esum , Double_t &deff,  Double_t &w0)
{
  //
  // Aug 31, 2001 
  // Applied for simulation data with threshold 3 adc
  // Calculate efective distance (deff) and weigh parameter (w0) 
  // for coordinate calculation; 0.5 GeV < esum <100 GeV.
  // Look to:  http://rhic.physics.wayne.edu/~pavlinov/ALICE/SHISHKEBAB/RES/CALIB/GEOMCORR/deffandW0VaEgamma_2.gif
  //
  Double_t e=0.0;
  const  Double_t kdp0=9.25147, kdp1=1.16700; // Hard coded now
  const  Double_t kwp0=4.83713, kwp1=-2.77970e-01, kwp2 = 4.41116;

  // No extrapolation here
  e = esum<0.5?0.5:esum;
  e = e>100.?100.:e;

  deff = kdp0 + kdp1*TMath::Log(e);
  w0   = kwp0 / (1. + TMath::Exp(kwp1*(e+kwp2)));
  //printf("<I> AliEMCALRecPoint::GetDeffW0 esum %5.2f : deff %5.2f : w0 %5.2f \n", esum, deff, w0); 
}

//______________________________________________________________________________
void AliEMCALRecPoint::EvalCoreEnergy(Float_t logWeight, TClonesArray * digits)
{
  // This function calculates energy in the core, 
  // i.e. within a radius rad = fCoreEnergy around the center. Beyond this radius
  // in accordance with shower profile the energy deposition 
  // should be less than 2%
  // Unfinished - Nov 15,2006
  // Distance is calculate in (phi,eta) units

  AliEMCALDigit * digit = 0 ;

  Int_t iDigit=0;

  if (!fLocPos.Mag()) {
    EvalLocalPosition(logWeight, digits);
  }
  
  Double_t phiPoint = fLocPos.Phi(), etaPoint = fLocPos.Eta();
  Double_t eta, phi, distance;
  for(iDigit=0; iDigit < fMulDigit; iDigit++) {
    digit = (AliEMCALDigit *) ( digits->At(fDigitsList[iDigit]) ) ;
    
    eta = phi = 0.0;
    fGeomPtr->EtaPhiFromIndex(digit->GetId(),eta, phi) ;
    phi = phi * TMath::DegToRad();
  
    distance = TMath::Sqrt((eta-etaPoint)*(eta-etaPoint)+(phi-phiPoint)*(phi-phiPoint));
    if(distance < fCoreRadius)
      fCoreEnergy += fEnergyList[iDigit] ;
  }
  
}
//____________________________________________________________________________
void  AliEMCALRecPoint::EvalElipsAxis(Float_t logWeight,TClonesArray * digits)
{
  // Calculates the axis of the shower ellipsoid in eta and phi
  // in cell units

  TString gn(fGeomPtr->GetName());

  Double_t wtot = 0.;
  Double_t x    = 0.;
  Double_t z    = 0.;
  Double_t dxx  = 0.;
  Double_t dzz  = 0.;
  Double_t dxz  = 0.;

  AliEMCALDigit * digit = 0;
	
  Double_t etai =0, phii=0, w=0; 
  int nSupMod=0, nModule=0, nIphi=0, nIeta=0;
  int iphi=0, ieta=0;
  for(Int_t iDigit=0; iDigit<fMulDigit; iDigit++) {
    digit = (AliEMCALDigit *) digits->At(fDigitsList[iDigit])  ;
    etai = phii = 0.; 
    // Nov 15,2006 - use cell numbers as coordinates
    // Copied for shish-kebab geometry, ieta,iphi is cast as double as eta,phi
    // We can use the eta,phi(or coordinates) of cell
    nSupMod = nModule = nIphi = nIeta = iphi = ieta = 0;

    fGeomPtr->GetCellIndex(digit->GetId(), nSupMod,nModule,nIphi,nIeta);
    fGeomPtr->GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, iphi,ieta);
	  
    // In case of a shared cluster, index of SM in C side, columns start at 48 and ends at 48*2
    // C Side impair SM, nSupMod%2=1; A side pair SM, nSupMod%2=0
    if(fSharedCluster && nSupMod%2) ieta+=AliEMCALGeoParams::fgkEMCALCols;
    
    etai=(Double_t)ieta;
    phii=(Double_t)iphi;
	  
    w = TMath::Max(0.,logWeight+TMath::Log(fEnergyList[iDigit]/fAmp ) ) ;
    // fAmp summed amplitude of digits, i.e. energy of recpoint 
    // Gives smaller value of lambda than log weight  
    // w = fEnergyList[iDigit] / fAmp; // Nov 16, 2006 - try just energy

    dxx  += w * etai * etai ;
    x    += w * etai ;
    dzz  += w * phii * phii ;
    z    += w * phii ; 

    dxz  += w * etai * phii ; 

    wtot += w ;
  }

  if ( wtot > 0 ) { 
    dxx /= wtot ;
    x   /= wtot ;
    dxx -= x * x ;
    dzz /= wtot ;
    z   /= wtot ;
    dzz -= z * z ;
    dxz /= wtot ;
    dxz -= x * z ;

    fLambda[0] =  0.5 * (dxx + dzz) + TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )  ;
    if(fLambda[0] > 0)
      fLambda[0] = TMath::Sqrt(fLambda[0]) ;
    else
      fLambda[0] = 0;
    
    fLambda[1] =  0.5 * (dxx + dzz) - TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )  ;

    if(fLambda[1] > 0) //To avoid exception if numerical errors lead to negative lambda.
      fLambda[1] = TMath::Sqrt(fLambda[1]) ;
    else
      fLambda[1]= 0. ;
  } else { 
    fLambda[0]= 0. ;
    fLambda[1]= 0. ;
  }

  //printf("AliEMCALRecPoint::EvalElipsAxis() lambdas  = %f,%f \n", fLambda[0],fLambda[1]) ; 

}

//______________________________________________________________________________
void  AliEMCALRecPoint::EvalPrimaries(TClonesArray * digits)
{
  // Constructs the list of primary particles (tracks) which 
  // have contributed to this RecPoint and calculate deposited energy 
  // for each track
    
  AliEMCALDigit * digit =0;
  Int_t * primArray = new Int_t[fMaxTrack] ;
  memset(primArray,-1,sizeof(Int_t)*fMaxTrack);
  Float_t * dEPrimArray = new Float_t[fMaxTrack] ;
  memset(dEPrimArray,-1,sizeof(Int_t)*fMaxTrack);
  
  Int_t index ;  
  for ( index = 0 ; index < GetDigitsMultiplicity() ; index++ ) { // all digits
    digit = dynamic_cast<AliEMCALDigit *>(digits->At( fDigitsList[index] )) ; 
    if(!digit) {
      AliError("No Digit!!");
      continue;
    }
    
    Int_t nprimaries = digit->GetNprimary() ;
    if ( nprimaries == 0 ) continue ;
    Int_t jndex ;
    for ( jndex = 0 ; jndex < nprimaries ; jndex++ ) { // all primaries in digit
      if ( fMulTrack > fMaxTrack ) {
        fMulTrack = fMaxTrack ;
        Error("EvalPrimaries", "increase fMaxTrack ")  ;
        break ;
      }
      Int_t newPrimary = digit->GetPrimary(jndex+1);
      Float_t dEPrimary = digit->GetDEPrimary(jndex+1);
      Int_t kndex ;
      Bool_t already = kFALSE ;
      for ( kndex = 0 ; kndex < fMulTrack ; kndex++ ) { //check if not already stored
        if ( newPrimary == primArray[kndex] ){
          already = kTRUE ;
          dEPrimArray[kndex] += dEPrimary; 
          break ;
        }
      } // end of check
      if ( !already && (fMulTrack < fMaxTrack)) { // store it
        primArray[fMulTrack] = newPrimary ; 
        dEPrimArray[fMulTrack] = dEPrimary ; 
        fMulTrack++ ;
      } // store it
    } // all primaries in digit
  } // all digits
  
  Int_t *sortIdx = new Int_t[fMulTrack];
  TMath::Sort(fMulTrack,dEPrimArray,sortIdx); 
  for(index = 0; index < fMulTrack; index++) {
    fTracksList[index] = primArray[sortIdx[index]] ;    
    fDETracksList[index] = dEPrimArray[sortIdx[index]] ;
  }
  delete [] sortIdx;
  delete [] primArray ;
  delete [] dEPrimArray ;
  
}

//______________________________________________________________________________
void  AliEMCALRecPoint::EvalParents(TClonesArray * digits)
{
  // Constructs the list of parent particles (tracks) which have contributed to this RecPoint

  AliEMCALDigit * digit=0 ;
  Int_t * parentArray = new Int_t[fMaxTrack] ;
  memset(parentArray,-1,sizeof(Int_t)*fMaxTrack);
  Float_t * dEParentArray = new Float_t[fMaxTrack] ;
  memset(dEParentArray,-1,sizeof(Int_t)*fMaxTrack);
  
  Int_t index ;  
  for ( index = 0 ; index < GetDigitsMultiplicity() ; index++ ) { // all digits
    if (fDigitsList[index] >= digits->GetEntries() || fDigitsList[index] < 0)
      AliError(Form("Trying to get invalid digit %d (idx in WriteRecPoint %d)",fDigitsList[index],index));
    digit = dynamic_cast<AliEMCALDigit *>(digits->At( fDigitsList[index] )) ; 
    if(!digit) {
      AliError("No Digit!!");
      continue;
    }
    
    Int_t nparents = digit->GetNiparent() ;
    if ( nparents == 0 ) continue ;
    
    Int_t jndex ;
    for ( jndex = 0 ; jndex < nparents ; jndex++ ) { // all primaries in digit
      if ( fMulParent > fMaxParent ) {
        fMulTrack = - 1 ;
        Error("EvalParents", "increase fMaxParent")  ;
        break ;
      }
      Int_t newParent = digit->GetIparent(jndex+1) ;
      Float_t newdEParent = digit->GetDEParent(jndex+1) ;
      Int_t kndex ;
      Bool_t already = kFALSE ;
      for ( kndex = 0 ; kndex < fMulParent ; kndex++ ) { //check if not already stored
        if ( newParent == parentArray[kndex] ){
          dEParentArray[kndex] += newdEParent;
          already = kTRUE ;
          break ;
        }
      } // end of check
      if ( !already && (fMulParent < fMaxParent)) { // store it
        parentArray[fMulParent] = newParent ; 
        dEParentArray[fMulParent] = newdEParent ; 
        fMulParent++ ;
      } // store it
    } // all parents in digit
  } // all digits
  
  if (fMulParent>0) {
    Int_t *sortIdx = new Int_t[fMulParent];
    TMath::Sort(fMulParent,dEParentArray,sortIdx); 
    for(index = 0; index < fMulParent; index++) {
      fParentsList[index] = parentArray[sortIdx[index]] ;      
      fDEParentsList[index] = dEParentArray[sortIdx[index]] ;
    }
    delete [] sortIdx;
  }
  
  delete [] parentArray;
  delete [] dEParentArray;
}

//____________________________________________________________________________
void AliEMCALRecPoint::GetLocalPosition(TVector3 & lpos) const
{
  // returns the position of the cluster in the local reference system
  // of the sub-detector
  
  lpos = fLocPos;
}

//____________________________________________________________________________
void AliEMCALRecPoint::GetGlobalPosition(TVector3 & gpos) const
{
  // returns the position of the cluster in the global reference system of ALICE
  // These are now the Cartesian X, Y and Z
  //  cout<<" geom "<<geom<<endl;
  // fGeomPtr->GetGlobal(fLocPos, gpos, fSuperModuleNumber);
  gpos = fGlobPos;
	
}

//____________________________________________________________________________
//void AliEMCALRecPoint::GetGlobalPosition(TVector3 & gpos, TMatrixF & gmat) const
//{
//  // returns the position of the cluster in the global reference system of ALICE
//  // These are now the Cartesian X, Y and Z
//  //  cout<<" geom "<<geom<<endl;
//
//  //To be implemented
//  fGeomPtr->GetGlobalEMCAL(this, gpos, gmat);
//
//}

//_____________________________________________________________________________
void AliEMCALRecPoint::EvalLocal2TrackingCSTransform()
{
  //Evaluates local to "tracking" c.s. transformation (B.P.).
  //All evaluations should be completed before calling for this
  //function.                           
  //See ALICE PPR Chapter 5 p.18 for "tracking" c.s. definition,
  //or just ask Jouri Belikov. :) 

  SetVolumeId(AliGeomManager::LayerToVolUID(AliGeomManager::kEMCAL,GetSuperModuleNumber()));

  const TGeoHMatrix* tr2loc = GetTracking2LocalMatrix();
  if(!tr2loc) AliFatal(Form("No Tracking2LocalMatrix found."));

  Double_t lxyz[3] = {fLocPos.X(),fLocPos.Y(),fLocPos.Z()};
  Double_t txyz[3] = {0,0,0};

  tr2loc->MasterToLocal(lxyz,txyz);
  SetX(txyz[0]); SetY(txyz[1]); SetZ(txyz[2]);

  if(AliLog::GetGlobalDebugLevel()>0) {
    TVector3 gpos; //TMatrixF gmat;
    //GetGlobalPosition(gpos,gmat); //Not doing anythin special, replace by next line.	
    fGeomPtr->GetGlobal(fLocPos, gpos, GetSuperModuleNumber());
    
    Float_t gxyz[3];
    GetGlobalXYZ(gxyz);
    AliInfo(Form("lCS-->(%.3f,%.3f,%.3f), tCS-->(%.3f,%.3f,%.3f), gCS-->(%.3f,%.3f,%.3f),  gCScalc-\
->(%.3f,%.3f,%.3f), supermodule %d",
                 fLocPos.X(),fLocPos.Y(),fLocPos.Z(),
                 GetX(),GetY(),GetZ(),
                 gpos.X(),gpos.Y(),gpos.Z(),
                 gxyz[0],gxyz[1],gxyz[2],GetSuperModuleNumber()));
  }

}

//____________________________________________________________________________
Float_t AliEMCALRecPoint::GetMaximalEnergy(void) const
{
  // Finds the maximum energy in the cluster
  
  Float_t menergy = 0. ;

  Int_t iDigit;
  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
 
    if(fEnergyList[iDigit] > menergy) 
      menergy = fEnergyList[iDigit] ;
  }
  return menergy ;
}

//____________________________________________________________________________
Int_t AliEMCALRecPoint::GetMaximalEnergyIndex(void) const
{
  // Finds the maximum energy in the cluster
  
  Float_t menergy = 0. ;
  Int_t mid       = 0  ;
  Int_t iDigit;
  
  for(iDigit=0; iDigit<fMulDigit; iDigit++) {
    
    if(fEnergyList[iDigit] > menergy){ 
      menergy = fEnergyList[iDigit] ;
      mid     = iDigit ;
    }
  }//loop on cluster digits
  
  return mid ;
}


//____________________________________________________________________________
Int_t AliEMCALRecPoint::GetMultiplicityAtLevel(Float_t H) const
{
  // Calculates the multiplicity of digits with energy larger than H*energy 
  
  Int_t multipl   = 0 ;
  Int_t iDigit ;
  for(iDigit=0; iDigit<fMulDigit; iDigit++) {

    if(fEnergyList[iDigit] > H * fAmp) 
      multipl++ ;
  }
  return multipl ;
}

//____________________________________________________________________________
Int_t  AliEMCALRecPoint::GetNumberOfLocalMax(AliEMCALDigit **  maxAt, Float_t * maxAtEnergy,
					       Float_t locMaxCut,TClonesArray * digits) const
{ 
  // Calculates the number of local maxima in the cluster using fLocalMaxCut as the minimum
  // energy difference between two local maxima

  AliEMCALDigit * digit  = 0;
  AliEMCALDigit * digitN = 0;
  
  Int_t iDigitN = 0 ;
  Int_t iDigit  = 0 ;

  for(iDigit = 0; iDigit < fMulDigit; iDigit++)
    maxAt[iDigit] = (AliEMCALDigit*) digits->At(fDigitsList[iDigit])  ;
  
  for(iDigit = 0 ; iDigit < fMulDigit; iDigit++) {   
    if(maxAt[iDigit]) {
      digit = maxAt[iDigit] ;
          
      for(iDigitN = 0; iDigitN < fMulDigit; iDigitN++) {	
        digitN = (AliEMCALDigit *) digits->At(fDigitsList[iDigitN]) ; 
	
	if ( AreNeighbours(digit, digitN) ) {
	  if (fEnergyList[iDigit] > fEnergyList[iDigitN] ) {    
	    maxAt[iDigitN] = 0 ;
	    // but may be digit too is not local max ?
	    if(fEnergyList[iDigit] < fEnergyList[iDigitN] + locMaxCut) 
	      maxAt[iDigit] = 0 ;
	  }
	  else {
	    maxAt[iDigit] = 0 ;
	    // but may be digitN too is not local max ?
	    if(fEnergyList[iDigit] > fEnergyList[iDigitN] - locMaxCut) 
	      maxAt[iDigitN] = 0 ; 
	  } 
	} // if Areneighbours
      } // while digitN
    } // slot not empty
  } // while digit
  
  iDigitN = 0 ;
  for(iDigit = 0; iDigit < fMulDigit; iDigit++) { 
    if(maxAt[iDigit] ){
      maxAt[iDigitN] = maxAt[iDigit] ;
      maxAtEnergy[iDigitN] = fEnergyList[iDigit] ;
      iDigitN++ ; 
    }
  }
  return iDigitN ;
}

//____________________________________________________________________________
Int_t AliEMCALRecPoint::GetPrimaryIndex() const  
{
  // Get the primary track index in TreeK which deposits the most energy 
  // in Digits which forms RecPoint. 

  if (fMulTrack)
    return fTracksList[0];
  return -12345;
}

//____________________________________________________________________________
void AliEMCALRecPoint::EvalTime(TClonesArray * digits){
  // time is set to the time of the digit with the maximum energy

  Float_t maxE  = 0;
  Int_t   maxAt = 0;
  for(Int_t idig=0; idig < fMulDigit; idig++){
    if(fEnergyList[idig] > maxE){
      maxE  = fEnergyList[idig] ;
      maxAt = idig;
    }
  }
  fTime = ((AliEMCALDigit*) digits->At(fDigitsList[maxAt]))->GetTime() ;
  
}

//______________________________________________________________________________
void AliEMCALRecPoint::Paint(Option_t *)
{
  // Paint this ALiRecPoint as a TMarker  with its current attributes
  
  TVector3 pos(0.,0.,0.)  ;
  GetLocalPosition(pos)   ;
  Coord_t x = pos.X()     ;
  Coord_t y = pos.Z()     ;
  Color_t markercolor = 1 ;
  Size_t  markersize  = 1.;
  Style_t markerstyle = 5 ;
  
  if (!gPad->IsBatch()) {
    gVirtualX->SetMarkerColor(markercolor) ;
    gVirtualX->SetMarkerSize (markersize)  ;
    gVirtualX->SetMarkerStyle(markerstyle) ;
  }
  gPad->SetAttMarkerPS(markercolor,markerstyle,markersize) ;
  gPad->PaintPolyMarker(1,&x,&y,"") ;
}

//_____________________________________________________________________
Double_t AliEMCALRecPoint::TmaxInCm(const Double_t e , const Int_t key)
{ 
  // e energy in GeV)
  // key  =  0(gamma, default)
  //     !=  0(electron)
  const Double_t ca   = 4.82;  // shower max parameter - first guess; ca=TMath::Log(1000./8.07)
  Double_t tmax = 0.;    // position of electromagnetic shower max in cm

  Double_t x0   = 1.31;  // radiation lenght (cm)
  //If old geometry in use
  if(!((fGeomPtr->GetEMCGeometry()->GetGeoName()).Contains("V1"))) x0 = 1.28;

  if(e>0.1) {
    tmax = TMath::Log(e) + ca;
    if      (key==0) tmax += 0.5; 
    else             tmax -= 0.5;
    tmax *= x0; // convert to cm
  }
  return tmax;
}

//______________________________________________________________________________
Float_t AliEMCALRecPoint::EtaToTheta(Float_t arg) const
{
  //Converts Theta (Radians) to Eta(Radians)
  return (2.*TMath::ATan(TMath::Exp(-arg)));
}

//______________________________________________________________________________
Float_t AliEMCALRecPoint::ThetaToEta(Float_t arg) const
{
  //Converts Eta (Radians) to Theta(Radians)
  return (-1 * TMath::Log(TMath::Tan(0.5 * arg)));
}

//____________________________________________________________________________
void AliEMCALRecPoint::Print(Option_t *opt) const
{
  // Print the list of digits belonging to the cluster
  if(strlen(opt)==0) return;
  TString message ; 
  message  = "AliEMCALRecPoint:\n" ;
  message +=  " digits # = " ; 
  AliInfo(message.Data()) ; 

  Int_t iDigit;
  for(iDigit=0; iDigit<fMulDigit; iDigit++)
    printf(" %d ", fDigitsList[iDigit] ) ;  
  printf("\n");

  AliInfo(" Energies = ") ;
  for(iDigit=0; iDigit<fMulDigit; iDigit++) 
    printf(" %f ", fEnergyList[iDigit] ) ;
  printf("\n");

  AliInfo("\n Abs Ids = ") ;
  for(iDigit=0; iDigit<fMulDigit; iDigit++) 
    printf(" %i ", fAbsIdList[iDigit] ) ;
  printf("\n");

  AliInfo(" Primaries  ") ;
  for(iDigit = 0;iDigit < fMulTrack; iDigit++)
    printf(" %d ", fTracksList[iDigit]) ;

  printf("\n Local x %6.2f y %7.2f z %7.1f \n", fLocPos[0], fLocPos[1], fLocPos[2]);

  message  = "       ClusterType     = %d" ;
  message += "       Multiplicity    = %d" ;
  message += "       Cluster Energy  = %f" ; 
  message += "       Core energy     = %f" ; 
  message += "       Core radius     = %f" ; 
  message += "       Number of primaries %d" ; 
  message += "       Stored at position %d" ; 
  AliInfo(Form(message.Data(), fClusterType, fMulDigit, fAmp, fCoreEnergy, fCoreRadius, fMulTrack, GetIndexInList()) ) ;  
}

//___________________________________________________________
Double_t  AliEMCALRecPoint::GetPointEnergy() const
{
  //Returns energy ....
  Double_t e=0.0;
  for(int ic=0; ic<GetMultiplicity(); ic++) e += double(fEnergyList[ic]);
  return e;
}
