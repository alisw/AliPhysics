
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
 
// Fill Unit Array 
// Called by ESD reader for jet analysis
// Author: Magali Estienne (magali.estienne@ires.in2p3.fr)

// --- Standard library ---
#include <Riostream.h>
#include <assert.h>

// --- ROOT system ---
#include <TSystem.h>
#include <TLorentzVector.h>
#include <TRefArray.h> 
#include "TTask.h"
#include <TGeoManager.h>
#include <TMath.h>
#include <TClonesArray.h>
#include <TArrayS.h>

// --- AliRoot header files ---
#include "AliJetFinder.h"
#include "AliJetReaderHeader.h"
#include "AliJetReader.h"
#include "AliJetESDReader.h"
#include "AliJetESDReaderHeader.h"
#include "AliESDEvent.h"
#include "AliESDVertex.h"
#include "AliJetDummyGeo.h"
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliJetUnitArray.h"
#include "AliJetFillUnitArrayEMCalDigits.h"
// Remove CDB dependence under construction
//#include "AliEMCALCalibData.h"
//#include "AliCDBManager.h"

//class AliCDBStorage;
//#include "AliCDBEntry.h"


ClassImp(AliJetFillUnitArrayEMCalDigits)

//_____________________________________________________________________________
AliJetFillUnitArrayEMCalDigits::AliJetFillUnitArrayEMCalDigits()
  : TTask("AliJetFillUnitArrayEMCalDigits","Fill Unit Array with tpc/its and emcal information"),
    fESD(0),
    fNumUnits(0), 
    fEtaMinCal(0.),
    fEtaMaxCal(0.),
    fPhiMinCal(0.),
    fPhiMaxCal(0.),
    fNIn(0),
    fOpt(0),
    fCluster(0),
    fDebug(0),
    fNCEMCAL(0),
    fNCPHOS(0),
    fNCCalo(0),
    fTPCGrid(0x0),
    fEMCalGrid(0x0),
    fECorrection(0),
    fReaderHeader(0x0),
    fMomentumArray(0x0),
    fUnitArray(0x0),
    fRefArray(0x0),
    fProcId(kFALSE),
    fGeom(0x0),
    fClus(0x0),
    fNDigitEmcal(0),
    fNDigitEmcalCut(0),
    fCalibData(0x0),
    fADCchannelECA(0),
    fADCpedestalECA(0)
{
  // constructor
}

//_____________________________________________________________________________
AliJetFillUnitArrayEMCalDigits::AliJetFillUnitArrayEMCalDigits(AliESDEvent *esd)
  : TTask("AliJetFillUnitArrayEMCalDigits","Fill Unit Array with tpc/its and emcal information"),
    fESD(esd),
    fNumUnits(0), 
    fEtaMinCal(0.),
    fEtaMaxCal(0.),
    fPhiMinCal(0.),
    fPhiMaxCal(0.),
    fNIn(0),
    fOpt(0),
    fCluster(0),
    fDebug(0),
    fNCEMCAL(0),
    fNCPHOS(0),
    fNCCalo(0),
    fTPCGrid(0x0),
    fEMCalGrid(0x0),
    fECorrection(0),
    fReaderHeader(0x0),
    fMomentumArray(0x0),
    fUnitArray(0x0),
    fRefArray(0x0),
    fProcId(kFALSE),
    fGeom(0x0),
    fClus(0x0),
    fNDigitEmcal(0),
    fNDigitEmcalCut(0),
    fCalibData(0x0),
    fADCchannelECA(0),
    fADCpedestalECA(0)
{
  // constructor 2
}

//_____________________________________________________________________________
AliJetFillUnitArrayEMCalDigits::AliJetFillUnitArrayEMCalDigits(const AliJetFillUnitArrayEMCalDigits &det)
  : TTask(det),//"AliJetFillUnitArrayEMCalDigits","Fill Unit Array with tpc/its and emcal information"),
    fESD(det.fESD),
    fNumUnits(det.fNumUnits), 
    fEtaMinCal(det.fEtaMinCal),
    fEtaMaxCal(det.fEtaMaxCal),
    fPhiMinCal(det.fPhiMinCal),
    fPhiMaxCal(det.fPhiMaxCal),
    fNIn(det.fNIn),
    fOpt(det.fOpt),
    fCluster(det.fCluster),
    fDebug(det.fDebug),
    fNCEMCAL(det.fNCEMCAL),
    fNCPHOS(det.fNCPHOS),
    fNCCalo(det.fNCCalo),
    fTPCGrid(det.fTPCGrid),
    fEMCalGrid(det.fEMCalGrid),
    fECorrection(det.fECorrection),
    fReaderHeader(det.fReaderHeader),
    fMomentumArray(det.fMomentumArray),
    fUnitArray(det.fUnitArray),
    fRefArray(det.fRefArray),
    fProcId(det.fProcId),
    fGeom(det.fGeom),
    fClus(det.fClus),
    fNDigitEmcal(det.fNDigitEmcal),
    fNDigitEmcalCut(det.fNDigitEmcalCut),
    fCalibData(det.fCalibData),
    fADCchannelECA(det.fADCchannelECA),
    fADCpedestalECA(det.fADCpedestalECA)
{
  // Copy constructor
}

//_____________________________________________________________________________
AliJetFillUnitArrayEMCalDigits& AliJetFillUnitArrayEMCalDigits::operator=(const AliJetFillUnitArrayEMCalDigits& other)
{
  // Assignment
  
  fESD = other.fESD;
  fNumUnits = other.fNumUnits; 
  fEtaMinCal = other.fEtaMinCal;
  fEtaMaxCal = other.fEtaMaxCal;
  fPhiMinCal = other.fPhiMinCal;
  fPhiMaxCal = other.fPhiMaxCal;
  fNIn = other.fNIn;
  fOpt = other.fOpt;
  fCluster = other.fCluster;
  fDebug = other.fDebug;
  fNCEMCAL = other.fNCEMCAL;
  fNCPHOS = other.fNCPHOS;
  fNCCalo = other.fNCCalo;
  fTPCGrid = other.fTPCGrid;
  fEMCalGrid = other.fEMCalGrid;
  fECorrection = other.fECorrection;
  fReaderHeader = other.fReaderHeader;
  fMomentumArray = other.fMomentumArray;
  fUnitArray = other.fUnitArray;
  fRefArray = other.fRefArray;
  fProcId = other.fProcId;
  fGeom = other.fGeom;
  fClus = other.fClus;
  fNDigitEmcal = other.fNDigitEmcal;
  fNDigitEmcalCut = other.fNDigitEmcalCut;
  fCalibData = other.fCalibData;
  fADCchannelECA = other.fADCchannelECA;
  fADCpedestalECA = other.fADCpedestalECA;

  return (*this);
  
}

//____________________________________________________________________________
void AliJetFillUnitArrayEMCalDigits::InitParameters()
{
  fNumUnits = fGeom->GetNCells();      // Number of towers in EMCAL

  fEtaMinCal = fGeom->GetArm1EtaMin();
  fEtaMaxCal = fGeom->GetArm1EtaMax();
  fPhiMinCal = fGeom->GetArm1PhiMin();
  fPhiMaxCal = fGeom->GetArm1PhiMax(); 
  fClus      = 0;

  // Get calibration parameters from file or digitizer default values.
  // Under construction
  //  GetCalibrationParameters() ;

  if(fDebug>1) printf("\n EMCAL parameters initiated ! \n");

}

//_____________________________________________________________________________
AliJetFillUnitArrayEMCalDigits::~AliJetFillUnitArrayEMCalDigits()
{
  // destructor
}

//_____________________________________________________________________________
void AliJetFillUnitArrayEMCalDigits::Exec(Option_t* /*option*/)
{
  //
  // Main method.
  // 
  
  fDebug   = fReaderHeader->GetDebug();
  fOpt     = fReaderHeader->GetDetector();
  fCluster = fReaderHeader->GetCluster();
  
  // Init parameters
  InitParameters();
  
  Int_t   goodDigit      = 0;
  Int_t   index          = 0;
  
  if(!fCluster) { // Keep all digit information
    // Loop over all cell information
    //------------------------------------------------------------------
    AliESDCaloCells &cells = *(fESD->GetEMCALCells());
    Int_t            ncell = cells.GetNumberOfCells() ;
    
    for (Int_t icell=  0; icell <  ncell; icell++) {
      Int_t     digitID   = cells.GetCellNumber(icell);
      Double_t  digitAmp  = cells.GetAmplitude(icell);
      Float_t digitEn   = digitAmp*0.0153; // Last correct
      //      Float_t digitEn = Calibrate((Int_t)digitAmp,digitID);
      
      Float_t etaD=-10., phiD=-10.;
      fGeom->EtaPhiFromIndex(digitID,etaD,phiD); 
      //  fEMCalGrid->GetEtaPhiFromIndex2(digitID,phiD,etaD);

      phiD = ((phiD < 0) ? phiD + 2.* TMath::Pi() : phiD);
      
      Float_t digitEt = digitEn*TMath::Abs(TMath::Sin(EtaToTheta(etaD)));

      AliJetUnitArray *uArray = (AliJetUnitArray*)fUnitArray->At(digitID);

      if(uArray->GetUnitEnergy() == 0.) goodDigit++;
      uArray->SetUnitTrackID(digitID);
      
      Float_t unitEnergy = 0.;
      Bool_t ok = kFALSE;
      unitEnergy = uArray->GetUnitEnergy();

      if(unitEnergy==0){
	if(!fProcId){
	  new(fRefArray) TRefArray(TProcessID::GetProcessWithUID(uArray));
	  fProcId = kTRUE;
	}
	fRefArray->Add(uArray);
	fNDigitEmcal++;
	ok = kTRUE;
      }
      
      // Detector flag
      if(unitEnergy>0.)
	uArray->SetUnitDetectorFlag(kAll);
      else uArray->SetUnitDetectorFlag(kEmcal);
      
      uArray->SetUnitEnergy(unitEnergy+digitEt);

      uArray->SetUnitCutFlag(kPtHigher);

      // To be modified !!!
      uArray->SetUnitSignalFlag(kGood);
	
      // This is for jet multiplicity
      uArray->SetUnitClusterID(index);
	
      if(fDebug > 1) printf("goodDigit : %d\n", goodDigit);

    } // End loop over cells

  } // end if !fCluster
  else { // Keep digit information from clusterization
    
    // Loop over calo clusters
    //------------------------------------------------------------------

    //select EMCAL clusters only
    TRefArray * caloClusters  = new TRefArray();
    fESD->GetEMCALClusters(caloClusters);

    // Total number of EMCAL cluster
    Int_t nclus = caloClusters->GetEntries() ;
    Int_t beg   = 0;
    Float_t pos[3] ;

    // Get reconstructed vertex position
    Double_t vertex_position[3] ;
    fESD->GetVertex()->GetXYZ(vertex_position) ;

    // Get CaloCells
    AliESDCaloCells &cells= *(fESD->GetEMCALCells());
    
    for(Int_t j = beg; j < nclus; j++) { // loop over clusters
      // Retrieve cluster from esd
      fClus = (AliESDCaloCluster *) caloClusters->At(j) ;

      // Get the cluster info
      Float_t energy       = fClus->E() ;
      Int_t iprim          = fClus->GetLabel();
      Int_t trackIndex     = fClus->GetTrackMatched();

      fClus->GetPosition(pos) ;
      TVector3 vpos(pos[0],pos[1],pos[2]) ;
      TLorentzVector p ;
      fClus->GetMomentum(p,vertex_position);

      Int_t     digMult = fClus->GetNCells() ;
      UShort_t *digID   = fClus->GetCellsAbsId() ;
      //Print cluster info
      if(fDebug>2) cout<<"Cluster "<< j <<"; digits mult "<<digMult<<"; type "<<(Int_t )fClus->GetClusterType()
		       <<"; Energy "<<energy<< "; transverse energy:" <<  energy*TMath::Abs(TMath::Sin(EtaToTheta(vpos.Eta())))
		       <<"; Phi "<<vpos.Phi()<<"; Eta "<<vpos.Eta() <<"; label "<<iprim<<endl;

      // Do double-counted electron correction 
      if (fECorrection != 0 && trackIndex !=-1 )
	{
	  // The electron correction go there
	  // Under construction !!!!

	} // End of Electron correction 

      // Get CaloCells of cluster and fill the unitArray
      for(Int_t i = 0; i < digMult ; i++)
	{
	  Int_t     digitID     = digID[i]; // or clus->GetCellNumber(i) ;
	  Double_t  digitAmp    = cells.GetCellAmplitude(digitID) ;
	  
	  // Calibration for an energy in GeV
	  Float_t digitEn = digitAmp*0.0153;
	  //	  Float_t digitEn = Calibrate((Int_t)digitAmp,digitID);	  

	  Float_t etaD=-10., phiD=-10.;
	  fGeom->EtaPhiFromIndex(digitID,etaD,phiD);
	  //          fEMCalGrid->GetEtaPhiFromIndex2(digitID,phiD,etaD);
	  
	  phiD = ((phiD < 0) ? phiD + 2.* TMath::Pi() : phiD);
	  
	  Float_t digitEt = digitEn*TMath::Abs(TMath::Sin(EtaToTheta(etaD)));
	  
	  AliJetUnitArray *uArray = (AliJetUnitArray*)fUnitArray->At(digitID);
	  if(uArray->GetUnitEnergy() == 0.) goodDigit++;
	  uArray->SetUnitTrackID(digitID);
	  
	  Float_t unitEnergy = 0.;
	  Bool_t ok = kFALSE;
	  unitEnergy = uArray->GetUnitEnergy();
	  
	  if(unitEnergy==0){
	    if(!fProcId){
	      new(fRefArray) TRefArray(TProcessID::GetProcessWithUID(uArray));
	      fProcId = kTRUE;
	    }
	    fRefArray->Add(uArray);
	    fNDigitEmcal++;
	    ok = kTRUE;
	  }
	  
	  // Detector flag
	  if(unitEnergy>0.)
	    uArray->SetUnitDetectorFlag(kAll);
	  else uArray->SetUnitDetectorFlag(kEmcal);
	  
	  uArray->SetUnitEnergy(unitEnergy+digitEt);
	  uArray->SetUnitCutFlag(kPtHigher);
	  
	  // To be modified !!!
	  uArray->SetUnitSignalFlag(kGood);
	  
	  // This is for jet multiplicity
	  uArray->SetUnitClusterID(index);
	  
	} // End loop over cells
    } // End loop over clusters
  } // end else
    
  fNIn += goodDigit;

  if(fDebug>1) 
    {
      printf("End of digits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
      printf("goodDigit : %d\n", goodDigit);
    }
  
}

// //____________________________________________________________________________
// void AliJetFillUnitArrayEMCalDigits::GetCalibrationParameters()
// {
//   // Set calibration parameters:
//   // if calibration database exists, they are read from database,
//   // otherwise, they are taken from digitizer.
//   //
//   // It is a user responsilibity to open CDB before reconstruction,
//   // for example:
//   // AliCDBStorage* storage = AliCDBManager::Instance()->GetStorage("local://CalibDB");

//   //Check if calibration is stored in data base

//   if(!fCalibData && (AliCDBManager::Instance()->IsDefaultStorageSet()))
//     {
//       AliCDBEntry *entry = (AliCDBEntry*)
//         AliCDBManager::Instance()->Get("EMCAL/Calib/Data");
//       if (entry) fCalibData =  (AliEMCALCalibData*) entry->GetObject();
//     }

//   if(!fCalibData)
//     printf("************* Calibration parameters not found in CDB! ****************");
// //    AliFatal("Calibration parameters not found in CDB!");


// }

// //____________________________________________________________________________
// Float_t  AliJetFillUnitArrayEMCalDigits::Calibrate(Int_t amp, Int_t AbsId)
// {

//   // Convert digitized amplitude into energy.
//   // Calibration parameters are taken from calibration data base for raw data,
//   // or from digitizer parameters for simulated data.

//   if(fCalibData){

//     if (fGeom==0)
//       printf("************* Did not get geometry from EMCALLoader ***************");
    
//     Int_t iSupMod  = -1;
//     Int_t nModule  = -1;
//     Int_t nIphi    = -1;
//     Int_t nIeta    = -1;
//     Int_t iphi     = -1;
//     Int_t ieta     = -1;

//     Bool_t bCell = fGeom->GetCellIndex(AbsId, iSupMod, nModule, nIphi, nIeta) ;
//     if(!bCell) {
//       // fGeom->PrintGeometry();
//       Error("Calibrate()"," Wrong cell id number : %i", AbsId);
//       assert(0);
//     }

//     fGeom->GetCellPhiEtaIndexInSModule(iSupMod,nModule,nIphi, nIeta,iphi,ieta);

//     fADCchannelECA  = fCalibData->GetADCchannel (iSupMod,ieta,iphi);
//     fADCpedestalECA = fCalibData->GetADCpedestal(iSupMod,ieta,iphi);

//    return -fADCpedestalECA + amp * fADCchannelECA ;

//   }
//   else //Return energy with default parameters if calibration is not available
//     return -fADCpedestalECA + amp * fADCchannelECA ;

// }


//_____________________________________________________________________________
Float_t  AliJetFillUnitArrayEMCalDigits::EtaToTheta(Float_t arg)
{
  //  return (180./TMath::Pi())*2.*atan(exp(-arg));
  return 2.*atan(exp(-arg));


}
