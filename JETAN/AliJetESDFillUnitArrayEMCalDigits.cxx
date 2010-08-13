
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
 
//---------------------------------------------------------------------
// Fill Unit Array 
// Called by ESD reader for jet analysis
// Author: Magali Estienne (magali.estienne@ires.in2p3.fr)
//---------------------------------------------------------------------

// --- ROOT system ---
#include <TLorentzVector.h>
#include <TMath.h>

// --- AliRoot header files ---
#include "AliESDCaloCluster.h"
#include "AliESDCaloCells.h"
#include "AliJetUnitArray.h"
#include "AliJetESDFillUnitArrayEMCalDigits.h"
// #include "AliEMCALCalibData.h"
// #include "AliCDBManager.h"
// // class AliCDBStorage;
// #include "AliCDBEntry.h"

// --- ROOT system ---
class TSystem;
class TGeoManager;
class TArrayS;

// --- AliRoot header files ---
class AliJetFinder;
class AliJetReader;
class AliJetESDReader;
// class AliEMCALCalibData;
// class AliCDBManager;
// class AliCDBStorage;
// class AliCDBEntry;

ClassImp(AliJetESDFillUnitArrayEMCalDigits)

//_____________________________________________________________________________
AliJetESDFillUnitArrayEMCalDigits::AliJetESDFillUnitArrayEMCalDigits()
  : AliJetFillUnitArray(),
    fESD(0),
    fNIn(0),
    fOpt(0),
    fCluster(0),
    fDebug(0),
    fNCEMCAL(0),
    fNCPHOS(0),
    fNCCalo(0),
    fApplyElectronCorrection(kFALSE),
    fApplyFractionHadronicCorrection(kFALSE),
    fFractionHadronicCorrection(0.3),
    fClus(0x0),
    fNDigitEmcal(0),
    fNDigitEmcalCut(0)
{
  // constructor
}

//_____________________________________________________________________________
AliJetESDFillUnitArrayEMCalDigits::AliJetESDFillUnitArrayEMCalDigits(AliESDEvent *esd)
  : AliJetFillUnitArray(),
    fESD(esd),
    fNIn(0),
    fOpt(0),
    fCluster(0),
    fDebug(0),
    fNCEMCAL(0),
    fNCPHOS(0),
    fNCCalo(0),
    fApplyElectronCorrection(kFALSE),
    fApplyFractionHadronicCorrection(kFALSE),
    fFractionHadronicCorrection(0.3),
    fClus(0x0),
    fNDigitEmcal(0),
    fNDigitEmcalCut(0)
{
  // constructor
}

//_____________________________________________________________________________
AliJetESDFillUnitArrayEMCalDigits::AliJetESDFillUnitArrayEMCalDigits(const AliJetESDFillUnitArrayEMCalDigits &det)
  : AliJetFillUnitArray(det),
    fESD(det.fESD),
    fNIn(det.fNIn),
    fOpt(det.fOpt),
    fCluster(det.fCluster),
    fDebug(det.fDebug),
    fNCEMCAL(det.fNCEMCAL),
    fNCPHOS(det.fNCPHOS),
    fNCCalo(det.fNCCalo),
    fApplyElectronCorrection(det.fApplyElectronCorrection),
    fApplyFractionHadronicCorrection(det.fApplyFractionHadronicCorrection),
    fFractionHadronicCorrection(det.fFractionHadronicCorrection),
    fClus(det.fClus),
    fNDigitEmcal(det.fNDigitEmcal),
    fNDigitEmcalCut(det.fNDigitEmcalCut)
{
  // Copy constructor
}

//_____________________________________________________________________________
AliJetESDFillUnitArrayEMCalDigits& AliJetESDFillUnitArrayEMCalDigits::operator=(const AliJetESDFillUnitArrayEMCalDigits& other)
{
  // Assignment

  fESD = other.fESD;
  fNIn = other.fNIn;
  fOpt = other.fOpt;
  fCluster = other.fCluster;
  fDebug = other.fDebug;
  fNCEMCAL = other.fNCEMCAL;
  fNCPHOS = other.fNCPHOS;
  fNCCalo = other.fNCCalo;
  fApplyElectronCorrection = other.fApplyElectronCorrection;
  fApplyFractionHadronicCorrection = other.fApplyFractionHadronicCorrection;
  fFractionHadronicCorrection = other.fFractionHadronicCorrection;
  fClus = other.fClus;
  fNDigitEmcal = other.fNDigitEmcal;
  fNDigitEmcalCut = other.fNDigitEmcalCut;

  return (*this);

}

//_____________________________________________________________________________
AliJetESDFillUnitArrayEMCalDigits::~AliJetESDFillUnitArrayEMCalDigits()
{
  // destructor
}

//_____________________________________________________________________________
void AliJetESDFillUnitArrayEMCalDigits::Exec(Option_t* const /*option*/)
{
  //
  // Main method.
  // Fill the unit array with the neutral particle information from the EMCal cells in ESD
  //

  fDebug = fReaderHeader->GetDebug();
  fOpt = fReaderHeader->GetDetector();
  fCluster = fReaderHeader->GetCluster();

  // Init parameters
  InitParameters();

//(not used ?)  Int_t   nDigitTot      = 0;
  Int_t   goodDigit      = 0;
//(not used ?)  Int_t   beg            = 0;
//(not used ?)  Int_t   end            = 0; 
  Int_t   index          = 0;
//(not used ?)  Int_t   count          = 0;

//(not used ?)  Int_t nRefEnt = fRefArray->GetEntries();

  if(!fCluster) { // Keep all digit information
    // Loop over all cell information
    //------------------------------------------------------------------
    AliESDCaloCells &cells= *(fESD->GetEMCALCells());
    Int_t ncell = cells.GetNumberOfCells() ;
//(not used ?)    Int_t type = cells.GetType();

    for (Int_t icell=  0; icell <  ncell; icell++) {
      Int_t   digitID   = cells.GetCellNumber(icell);
      Float_t digitAmp  = cells.GetAmplitude(icell);
//(not used ?)      Float_t digitTime = cells.GetTime(icell);
      Float_t digitEn   = digitAmp*0.0153; // Last correct
      //  Float_t digitEn = Calibrate(digitAmp,digitID);

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
    Double_t vertexPosition[3] ;
    fESD->GetVertex()->GetXYZ(vertexPosition) ;

    // Get CaloCells
    AliESDCaloCells &cells= *(fESD->GetEMCALCells());
//(not used ?)    Int_t ncell = cells.GetNumberOfCells() ;
//(not used ?)    Int_t type = cells.GetType();

    for(Int_t j = beg; j < nclus; j++) { // loop over clusters
      // Retrieve cluster from esd
      fClus = (AliESDCaloCluster *) caloClusters->At(j) ;

      // Get the cluster info
//(not used ?)      Float_t energy         = fClus->E() ;
//(not used ?)      Int_t   iprim          = fClus->GetLabel();

      fClus->GetPosition(pos) ;
      TVector3 vpos(pos[0],pos[1],pos[2]) ;
      TLorentzVector p ;
      fClus->GetMomentum(p,vertexPosition);

      Int_t     digMult = fClus->GetNCells() ;
      UShort_t *digID   = fClus->GetCellsAbsId() ;
//(not used ?)      Double_t *digAmpFrac = fClus->GetCellsAmplitudeFraction() ;
      Int_t     trackIndex = fClus->GetTrackMatchedIndex();

      // Do double-counted electron correction 
      if (fApplyElectronCorrection != 0 && trackIndex !=-1 )
	{
          // The electron correction go there
          // Under construction !!!!
	}  // End of Electron correction

      // Get CaloCells of cluster and fill the unitArray
      for(Int_t i = 0; i < digMult ; i++){
        Int_t    digitID   = digID[i]; // or clus->GetCellNumber(i) ;
//(not used ?)        Double_t digitAmpFrac = digAmpFrac[i];
        Float_t  digitAmp     = cells.GetCellAmplitude(digitID) ;
//(not used ?)        Float_t  digitTime    = cells.GetCellTime(digitID);

        // Calibration for an energy in GeV
        Float_t digitEn = digitAmp*0.0153;

        Float_t etaD=-10., phiD=-10.;
        fGeom->EtaPhiFromIndex(digitID,etaD,phiD);
        //  fEMCalGrid->GetEtaPhiFromIndex2(digitID,phiD,etaD);

        phiD = ((phiD < 0) ? phiD + 2.* TMath::Pi() : phiD);

        Float_t digitEt = digitEn*TMath::Abs(TMath::Sin(EtaToTheta(etaD)));

	cout << "Digit " << i << ", eta: " << etaD << ", phi: " << phiD << endl;

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
	
	if(fDebug > 12) printf("goodDigit : %d\n", goodDigit);
      
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

/*
//____________________________________________________________________________
void AliJetESDFillUnitArrayEMCalDigits::GetCalibrationParameters()
{
  // Set calibration parameters:
  // if calibration database exists, they are read from database,
  // otherwise, they are taken from digitizer.
  //
  // It is a user responsilibity to open CDB before reconstruction,
  // for example:
  // AliCDBStorage* storage = AliCDBManager::Instance()->GetStorage("local://CalibDB");

  //Check if calibration is stored in data base

  if(!fCalibData && (AliCDBManager::Instance()->IsDefaultStorageSet()))
    {
      AliCDBEntry *entry = (AliCDBEntry*)
        AliCDBManager::Instance()->Get("EMCAL/Calib/Data");
      if (entry) fCalibData =  (AliEMCALCalibData*) entry->GetObject();
    }

  if(!fCalibData)
    printf("************* Calibration parameters not found in CDB! ****************");
//    AliFatal("Calibration parameters not found in CDB!");


}

//____________________________________________________________________________
Float_t  AliJetESDFillUnitArrayEMCalDigits::Calibrate(Int_t amp, Int_t AbsId)
{

  // Convert digitized amplitude into energy.
  // Calibration parameters are taken from calibration data base for raw data,
  // or from digitizer parameters for simulated data.

  if(fCalibData){

    if (fGeom==0)
      printf("************* Did not get geometry from EMCALLoader ***************");
//      AliFatal("Did not get geometry from EMCALLoader") ;

    Int_t iSupMod = -1;
    Int_t nModule  = -1;
    Int_t nIphi   = -1;
    Int_t nIeta   = -1;
    Int_t iphi    = -1;
    Int_t ieta    = -1;

    Bool_t bCell = fGeom->GetCellIndex(AbsId, iSupMod, nModule, nIphi, nIeta) ;
    if(!bCell) {
      // fGeom->PrintGeometry();
      Error("Calibrate()"," Wrong cell id number : %i", AbsId);
      assert(0);
    }

    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,nModule,nIphi, nIeta,iphi,ieta);

    fADCchannelECA  = fCalibData->GetADCchannel (iSupMod,ieta,iphi);
    fADCpedestalECA = fCalibData->GetADCpedestal(iSupMod,ieta,iphi);

   return -fADCpedestalECA + amp * fADCchannelECA ;

  }
  else //Return energy with default parameters if calibration is not available
    return -fADCpedestalECA + amp * fADCchannelECA ;

}
*/

