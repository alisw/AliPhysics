
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

//------------------------------------------------------------- 
// Fill Unit Array with EMCal information 
// Called by ESD reader for jet analysis
// Author: Magali Estienne (magali.estienne@ires.in2p3.fr)
//------------------------------------------------------------- 

// --- Standard library ---
#include <Riostream.h>

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
//#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliJetDummyGeo.h"
#include "AliESDCaloCluster.h"
#include "AliJetUnitArray.h"
#include "AliJetFillUnitArrayEMCalDigits.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"

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
    fDebug(0),
    fNCEMCAL(0),
    fNCPHOS(0),
    fNCCalo(0),
    fTPCGrid(0x0),
    fEMCalGrid(0x0),
    fReaderHeader(0x0),
    fMomentumArray(0x0),
    fUnitArray(0x0),
    fRefArray(0x0),
    fGeom(0x0),
    fClus(0x0),
    fNDigitEmcal(0),
    fNDigitEmcalCut(0)
{
  // constructor
}

//_____________________________________________________________________________
AliJetFillUnitArrayEMCalDigits::AliJetFillUnitArrayEMCalDigits(AliESDEvent */*esd*/)
  : TTask("AliJetFillUnitArrayEMCalDigits","Fill Unit Array with tpc/its and emcal information"),
    fESD(0),    
    fNumUnits(0), 
    fEtaMinCal(0.),
    fEtaMaxCal(0.),
    fPhiMinCal(0.),
    fPhiMaxCal(0.),
    fNIn(0),
    fOpt(0),
    fDebug(0),
    fNCEMCAL(0),
    fNCPHOS(0),
    fNCCalo(0),
    fTPCGrid(0x0),
    fEMCalGrid(0x0),
    fReaderHeader(0x0),
    fMomentumArray(0x0),
    fUnitArray(0x0),
    fRefArray(0x0),
    fGeom(0x0),
    fClus(0x0),
    fNDigitEmcal(0),
    fNDigitEmcalCut(0)
{
  // constructor
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
  // Explain

  fDebug = fReaderHeader->GetDebug();
  fOpt = fReaderHeader->GetDetector();

  // Init parameters
  InitParameters();

  // Get number of clusters from EMCAL
  
  Int_t   nDigitTot      = 0;
  Int_t   goodDigit      = 0;
  Int_t   beg            = 0;
  Int_t   end            = 0;
  Float_t ptMin          = fReaderHeader->GetPtCut();

  // Loop over calo clusters 
  //------------------------------------------------------------------
  Int_t type = 0; 
  Int_t index = 0;

  // Total number of EMCAL cluster
  end =  fESD->GetNumberOfCaloClusters();

  for(Int_t j = beg; j < end; j++) {
      fClus = fESD->GetCaloCluster(j);
      if(!fClus->IsEMCAL()) continue;

      type = fClus->GetClusterType(); 
      index = fClus->GetID();
      nDigitTot = fClus->GetNumberOfDigits();
      
      // Keep clusters or pseudo clusters
      if (type != AliESDCaloCluster::kEMCALClusterv1) continue;
      //      if (type != AliESDCaloCluster::kPseudoCluster) continue;

      // Get the digit index and the digit information
      //============================================================
	  
      // Get number of digits in a cluster
      Int_t nD = fClus->GetNumberOfDigits();

      TArrayS *digID = fClus->GetDigitIndex();      
      TArrayS *digEnergy = fClus->GetDigitAmplitude();
      Float_t *digitEnergy = new Float_t[nD];      
      //      Float_t digitEn = 0.;
	  
      // Loop over digits
      for(Int_t k=0; k<nD; k++) {
	
	// Convert energy in GeV
	Int_t idF = (Int_t)digID->At(k);
	// Calibration for an energy in GeV
	digitEnergy[k] = (Float_t)digEnergy->At(k)/500.; 

	// Second method to extract eta, phi positions of a digit
	//=================================================================

	Float_t etaD=-10., phiD=-10.;
	fGeom->EtaPhiFromIndex(idF,etaD,phiD); 
	//	    fEMCalGrid->GetEtaPhiFromIndex2(idF,phiD,etaD);
	phiD = ((phiD < 0) ? phiD + 2.* TMath::Pi() : phiD);

	Float_t etDigit = digitEnergy[k]*TMath::Abs(TMath::Sin(EtaToTheta(etaD)));

	AliJetUnitArray *uArray = (AliJetUnitArray*)fUnitArray->At(idF);
	if(uArray->GetUnitEnergy() == 0.) goodDigit++;

	Float_t unitEnergy = 0.;
	Bool_t ok = kFALSE;
	unitEnergy = uArray->GetUnitEnergy();
	if(unitEnergy==0){
	  fRefArray->Add(uArray);
	  fNDigitEmcal++;
	  ok = kTRUE;
	}
	uArray->SetUnitEnergy(unitEnergy+etDigit);
	// Put a cut flag
	if(uArray->GetUnitEnergy()<ptMin)
	  uArray->SetUnitCutFlag(kPtSmaller);
	else {
	  uArray->SetUnitCutFlag(kPtHigher);
	  if(ok) fNDigitEmcalCut++;
	}
	// Detector flag
	if(unitEnergy>0.) 
	  uArray->SetUnitDetectorFlag(kAll);
	else uArray->SetUnitDetectorFlag(kEmcal);
	
	// This is for jet multiplicity
	uArray->SetUnitClusterID(index);
	
	if(fDebug > 12) printf("goodDigit : %d\n", goodDigit);
	
      } // End loop over digits
      
  } // End loop over clusters
  
  fNIn += goodDigit;

}

//_____________________________________________________________________________
Float_t  AliJetFillUnitArrayEMCalDigits::EtaToTheta(Float_t arg)
{
  //  return (180./TMath::Pi())*2.*atan(exp(-arg));
  return 2.*atan(exp(-arg));


}
