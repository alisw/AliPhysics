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


//======================================================================
// ***July 2009
// Fill Unit Array class 
// Base class used by AliJetESDReader to fill a UnitArray from the information extracted 
// from the particle tracks or emcal cells
// Author: magali.estienne@subatech.in2p3.fr
//======================================================================


#include "AliJetFillUnitArray.h"

// --- ROOT system ---
class TSystem;
class TLorentzVector;
class TVector3;
class TGeoManager;
class TProcessID;

// --- AliRoot header files ---
class AliJetFinder;
class AliJetReader;
class AliJetESDReader;
class AliJetESDReaderHeader;
class AliJetUnitArray;

ClassImp(AliJetFillUnitArray)

//_____________________________________________________________________________
AliJetFillUnitArray::AliJetFillUnitArray()
  : TTask("AliJetFillUnitArray","Fill Unit Array with tpc/its and emcal information"),
    fNTracks(0),
    fNTracksCut(0),
    fOpt(0),
    fDZ(0),
    fDebug(0),
    fReaderHeader(0x0),
    fMomentumArray(0x0),
    fUnitArray(0x0),
    fRefArray(0x0),
    fRef(0x0),
    fSignalFlag(0),
    fCutFlag(0),
    fProcId(kFALSE),
    fTPCGrid(0x0),
    fEMCalGrid(0x0),
    fGeom(0x0),
    fNphi(0),
    fNeta(0),
    fGrid(0),
    fPhi2(0),
    fEta2(0),
    fIndex(0x0),
    fParams(0x0),
    fPhiMin(0),
    fPhiMax(0),
    fEtaMin(0),
    fEtaMax(0),
    fEtaBinInTPCAcc(0),
    fPhiBinInTPCAcc(0),
    fEtaBinInEMCalAcc(0),
    fPhiBinInEMCalAcc(0),
    fNbinPhi(0)
{
  // constructor
}

AliJetFillUnitArray::AliJetFillUnitArray(const AliJetFillUnitArray& cpfrom)
  : TTask("AliJetFillUnitArray","Fill Unit Array with tpc/its and emcal information"),
    fNTracks(0),
    fNTracksCut(0),
    fOpt(0),
    fDZ(0),
    fDebug(0),
    fReaderHeader(0x0),
    fMomentumArray(0x0),
    fUnitArray(0x0),
    fRefArray(0x0),
    fRef(0x0),
    fSignalFlag(0),
    fCutFlag(0),
    fProcId(kFALSE),
    fTPCGrid(0x0),
    fEMCalGrid(0x0),
    fGeom(0x0),
    fNphi(0),
    fNeta(0),
    fGrid(0),
    fPhi2(0),
    fEta2(0),
    fIndex(0x0),
    fParams(0x0),
    fPhiMin(0),
    fPhiMax(0),
    fEtaMin(0),
    fEtaMax(0),
    fEtaBinInTPCAcc(0),
    fPhiBinInTPCAcc(0),
    fEtaBinInEMCalAcc(0),
    fPhiBinInEMCalAcc(0),
    fNbinPhi(0)
{
    //
    // Copy constructor
    //
    fNTracks          = cpfrom.fNTracks;
    fNTracksCut       = cpfrom.fNTracksCut;
    fOpt              = cpfrom.fOpt;
    fDZ               = cpfrom.fDZ;
    fDebug            = cpfrom.fDebug;
    fReaderHeader     = cpfrom.fReaderHeader;
    fMomentumArray    = cpfrom.fMomentumArray;
    fUnitArray        = cpfrom.fUnitArray;
    fRefArray         = cpfrom.fRefArray;
    fRef              = cpfrom.fRef;
    fSignalFlag       = cpfrom.fSignalFlag;
    fCutFlag          = cpfrom.fCutFlag;  
    fProcId           = cpfrom.fProcId;
    fTPCGrid          = cpfrom.fTPCGrid;
    fEMCalGrid        = cpfrom.fEMCalGrid;
    fGeom             = cpfrom.fGeom;
    fNphi             = cpfrom.fNphi;   
    fNeta             = cpfrom.fNeta;  
    fGrid             = cpfrom.fGrid;
    fPhi2             = cpfrom.fPhi2; 
    fEta2             = cpfrom.fEta2; 
    fIndex            = cpfrom.fIndex;
    fParams           = cpfrom.fParams;
    fPhiMin           = cpfrom.fPhiMin;
    fPhiMax           = cpfrom.fPhiMax;
    fEtaMin           = cpfrom.fEtaMin;
    fEtaMax           = cpfrom.fEtaMax;
    fEtaBinInTPCAcc   = cpfrom.fEtaBinInTPCAcc;
    fPhiBinInTPCAcc   = cpfrom.fPhiBinInTPCAcc;
    fEtaBinInEMCalAcc = cpfrom.fEtaBinInEMCalAcc;
    fPhiBinInEMCalAcc = cpfrom.fPhiBinInEMCalAcc;
    fNbinPhi          = cpfrom.fNbinPhi;
}

AliJetFillUnitArray& AliJetFillUnitArray::operator=(const AliJetFillUnitArray& rhs)
{
    //
    // Assignment operator
    //
    if (this != &rhs) {
	fNTracks          = rhs.fNTracks;
	fNTracksCut       = rhs.fNTracksCut;
	fOpt              = rhs.fOpt;
	fDZ               = rhs.fDZ;
	fDebug            = rhs.fDebug;
	fReaderHeader     = rhs.fReaderHeader;
	fMomentumArray    = rhs.fMomentumArray;
	fUnitArray        = rhs.fUnitArray;
	fRefArray         = rhs.fRefArray;
	fRef              = rhs.fRef;
	fSignalFlag       = rhs.fSignalFlag;
	fCutFlag          = rhs.fCutFlag;  
	fProcId           = rhs.fProcId;
	fTPCGrid          = rhs.fTPCGrid;
	fEMCalGrid        = rhs.fEMCalGrid;
	fGeom             = rhs.fGeom;
	fNphi             = rhs.fNphi;   
	fNeta             = rhs.fNeta;  
	fGrid             = rhs.fGrid;
	fPhi2             = rhs.fPhi2; 
	fEta2             = rhs.fEta2; 
	fIndex            = rhs.fIndex;
	fParams           = rhs.fParams;
	fPhiMin           = rhs.fPhiMin;
	fPhiMax           = rhs.fPhiMax;
	fEtaMin           = rhs.fEtaMin;
	fEtaMax           = rhs.fEtaMax;
	fEtaBinInTPCAcc   = rhs.fEtaBinInTPCAcc;
	fPhiBinInTPCAcc   = rhs.fPhiBinInTPCAcc;
	fEtaBinInEMCalAcc = rhs.fEtaBinInEMCalAcc;
	fPhiBinInEMCalAcc = rhs.fPhiBinInEMCalAcc;
	fNbinPhi          = rhs.fNbinPhi;
    }
	return *this;
}

//_____________________________________________________________________________
AliJetFillUnitArray::~AliJetFillUnitArray()
{
  // destructor
}

//_____________________________________________________________________________
void AliJetFillUnitArray::GetEtaPhiFromIndex(Int_t index, Float_t &eta, Float_t &phi)
{
  // Get the eta,phi position from the index

  for(Int_t j=0; j<fNphi+1; j++) {
    for(Int_t i=0; i<fNeta+1; i++) {

      // TPC grid only 
      //-------------------------------------
      if(fGrid==0) {	
	if(j*(fNeta+1)+i == index) {
	  eta = fEta2->At(i); 
	  phi = fPhi2->At(j);
	}
      }

      // TPC-EMCAL grid
      //-------------------------------------
      Int_t ii = 0;
      if(i==0) ii = 0;
      if(i>0 && i<(fEtaBinInTPCAcc-fEtaBinInEMCalAcc)/2) ii = i; 
      if(i>=(fEtaBinInTPCAcc+fEtaBinInEMCalAcc)/2 && i<fNeta+1) ii = i-fEtaBinInEMCalAcc;

      if(fGrid==1) {
	if(j<(fNbinPhi+1) && j*(fNeta+1)+i == index) {
	  eta = fEta2->At(i);
	  phi = fPhi2->At(j);
	}  

	if((j>=(fNbinPhi+1) && j<(fNbinPhi+1+fPhiBinInEMCalAcc)) && 
	   ((fNbinPhi+1)*(fNeta+1) + (j-fNbinPhi-1)*(fEtaBinInTPCAcc-fEtaBinInEMCalAcc) + ii)== index ) {
	  if(ii==0) {Int_t ind = 0; eta = fEta2->At(ind);}
	  else eta = fEta2->At(i);
	  phi = fPhi2->At(j);
	}

	if(j>=(fNbinPhi+1+fPhiBinInEMCalAcc) && ((fNbinPhi+1)*(fNeta+1)+fPhiBinInEMCalAcc*((fEtaBinInTPCAcc-fEtaBinInEMCalAcc))+(j-(fNbinPhi+1+fPhiBinInEMCalAcc))*(fNeta+1)+i == index)) {
	  eta = fEta2->At(i);
	  phi = fPhi2->At(j);
	}
      }
    }
  }
}

//_____________________________________________________________________________
Float_t  AliJetFillUnitArray::EtaToTheta(Float_t arg)
{
  return 2.*atan(exp(-arg));
}







