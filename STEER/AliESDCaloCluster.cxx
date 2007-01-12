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
/* $Log $ */

//-----------------------------------------------------------------
//           Implementation of the ESD Calorimeter cluster class
//   ESD = Event Summary Data
//   This is the class to deal with during the phisics analysis of data
//
//   J.L. Klay (LLNL)
//-----------------------------------------------------------------

#include <TLorentzVector.h>
#include "AliESDCaloCluster.h"

ClassImp(AliESDCaloCluster)

//_______________________________________________________________________
AliESDCaloCluster::AliESDCaloCluster() : 
  fID(0),
  fClusterType(-1),
  fEMCALCluster(kFALSE),
  fPHOSCluster(kFALSE),
  fEnergy(-1),
  fDispersion(-1),
  fChi2(-1),
  fPrimaryIndex(-1),
  fM20(0),
  fM02(0),
  fM11(0),
  fNExMax(0),
  fEmcCpvDistance(9999),
  fNumberOfDigits(0),
  fDigitAmplitude(0x0),
  fDigitTime(0x0),
  fDigitIndex(0x0)
{
  //
  // The default ESD constructor 
  //
  fGlobalPos[0] = fGlobalPos[1] = fGlobalPos[2] = 0.;
  for(Int_t i=0; i<AliPID::kSPECIESN; i++) fPID[i] = 0.;
}

//_______________________________________________________________________
AliESDCaloCluster::AliESDCaloCluster(const AliESDCaloCluster& clus) : 
  TObject(clus),
  fID(clus.fID),
  fClusterType(clus.fClusterType),
  fEMCALCluster(clus.fEMCALCluster),
  fPHOSCluster(clus.fPHOSCluster),
  fEnergy(clus.fEnergy),
  fDispersion(clus.fDispersion),
  fChi2(clus.fChi2),
  fPrimaryIndex(clus.fPrimaryIndex),
  fM20(clus.fM20),
  fM02(clus.fM02),
  fM11(clus.fM11),
  fNExMax(clus.fNExMax),
  fEmcCpvDistance(clus.fEmcCpvDistance),
  fNumberOfDigits(clus.fNumberOfDigits),
  fDigitAmplitude(0x0),
  fDigitTime(0x0),
  fDigitIndex(0x0)
{
  //
  // The copy constructor 
  //
  fGlobalPos[0] = clus.fGlobalPos[0];
  fGlobalPos[1] = clus.fGlobalPos[1];
  fGlobalPos[2] = clus.fGlobalPos[2];

  for(Int_t i=0; i<AliPID::kSPECIESN; i++) fPID[i] = clus.fPID[i];

  if (clus.fNumberOfDigits > 0) {
    if (clus.fDigitAmplitude) {
      fDigitAmplitude = new UShort_t[clus.fNumberOfDigits];
      for (Int_t i=0; i<clus.fNumberOfDigits; i++)
	fDigitAmplitude[i]=clus.fDigitAmplitude[i];
    }
    if (clus.fDigitTime) {
      fDigitTime = new UShort_t[clus.fNumberOfDigits];
      for (Int_t i=0; i<clus.fNumberOfDigits; i++)
	fDigitTime[i]=clus.fDigitTime[i];
    }
    if (clus.fDigitIndex) {
      fDigitIndex = new UShort_t[clus.fNumberOfDigits];
      for (Int_t i=0; i<clus.fNumberOfDigits; i++)
	fDigitIndex[i]=clus.fDigitIndex[i];
    }
  }
}

//_______________________________________________________________________
AliESDCaloCluster &AliESDCaloCluster::operator=(const AliESDCaloCluster& source)
{
  // assignment operator

  if(&source == this) return *this;

  fID = source.fID;
  fClusterType = source.fClusterType;
  fEMCALCluster = source.fEMCALCluster;
  fPHOSCluster = source.fPHOSCluster;
  fEnergy = source.fEnergy;
  fDispersion = source.fDispersion;
  fChi2 = source.fChi2;
  fPrimaryIndex = source.fPrimaryIndex;
  fM20 = source.fM20;
  fM02 = source.fM02;
  fM11 = source.fM11;
  fNExMax = source.fNExMax;
  fEmcCpvDistance = source.fEmcCpvDistance;
  fNumberOfDigits = source.fNumberOfDigits;
  delete fDigitAmplitude; fDigitAmplitude=0x0;
  delete fDigitTime; fDigitTime = 0x0;
  delete fDigitIndex; fDigitIndex = 0x0;

  fGlobalPos[0] = source.fGlobalPos[0];
  fGlobalPos[1] = source.fGlobalPos[1];
  fGlobalPos[2] = source.fGlobalPos[2];

  for(Int_t i=0; i<AliPID::kSPECIESN; i++) fPID[i] = source.fPID[i];

  if (source.fNumberOfDigits > 0) {
    if (source.fDigitAmplitude) {
      fDigitAmplitude = new UShort_t[source.fNumberOfDigits];
      for (Int_t i=0; i<source.fNumberOfDigits; i++)
	fDigitAmplitude[i]=source.fDigitAmplitude[i];
    }
    if (source.fDigitTime) {
      fDigitTime = new UShort_t[source.fNumberOfDigits];
      for (Int_t i=0; i<source.fNumberOfDigits; i++)
	fDigitTime[i]=source.fDigitTime[i];
    }
    if (source.fDigitIndex) {
      fDigitIndex = new UShort_t[source.fNumberOfDigits];
      for (Int_t i=0; i<source.fNumberOfDigits; i++)
	fDigitIndex[i]=source.fDigitIndex[i];
    }
  }

  return *this;

}


//_______________________________________________________________________
AliESDCaloCluster::~AliESDCaloCluster(){ 
  //
  // This is destructor according Coding Conventrions 
  //
  // AliESDCaloCluster is the owner of the arrays
  // even if they are created outside

  delete[] fDigitAmplitude;
  delete[] fDigitTime;
  delete[] fDigitIndex;

}

//_______________________________________________________________________
void AliESDCaloCluster::SetPid(const Float_t *p) {
  // Sets the probability of each particle type
  // Copied from AliESDtrack SetPIDValues
  // This function copies "n" PID weights from "scr" to "dest"
  // and normalizes their sum to 1 thus producing conditional
  // probabilities.
  // The negative weights are set to 0.
  // In case all the weights are non-positive they are replaced by
  // uniform probabilities

  Int_t n = AliPID::kSPECIESN;

  Float_t uniform = 1./(Float_t)n;

  Float_t sum = 0;
  for (Int_t i=0; i<n; i++)
    if (p[i]>=0) {
      sum+=p[i];
      fPID[i] = p[i];
    }
    else {
      fPID[i] = 0;
    }

  if(sum>0)
    for (Int_t i=0; i<n; i++) fPID[i] /= sum;
  else
    for (Int_t i=0; i<n; i++) fPID[i] = uniform;

}

//_______________________________________________________________________
void AliESDCaloCluster::GetMomentum(TLorentzVector& p) {
  // Returns TLorentzVector with momentum of the cluster. Only valid for clusters 
  // identified as photons or pi0 (overlapped gamma) produced on the vertex
  
  Double_t r = TMath::Sqrt(fGlobalPos[0]*fGlobalPos[0]+
		            fGlobalPos[1]*fGlobalPos[1]+
		            fGlobalPos[2]*fGlobalPos[2]   ) ; 

  p.SetPxPyPzE( fEnergy*fGlobalPos[0]/r,  fEnergy*fGlobalPos[1]/r,  fEnergy*fGlobalPos[2]/r,  fEnergy) ; 
  
}
