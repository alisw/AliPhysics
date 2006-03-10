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

#include "AliESDCaloCluster.h"

ClassImp(AliESDCaloCluster)

//_______________________________________________________________________
AliESDCaloCluster::AliESDCaloCluster() : 
  fID(0),
  fClusterType(-1),
  fEMCALCluster(kFALSE),
  fEnergy(-1),
  fDispersion(-1),
  fChi2(-1),
  fNumberOfDigits(0),
  fDigitAmplitude(0),
  fDigitTime(0),
  fDigitIndex(0)
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
  fEnergy(clus.fEnergy),
  fDispersion(clus.fDispersion),
  fChi2(clus.fChi2),
  fNumberOfDigits(clus.fNumberOfDigits)
{
  //
  // The copy constructor 
  //
  fGlobalPos[0] = clus.fGlobalPos[0];
  fGlobalPos[1] = clus.fGlobalPos[1];
  fGlobalPos[2] = clus.fGlobalPos[2];

  for(Int_t i=0; i<AliPID::kSPECIESN; i++) fPID[i] = clus.fPID[i];

  fDigitAmplitude = clus.fDigitAmplitude;
  fDigitTime = clus.fDigitTime;
  fDigitIndex = clus.fDigitIndex;

}


//_______________________________________________________________________
AliESDCaloCluster::~AliESDCaloCluster(){ 
  //
  // This is destructor according Coding Conventrions 
  //
  //printf("Delete cluster\n");

  //Not sure why but it won't let me delete these in the dtor here.
  //The Reconstruction gives me the error
  //*** glibc detected *** double free or corruption (!prev):
  //0x0c1550b0 ***
  /*
  if(fDigitAmplitude)
    delete[] fDigitAmplitude;
  if(fDigitTime)
    delete[] fDigitTime;
  if(fDigitIndex)
    delete[] fDigitIndex;
  */
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
