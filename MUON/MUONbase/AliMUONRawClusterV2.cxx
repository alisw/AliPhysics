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

//-----------------------------------------------------------------------------
/// \class AliMUONRawClusterV2
///
/// Class for the MUON RecPoint
///
/// \author Philippe Pillot, Subatech
//-----------------------------------------------------------------------------


#include "AliMUONRawClusterV2.h"

#include "AliLog.h"

#include <TClonesArray.h>
#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMUONRawClusterV2)
/// \endcond


//____________________________________________________
AliMUONRawClusterV2::AliMUONRawClusterV2() 
  : AliMUONVCluster(),
    fX(FLT_MAX),
    fY(FLT_MAX),
    fZ(FLT_MAX),
    fErrX2(FLT_MAX),
    fErrY2(FLT_MAX),
    fQ(0.),
    fChi2(0.),
    fNDigits(0),
    fDigitsId(0x0),
    fMCLabel(-1)
{
  /// Default Constructor
}

//_____________________________________________________________________________
AliMUONRawClusterV2::AliMUONRawClusterV2(Int_t chamberId, Int_t detElemId, Int_t clusterIndex)
  : AliMUONVCluster(chamberId, detElemId, clusterIndex),
    fX(FLT_MAX),
    fY(FLT_MAX),
    fZ(FLT_MAX),
    fErrX2(FLT_MAX),
    fErrY2(FLT_MAX),
    fQ(0.),
    fChi2(0.),
    fNDigits(0),
    fDigitsId(0x0),
    fMCLabel(-1)
{
  /// Constructor
}

//____________________________________________________
AliMUONRawClusterV2::~AliMUONRawClusterV2() 
{
  /// Destructor
  delete [] fDigitsId;
}

//____________________________________________________
AliMUONRawClusterV2::AliMUONRawClusterV2(const AliMUONRawClusterV2& cluster)
  : AliMUONVCluster(cluster),
    fX(cluster.fX),
    fY(cluster.fY),
    fZ(cluster.fZ),
    fErrX2(cluster.fErrX2),
    fErrY2(cluster.fErrY2),
    fQ(cluster.fQ),
    fChi2(cluster.fChi2),
    fNDigits(cluster.fNDigits),
    fDigitsId(0x0),
    fMCLabel(cluster.fMCLabel)

{
  /// Copy constructor
  
  if (cluster.fDigitsId) {
    fDigitsId = new UInt_t[fNDigits];
    memcpy(fDigitsId,cluster.fDigitsId, fNDigits*sizeof(UInt_t));
  }
}

  //__________________________________________________________________________
AliMUONRawClusterV2 & AliMUONRawClusterV2::operator=(const AliMUONRawClusterV2& cluster)
{
  /// Asignment operator
  
  // check assignement to self
  if (this == &cluster)
    return *this;

  // base class assignement
  AliMUONVCluster::operator=(cluster);
  
  fX = cluster.fX;
  fY = cluster.fY;
  fZ = cluster.fZ;
  fErrX2 = cluster.fErrX2;
  fErrY2 = cluster.fErrY2;
  fQ = cluster.fQ;
  fChi2 = cluster.fChi2;
  SetDigitsId(cluster.fNDigits,cluster.fDigitsId);
  fMCLabel = cluster.fMCLabel;

  return *this;
}

//____________________________________________________
void AliMUONRawClusterV2::Clear(Option_t*)
{
  /// clear memory
  delete [] fDigitsId;
  fDigitsId = 0x0;
  fNDigits = 0;
}

//____________________________________________________
void AliMUONRawClusterV2::SetDigitsId(Int_t nDigits, const UInt_t *digitsId)
{
  /// Set size of array of digits Id to n ints and set the content
  /// if digitsId is not given the array is filled with id=0
  
  if (fDigitsId && fNDigits != nDigits) {
    delete [] fDigitsId;
    fDigitsId = 0;
  }
  fNDigits = nDigits;
  if (fNDigits == 0) return;
  if (!fDigitsId) fDigitsId = new UInt_t[fNDigits];
  if (digitsId == 0)
    for (Int_t i=0; i<fNDigits; i++) fDigitsId[i] = 0;
  else
    memcpy(fDigitsId,digitsId, fNDigits*sizeof(UInt_t));
}

//____________________________________________________
void AliMUONRawClusterV2::AddDigitId(UInt_t id)
{
  /// Reset size of array of digits Id and add the new id to its content
  
  UInt_t *digitsIdNew = new UInt_t[fNDigits+1];
  memcpy(digitsIdNew,fDigitsId, fNDigits*sizeof(UInt_t));
  digitsIdNew[fNDigits++] = id;
  delete[] fDigitsId;
  fDigitsId = digitsIdNew;
}

//____________________________________________________
Int_t AliMUONRawClusterV2::Compare(const TObject *obj) const
{
/// Compare

  const AliMUONRawClusterV2* raw = static_cast<const AliMUONRawClusterV2*>(obj);
  if ( GetCharge() > raw->GetCharge() ) 
  {
    return 1;
  }
  else if ( GetCharge() < raw->GetCharge() ) 
  {
    return -1;
  }
  return 0;
}
