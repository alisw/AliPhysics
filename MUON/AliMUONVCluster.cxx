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

// $Id$

//-----------------------------------------------------------------------------
/// \class AliMUONVCluster
///
/// An abstract base class for clusters
///
/// \author Philippe Pillot, Subatech
//-----------------------------------------------------------------------------

#include "AliMUONVCluster.h"

#include "AliLog.h"

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMUONVCluster)
/// \endcond

//_____________________________________________________________________________
AliMUONVCluster::AliMUONVCluster()
{
  /// default constructor
}

//_____________________________________________________________________________
AliMUONVCluster::AliMUONVCluster(Int_t chamberId, Int_t detElemId, Int_t clusterIndex)
  : TObject() 
{
  /// constructor
  SetUniqueID(BuildUniqueID(chamberId, detElemId, clusterIndex));
}

//_____________________________________________________________________________
AliMUONVCluster::~AliMUONVCluster()
{
  /// destructor
}

//_____________________________________________________________________________
void AliMUONVCluster::Print(Option_t *option) const
{
  /// print cluster content
  /// if option=FULL print also all Digit ID
  UInt_t cId = GetUniqueID();
  Int_t nDigits = GetNDigits();
  
  cout<<Form("clusterID=%u (ch=%d, det=%d, index=%d)",
	     cId,GetChamberId(),GetDetElemId(),GetClusterIndex(cId))<<endl;
  
  cout<<Form("position=(%5.2f, %5.2f, %5.2f), sigma=(%5.2f, %5.2f, 0.0), charge=%5.2f, chi2=%5.2f, MClabel=%d",
	     GetX(),GetY(),GetZ(),GetErrX(),GetErrY(),GetCharge(),GetChi2(),GetMCLabel())<<endl;
  
  if (strcmp(option,"FULL") == 0) {
    cout<<"nDigits="<<nDigits<<" digitID=(";
    for (Int_t i=0; i<nDigits; i++) cout<<GetDigitId(i)<<", ";
    cout<<")"<<endl;
  }
  
}
