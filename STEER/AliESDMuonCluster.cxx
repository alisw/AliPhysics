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
/// \class AliESDMuonCluster
///
/// Class to describe the MUON clusters in the Event Summary Data
///
/// \author Philippe Pillot, Subatech
//-----------------------------------------------------------------------------

#include "AliESDMuonCluster.h"

#include "AliLog.h"

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliESDMuonCluster)
/// \endcond

//_____________________________________________________________________________
AliESDMuonCluster::AliESDMuonCluster()
: TObject()
{
  /// default constructor
  fXYZ[0] = fXYZ[1] = fXYZ[2] = 0.;
  fErrXY[0] = fErrXY[1] = 0.;
}

//_____________________________________________________________________________
AliESDMuonCluster::AliESDMuonCluster (const AliESDMuonCluster& cluster)
: TObject(cluster)
{
  /// Copy constructor
  fXYZ[0] = cluster.fXYZ[0];
  fXYZ[1] = cluster.fXYZ[1];
  fXYZ[2] = cluster.fXYZ[2];
  fErrXY[0] = cluster.fErrXY[0];
  fErrXY[1] = cluster.fErrXY[1];
}

//_____________________________________________________________________________
AliESDMuonCluster& AliESDMuonCluster::operator=(const AliESDMuonCluster& cluster)
{
  /// Equal operator
  if (this == &cluster) return *this;
  
  TObject::operator=(cluster); // don't forget to invoke the base class' assignment operator
  
  fXYZ[0] = cluster.fXYZ[0];
  fXYZ[1] = cluster.fXYZ[1];
  fXYZ[2] = cluster.fXYZ[2];
  fErrXY[0] = cluster.fErrXY[0];
  fErrXY[1] = cluster.fErrXY[1];
  
  return *this;
}

//_____________________________________________________________________________
void AliESDMuonCluster::Print(Option_t */*option*/) const
{
  /// print cluster content
  UInt_t cId = GetUniqueID();
  
  cout<<Form("clusterID=%u (ch=%d, det=%d, index=%d)",
	     cId,GetChamberId(),GetDetElemId(),GetClusterIndex())<<endl;
  
  cout<<Form("position=(%5.2f, %5.2f, %5.2f), sigma=(%5.2f, %5.2f, 0.0)",
	     GetX(),GetY(),GetZ(),GetErrX(),GetErrY())<<endl;
}

