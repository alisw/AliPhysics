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

/// \class AliMUONVClusterServer
/// 
/// Interface of a cluster finder for combined tracking
/// 
/// The tracking, when in needs for clusters, will ask the cluster server
/// to add clusters, from a given region in a given chamber, to an existing 
/// clusterstore, using the 
///
/// Clusterize(Int_t chamberId, AliMUONVClusterStore& clusterStore, const AliMpArea& area)
///
/// method
///
/// Cluster server must be instructed (at the beginning of each event) 
/// about which digits to use, using the 
///
/// UseDigitStore(const AliMUONVDigitStore& digitStore)
///
/// method.
///
/// \author Laurent Aphecetche, Subatech

#include "AliMUONVClusterServer.h"

/// \cond CLASSIMP  
ClassImp(AliMUONVClusterServer)
/// \endcond

//_____________________________________________________________________________
AliMUONVClusterServer::AliMUONVClusterServer()
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONVClusterServer::~AliMUONVClusterServer()
{
  /// dtor
}
