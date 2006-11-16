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

#include "AliMUONVClusterFinder.h"
#include "AliLog.h"

/// \class AliMUONVClusterFinder
///
/// Defines an interface for a cluster finder.
/// 
/// A cluster finder is supposed to work on a single detection element at a
/// time, thus the Prepare function (which sets up the cluster finder for a
/// particular DE)
///
/// \author Laurent Aphecetche

/// \cond CLASSIMP
ClassImp(AliMUONVClusterFinder)
/// \endcond

//_____________________________________________________________________________
AliMUONVClusterFinder::AliMUONVClusterFinder() : TObject()
{
}

//_____________________________________________________________________________
AliMUONVClusterFinder::~AliMUONVClusterFinder()
{
}

//_____________________________________________________________________________
Bool_t
AliMUONVClusterFinder::UsePad(const AliMUONPad&)
{
  AliError("Not implemented");
  return kFALSE;
}
