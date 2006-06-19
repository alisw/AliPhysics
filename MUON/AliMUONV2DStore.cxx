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

#include "AliMUONV2DStore.h"

/// \class AliMUONV2DStore
/// Defines an interface equivalent to a 2D array of TObject, indexed
/// by a pair of integers (somehow a matrix, 
/// except that indices are not necessarily sequential).
/// 
/// It's extremely simple and hopefully allow many implementations.
/// It also makes the object ownership self-evident.
///
/// \author Laurent Aphecetche

/// \cond CLASSIMP
ClassImp(AliMUONV2DStore)
/// \endcond

//_____________________________________________________________________________
AliMUONV2DStore::AliMUONV2DStore()
{
/// Default constructor
}

//_____________________________________________________________________________
AliMUONV2DStore::~AliMUONV2DStore()
{
/// Destructor
}




