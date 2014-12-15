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

/// \class AliMUONTrackerDataWrapper
///
/// A simple wrapper to convert an AliMUONVTrackerData object into
/// an AliMUONVTrackerDataMaker object.
///
/// This is mainly to offer backward compatibility : the mchview program
/// used to save AliMUONVTrackerData objects, while it now saves 
/// AliMUONVTrackerDataMaker ones.
/// So to read back old files, we need to be able to do the "conversion".
///
/// \author Laurent Aphecetche, Subatech

#include "AliMUONTrackerDataWrapper.h"

#include "AliLog.h"
#include "AliMUONVTrackerData.h"

/// \cond CLASSIMP
ClassImp(AliMUONTrackerDataWrapper)
/// \endcond

//_____________________________________________________________________________
AliMUONTrackerDataWrapper::AliMUONTrackerDataWrapper(AliMUONVTrackerData* data)
: AliMUONVTrackerDataMaker(), fData(data)
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONTrackerDataWrapper::~AliMUONTrackerDataWrapper()
{
  /// dtor
  delete fData;
}

//_____________________________________________________________________________
Long64_t 
AliMUONTrackerDataWrapper::Merge(TCollection*)
{
  /// Merge
  AliError("Not implemented yet");
  return 0;
}

//_____________________________________________________________________________
Int_t
AliMUONTrackerDataWrapper::NumberOfEvents() const
{
  /// Get the number of events the data has seen
  if ( Data() ) 
  {
    return Data()->NumberOfEvents(-1);
  }
  return 0;
}
