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

#include "AliMUONVTrackerData.h"

/// \class AliMUONVTrackerData
///
/// Base class for MUON data that can be presented at different levels
/// in the hierarchy of the MUON system.
/// 
/// We always feed the AliMUONVTrackerData with data at the channel level,
/// and it then computes the same data at upper levels, such as manu, pcb,
/// bus patches, detection elements, and even chamber wise.
///
/// The dimension (or dim) parameter that appears in many methods is the 
/// "number of parameters" per channel.
///
/// \author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONVTrackerData)
///\endcond

//_____________________________________________________________________________
AliMUONVTrackerData::AliMUONVTrackerData(const char* name, const char* title,
                                         Bool_t)
: TNamed(name,title), TQObject()
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONVTrackerData::~AliMUONVTrackerData()
{
  /// dtor
}

//_____________________________________________________________________________
Bool_t 
AliMUONVTrackerData::HasChannel(Int_t detElemId, Int_t manuId, Int_t manuChannel) const
{
  /// Whether we have data for a given channel
  
  return (Count(detElemId,manuId,manuChannel) > 0.0);
}

//_____________________________________________________________________________
void 
AliMUONVTrackerData::NumberOfEventsChanged()
{
  /// Signal that our number of events changed
  Emit("NumberOfEventsChanged()");
}

//_____________________________________________________________________________
void
AliMUONVTrackerData::Print(Option_t* wildcard) const
{
  /// Printout
  Print(wildcard,"summary");
}

//_____________________________________________________________________________
Bool_t
AliMUONVTrackerData::Replace(const AliMUONVStore& /*store*/)
{
  Emit("Replace(const AliMUONVStore&)");
  return kTRUE;
}
