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

#include "AliMpManuIterator.h"

#include "AliMpBusPatch.h"
#include "AliMpDDLStore.h"
#include "TExMap.h"
#include "AliLog.h"

/// \class AliMpManuIterator
///
/// Class to loop over all manus of MUON Tracker
/// 
/// \author Laurent Aphecetche, Subatech

/// \cond CLASSIMP
ClassImp(AliMpManuIterator)
/// \endcond

//_____________________________________________________________________________
AliMpManuIterator::AliMpManuIterator()
: TObject(), 
fIterator(AliMpDDLStore::Instance()->CreateBusPatchIterator()),
fCurrentBusPatch(0x0),
fCurrentManuIndex(-1)
{
  /// ctor
  Reset();
}

//_____________________________________________________________________________
AliMpManuIterator::~AliMpManuIterator()
{
  /// dtor
  delete fIterator;
}

//_____________________________________________________________________________
Bool_t
AliMpManuIterator::Next(Int_t& detElemId, Int_t& manuId)
{
  /// Set the next (de,manu) pair and return kTRUE, or kFALSE if ended.
  
  ++fCurrentManuIndex;
  
  if ( fCurrentManuIndex < fCurrentBusPatch->GetNofManus() ) 
  {
    detElemId = fCurrentBusPatch->GetDEId();
    manuId = fCurrentBusPatch->GetManuId(fCurrentManuIndex);
    return kTRUE;
  }
  else
  {
    fCurrentBusPatch = static_cast<AliMpBusPatch*>(fIterator->Next());
    if (!fCurrentBusPatch ) 
    {
      return kFALSE;
    }
    fCurrentManuIndex = -1;
    return Next(detElemId,manuId);
  }
}

//_____________________________________________________________________________
void
AliMpManuIterator::Reset()
{
  /// Rewind the iterator
  fIterator->Reset();
  
  fCurrentBusPatch = static_cast<AliMpBusPatch*>(fIterator->Next());
  
  fCurrentManuIndex = -1;
}
