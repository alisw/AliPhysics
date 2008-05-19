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
// $MpId: AliMpTrigger.cxx,v 1.4 2006/05/24 13:58:52 ivana Exp $

//-----------------------------------------------------------------------------
// Class AliMUONTriggerCrateConfig
// --------------------
// The class defines the configuration of trigger crate
// Author: Ch. Finck, Subatech Nantes
//-----------------------------------------------------------------------------

#include "AliMUONTriggerCrateConfig.h"

#include "AliLog.h"

#include <Riostream.h>

/// \cond CLASSIMP
ClassImp(AliMUONTriggerCrateConfig)
/// \endcond


 //______________________________________________________________________________
AliMUONTriggerCrateConfig::AliMUONTriggerCrateConfig()
  : TNamed("Trigger Crate","configuration trigger crate"),
    fId(0),
    fMask(0),
    fMode(0),
    fCoinc(0),
    fLocalBoard(false)
{
/// Standard constructor for Shuttle + DA
}


 //______________________________________________________________________________
AliMUONTriggerCrateConfig::AliMUONTriggerCrateConfig(const Char_t* name, UShort_t id, UShort_t mask, UShort_t mode, UShort_t coinc)
  : TNamed(name, "configuration trigger crate"),
    fId(id),
    fMask(mask),
    fMode(mode),
    fCoinc(coinc),
    fLocalBoard(false)
{
/// Standard constructor for Shuttle + DA
}


//______________________________________________________________________________
AliMUONTriggerCrateConfig::~AliMUONTriggerCrateConfig()
{
/// Destructor
}

//
// public methods
//

//______________________________________________________________________________
Bool_t AliMUONTriggerCrateConfig::AddLocalBoard(Int_t localBoardId)
{
/// Add detection element with given detElemId.
/// Return true if the detection element was added

  if ( HasLocalBoard(localBoardId) ) {
    AliWarningStream() 
      << "Local board with Id=" << localBoardId << " already present."
      << endl;
    return false;
  }    

  fLocalBoard.Add(localBoardId);
  return true;
}   


//______________________________________________________________________________
Int_t AliMUONTriggerCrateConfig::GetNofLocalBoards() const
{  
/// Return the number of local board in this crate

  return fLocalBoard.GetSize(); 
}

//______________________________________________________________________________
Int_t  AliMUONTriggerCrateConfig::GetLocalBoardId(Int_t index) const
{  
/// Return the local board by index (in loop)

  if (index >= 0 && index < fLocalBoard.GetSize())
      return fLocalBoard.GetValue(index); 
  else 
      return 0; // begin at 1
}

//______________________________________________________________________________
Bool_t  AliMUONTriggerCrateConfig::HasLocalBoard(Int_t localBoardId) const
{  
/// Return true if crate has local boardwith given localBoardId

  return fLocalBoard.HasValue(localBoardId); 
}

