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
///  \class AliMUONTriggerBoard
///
///  Trigger board super class implementation.
///  Can be a local, regional, or global board
///  Regional board is per convention always in the slot 0
///
///  \author Rachid Guernane (LPCCFd)
//-----------------------------------------------------------------------------

#include "AliMUONTriggerBoard.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliMUONTriggerBoard)
/// \endcond

//___________________________________________
AliMUONTriggerBoard::AliMUONTriggerBoard()
    : TNamed(),
      fSlot(0),
      fResponse(0)
{
/// Default Ctor

}

//___________________________________________
AliMUONTriggerBoard::AliMUONTriggerBoard(const char *Name, Int_t islot) 
    : TNamed(Name,"Trigger board"),
      fSlot(islot),
      fResponse(0)
{
/// Standard Ctor

}


//___________________________________________
AliMUONTriggerBoard::AliMUONTriggerBoard(const AliMUONTriggerBoard &rhs) 
    : TNamed(rhs),
      fSlot(rhs.fSlot),
      fResponse(rhs.fResponse)
{
  //
  /// Copy constructor
  //
}


//___________________________________________
AliMUONTriggerBoard& AliMUONTriggerBoard::operator=(const AliMUONTriggerBoard &rhs)
{
/// Assigment operator;

  if (this == &rhs)
    return *this;

  // base class assignement
  TNamed::operator=(rhs);

  fSlot = rhs.fSlot;
  fResponse = rhs.fResponse;

  return *this;
}


//___________________________________________
AliMUONTriggerBoard::~AliMUONTriggerBoard()
{
/// Destructor
}  




