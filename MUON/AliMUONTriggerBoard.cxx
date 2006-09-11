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

//*-- Author: Rachid Guernane (LPCCFd)
//    TRIGGER BOARD SUPER CLASS IMPLEMENTATION
//    CAN BE A LOCAL, REGIONAL, OR GLOBAL BOARD
//    REGIONAL BOARD IS PER CONVENTION
//    ALWAYS IN THE SLOT 0

#include "AliMUONTriggerBoard.h"
#include "AliLog.h"

//___________________________________________
AliMUONTriggerBoard::AliMUONTriggerBoard()
    : TNamed(),
      fSlot(0),
      fResponse(0)
{
// Default Ctor

}

//___________________________________________
AliMUONTriggerBoard::AliMUONTriggerBoard(const char *Name, Int_t islot) 
    : TNamed(Name,"Trigger board"),
      fSlot(islot),
      fResponse(0)
{
// Standard Ctor

}

//___________________________________________
AliMUONTriggerBoard::AliMUONTriggerBoard(const AliMUONTriggerBoard &board)
    : TNamed(board),
      fSlot(board.fSlot),
      fResponse(board.fResponse)
{
// Dummy Copy Ctor
//   board.Copy(*this);
}

//___________________________________________
AliMUONTriggerBoard& AliMUONTriggerBoard::operator=(const AliMUONTriggerBoard&)
{
   AliFatal("Assignment operator not implemented");
   return *this;
}

//___________________________________________
void AliMUONTriggerBoard::Copy(TObject&) const
{
   Fatal("Copy","Not implemented!\n");
}

ClassImp(AliMUONTriggerBoard)


