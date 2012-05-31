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
/// \class AliMUONTriggerCrate
///  Collection of trigger boards
///  - one regional
///  - sixteen local
///  slot 0 holds the regional board
/// \author Rachid Guernane (LPCCFd)
//-----------------------------------------------------------------------------

#include <TObjArray.h>

#include "AliMUONTriggerBoard.h"
#include "AliMUONTriggerCrate.h"

ClassImp(AliMUONTriggerCrate)

//___________________________________________
AliMUONTriggerCrate::AliMUONTriggerCrate()
    : fNslots(0),
      fNboards(0),
      fBoards(0x0),
      fSourceFileName(0)    
{
/// Default constructor
}

//___________________________________________
AliMUONTriggerCrate::~AliMUONTriggerCrate()
{
/// Destructor

  delete fBoards;
}

//___________________________________________
AliMUONTriggerCrate::AliMUONTriggerCrate(const char *name, Int_t n) : 
    TNamed(name,"Regional trigger crate"),
    fNslots(n),
    fNboards(0),
    fBoards(new TObjArray(fNslots)),
    fSourceFileName(0)
{
  /// Standard constructor
  fBoards->SetOwner(kTRUE);
}

//___________________________________________
void AliMUONTriggerCrate::AddBoard(AliMUONTriggerBoard *board, Int_t i)
{
/// Add board in crate container
   fBoards->AddAt(board,i);
   fNboards++;
}

