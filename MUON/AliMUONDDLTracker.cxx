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
 
#include "AliMUONDDLTracker.h"
#include "AliRawDataHeader.h"

ClassImp(AliMUONDDLTracker)
 const Int_t AliMUONDDLTracker::fgkBlkHeaderLength = 8;
 const Int_t AliMUONDDLTracker::fgkDspHeaderLength = 8;
 const Int_t AliMUONDDLTracker::fgkEndOfDDL = 0x0FFFFFFFF;

//___________________________________________
AliMUONDDLTracker::AliMUONDDLTracker()
  :  TObject(),
     fTotalBlkLength(0),
     fBlkLength(0),
     fDSPId(0),
     fPadding(0x0DEADDEAD),
     fTotalDspLength(0),
     fDspLength(0),
     fDSPId1(0),
     fEventWord(0) 
{
  //ctor
  for (Int_t i = 0; i < 4; i++)
    fBlkTriggerWord[i] = fDspTriggerWord[i] = 0;
}
