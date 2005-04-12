//  **************************************************************************
//  * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
//  *                                                                        *
//  * Author: The ALICE Off-line Project.                                    *
//  * Contributors are mentioned in the code where appropriate.              *
//  *                                                                        *
//  * Permission to use, copy, modify and distribute this software and its   *
//  * documentation strictly for non-commercial purposes is hereby granted   *
//  * without fee, provided that the above copyright notice appears in all   *
//  * copies and that both the copyright notice and this permission notice   *
//  * appear in the supporting documentation. The authors make no claims     *
//  * about the suitability of this software for any purpose. It is          *
//  * provided "as is" without express or implied warranty.                  *
//  **************************************************************************

#include "AliRICHDigit.h"
#include <AliLog.h>

ClassImp(AliRICHDigit)

//__________________________________________________________________________________________________
void AliRICHDigit::Print(Option_t*)const
{
//Print current digit  
  AliInfo(Form("cfm=%9i, cs=%2i, x=%3i, y=%3i, q=%8.3f, TID1=%5i, TID2=%5i, TID3=%5i",
                  fCFM,fChamber,fPadX,fPadY,fQdc,fTracks[0],fTracks[1],fTracks[2]));
}
//__________________________________________________________________________________________________
