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

#include "AliRICHDigit.h"//class header

ClassImp(AliRICHDigit)

//__________________________________________________________________________________________________
void AliRICHDigit::Print(Option_t*)const
{
//Print current digit  
//Arguments: option string not used
//  Returns: none    
  Printf("pad=(%2i,%2i,%3i,%3i), q=%8.3f, cfm=%9i, TID=(%5i,%5i,%5i)",
               Chamber(),Sector(),PadX(),PadY() ,  Qdc(),    Cfm()  , fTracks[0],fTracks[1],fTracks[2]);
}
//__________________________________________________________________________________________________
void AliRICHDigit::Test()
{
}
