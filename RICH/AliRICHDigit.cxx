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
  Printf("pad=(%2i,%2i,%3i,%3i), QDC=%8.3f, cfm=%9i, TID=(%5i,%5i,%5i)",
               C(),Pad2Sec(PadX(),PadY()),PadX(),PadY() ,  Qdc(),    Cfm()  , fTracks[0],fTracks[1],fTracks[2]);
}
//__________________________________________________________________________________________________
void AliRICHDigit::Test()
{
  Printf("Test of Pad2Sec:");
  Int_t x1=kFirstPad, x2=kPadsSecX, x3=kPadsSecX+1, x4=kPadsChamX;//all possible padx for corners
  Int_t y;
  
  y=kPadsChamY;
  Printf("Sector5:(%3i,%3i)->%i (%3i,%3i)->%i  Sector6:(%3i,%3i)->%i (%3i,%3i)->%i",x1,y,Pad2Sec(x1,y),x2,y,Pad2Sec(x2,y),x3,y,Pad2Sec(x3,y),x4,y,Pad2Sec(x4,y));
  y=2*kPadsSecY+1;
  Printf("Sector5:(%3i,%3i)->%i (%3i,%3i)->%i  Sector6:(%3i,%3i)->%i (%3i,%3i)->%i",x1,y,Pad2Sec(x1,y),x2,y,Pad2Sec(x2,y),x3,y,Pad2Sec(x3,y),x4,y,Pad2Sec(x4,y));
  Printf("");
  y=2*kPadsSecY;
  Printf("Sector3:(%3i,%3i)->%i (%3i,%3i)->%i  Sector4:(%3i,%3i)->%i (%3i,%3i)->%i",x1,y,Pad2Sec(x1,y),x2,y,Pad2Sec(x2,y),x3,y,Pad2Sec(x3,y),x4,y,Pad2Sec(x4,y));
  y=kPadsSecY+1;
  Printf("Sector3:(%3i,%3i)->%i (%3i,%3i)->%i  Sector4:(%3i,%3i)->%i (%3i,%3i)->%i",x1,y,Pad2Sec(x1,y),x2,y,Pad2Sec(x2,y),x3,y,Pad2Sec(x3,y),x4,y,Pad2Sec(x4,y));
  Printf("");
  y=kPadsSecY;
  Printf("Sector1:(%3i,%3i)->%i (%3i,%3i)->%i  Sector2:(%3i,%3i)->%i (%3i,%3i)->%i",x1,y,Pad2Sec(x1,y),x2,y,Pad2Sec(x2,y),x3,y,Pad2Sec(x3,y),x4,y,Pad2Sec(x4,y));
  y=kFirstPad;
  Printf("Sector1:(%3i,%3i)->%i (%3i,%3i)->%i  Sector2:(%3i,%3i)->%i (%3i,%3i)->%i",x1,y,Pad2Sec(x1,y),x2,y,Pad2Sec(x2,y),x3,y,Pad2Sec(x3,y),x4,y,Pad2Sec(x4,y));
}
