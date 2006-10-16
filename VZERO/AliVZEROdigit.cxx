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

#include "AliVZEROdigit.h"

ClassImp(AliVZEROdigit)

AliVZEROdigit::AliVZEROdigit()
{
  fPMNumber   = 0;
  fADC        = 0;
  fTime       = 0;
}

AliVZEROdigit::AliVZEROdigit(Int_t* tracks, Int_t *digits):
  AliDigit(tracks){
  
  //
  // Creates VZERO digits  
  // Fills the AliVZEROdigit data members from the array digits. 
  //
  
  fTrack      =  tracks[0];
  fEvent      =  digits[0];  
}

AliVZEROdigit::AliVZEROdigit(Int_t PMnumber, Int_t adc, Int_t time)
{  
  fPMNumber   = PMnumber;
  fADC        = adc;
  fTime       = time;
}



