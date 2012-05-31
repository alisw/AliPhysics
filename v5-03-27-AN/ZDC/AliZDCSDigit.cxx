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

///_________________________________________________________________________
///
///
///   ZDC summable digit = Cerenkov light for each PM 
///
///_________________________________________________________________________

#include "AliZDCSDigit.h"


ClassImp(AliZDCSDigit)

//____________________________________________________________________________
AliZDCSDigit::AliZDCSDigit() :
  fLightPM(0),
  fTrackTime(0)
{
  // Default constructor 
  
  fSector[0]   = 0;
  fSector[1]   = 0;
}

//____________________________________________________________________________
AliZDCSDigit::AliZDCSDigit(Int_t* sector, Float_t lightPM, Float_t trackTime) :
  fLightPM(lightPM),
  fTrackTime(trackTime)
{  
  // Constructor 
 
  fSector[0] = sector[0];
  fSector[1] = sector[1];
}
