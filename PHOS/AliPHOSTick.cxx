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

//_________________________________________________________________________
// Auxiliary class to help calculate the time of crossing 
// of the threshold by the front edge of the time signal
//
//*-- Author :  Dmitri Peressounko (SUBATECH) 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

// --- Standard library ---

// --- AliRoot header files ---
#include "AliPHOSTick.h"

ClassImp(AliPHOSTick)


//____________________________________________________________________________ 
AliPHOSTick::AliPHOSTick():
  fTime(0),
  fA(0),
  fB(0)
{
}

//____________________________________________________________________________ 
AliPHOSTick::AliPHOSTick(Float_t time, Float_t a, Float_t slope):
  fTime(time),
  fA(a),
  fB(slope)
{
}

//____________________________________________________________________________ 
Int_t AliPHOSTick::Compare(const TObject * obj) const {
  if(obj->InheritsFrom("AliPHOSTick")){
    AliPHOSTick * tick = (AliPHOSTick *) obj ;
    if(fTime < tick->fTime)
      return -1 ;
    else
      if(fTime == tick->fTime)
	return 0 ;
      else
	return 1 ;
  }
  else
    return 1 ;
} 
//____________________________________________________________________________
void AliPHOSTick::operator+=(AliPHOSTick const & tick) 
{
  // Adds the amplitude of digits and completes the list of primary particles
  // if amplitude is larger than 
    
  fA = fA + fB*(tick.fTime - fTime) + tick.fA ;
  fB = fB + tick.fB ;
  if(tick.fTime > fTime) 
    fTime = tick.fTime ;
  
}
