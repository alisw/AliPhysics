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

/*
  $Log$
*/

////////////////////////////////////////////////
//  Digit class for set CPV                   //
//                                            //
//  Author: Yuri Kharlov, IHEP, Protvino      //
//  e-mail: Yuri.Kharlov@cern.ch              //
//  Last modified: 2 November 2000            //
////////////////////////////////////////////////
 
// --- aliroot header files ---
#include "AliPHOSCPVDigit.h"

ClassImp(AliPHOSCPVDigit)

//______________________________________________________________________________

AliPHOSCPVDigit::AliPHOSCPVDigit(Int_t x, Int_t y, Float_t q)
{
  //
  // Create a CPV digit object
  //

  fXpad = x;
  fYpad = y;
  fQpad = q;
}

//______________________________________________________________________________
