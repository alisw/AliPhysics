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

///////////////////////////////////////////////////////////////////////////////
///                                                                          //
///  Time Projection Chamber Digit                                           //
///                                                                          //
///////////////////////////////////////////////////////////////////////////////

#include "AliTPCdigit.h"


ClassImp(AliTPCdigit)
  //_____________________________________________________________________
  AliTPCdigit::AliTPCdigit()
              :AliDigit(),
	       fSector(0),
               fPadRow(0),
               fPad(0),
               fTime(0),
               fSignal(0)
{
  //
  //   default constructor
  //
} 
//_____________________________________________________________________________
AliTPCdigit::AliTPCdigit(Int_t *tracks, Int_t *digits)
            :AliDigit(tracks),
	       fSector(0),
               fPadRow(0),
               fPad(0),
	       fTime(0),
               fSignal(0)
{
  //
  // Creates a TPC digit object
  //
  fSector     = digits[0];
  fPadRow     = digits[1];
  fPad        = digits[2];
  fTime       = digits[3];
  fSignal     = digits[4];
}

