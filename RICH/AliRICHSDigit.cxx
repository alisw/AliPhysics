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
  Revision 1.2  2000/12/18 17:45:17  jbarbosa
  Cleaned up PadHits object.

  Revision 1.1  2000/06/12 15:25:31  jbarbosa
  Cleaned up version.

*/


#include "AliRICHSDigit.h"

ClassImp(AliRICHSDigit)
//______________________________________________________________
AliRICHSDigit::AliRICHSDigit(Int_t *clhits)
{

// Default constructor for AliRICHSDigits

  fHitNumber=clhits[0];
  fQpad=clhits[1];
  fPadX=clhits[2];
  fPadY=clhits[3];
  fRSec=clhits[4];
}

