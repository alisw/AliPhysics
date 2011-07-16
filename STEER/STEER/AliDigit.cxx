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

//------------------------------------------------------------------------
//             Main digit class
//  Base class for Alice digits
//  Author:
//------------------------------------------------------------------------

#include "AliDigit.h"
 
ClassImp(AliDigit)

AliDigit::AliDigit()
{
  //
  // Default constructor
  //
  fTracks[0]=fTracks[1]=fTracks[2]=-1;  
}

AliDigit::AliDigit(Int_t *tracks)
{
  //
  // Standard constructor
  //
  fTracks[0] = tracks[0];
  fTracks[1] = tracks[1];
  fTracks[2] = tracks[2];
}

	 
