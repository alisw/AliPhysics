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

/* $Id: $ */

#include "AliCaloFitSubarray.h"


// Container class to hold info from bunches/samples
// selected for signal fitting.
// Variables are:
//  Int_t   fBunchIndex;  // Index for selected bunch
//  Int_t   fMaxRev;      // Max index in reversed array
//  Int_t   fFirst;   // first index in array used for fit
//  Int_t   fLast;    // last index in array used for fit

AliCaloFitSubarray::AliCaloFitSubarray(const Int_t bunchIndex, 
				       const Int_t maxrev, 
				       const Int_t first, 
				       const Int_t last ) : 
  fBunchIndex(bunchIndex),
  fMaxRev(maxrev), 
  fFirst(first), 
  fLast(last)   
{
}

AliCaloFitSubarray::AliCaloFitSubarray(const Int_t init) : 
  fBunchIndex(init),
  fMaxRev(init), 
  fFirst(init), 
  fLast(init)   
{
}

AliCaloFitSubarray::AliCaloFitSubarray(const AliCaloFitSubarray & fitS) :
  fBunchIndex( fitS.fBunchIndex ),
  fMaxRev( fitS.fMaxRev ), 
  fFirst( fitS.fFirst ), 
  fLast( fitS.fLast )   
{
}

AliCaloFitSubarray::~AliCaloFitSubarray()
{
}

