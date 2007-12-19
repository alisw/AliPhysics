/* Copyright (C) 2007 Christian Holm Christensen <cholm@nbi.dk>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
 * USA
 */
//____________________________________________________________________
//
// Object used by AliFMDFlowBinned1D to split an event into sub-events.
// Default is to split randomly.  
// User defined derived classes can do other stuff.
//
#include "flow/AliFMDFlowSplitter.h"
#include "flow/AliFMDFlowAxis.h"
#include <TRandom.h>
#include <cmath>
#include <iostream>

//____________________________________________________________________
Bool_t
AliFMDFlowSplitter::Select(ULong_t) const
{
  // Decide whether entry should go in A or B sub-event. 
  // Parameters: 
  //   entry	The entry number 
  // Return 
  //   true if this should go in sub-event A 
  return (Float_t(rand()) / RAND_MAX > 0.5);
}


//____________________________________________________________________
void
AliFMDFlowShuffle::Event(Double_t*, Double_t*, ULong_t n)
{
  // Prepare for an event 
  // Parameters 
  //   phis List of phis. 
  //   xs   List of bin variable 
  //   n    Number of entries in @a phis and @a n 
  // Return 
  fN = n;
  Shuffle();
}

//____________________________________________________________________
Bool_t
AliFMDFlowShuffle::Select(ULong_t entry) const
{
  // Decide whether entry should go in A or B sub-event. 
  // Parameters: 
  //   entry	The entry number 
  // Return 
  //   true if this should go in sub-event A 
  if (entry >= fN || Int_t(entry) >= fIdx.fN) return false;
  ULong_t n = fIdx[entry];
  return (n < fN/2);
}

//____________________________________________________________________
void
AliFMDFlowShuffle::Shuffle()
{
  // Suffle index 
  if (fIdx.fN < Int_t(fN)) fIdx.Set(fN);
  for (ULong_t i = 0; i < fN; i++) fIdx[i] = i;
  for (ULong_t i = 0; i < fN; i++) { 
    // Swap 2 random locations 
    ULong_t j = ULong_t(gRandom->Rndm()*(fN-1));
    ULong_t k = fIdx[j];
    fIdx[j]   = fIdx[i];
    fIdx[i]   = k;
  }
}

//____________________________________________________________________
//
// EOF
//
