#ifndef ALIMPENCODEPAIR_H
#define ALIMPENCODEPAIR_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

#include <Rtypes.h>
#include <Riosfwd.h>

typedef Int_t MpPair_t;

namespace AliMp
{
  /// Encode the pair of integers to another integer.
  MpPair_t Pair(Int_t first, Int_t second);
  
  /// Decode the first integer from encoded pair
  Int_t  PairFirst(MpPair_t pair);

  /// Decode the second integer from encoded pair
  Int_t  PairSecond(MpPair_t pair);

  /// A special printing for encoded pair.
  ostream& PairPut(ostream& s, MpPair_t pair);
}

#endif //ALIMPENCODEPAIR_H
