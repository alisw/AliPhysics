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
///
/// This is a service class for packing and unpacking bits in a 32 bit word.
///
///////////////////////////////////////////////////////////////////////////////

#include <TError.h>

#include "AliBitPacking.h"


ClassImp(AliBitPacking)


//_____________________________________________________________________________
Bool_t AliBitPacking::PackWord(UInt_t data, UInt_t &word, 
			       Int_t startBit, Int_t stopBit)
{
// Packs data into the word buffer from startBit bit up to stopBit bit

  // check that startBit and stopBit are reasonable
  if (startBit > stopBit) {
    ::Error("AliBitPacking::PackWord", 
	    "startBit is larger than stopBit");
    return kFALSE;
  }
  if (stopBit > 31) {
    ::Error("AliBitPacking::PackWord", 
	    "stopBit exceeds valid range of 32 bits");
    return kFALSE;
  }

  // create a word with the bits 0 to (stopBit-startBit) set
  UInt_t bits = 0xFFFFFFFF;
  if (stopBit-startBit < 31) bits = (1 << (stopBit-startBit+1)) - 1;

  // check that the data fits into the given bit range
  if (data > bits){
    ::Error("AliBitPacking::PackWord", 
	    "Word to be filled is not within desired length");
    return kFALSE;
  }

  // clear the bits from startBit to stopBit
  word &= (0xFFFFFFFF ^ (bits << startBit));

  // fill in the data bits
  word |= (data << startBit);

  return kTRUE;
}

//_____________________________________________________________________________
UInt_t AliBitPacking::UnpackWord(UInt_t word, Int_t startBit, Int_t stopBit)
{
// Unpacks data of stopBit-startBit+1 bits from word buffer starting from 
// the position indicated by startBit

  // check that startBit and stopBit are reasonable
  if (startBit > stopBit) {
    ::Error("AliBitPacking::UnpackWord", 
	    "startBit is larger than stopBit");
    return 0xFFFFFFFF;
  }
  if (stopBit > 31) {
    ::Error("AliBitPacking::UnpackWord", 
	    "stopBit exceeds valid range of 32 bits");
    return 0xFFFFFFFF;
  }

  // create a word with the bits 0 to (stopBit-startBit) set
  UInt_t bits = 0xFFFFFFFF;
  if (stopBit-startBit < 31) bits = (1 << (stopBit-startBit+1)) - 1;

  // pick out the requested bits
  return ((word >> startBit) & bits);
}
