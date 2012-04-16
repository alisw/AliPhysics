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

/* 
 * Class for reading TRU data from a bunch from a raw datastream.
 * Author: Henrik Qvigstad <henrik.qvigstad@cern.ch>
 */

#include "AliPHOSTRURawReader.h"

#include "AliCaloRawStreamV3.h"
#include <bitset> // std::bitset

ClassImp(AliPHOSTRURawReader)


//________________________________________________________________
AliPHOSTRURawReader::AliPHOSTRURawReader()
  : TObject(),
    fSignals(),
    fFlags(),
    fActive(0),
    fActiveTime(),
    fHasSignal(false),
    fHasSignalTime()

{
  // default constructor
  
  // fSignals Initialization:
  for(Int_t row = 0; row < kN2x2XPrTRURow; ++row) {
    for(Int_t branch = 0; branch < kN2x2ZPrBranch; ++branch) {
      for(Int_t timeBin = 0; timeBin < kNTimeBins; ++timeBin) {
	fSignals[row][branch][timeBin] = kDefaultSignalValue;
      }
    }
  }
      
  // fFlags Initialization
  for(Int_t row = 0; row < kN4x4XPrTRURow; ++row){
    for(Int_t branch = 0; branch < kN4x4ZPrBranch; ++branch){
      for(Int_t timeBin = 0; timeBin < kNTimeBins; ++timeBin){
	fFlags[row][branch][timeBin] = kFALSE;
      }
    }
  }
  
  // fActiveTime Initialization
  for(Int_t timeBin = 0; timeBin < kNTimeBins; ++timeBin){
    fActiveTime[timeBin] = kFALSE;
    fHasSignalTime[timeBin] = kFALSE;
  }
}


//________________________________________________________________
AliPHOSTRURawReader::~AliPHOSTRURawReader()
{
  // destructor
}

//________________________________________________________________
void AliPHOSTRURawReader::ReadFromStream(AliCaloRawStreamV3* rawStream)
{
 // reads the trigger signal amplitudes and trigger flags from the rawStream
  
  const UShort_t * const signal = rawStream->GetSignals(); // stream of 10-bit words, buffered as 16-bit words
  const Int_t signalLength = rawStream->GetBunchLength();
  const Int_t timeBin = rawStream->GetHWAddress() & 0x7f; // , i.e. & '01111111', strips the first 7 bits ([0:6]) and interprits it as a int

  // TODO: this may need to be an error
  if(signalLength != 16 && signalLength != 16+112)
    Warning("ReadFromStream", " signalLength: !16 && !16+112, reader logic may not be valid.");
  

  fActive = kTRUE;
  fActiveTime[timeBin] = kTRUE;

  /* There are 91 4x4 Sliding Window signal sum trigger flags.
   * We read the trigger location information from first 12 10-bit words of the signal
   * (these are buffered to 12 16-bit words in the rawstream object)
   * , but only 10 of the 10-bit words are used for trigger location within the TRU.
   * The flags are found in word_2[0] to word_11[9].
   * We read the information as following. */

  for(Int_t wIdx = 2; wIdx < 12; ++wIdx) { // words
    const std::bitset<10> word(signal[wIdx]); // @param word represent a 10-bit word
    for(Int_t bIdx = 0; bIdx < 10; ++bIdx) { // bits
      //const Int_t index = 119 - (wIdx*10 + (9-bi)); // index used in TRU documentation,
      const Int_t index = (110 - wIdx*10) + bIdx; // equivalent
	if( index < 91 ) { // we are only interrested in these words/bits
	  const Int_t xIdx = index % 7; // x index in TRU internal 2x2 coordinate system
	  const Int_t zIdx = index / 7; // z index in TRU internal 2x2 coordinate system
	  // fFlags[xIdx][zIdx][time] = (signal[wIdx] & (0x01 * 2**i)) != 1;
	  fFlags[xIdx][zIdx][timeBin] = (bool) word[bIdx];
	}
    } // end bit loop
  }// end word loop



  /* The 2x2 trigger signal sum may follow.
   * If so, it is found in the following 112 10-bit words
   * (, 16-bit words in buffer.)
   * We read the signal as following.
   */

  if( 16+112 == signalLength) {
    fHasSignal = kTRUE;
    fHasSignalTime[timeBin] = kTRUE;
    for (Int_t idx = 0; idx < 112; idx++)
      {
	const Int_t xIdx = 7 - idx % 8;  // x index in TRU
	const Int_t zIdx = 13 - idx / 8; // z index in TRU
	const Int_t wIdx = idx + 16;     // word index in signal array
	fSignals[xIdx][zIdx][timeBin] = signal[wIdx];
      }
  }
}


//________________________________________________________________
void AliPHOSTRURawReader::Reset()
{
  // Reset to default values
  
  if( ! fActive )
    return;

  for(Int_t timeBin = 0; timeBin < kNTimeBins; ++timeBin) { // loop timeBins
    if( fActiveTime[timeBin] ) {
      for(Int_t xIdx = 0; xIdx < kN2x2XPrTRURow; ++xIdx) { // loop 2x2
	for(Int_t zIdx = 0; zIdx < kN2x2ZPrBranch; ++zIdx) {
	  fSignals[xIdx][zIdx][timeBin] = kDefaultSignalValue;
	} // zIdx
      } // xIdx
      for(Int_t xIdx = 0; xIdx < kN4x4XPrTRURow; ++xIdx) { // loop 4x4
	for(Int_t zIdx = 0; zIdx < kN4x4ZPrBranch; ++zIdx) {
	  fFlags[xIdx][zIdx][timeBin] = false;
	} // zIdx
      } // xIdx
    }// end if fActiveTime
    fActiveTime[timeBin] = false;
    fHasSignalTime[timeBin] = false;
  } // timeBin

  fActive = false;
  fHasSignal = false;
}
