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
 * Author: Jussi Viinikainen <jussi.viinikainen@cern.ch> (adaptation to Run2 format)
 */

#include "AliPHOSTRURawReader.h"
#include "AliCaloRawStreamV3.h"
#include "AliLog.h"

ClassImp(AliPHOSTRURawReader)


//________________________________________________________________
AliPHOSTRURawReader::AliPHOSTRURawReader()
: TObject(),
fSignals(),
fFlags(),
fFlags2x2(),
fActive(false),
fHasSignal(false),
fActiveTime(),
fHasSignalTime(),
fUse4x4Flags(false)
{
  // default constructor
  
  // fSignals Initialization:
  for(Int_t row = 0; row < fgkN2x2XPrTRURow; ++row) {
    for(Int_t branch = 0; branch < fgkN2x2ZPrBranch; ++branch) {
      for(Int_t timeBin = 0; timeBin < fgkNTimeBins; ++timeBin) {
        fSignals[row][branch][timeBin] = fgkDefaultSignalValue;
        fFlags2x2[row][branch][timeBin] = kFALSE;
      }
    }
  }
  
  // fFlags Initialization
  for(Int_t row = 0; row < fgkN4x4XPrTRURow; ++row){
    for(Int_t branch = 0; branch < fgkN4x4ZPrBranch; ++branch){
      for(Int_t timeBin = 0; timeBin < fgkNTimeBins; ++timeBin){
        fFlags[row][branch][timeBin] = kFALSE;
      }
    }
  }
  
  // fActiveTime Initialization
  for(Int_t timeBin = 0; timeBin < fgkNTimeBins; ++timeBin){
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
  const Int_t signalLength = rawStream->GetBunchLength();  // The length of the signal in time steps
  const Int_t index = rawStream->GetColumn();  // For some reason the index of the readout channel is given by GetColumn function
  Int_t timeBin = rawStream->GetStartTimeBin(); // Find the time bin of the first time step
  if (timeBin <0 || timeBin >= fgkNTimeBins) {
    AliError(Form("Wrong number of time bins: %d (<0 or >%d!)\n",timeBin,fgkNTimeBins));
    return;
  }
  
  Int_t channelIndex = index;

  if (channelIndex >= 2048)
    channelIndex-=2048; // branch 1: 0 <= z < 28
  
  fActive = kTRUE; // Set the TRU active
  
  /* Channels in TRU:
   *
   * There are 112 readout channels and 12 channels reserved for production flags
   * The channels are indexed as follows:
   *
   *  Channels 0-111: channel data readout
   *  Channels 112-123: production flags
   */
  
      
  if(channelIndex < fgkNReadoutChannels){  // Channel data
    
    /* Channel data:
     *
     * The channel data is read one channel at a time
     * The x and z indices corresponding to the channel are calculated from the channel index
     */
    
    fHasSignal = kTRUE;
    const Int_t xBin = 7 - channelIndex % 8;  // x index in TRU internal 2x2 coordinate system
    const Int_t zBin = 13 - channelIndex / 8; // z index in TRU internal 2x2 coordinate system
    
    // Loop over all the time steps in the signal
    for(Int_t i = 0; i < signalLength; i++){
      if (timeBin > -1 && timeBin < fgkNTimeBins) {
	fSignals[xBin][zBin][timeBin] = signal[i];
	fActiveTime[timeBin] = kTRUE;
	fHasSignalTime[timeBin] = kTRUE;
      }
      timeBin--; // The time bins come in reverse order from raw data
    }
  } else { // Production flags

    /* Production flags:
     *
     * Production flags are supplied in channels 112 - 123
     * Each of the channels is 10 bit wide
     * The bits inside the channel (indexing starting from the first bit of channel 112) is as follows:
     *
     *  Bits 0-111: Trigger flags for corresponding channel index
     *              If using 4x4 algorithm, only 91 first bits are used of these
     *  Bit 112: Marker for 4x4 algorithm (1 active, 0 not active)
     *  Bit 113: Marker for 2x2 algorithm (1 active, 0 not active)
     *  Bit 114: Global L0 OR of all patches in the TRU
     */
    
    Int_t channel, xBin, zBin;
    for(Int_t i = 0; i < signalLength; i++){
      fActiveTime[timeBin] = kTRUE;
      
      // If bit 112 is 1, we are considering 4x4 algorithm
      if(channelIndex == fgkFinalProductionChannel){
        if( (signal[i] & ( 1 << 2 )) > 0 ){ // Check the bit number 112
          fUse4x4Flags = kTRUE;
        }
      }
      
      // Assign the bits in the words to corresponding channels
      for(Int_t bitIndex = 0; bitIndex < fgkWordLength; bitIndex++){
        
        /* Find the correct channel number assuming that
         * channelIndex 112 = bits 0-9 corresponding trigger flags in channels 0-9
         * channelIndex 113 = bits 10-19 corresponding trigger flags in channels 10-19
         * and so on
         */
        channel = (channelIndex - fgkNReadoutChannels) * fgkWordLength + bitIndex;
        
        /*
         * Note: flags corresponding to 2x2 sums and 4x4 sums need to be filled
         *       in different loops since the internal xz coordinates are different
         *       for 2x2 case and 4x4 case. We do not know which data we are filling
         *       at the time we do the filling. This information comes only in the 
          *      channel 123 in bit 112.
         */
        
        if(channel < fgkNReadoutChannels){ // Fill histogram for 2x2 trigger flags
          xBin = 7 - channel % 8;  // x index in TRU internal 2x2 coordinate system
          zBin = 13 - channel / 8; // z index in TRU internal 2x2 coordinate system
          
          // check if the bit bitIndex is 1
          if( (signal[i] & ( 1 << bitIndex )) > 0 ){
            fFlags2x2[xBin][zBin][timeBin] = kTRUE;
          }
        }
        
        if(channel < fgkN4x4TriggerFlags){ // Fill histogram for 4x4 trigger flags
          xBin = channel % 7;  // x index in TRU internal 4x4 coordinate system
          zBin = 12 - channel / 7; // z index in TRU internal 4x4 coordinate system
          
          // check if the bit bitIndex is 1
          if( (signal[i] & ( 1 << bitIndex )) > 0 ){
            fFlags[xBin][zBin][timeBin] = kTRUE;
          }
        }
      } // Bits in one word
      timeBin--; // The time bins come in reverse order from raw data
    } // Length of signal
  } // Production flags
  
}


//________________________________________________________________
void AliPHOSTRURawReader::Reset()
{
  // Reset to default values
  
  if( ! fActive )
    return;
  
  for(Int_t timeBin = 0; timeBin < fgkNTimeBins; ++timeBin) { // loop timeBins
    if( fActiveTime[timeBin] ) {
      for(Int_t xIdx = 0; xIdx < fgkN2x2XPrTRURow; ++xIdx) { // loop 2x2
        for(Int_t zIdx = 0; zIdx < fgkN2x2ZPrBranch; ++zIdx) {
          fSignals[xIdx][zIdx][timeBin] = fgkDefaultSignalValue;
          fFlags2x2[xIdx][zIdx][timeBin] = kFALSE;
        } // zIdx
      } // xIdx
      for(Int_t xIdx = 0; xIdx < fgkN4x4XPrTRURow; ++xIdx) { // loop 4x4
        for(Int_t zIdx = 0; zIdx < fgkN4x4ZPrBranch; ++zIdx) {
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

//___________________________________________________________________
Bool_t AliPHOSTRURawReader::GetTriggerFlag(Int_t xIdx, Int_t zIdx, Int_t timeBin) const {
  // Getter for trigger flags
  if(fUse4x4Flags){
    return fFlags[xIdx][zIdx][timeBin];
  }
  return fFlags2x2[xIdx][zIdx][timeBin];
}
