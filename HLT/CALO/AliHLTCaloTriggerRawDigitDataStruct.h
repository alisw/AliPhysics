/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Markus Fasel                                          *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#ifndef ALIHLTEMCALTRIGGERRAWDIGITDATASTRUCT_H
#define ALIHLTEMCALTRIGGERRAWDIGITDATASTRUCT_H
/**
 * Trigger raw digit data struct for EMCAL HLT
 * Based on AliEMCALTriggerRawDigit
 *
 * @file   AliHLTEMCALTriggerRawDigitDataStruct.h
 * @author Markus Fasel
 * @date
 * @brief  Trigger raw digit data struct for EMCAL HLT
 */

#include <Rtypes.h>
#include "AliEMCALTriggerTypes.h"

/**
 * @struct AliHLTEMCALTriggerRawDigitDataStruct
 * Trigger raw digit struct for EMCAL HLT
 *
 * @ingroup alihlt_calo
 */
struct AliHLTCaloTriggerRawDigitDataStruct {
  /**  Unique ID */
  Int_t       fID;
  /** Trigger bits */
  Int_t       fTriggerBits;
  /** N L0 times */
  UChar_t     fNL0Times;
  /** L0 times */
  UChar_t     fL0Times[10];
  /** L1 time sum */
  Int_t       fL1TimeSum;
  /** Number of time samples */
  UChar_t     fNTimeSamples;
  /** Time samples */
  Int_t       fTimeSamples[15];
  /** Subregion */
  Int_t       fL1SubRegion;
};

/**
 * Initialize the raw digit
 * @param rawdigit Input digit
 */
void InitializeRawDigit(AliHLTCaloTriggerRawDigitDataStruct &rawdigit);

/**
 * Set the unique ID of the trigger digit (fastor index)
 * @param digit Input digit
 * @param id ID
 */
void SetRawDigitID(AliHLTCaloTriggerRawDigitDataStruct &digit, Int_t id);

/**
 * Set trigger bit to the digit
 * @param dig Input digit
 * @param bit Bit to be set
 * @param mode Mode (Data or MC)
 */
void SetTriggerBit(AliHLTCaloTriggerRawDigitDataStruct &dig, Int_t bit, Int_t mode);

/**
 * Set L0 time to the digit
 * @param dig Input digit
 * @param itime
 */
void SetL0Time(AliHLTCaloTriggerRawDigitDataStruct &dig, Int_t i);

/**
 * Set L1 time sum
 * @param dig Input digit
 * @param l1timeSum
 */
void SetL1TimeSum(AliHLTCaloTriggerRawDigitDataStruct &dig, Int_t l1timeSum);

/**
 * Set the L1 subregion
 * @param dig Input digit
 * @param l1subregion l1 subregion
 */
void SetL1SubRegion(AliHLTCaloTriggerRawDigitDataStruct &dig, Int_t l1subregion);

/**
 * Set time samples to the digit
 * @param dig Input digit
 * @param nsamples Number of time samples
 * @param samples Array with time samples
 */
void SetTimeSamples(AliHLTCaloTriggerRawDigitDataStruct &dig, Int_t nsamples, Int_t *samples);

/**
 * Get the raw digit ID
 * @param dig Input digit
 * @return ID of the digit
 */
Int_t GetRawDigitID(const AliHLTCaloTriggerRawDigitDataStruct &dig);

/**
 * Get L0 times
 * @param dig Input digit
 * @param i Index
 * @param time L0 time
 * @return False in case of error, true in case of success
 */
Bool_t GetL0Time(const AliHLTCaloTriggerRawDigitDataStruct &dig, Int_t i, Int_t& time);

/**
 * Get L0 times
 * @param dig Input digit
 * @param times Output buffer for the L0 times
 */
void GetL0Times(const AliHLTCaloTriggerRawDigitDataStruct &dig, Int_t times[]);

/**
 * returns the time and amplitude of a given time sample and if the sample was ok
 * @param dig Input digit
 * @param iSample sample index
 * @param timeBin Output container for time bin
 * @param amp Amp at time bin
 * @return False in case of error, true in case of success
 */
Bool_t GetTimeSample(const AliHLTCaloTriggerRawDigitDataStruct &dig, Int_t iSample, Int_t& timeBin, Int_t& amp);

/**
 * Get L0 time sum
 * @param dig Input digit
 * @param time Time bin
 * @return Level0 time sum
 */
Int_t GetL0TimeSum(const AliHLTCaloTriggerRawDigitDataStruct &dig, Int_t time);

/**
 * Get the L1 subregion
 * @param dig Input digit
 * @return Level1 subregion
 */
Int_t GetL1SubRegion(const AliHLTCaloTriggerRawDigitDataStruct &dig);

/**
 * Get trigger bit
 * @param dig Input digit
 * @param type Trigger bit type
 * @param mode Mode (Data or MC)
 * @return Trigger bit (0 or 1)
 */
Int_t GetTriggerBit(const AliHLTCaloTriggerRawDigitDataStruct &dig, const TriggerType_t type, const Int_t mode);

/**
 * Checks the maximum amplitude in the time sample
 * @param dig
 * @param amplitude
 * @param time
 * @return
 */
Bool_t GetRawDigitMaximumAmplitude(const AliHLTCaloTriggerRawDigitDataStruct &dig, Int_t& amplitude, Int_t& time);

/**
 * Print raw digit information
 * @param dig Input digit
 */
void PrintRawDigit(const AliHLTCaloTriggerRawDigitDataStruct &dig);

/**
 * Print raw digit information into an output buffer
 * @param dig Input digit
 * @param outputbuffer Buffer the string is written to
 */
void SPrintRawDigit(const AliHLTCaloTriggerRawDigitDataStruct &dig, char *outputbuffer);

#endif
