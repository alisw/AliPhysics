
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#ifndef ALIHLTPHOSDIGITMAKER_H
#define ALIHLTPHOSDIGITMAKER_H

/**
 * Class makes digits from information from raw data
 *
 * @file   AliHLTPHOSDigitMaker.h
 * @author Oystein Djuvsland
 * @date
 * @brief  Digit maker for PHOS HLT
 */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSBase.h"
#include "AliHLTPHOSConstants.h"

/**
 * @class AliHLTPHOSDigitMaker
 * Digit maker for PHOS HLT. Takes input from AliHLTPHOSRawAnalyzer, and 
 * outputs an AliHLTPHOSDigitDataStruct container, or a TClonesArray of
 * AliHLTPHOSDigit. Can do software zero suppression
 * @ingroup alihlt_phos
 */

class AliHLTPHOSDigit;
class TClonesArray;
class TTree;
class AliHLTPHOSValidCellDataStruct;
class AliHLTPHOSRcuCellEnergyDataStruct;
class AliHLTPHOSDigitContainerDataStruct;
class AliHLTPHOSDigitDataStruct;
       
using namespace PhosHLTConst;

class AliHLTPHOSDigitMaker : public AliHLTPHOSBase
{
public:

  /** Constructor */
  AliHLTPHOSDigitMaker();

  /** Destructor */
  virtual ~AliHLTPHOSDigitMaker();
 
  // void SetValidCellData(AliHLTPHOSValidCellDataStruct *data) { fCellDataPtr = data; }
  //  void SetDigitContainerStruct(AliHLTPHOSDigitContainerStruct *container) 
  //{ fDigitContainerStructPtr = container; }
   
  /**
   * Sets the AliHLTPHOSDigitDataStruct container 
   * @param container the digit container
   */
  void SetDigitContainerStruct(AliHLTPHOSDigitContainerDataStruct *container) 
  { fDigitContainerStructPtr = container; }

  /** 
   * Sets the TClonesArray of AliHLTPHOSDigit 
   * @param array the array
   */
  void SetDigitArray(TClonesArray *array) { fDigitArrayPtr = array; }

  /** 
   * Sets the digit threshold 
   * @param threshold the threshold
   */
  void SetDigitThreshold(Int_t threshold) { fDigitThreshold = threshold; }

  /**
   * Sets the number of pre samples
   * @param n the number of pre samples
   */
  void SetNrPresamples(Int_t n) { fNrPresamples = n; }

  /**
   * Make the digits for one event.
   * @param rcuCellEnergies is the data from the AliHLTPHOSRawAnalyzer
   * @return the number of digits found
   */
  Int_t MakeDigits(AliHLTPHOSRcuCellEnergyDataStruct* rcuCellEnergies);

  /** Reset the digit maker */
  void Reset();

private:

  /** Pointer to valid cell list */
  AliHLTPHOSValidCellDataStruct *fCellDataPtr;                   //! transient

  /** Pointer to the digit container */
  AliHLTPHOSDigitContainerDataStruct *fDigitContainerStructPtr;  //! transient

  /** Pointer to the digit TClonesArray */
  TClonesArray *fDigitArrayPtr;                                  //! transient

  /** Pointer to a AliHLTPHOSDigit */
  AliHLTPHOSDigit *fDigitPtr;                                    //! transient

  /** Pointer to a AliHLTPHOSDigitDataStruct */
  AliHLTPHOSDigitDataStruct *fDigitStructPtr;                    //! transient

  /** Digit count */
  Int_t fDigitCount; 

  /** Number of presamples */
  Int_t fNrPresamples; 

  /** Threshold for making digit ( zero suppression threshold) */
  Float_t fDigitThreshold; 

  ClassDef(AliHLTPHOSDigitMaker, 1); 
};


#endif
 
