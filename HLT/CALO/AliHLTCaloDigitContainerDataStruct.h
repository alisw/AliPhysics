//-*- Mode: C++ -*-
// $Id: AliHLTCaloDigitContainerDataStruct.h 35107 2009-09-30 01:45:06Z phille $

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                                      *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
x * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#ifndef ALIHLTCALODIGITCONTAINERSTRUCT_H
#define ALIHLTCALODIGITCONTAINERSTRUCT_H

/**
 * Digit data struct for CALO HLT
 *
 * @file   AliHLTCaloDigitContainerDataStruct.h
 * @author Oystein Djuvsland
 * @date
 * @brief  Digit container data struct for CALO HLT
 */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "Rtypes.h"
#include "AliHLTCaloDigitDataStruct.h"



// using namespace CaloHLTConst;

/**
 * @struct AliHLTCaloDigitContainerDataStruct
 * Digit container data struct for Calo HLT
 *
 * @ingroup alihlt_calo
 */
struct AliHLTCaloDigitContainerDataStruct
{

  /** Number of digits in container */
  UInt_t fNDigits;                                     //COMMENT

  /** Module number */
  UInt_t fCaloModule;                                  //COMMENT

  /** Array of digits in container */
  AliHLTCaloDigitDataStruct fDigitDataStruct[10000];

};


#endif
