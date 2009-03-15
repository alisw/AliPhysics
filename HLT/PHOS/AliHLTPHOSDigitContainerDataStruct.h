//-*- Mode: C++ -*-
// $Id$

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


#ifndef ALIHLTPHOSDIGITCONTAINERSTRUCT_H
#define ALIHLTPHOSDIGITCONTAINERSTRUCT_H

/**
 * Digit data struct for PHOS HLT
 *
 * @file   AliHLTPHOSDigitContainerDataStruct.h
 * @author Oystein Djuvsland
 * @date
 * @brief  Digit container data struct for PHOS HLT
 */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "Rtypes.h"
#include "AliHLTPHOSDigitDataStruct.h"

/**
 * @struct AliHLTPHOSDigitContainerDataStruct
 * Digit container data struct for PHOS HLT
 *
 * @ingroup alihlt_phos
 */
struct AliHLTPHOSDigitContainerDataStruct
{

  /** Number of digits in container */
  UInt_t fNDigits;                                     //COMMENT

  /** Module number */
  UInt_t fPHOSModule;                                  //COMMENT

  /** Array of digits in container */
  AliHLTPHOSDigitDataStruct fDigitDataStruct[NXCOLUMNSRCU*NZROWSRCU*NGAINS];    //COMMENT

};


#endif
