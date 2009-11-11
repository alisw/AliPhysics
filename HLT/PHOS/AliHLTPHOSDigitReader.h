/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Oystein Djuvsland <oysteind@ift.uib.no>                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIHLTPHOSDIGITREADER_H
#define ALIHLTPHOSDIGITREADER_H


/** @file   AliHLTPHOSClusterizerComponent.cxx
    @author Oystein Djuvsland
    @date   
    @brief  A clusterizer component for PHOS HLT
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


#include "AliHLTPHOSDigitDataStruct.h"

/** 
 * @class AliHLTPHOSDigitReader
 * Class takes as input a AliHLTPHOSDigitHeaderStruct and iterates through 
 * the list of digits following the header
 * @ingroup alihlt_phos
 */
class AliHLTPHOSDigitReader
{

public:
  AliHLTPHOSDigitReader();			
  virtual ~AliHLTPHOSDigitReader();

  void SetDigitHeader(AliHLTPHOSDigitHeaderStruct *digitHeader) 
  { 
    fDigitHeader = digitHeader; 
    if(fDigitHeader->fNDigits != 0)
      {
	fFirstDigit = reinterpret_cast<AliHLTPHOSDigitDataStruct*>(reinterpret_cast<Long_t>(fDigitHeader) + sizeof(AliHLTPHOSDigitHeaderStruct) + fDigitHeader->fFirstDigitOffset);
      }
    else 
      {
	fFirstDigit = 0;
      }
    fNextDigit = fFirstDigit;
  }
    
  void SetCurrentDigit(AliHLTPHOSDigitDataStruct *currentDigit) { fCurrentDigit = currentDigit; }

  AliHLTPHOSDigitDataStruct* NextDigit();

  void DropDigit();

  void Rewind() { fCurrentDigit = fFirstDigit; }

  Int_t GetCurrentDigitOffset() { return reinterpret_cast<Long_t>(fCurrentDigit) - reinterpret_cast<Long_t>(fDigitHeader) + sizeof(AliHLTPHOSDigitHeaderStruct); }
private:
  
  /** Pointer to the digit header */
  AliHLTPHOSDigitHeaderStruct *fDigitHeader;    //COMMENT

  /** Pointer to the current digit */
  AliHLTPHOSDigitDataStruct *fCurrentDigit;    //COMMENT

  /** Pointer to the next digit */
  AliHLTPHOSDigitDataStruct *fNextDigit;    //COMMENT

  /** Pointer to the current digit */
  AliHLTPHOSDigitDataStruct *fPrevDigit;    //COMMENT

  /** Pointer to the first digit */
  AliHLTPHOSDigitDataStruct *fFirstDigit;    //COMMENT

  AliHLTPHOSDigitReader (const AliHLTPHOSDigitReader  & );
  AliHLTPHOSDigitReader & operator = (const AliHLTPHOSDigitReader &);


};


#endif
