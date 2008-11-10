//-*- Mode: C++ -*-
// $Id$


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

#ifndef ALIHLTPHOSDIGITDATASTRUCT_H
#define ALIHLTPHOSDIGITDATASTRUCT_H

/**
 * Digit struct for PHOS HLT
 *
 * @file   AliHLTPHOSDigitDataStruct.h
 * @author Oystein Djuvsland
 * @date
 * @brief  Digit struct for PHOS HLT
 */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

/**
 * @struct AliHLTPHOSDigitDataStruct
 * Digit struct for PHOS HLT
 *
 * @ingroup alihlt_phos
 */
struct AliHLTPHOSDigitDataStruct
{
  /** The x coordinate */
  Int_t fX;

  /** The x coordinate */
  Int_t fZ;

  /** The module number */
  Int_t fModule;

  /** The amplitude in ADC counts */
  Float_t fAmplitude;

  /** The time in sample count */ 
  Float_t fTime;

  /* The energy in GeV */
  Float_t fEnergy;

  /** The gain */
  Int_t fGain;
  
  /** The crazyness */
  Int_t fCrazyness; 

  /** The baseline */
  Float_t fBaseline;

};

#endif

