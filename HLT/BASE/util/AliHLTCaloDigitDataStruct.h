//-*- Mode: C++ -*-
// $Id: AliHLTCaloDigitDataStruct.h 35319 2009-10-07 15:27:09Z odjuvsla $


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

#ifndef ALIHLTCALODIGITDATASTRUCT_H
#define ALIHLTCALODIGITDATASTRUCT_H

#include "Rtypes.h"

/**
 * Digit struct for Calo HLT
 *
 * @file   AliHLTCaloDigitDataStruct.h
 * @author Oystein Djuvsland
 * @date
 * @brief  Digit struct for calo HLT
 */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

/**
 * @struct AliHLTCaloDigitDataStruct
 * Digit struct for Calo HLT
 *
 * @ingroup alihlt_calo
 */
struct AliHLTCaloDigitDataStruct
{

  /** Unique ID (in the module) for this digit */
  UInt_t fID;
  
  /** The x coordinate */
  UShort_t fX;

  /** The z coordinate */
  UShort_t fZ;

  /** The module number */
  Int_t fModule;

  /** The amplitude in ADC counts */
  Float_t fAmplitude;

  /** The time in sample count */ 
  Float_t fTime;

  /** The energy in GeV */
  Float_t fEnergy;

  /** The gain */
  Int_t fGain;
  
  /** The crazyness */
  Int_t fCrazyness; 

  /** The baseline */
  Float_t fBaseline;

  /** Energy from overflow in channel? */
  Bool_t fOverflow;

  /** ID of associated  cluster (-1 for no cluster) */
  Short_t fAssociatedCluster;
  
};

#endif

