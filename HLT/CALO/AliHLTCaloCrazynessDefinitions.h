#ifndef ALIHLTCALOCRAZYDEFINITIONS_H
#define ALIHLTCALOCRAZYDEFINITIONS_H

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

/**
 * @file   AliHLTCaloCrazynessDefinitions.h
 * @author Oystein Djuvsland
 * @date
 * @brief  Defines the crazyness variable
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


#include "Rtypes.h"

class AliHLTCaloCrazynessDefinitions
  {

  public:

    /** Always set if crazy */
    static const Short_t fgkCrazyBit   = 0x0001;

    /** Specifies if raw data is included */
    static const Short_t fgkRawDataBit = 0x8000;

    /**
     * These bits are not taken at the moment,
     * feel free to take one!
     */
    static const Short_t fkDummyBits  =
      0x4000 |
      0x2000 |
      0x1000 |
      0x0800 |
      0x0400 |
      0x0200 |
      0x0100 |
      0x0080 |
      0x0040 |
      0x0020 ;

    /** If the channel is completely crazy */
    static const Short_t fgkBadShitBit = 0x0010|fgkCrazyBit;

    /** If the channel has been healed from crazyness */
    static const Short_t fgkHealedBit  = 0x0008|fgkCrazyBit;

    /** If the crazyness is due to a spike */
    static const Short_t fgkSpikeBit   = 0x0004|fgkCrazyBit;

    /** If the crazyness is due to a bad E/T estimate */
    static const Short_t fgkBadEstBit  = 0x0002|fgkCrazyBit;


  };

#endif
