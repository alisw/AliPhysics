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
#ifndef ALIHLTEMCALTRIGGERPATCHSTRUCT_H
#define ALIHLTEMCALTRIGGERPATCHSTRUCT_H
/**
 * Trigger patch data struct for CALO HLT
 *
 * @file   AliHLTCaloTriggerPatchDataStruct.h
 * @author Markus Fasel
 * @date
 * @brief  Trigger patch data struct for CALO HLT
 */

#include <Rtypes.h>

/**
 * @struct AliHLTCaloTriggerPatchDataStruct
 * Trigger patch data struct for Calo HLT
 *
 * @ingroup alihlt_calo
 */
struct AliHLTCaloTriggerPatchDataStruct {
  /** Starting column of the trigger patch */
  UChar_t                   fCol;
  /** Starting row of the trigger patch */
  UChar_t                   fRow;
  /** Size (in number of FASTORS) of the trigger patch */
  UChar_t                   fSize;
  /** Accumulated ADC counts in the trigger patch */
  UInt_t                    fADC;
  /** Bit mask of the trigger patch */
  UInt_t                    fBitMask;
  /** Data origin (STU, RECAL, Digits) */
  UChar_t                   fOrigin;
};

#endif
