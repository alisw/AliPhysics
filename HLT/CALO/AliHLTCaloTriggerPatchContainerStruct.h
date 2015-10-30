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
#ifndef ALIHLTCALOTRIGGERPATCHCONTAINERSTRUCT_H
#define ALIHLTCALOTRIGGERPATCHCONTAINERSTRUCT_H
/**
 * Rec point container struct for CALO HLT
 *
 * @file   AliHLTCaloTriggerPatchContainerStruct.h
 * @author Markus Fasel
 * @date
 * @brief  Trigger patch container struct for CALO HLT
 */

#include "AliHLTCaloTriggerPatchDataStruct.h"

/**
 * @struct AliHLTCaloTriggerPatchContainerStruct
 * Trigger patch container struct for CALO HLT
 *
 * @ingroup alihlt_calo
 */
struct AliHLTCaloTriggerPatchContainerStruct {

  /** Number of trigger patches in the container */
  UInt_t                                        fNTriggerPatches;
  /** Array of reconstructed trigger patches */
  AliHLTCaloTriggerPatchDataStruct              fTriggerPatches[100];
};


#endif
