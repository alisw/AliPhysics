//-*- Mode: C++ -*-
// $Id: AliHLTCaloRecPointContainerStruct.h 29824 2008-11-10 13:43:55Z richterm $

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                                      *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIHLTCALORECPOINTCOINTAINERSTRUCT_H
#define ALIHLTCALORECPOINTCOINTAINERSTRUCT_H
/**
 * Rec point container struct for CALO HLT
 *
 * @file   AliHLTCaloRecPointContainerStruct.h
 * @author Oystein Djuvsland
 * @date
 * @brief  Rec point container struct for CALO HLT
 */

#include "AliHLTCaloRecPointDataStruct.h"

/**
 * @struct AliHLTCaloRecPointContainerStruct
 * Rec point container struct for CALO HLT
 *
 * @ingroup alihlt_calo
 */
struct AliHLTCaloRecPointContainerStruct
{

  /** The CALO module number */
  UInt_t fCaloModule;

  /** The number of rec points */
  UInt_t fNRecPoints; 

  /** Array of rec points in the container */
  AliHLTCaloRecPointDataStruct fRecPointArray[100];

};

#endif
