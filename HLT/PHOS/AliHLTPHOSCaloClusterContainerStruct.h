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

#ifndef ALIHLTPHOSCALOCLUSTERCOINTAINERSTRUCT_H
#define ALIHLTPHOSCALOCLUSTERCOINTAINERSTRUCT_H
/**
 * Calo cluster container struct for PHOS HLT
 *
 * @file   AliHLTPHOSCaloClusterContainerStruct.h
 * @author Oystein Djuvsland
 * @date
 * @brief  Rec point container struct for PHOS HLT
 */

#include "AliHLTPHOSCaloClusterDataStruct.h"

/**
 * @struct AliHLTPHOSCaloClusterContainerStruct
 * Calo cluster container struct for PHOS HLT
 *
 * @ingroup alihlt_phos
 */
struct AliHLTPHOSCaloClusterContainerStruct
{

  /** The PHOS module number */
  UInt_t fPHOSModule;

  /** The number of clusters */
  UInt_t fNCaloClusters; 

  /** Array of rec points in the container */
  AliHLTPHOSCaloClusterDataStruct fCaloClusterArray[100];

};

#endif
