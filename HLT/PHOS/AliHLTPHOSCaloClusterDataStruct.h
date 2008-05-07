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

#ifndef ALIHLTPHOSCALOCLUSTERDATASTRUCT_H
#define ALIHLTPHOSCALOCLUSTERDATASTRUCT_H

/**
 * Rec point data struct for PHOS HLT
 *
 * @file   AliHLTPHOSClusterDataStruct.h
 * @author Oystein Djuvsland
 * @date
 * @brief  Rec point data struct for PHOS HLT
 */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliPID.h"

/**
 * @struct AliHLTPHOSCaloClusterDataStruct
 * Calorimeter cluster data struct for PHOS HLT
 *
 * @ingroup alihlt_phos
 */

struct AliHLTPHOSCaloClusterDataStruct
{

  /** Number of cells in the cluster */
  UInt_t fNCells;                                //COMMENT

  /** */
  UShort_t fCellsAbsId[64];                      //COMMENT

  /** */
  Double32_t fCellsAmpFraction[64];              //[0.,1.,16]

  /** Global position */
  Double32_t fGlobalPos[3];                      //COMMENT

  /** The total energy of the cell */
  Double32_t fEnergy;                            //COMMENT

  /** Dispersion */
  Double32_t fDispersion;                        //COMMENT

  /** Quality of cluster fit */
  Double32_t fFitQuality;                        //COMMENT

  /** Second moment along the main eigen axis */
  Double32_t fM20;                               //COMMENT

  /** Second moment along the second eigen axis */ 
  Double32_t fM02;                               //COMMENT

  /** Second mixed  moment Mxy */
  //  Double32_t fM11;                               //COMMENT
  
  /** Distance to closest CPV rec point */
  Double32_t fEmcCpvDistance;                    //COMMENT

  /** Distance to nearest bad channel */
  Double32_t fDistToBadChannel;                  //COMMENT

  /** PID */
  Double32_t fPID[AliPID::kSPECIESN];            //[0.,1.,8]

  /** Unique ID of the cluster*/
  Int_t fID;                                     //COMMENT

  /** Number of (Ex) Maxima */
  UChar_t fNExMax;                               //COMMENT 

  /** Flag for differtent cluster type/versions */
  Char_t fClusterType;                           //COMMENT

  /** Distance to nearest bad channel */
  Double32_t fDistanceToBadChannel;              //COMMENT
  
};

#endif
