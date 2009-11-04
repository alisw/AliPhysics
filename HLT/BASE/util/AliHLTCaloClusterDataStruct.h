
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

#ifndef ALIHLTCALOCLUSTERDATASTRUCT_H
#define ALIHLTCALOCLUSTERDATASTRUCT_H

/**
 * Calo cluster struct for  HLT
 *
 * @file   AliHLTCaloClusterDataStruct.h
 * @author Oystein Djuvsland
 * @date
 * @brief  Calo cluster struct for HLT
 */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliPID.h"
#include "TArrayI.h"
/**
 * @struct AliHLTCaloClusterHeaderStruct
 * Calorimeter cluster header describing the number of 
 * clusters in the following block
 *
 * @ingroup alihlt_phos
 */

struct AliHLTCaloClusterHeaderStruct
{
  Short_t fNClusters;
};

/**
 * @struct AliHLTCaloClusterDataStruct
 * Calorimeter cluster data struct for  HLT
 * Similar to the AliESDCaloCluster class
 * @ingroup alihlt_phos
 */

struct AliHLTCaloClusterDataStruct
{

  /** Number of cells in the cluster */
  UInt_t fNCells;                                //COMMENT

  /** Global position */
  Float_t fGlobalPos[3];                      //COMMENT

  /** The total energy of the cell */
  Float_t fEnergy;                            //COMMENT

  /** The time of flight */
  Float_t fTOF;                               //COMMENT

  /** Dispersion */
  Float_t fDispersion;                        //COMMENT

  /** Quality of cluster fit */
  Float_t fFitQuality;                        //COMMENT

  /** Second moment along the main eigen axis */
  Float_t fM20;                               //COMMENT

  /** Second moment along the second eigen axis */ 
  Float_t fM02;                               //COMMENT

  /** Second mixed  moment Mxy */
  //  Float_t fM11;                               //COMMENT
  
  /** Distance to closest CPV rec point */
  Float_t fEmcCpvDistance;                    //COMMENT

  /** Distance to nearest bad channel */
  Float_t fDistToBadChannel;                  //COMMENT

  /** PID */
  Float_t fPID[AliPID::kSPECIESN];            //COMMENT

  /** Unique ID of the cluster*/
  Int_t fID;                                     //COMMENT

  /** Number of (Ex) Maxima */
  UChar_t fNExMax;                               //COMMENT 

  /** Flag for differtent cluster type/versions */
  Char_t fClusterType;                           //COMMENT

  /** Distance to nearest bad channel */
  Float_t fDistanceToBadChannel;              //COMMENT

  /** The absolute IDs of the cells*/
  UShort_t fCellsAbsId;                      //COMMENT

  /** */
  Float_t fCellsAmpFraction;              //COMMENT

  /** Labels of matching tracks. First entry is best match */
  TArrayI* fTracksMatched;  //COMMENT
  
};


#endif
