//-*- Mode: C++ -*-
// $Id: AliHLTCaloRecPointDataStruct.h 29824 2008-11-10 13:43:55Z richterm $


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

#ifndef ALIHLTCALORECPOINTDATASTRUCT_H
#define ALIHLTCALORECPOINTDATASTRUCT_H

/**
 * Rec point data struct for CALO HLT
 *
 * @file   AliHLTCaloRecPointDataStruct.h
 * @author Oystein Djuvsland
 * @date
 * @brief  Rec point data struct for CALO HLT
 */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTCaloDigitDataStruct.h"

/**
 * @struct AliHLTCaloRecPointDataStruct
 * Rec point data struct for Calo HLT
 *
 * @ingroup alihlt_calo
 */
struct AliHLTCaloRecPointDataStruct
{

  //AliHLTUInt8_t fMultiplicity; 
  /** Multiplicity of digits in the rec point */
  UInt_t fMultiplicity;                       //COMMENT

  /** x coordinate */
  Float_t fX;                                 //COMMENT
  
  /** y coordinate */
  Float_t fY;									//added, federico
  
  /** z coordinate */ 
  Float_t fZ;                                 //COMMENT

  /** Module number */
  Int_t fModule;                              //COMMENT
   
   /** Particle type */
  Int_t fParticle;                              //COMMENT
  
  /** The total energy of the rec point */
  Float_t fAmp;                               //COMMENT

  /** Second moment along x axis */
  Float_t fM2x;                               //COMMENT

  /** Second moment along z axis */ 
  Float_t fM2z;                               //COMMENT

  /** Third moment along x axis */
  Float_t fM3x;                               //COMMENT

  /** Fourth moment along z axis */ 
  Float_t fM4z;                               //COMMENT

  /** Angle between cog vector and eigen vector */
  Float_t fPhixe;                             //COMMENT

  /** Distance to nearest bad channel */
  Float_t fDistanceToBadChannel;              //COMMENT

  /** Index of the digits in the rec point */
  Int_t fDigits;                              //COMMENT

};

#endif
