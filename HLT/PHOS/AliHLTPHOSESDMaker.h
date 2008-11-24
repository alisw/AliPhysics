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

#ifndef ALIHLTPHOSESDMAKER_H
#define ALIHLTPHOSESDMAKER_H

/**
 * Class writes ESDs
 *
 * @file   AliHLTPHOSESDMaker.h
 * @author Oystein Djuvsland
 * @date
 * @brief  ESD writer for PHOS HLT
 */

// see header file for class documentation
// or
// refer to README to build package 
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSBase.h"

class AliHLTPHOSCaloClusterContainerStruct;
class TClonesArray;
class AliESDEvent;

/** 
 * @class AliHLTPHOSESDMaker
 * Makes ESDs out of reconstruction points, first it creates
 * the AliESDCaloClusters, then puts them into an AliESD object
 *
 * @ingroup alihlt_phos
 */
class AliHLTPHOSESDMaker : public AliHLTPHOSBase
{

public:
  
  /** Constructor */
  AliHLTPHOSESDMaker();
  
  /** Destructor */
  virtual ~AliHLTPHOSESDMaker();

  /** Copy constructor */  
  AliHLTPHOSESDMaker(const AliHLTPHOSESDMaker &) : 
    AliHLTPHOSBase(),
    fNCaloClusters(0),
    fCaloClustersPtr(0),
    fESDEventPtr(0),
    fCaloClusterContainerPtr(0)
  {
    //Copy constructor not implemented
  }
  
  /** Assignment */
  AliHLTPHOSESDMaker & operator = (const AliHLTPHOSESDMaker)
  {
    //Assignment
    return *this; 
  }

  /** 
   * Set the AliESDEvent object to be filled
   * @param esdEventPtr is a pointer to the AliESDEvent
   */
  void SetESDEvent(AliESDEvent* esdEventPtr) { fESDEventPtr = esdEventPtr; }

  /**
   * Set the AliHLTPHOSCaloClusterContainerStruct
   * @param clusterContainerPtr is a pointer to the cluster container
   */
  void SetCaloClusterContainer(AliHLTPHOSCaloClusterContainerStruct* clusterContainerPtr) 
  { fCaloClusterContainerPtr = clusterContainerPtr; }
  /** 
   * Create AliESDCaloClusters from the AliHLTPHOSCaloClusterContainerStruct
   * @return
   */
  Int_t FillESDCaloClusters();
  
  /** 
   * Fill the AliESDEvent object
   * @return
   */
  Int_t FillESDEvent();

  /** 
   * Fill the AliESDEvent object with clusters from a calo cluster container
   * @param caloClusterContainerPtr is a pointer to a cluster container
   * @return
   */
  Int_t FillESDEvent(AliHLTPHOSCaloClusterContainerStruct* caloClusterContainerPtr);

  /**
   * Reset the ESD and ESDCaloCluster array
   */
  void ResetESD();

  /** 
   * Get the AliESDCaloClusters
   * @return a pointer to a TClonesArray of AliESDCaloClusters
   */
  TClonesArray* GetESDCaloClusters() { return fCaloClustersPtr; }
    
private:

  /** Number of calo clusters */
  Int_t fNCaloClusters;                                             //COMMENT

  /** Array of AliESDCaloClusters */
  TClonesArray* fCaloClustersPtr;                                   //! transient

  /** The AliESDEvent object to fill */
  AliESDEvent* fESDEventPtr;                                        //! transient

  /** The AliHLTPHOSCaloClusterContainerStruct */               
  AliHLTPHOSCaloClusterContainerStruct* fCaloClusterContainerPtr;   //! transient

  ClassDef(AliHLTPHOSESDMaker, 1);
};

#endif
