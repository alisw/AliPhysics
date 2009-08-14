
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

#ifndef ALIHLTPHOSESDCALOCLUSTERMAKER_H
#define ALIHLTPHOSESDCALOCLUSTERMAKER_H

/**
 * Class writes ESDs
 *
 * @file   AliHLTPHOSESDCaloClusterMaker.h
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

class AliHLTPHOSCaloClusterHeaderStruct;
class AliHLTPHOSCaloClusterReader;
class TClonesArray;
class AliESDCaloCluster;

/** 
 * @class AliHLTPHOSESDCaloClusterMaker
 * Makes ESD Clusters out of AliHLTPHOSCaloClusterContainerStructs
 * @ingroup alihlt_phos
 */

class AliHLTPHOSESDCaloClusterMaker : public AliHLTPHOSBase
{

 public:
  
  
  /** Constructor */
  AliHLTPHOSESDCaloClusterMaker();
  
  /** Destructor */
  virtual ~AliHLTPHOSESDCaloClusterMaker();

  /** Copy constructor */  
  AliHLTPHOSESDCaloClusterMaker(const AliHLTPHOSESDCaloClusterMaker &) : 
    AliHLTPHOSBase(),
    fNCaloClusters(0),
    fClusterReaderPtr(0),
    fESDCaloClusterPtr(0)
    {
      //Copy constructor not implemented
    }
  
  /** Assignment */
  AliHLTPHOSESDCaloClusterMaker & operator = (const AliHLTPHOSESDCaloClusterMaker)
    {
      //Assignment
      return *this; 
    }

  /** 
   * Create AliESDCaloClusters from the AliHLTPHOSCaloClusterContainerStruct
   * @return number of clusters created
   */
  Int_t FillESDCaloClusters(TClonesArray* esdClusters, AliHLTPHOSCaloClusterHeaderStruct* clusterContainer);

 private: 
  
  /* Number of clusters */
  Int_t fNCaloClusters;  // Number of clusters

  /* The reader */
  AliHLTPHOSCaloClusterReader* fClusterReaderPtr; // !transient The reader

  /* An ESD Calo Cluster */
  AliESDCaloCluster* fESDCaloClusterPtr; // !transient An ESD Calo Cluster

  ClassDef(AliHLTPHOSESDCaloClusterMaker, 0);
};

#endif
