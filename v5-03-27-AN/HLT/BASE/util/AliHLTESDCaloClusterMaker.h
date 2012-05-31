// -*- Mode: C++ -*-
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

#ifndef ALIHLTESDCALOCLUSTERMAKER_H
#define ALIHLTESDCALOCLUSTERMAKER_H

/**
 * Class writes ESDs
 *
 * @file   AliHLTESDCaloClusterMaker.h
 * @author Oystein Djuvsland
 * @date
 * @brief  ESD writer for  HLT
 */

// see header file for class documentation
// or
// refer to README to build package 
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

class AliESDEvent;
class AliHLTCaloClusterHeaderStruct;
class AliHLTCaloClusterReader;
class TClonesArray;
class AliESDCaloCluster;

/** 
 * @class AliHLTESDCaloClusterMaker
 * Makes ESD Clusters out of AliHLTCaloClusterDataStructs
 * @ingroup alihlt_phos
 */

class AliHLTESDCaloClusterMaker 
{

 public:
  
  
  /** Constructor */
  AliHLTESDCaloClusterMaker();
  
  /** Destructor */
  virtual ~AliHLTESDCaloClusterMaker();

  /**
   * Convert AliHLTCaloClusterDataStruct clusters and fill an ESDEvent object with 
   * AliESDCaloCluster clusters
   * @return number of clusters converted and filled
   */
  Int_t FillESD(AliESDEvent *esdPtr, const AliHLTCaloClusterHeaderStruct *clusterHeader);

 private: 
  /** Copy constructor prohibited */  
  AliHLTESDCaloClusterMaker(const AliHLTESDCaloClusterMaker &);
  /** Assignment operator prohibited*/
  AliHLTESDCaloClusterMaker & operator = (const AliHLTESDCaloClusterMaker);

  /* The cluster struct reader */
  AliHLTCaloClusterReader* fClusterReaderPtr; // !transient The reader

  ClassDef(AliHLTESDCaloClusterMaker, 0);
};

#endif
