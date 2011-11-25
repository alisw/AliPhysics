
#ifndef ALIHLTCALOCLUSTERREADER_H
#define ALIHLTCALOCLUSTERREADER_H

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

#include "Rtypes.h"

struct AliHLTCaloClusterDataStruct;
class AliHLTCaloClusterHeaderStruct;
class AliHLTCaloDigitDataStruct;
/** 
 * @class AliHLTCaloClusterReader
 * Reads clusters from a memory block of AliHLTCaloClusters
 *
 * @ingroup alihlt_phos
 */
class  AliHLTCaloClusterReader
{
public:
  /** Default constructor */
  AliHLTCaloClusterReader();

  /** Destructor */
  virtual ~AliHLTCaloClusterReader();

  /** 
   * Get the next cluster from the memory buffer
   */
  AliHLTCaloClusterDataStruct*   NextCluster();

  
  /**
   * Get cell properties from a cluster
   * @param clusterPtr is a pointer to the cluster
   * @param cellID is a reference to the variable containing the absolute ID
   * @param cellAmp is a reference to the variable containing the cell amplitude fraction
   * @param index is the index of the cell
   * @return true if the cell exists
   */
  bool GetCell(AliHLTCaloClusterDataStruct *clusterPtr, UShort_t &cellID, Double32_t &cellAmp, UInt_t index);


  /** 
   * Set the memory buffer containing the clusters
   * @param clusterHeaderPtr pointer the the first entry in the buffer, the header
   */
  void SetMemory(const AliHLTCaloClusterHeaderStruct* clusterHeaderPtr);

  /**
   * Reset the memory pointer and number of counts
   */
  void Reset();

  /**
  * Get the digits shipped with the clusters
  */
  AliHLTCaloDigitDataStruct* GetDigits() {return fDigitsPointer; } 
  
  /** 
  * Get the number of digits shipped with the clusters
  */
  Int_t GetNDigits() {return fNDigits; }
  
 private:
  
  /**
   * Copy constructor
   */
  AliHLTCaloClusterReader(const  AliHLTCaloClusterReader & );
  
  /**
   * Assignment operator
   */
  AliHLTCaloClusterReader & operator = (const  AliHLTCaloClusterReader &);
  
  /* Pointer to the current cluster to be read */
  AliHLTCaloClusterDataStruct* fCurrentClusterPtr; // !transient Pointer to the current cluster to be read

  /* Check if the memory has been set */
  bool fIsSetMemory; //Check if the memory has been set

  /* Max number of clusters */
  int fMaxCnt;  // Max number of clusters

  /* The current number of clusters */
  int fCurrentCnt; // The current number of clusters

  /** Pointer to the digits in the data block */
  AliHLTCaloDigitDataStruct* fDigitsPointer;    // Pointer to the digits in the data block 
    
   /** Number of digits in the data block */
   Int_t fNDigits; // Number of digits in the data block 
   
};

#endif
