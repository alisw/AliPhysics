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

#ifndef ALIHLTCALOCLUSTERIZER_H
#define ALIHLTCALOCLUSTERIZER_H


/**
 * Class does clusterization in for Calorimeters on an event basis. It is intended 
 * for use in HLT, but can also be used offline
 *
 * @file   AliHLTCaloClusterizer.h
 * @author Oystein Djuvsland
 * @date
 * @brief  Clusterizer for CALO HLT
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

//#include "AliHLTCaloBase.h"

#include "AliHLTCaloRecPointContainerStruct.h"
#include "AliHLTCaloRecPointDataStruct.h"
#include "AliHLTCaloDigitContainerDataStruct.h"
#include "AliHLTCaloDigitDataStruct.h"
#include "TString.h"
#include "AliHLTCaloConstantsHandler.h"

//#include "AliPHOSGeometry.h"
#include "AliHLTLogging.h"

class TClonesArray;
class TString;
//class AliPHOSDigit;
//class AliPHOSRecoParamEmc;
//class AliPHOSRecoParam;

/** 
 * @class AliHLTCaloClusterizer
 * Clusterizer for CALO HLT. The clusterizer takes digits as input, either
 * in the form of a container of AliHLTCaloDigitDataStruct or a
 * TClonesArray of AliPHOSDigit through an instance of a AliPHOSLoader
 *
 * @ingroup alihlt_calo
 */


class AliHLTCaloClusterizer : public AliHLTCaloConstantsHandler, public AliHLTLogging
{
  
public:
  
  /** Constructor */
  AliHLTCaloClusterizer(TString det);    

  /** Destructor */
  virtual ~AliHLTCaloClusterizer();
  
  /** Set digit container */
  void SetDigitContainer(AliHLTCaloDigitContainerDataStruct* digitContainerPtr)
  { fDigitContainerPtr = digitContainerPtr; }

  /** Set array with digits */
  void SetDigitArray(AliHLTCaloDigitDataStruct **digitPointerArr)
  { fDigitsPointerArray = digitPointerArr; } 

  /** Set rec point data buffer */
  void SetRecPointDataPtr(AliHLTCaloRecPointDataStruct* recPointDataPtr);

  /** Set reco parameters */
  //  void SetRecoParameters(AliPHOSRecoParam* recoPars);

  /** Set emc clustering threshold */
  void SetEmcClusteringThreshold(Float_t threshold) { fEmcClusteringThreshold = threshold; }

  /** Set emc min energy threshold */
  void SetEmcMinEnergyThreshold(Float_t threshold) { fEmcMinEnergyThreshold = threshold; }

  /** Set emc time gate */
  void SetEmcTimeGate(Float_t gate) { fEmcTimeGate = gate; }
  
  /** Starts clusterization of the event */ 
  virtual Int_t ClusterizeEvent(Int_t nDigits);
  
  /**
   * For a given digit this digit scans for neighbouring digits which 
   * passes the threshold for inclusion in a rec point. If one is found 
   * it is added to the current rec point
   * @param digIndex index of the digit in the digit container
   * @param recPoint pointer to the current rec point
   */
  virtual Int_t ScanForNeighbourDigits(Int_t digIndex, AliHLTCaloRecPointDataStruct* recPoint);

  /**
   * Checks if two digits are neighbours
   * @param d1 first digit
   * @param d2 second digit
   */
  virtual Int_t AreNeighbours(AliHLTCaloDigitDataStruct* d1, AliHLTCaloDigitDataStruct* d2);

  /**
  * Get pointer to the rec points array
  */
  AliHLTCaloRecPointDataStruct** GetRecPoints() const { return fRecPointArray; }

  /** 
  * Sort the digits by energy
  */
  void SetSortDigitsByEnergy();
  
  /** 
  * Sort the digits by position
  */
  void SetSortDigitsByPosition();
  
  /** 
  * Set the sorting function (as required by stdlib's qsort) if you don't want to use the provided ones 
  */
  void SetSortingFunction(Int_t (*compare)(const void*, const void*)) { fCompareFunction = compare; }
  
  
  
protected:

   /** 
   * Check the rec point buffer size and resize the buffer if necessary
   */
  virtual Int_t CheckBuffer(); //COMMENT
  
   /** 
   * Check the rec point array size and resize the array if necessary
   */
  virtual Int_t CheckArray(); //COMMENT
  
  /** 
  * Sort the digits
  */
  void SortDigits();

  /** 
  * Compare digits by position
  */
  static Int_t CompareDigitsByPosition(const void *dig0, const void *dig);
  
  /** 
  * Compare digits by energy
  */
  static Int_t CompareDigitsByEnergy(const void *dig0, const void *dig);
  
  /** 
  * Pointer to the compare function for the sorting of digits
  */
  //Int_t (AliHLTCaloClusterizer::*fCompareFunction)(const void*, const void*);
  Int_t (*fCompareFunction)(const void*, const void*);
  
  /** Array of pointers to the rec point output */
  AliHLTCaloRecPointDataStruct **fRecPointArray; //COMMENT

   /** Pointer to the rec point output */
  AliHLTCaloRecPointDataStruct *fRecPointDataPtr; //COMMENT

  /** The first rec point in the list */
  AliHLTCaloRecPointDataStruct *fFirstRecPointPtr; //COMMENT

  /** Size of the rec point array */
  Int_t fArraySize;
  
  /** Available size for the rec point output */
  Int_t fAvailableSize;

  /** The used size for the rec point output */
  Int_t fUsedSize;
  
  /** Number of rec points created so far */
  Int_t fNRecPoints;
  
  /** Pointer to the digit index array in the rec point */
  Int_t* fDigitIndexPtr;                                       //! transient

  /** Energy threshold for starting a cluster for the calorimeter */
  Float_t fEmcClusteringThreshold;                             //COMMENT

  /** Energy threshold for including a crystal in a cluster */
  Float_t fEmcMinEnergyThreshold;                              //COMMENT

  /** Maximum time difference for inclusion in a rec point */
  Float_t fEmcTimeGate;                                        //COMMENT

  /** Counts the digits in a rec point */
  Int_t fDigitsInCluster;                                      //COMMENT

  /** Array of our digits */
  AliHLTCaloDigitDataStruct **fDigitsPointerArray;             //! transient

  /** Contains the digits from one event */
  AliHLTCaloDigitContainerDataStruct *fDigitContainerPtr;      //! transient

  /** Maximum difference in index to be a neighbour */
  Int_t fMaxDigitIndexDiff;                                    //COMMENT

  /** Number of digits in event */
  Int_t fNDigits;                                              //COMMENT
  
  /** Are we sorting digits by position? */
  Bool_t fSortedByPosition; //COMMENT

  /** Are we sorting digits by energy? */
  Bool_t fSortedByEnergy; //COMMENT

   /** Are we sorting at all? */
   Bool_t fSortDigits; //COMMENT

private:

  /** Default constructor, prohibited */
  AliHLTCaloClusterizer();                          // COMMENT
  
  /** Copy constructor, prohibited */
  AliHLTCaloClusterizer (const AliHLTCaloClusterizer &); //COMMENT
  
  /** Assignment operator, prohibited */
  AliHLTCaloClusterizer & operator = (const AliHLTCaloClusterizer &); //COMMENT

  ClassDef(AliHLTCaloClusterizer, 0);

};

#endif
