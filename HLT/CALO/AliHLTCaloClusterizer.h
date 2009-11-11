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
#include "AliHLTCaloConstantsHandler.h"
#include "TString.h"

//#include "AliPHOSGeometry.h"

class TClonesArray;
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
//class AliHLTCaloClusterizer : public AliHLTCaloBase



class AliHLTCaloClusterizer : public AliHLTCaloConstantsHandler
{
  
public:
  
  /** Constructor */
  AliHLTCaloClusterizer(TString det);    
  
  /** Destructor */
  virtual ~AliHLTCaloClusterizer();

//   /** Copy constructor */  
//   AliHLTCaloClusterizer(const AliHLTCaloClusterizer &) : 
//     //    AliHLTCaloBase(),
//     AliHLTCaloConstantsHandler(new TString("BALLE")),
//     fRecPointDataPtr(0),
//     fDigitDataPtr(0),
//     fEmcClusteringThreshold(0),
//     fEmcMinEnergyThreshold(0),
//     fEmcTimeGate(0),
//     fDigitsInCluster(0),
//     fDigitContainerPtr(0),
//     fMaxDigitIndexDiff(2*NZROWSMOD)
//   {
//     //Copy constructor not implemented
//   }
  
  /** Assignment */
  AliHLTCaloClusterizer & operator = (const AliHLTCaloClusterizer)
  {
    //Assignment
    return *this; 
  }
  
  /** Set digit container */
  void SetDigitContainer(AliHLTCaloDigitContainerDataStruct* digitContainerPtr)
  { fDigitContainerPtr = digitContainerPtr; }

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
  virtual Int_t ClusterizeEvent(UInt_t availableSize, UInt_t& totSize);
  
  /**
   * For a given digit this digit scans for neighbouring digits which 
   * passes the threshold for inclusion in a rec point. If one is found 
   * it is added to the current rec point
   * @param digIndex index of the digit in the digit container
   * @param recPoint pointer to the current rec point
   */
  virtual void ScanForNeighbourDigits(Int_t digIndex, AliHLTCaloRecPointDataStruct* recPoint);

  /**
   * Checks if two digits are neighbours
   * @param d1 first digit
   * @param d2 second digit
   */
  virtual Int_t AreNeighbours(AliHLTCaloDigitDataStruct* d1, AliHLTCaloDigitDataStruct* d2);


protected:

  /** Pointer to the rec point output */
  AliHLTCaloRecPointDataStruct* fRecPointDataPtr;              //! transient

  /** Pointer to the digit output */
  AliHLTCaloDigitDataStruct* fDigitDataPtr;                    //! transient

  /** Energy threshold for starting a cluster for the calorimeter */
  Float_t fEmcClusteringThreshold;                             //COMMENT

  /** Energy threshold for including a crystal in a cluster */
  Float_t fEmcMinEnergyThreshold;                              //COMMENT

  /** Maximum time difference for inclusion in a rec point */
  Float_t fEmcTimeGate;                                        //COMMENT

  /** Counts the digits in a rec point */
  Int_t fDigitsInCluster;                                      //COMMENT

  /** Contains the digits from one event */
  AliHLTCaloDigitContainerDataStruct *fDigitContainerPtr;      //! transient

  /** Maximum difference in index to be a neighbour */
  Int_t fMaxDigitIndexDiff;                                    //COMMENT

private:
  AliHLTCaloClusterizer();

  AliHLTCaloClusterizer (const AliHLTCaloClusterizer  & );
  //  AliHLTCaloClusterizer & operator = (const AliHLTCaloClusterizer &);


  ClassDef(AliHLTCaloClusterizer, 0);


};

#endif
