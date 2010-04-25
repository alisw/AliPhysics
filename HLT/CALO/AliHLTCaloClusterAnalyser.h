//-*- Mode: C++ -*-
// $Id: AliHLTCaloClusterAnalyser.h 35107 2009-09-30 01:45:06Z phille $

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

#ifndef ALIHLTCALOCLUSTERANALYSER_H
#define ALIHLTCALOCLUSTERANALYSER_H

#include "AliHLTLogging.h"

/**
 * Class calculates properties of rec points
 *
 * @file   AliHLTCaloClusterAnalyser.h
 * @author Oystein Djuvsland
 * @date
 * @brief  Cluster analyser for CALO HLT
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

//#include "AliHLTCaloBase.h"

//class AliHLTCaloPhysicsAnalyzer;
class AliHLTCaloRecPointHeaderStruct;
class AliHLTCaloRecPointDataStruct;
class AliHLTCaloDigitDataStruct;
class AliHLTCaloClusterHeaderStruct;
class AliHLTCaloClusterDataStruct;
class AliHLTCaloGeometry;
class AliHLTCaloRecoParamHandler;

class TH1F; //DEBUG


/** 
 * @class AliHLTCaloClusterAnalyser
 * ClusterAnalyser for CALO HLT. Algorithms for center of gravity
 * and moment calculations are based on classes from the CALO
 * offline analysis directory
 *
 * @ingroup alihlt_calo
 */
//class AliHLTCaloClusterAnalyser : public AliHLTCaloBase
class AliHLTCaloClusterAnalyser : public AliHLTLogging
{
public:

  /** Constructor */
  AliHLTCaloClusterAnalyser();

  /** Destructor */
  virtual ~AliHLTCaloClusterAnalyser();
  
  /**
   * Set the rec point data buffer
   * @param recPointDataPtr is a pointer to the rec points
   */
  void SetRecPointArray(AliHLTCaloRecPointDataStruct **recPointDataPtr, Int_t nRecPoints);

  /** 
   * Set the calo cluster output buffer
   * @param caloClusterDataPtr is a pointer to the calo cluster buffer
   */
  void SetCaloClusterData(AliHLTCaloClusterDataStruct *caloClusterDataPtr);

  /** 
   * Calculates the center of gravity for the reconstruction points in the container
   * @return
   */
  Int_t CalculateCenterOfGravity();

  /** 
   * Calculates the moments for the reconstruction points in the container
   * @return 
   */
  Int_t CalculateRecPointMoments();

  /** 
   * Calculates the moments for a certain cluster
   * @return 
   */
  Int_t CalculateClusterMoments(AliHLTCaloRecPointDataStruct *recPointPtr, AliHLTCaloClusterDataStruct* clusterPtr);

  /** 
   * Deconvolute the clusters in an AliHLTCaloRecPointContainerStruct
   * @return
   */
  Int_t DeconvoluteClusters();

  /**
   * Convert the rec points into calo clusters
   * @return
   */
  Int_t CreateClusters(Int_t nRecPoints, UInt_t availableSize, UInt_t& totSize);

  /**
   * Fit a cluster
   * param recPointPtr is a pointer to the rec point to fit
   * @return 
   */
  Int_t FitCluster(AliHLTCaloRecPointDataStruct* /*recPointPtr*/) { return 0; }

  /**
   * Get the distance to the nearest CPV rec point
   * param recPointPtr is the pointer to the emc rec point
   * @return the distance
   */
  Float_t GetCPVDistance(AliHLTCaloRecPointDataStruct* /*recPointPtr*/) { return 0; };

  /**
   * Do partice identification
   * param clusterPtr is the pointer to the emc cluster
   * @return 
   */
  Int_t DoParticleIdentification(AliHLTCaloClusterDataStruct* /*clusterPtr*/) { return 0; }
  
  /**
   * Get the distance to the neares bad channel
   * param clusterPtr is a pointer to the calo cluster
   * @return the distance
   */
  Float_t GetDistanceToBadChannel(AliHLTCaloClusterDataStruct* /*clusterPtr*/) { return 0; }

  /**
   * Set do cluster fit
   */
  void SetDoClusterFit() { fDoClusterFit = true; }
  
  /**
   * Set have cpv info
   */
  void SetHaveCPVInfo() { fHaveCPVInfo = true; }

  /** 
   * Set do PID
   */
  void SetDoPID() { fDoPID = true; }

  /**
   * Set have distance to bad channel
   */
  void SetHaveDistanceToBadChannel() { fHaveDistanceToBadChannel = true; }

  /**
  * Set the geometry object (different for EMCAL and PHOS)
  */
  void SetGeometry(AliHLTCaloGeometry *geometry) { fGeometry = geometry; }

  /** 
  * Set pointer to the digits
  */
  void SetDigitDataArray(AliHLTCaloDigitDataStruct *digits);

  /**
  * Set the cluster type 
  */
  void SetClusterType(Char_t clusterType) { fClusterType = clusterType; }
  
  /** 
  * Set the reconstruction parameters handler 
  */
  void SetRecoParamHandler(AliHLTCaloRecoParamHandler *recoParams) { fRecoParamsPtr = recoParams; }
  
  /** Set cut on single cell clusters */
  void SetCutOnSingleCellClusters(Bool_t doCut, Float_t energyCut) {fCutOnSingleCellClusters = doCut; fSingleCellEnergyCut = energyCut; }
  
 
private:
  
  /** Used for calculation of center of gravity */
  Float_t fLogWeight;                                       //COMMENT
  
  /** Pointer to the rec points */
  AliHLTCaloRecPointDataStruct **fRecPointArray;         //! transient

  /** Pointer to the digits */
  AliHLTCaloDigitDataStruct *fDigitDataArray;         //! transient

  /** Number of rec points */
  Int_t fNRecPoints;                                      //COMMENT

  /** Pointer to the cluster buffer */
  AliHLTCaloClusterDataStruct *fCaloClusterDataPtr;   //! transient

  /** Pointer to the cluster header */
  AliHLTCaloClusterHeaderStruct *fCaloClusterHeaderPtr;   //! transient

  /** Should we do cluster fitting? */
  Bool_t fDoClusterFit;                                     //COMMENT
  
  /** Do we have CPV info? */
  Bool_t fHaveCPVInfo;                                      //COMMENT

  /** Should we do PID? */ 
  Bool_t fDoPID;                                            //COMMENT

  /** Do we have distance to bad channel? */
  Bool_t fHaveDistanceToBadChannel;                         //COMMENT
  
  /** The geometry object */
  AliHLTCaloGeometry* fGeometry;                                   //! transient

  /** The cluster type */
  Char_t fClusterType;                   //COMMENT
  
  /** Handler to get the hold of reconstruction parameters */
  AliHLTCaloRecoParamHandler *fRecoParamsPtr; //COMMENT
  
  /** Should we cut out single celled high energy clusters? */
  Bool_t fCutOnSingleCellClusters;       //COMMENT
  
  /** If we cut on single celled clusters, what is our energy cut? */
  Float_t fSingleCellEnergyCut; //COMMENT
 
 /** Copy constructor  not implemented */
 AliHLTCaloClusterAnalyser ( const AliHLTCaloClusterAnalyser &); // not implemented
    
 /** Assignment */
AliHLTCaloClusterAnalyser & operator = ( const AliHLTCaloClusterAnalyser &); // not implemented
    
};

#endif
