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

#include "Rtypes.h"

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
class AliHLTCaloClusterHeaderStruct;
class AliHLTCaloClusterDataStruct;
class AliPHOSGeoUtils;

/** 
 * @class AliHLTCaloClusterAnalyser
 * ClusterAnalyser for CALO HLT. Algorithms for center of gravity
 * and moment calculations are based on classes from the CALO
 * offline analysis directory
 *
 * @ingroup alihlt_calo
 */
//class AliHLTCaloClusterAnalyser : public AliHLTCaloBase
class AliHLTCaloClusterAnalyser 
{
public:

  /** Constructor */
  AliHLTCaloClusterAnalyser();

  /** Destructor */
  virtual ~AliHLTCaloClusterAnalyser();
  
  /** Copy constructor */
  AliHLTCaloClusterAnalyser(const AliHLTCaloClusterAnalyser &) : 
    //   AliHLTCaloBase(),
    fLogWeight(0),
    fRecPointDataPtr(0),
    fNRecPoints(0),
    fCaloClusterDataPtr(0),
    fCaloClusterHeaderPtr(0),
    fPHOSGeometry(0),
    //fAnalyzerPtr(0),
    fDoClusterFit(false),
    fHaveCPVInfo(false),
    fDoPID(false),
    fHaveDistanceToBadChannel(false)
    
  {
    //Copy constructor not implemented
  }
  
  /** Assignment */
  AliHLTCaloClusterAnalyser & operator = (const AliHLTCaloClusterAnalyser)
    {
      //Assignment
      return *this; 
    }
  
  /**
   * Set the rec point data buffer
   * @param recPointDataPtr is a pointer to the rec points
   */
  void SetRecPointDataPtr(AliHLTCaloRecPointHeaderStruct *recPointDataPtr);

  /** 
   * Set the calo cluster output buffer
   * @param caloClusterDataPtr is a pointer to the calo cluster buffer
   */
  void SetCaloClusterDataPtr(AliHLTCaloClusterDataStruct *caloClusterDataPtr);

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
  Int_t CreateClusters(UInt_t availableSize, UInt_t& totSize);

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

private:
  
  /** Used for calculation of center of gravity */
  Float_t fLogWeight;                                       //COMMENT
  
  /** Pointer to the rec points */
  AliHLTCaloRecPointDataStruct *fRecPointDataPtr;         //! transient

  /** Number of rec points */
  Int_t fNRecPoints;                                      //COMMENT

  /** Pointer to the cluster buffer */
  AliHLTCaloClusterDataStruct *fCaloClusterDataPtr;   //! transient

  /** Pointer to the cluster header */
  AliHLTCaloClusterHeaderStruct *fCaloClusterHeaderPtr;   //! transient

  /** Instance of the PHOS geometry */
  AliPHOSGeoUtils *fPHOSGeometry;                           //! transient

  //TODO: should not use PhysicsAnalyzer for global coord!
  /** */
  //  AliHLTCaloPhysicsAnalyzer *fAnalyzerPtr;                  //! transient

  /** Should we do cluster fitting? */
  Bool_t fDoClusterFit;                                     //COMMENT
  
  /** Do we have CPV info? */
  Bool_t fHaveCPVInfo;                                      //COMMENT

  /** Should we do PID? */ 
  Bool_t fDoPID;                                            //COMMENT

  /** Do we have distance to bad channel? */
  Bool_t fHaveDistanceToBadChannel;                         //COMMENT
  
  
};

#endif
