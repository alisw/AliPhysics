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

#ifndef ALIHLTPHOSCLUSTERANALYSER_H
#define ALIHLTPHOSCLUSTERANALYSER_H

/**
 * Class calculates properties of rec points
 *
 * @file   AliHLTPHOSClusterAnalyser.h
 * @author Oystein Djuvsland
 * @date
 * @brief  Cluster analyser for PHOS HLT
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTPHOSBase.h"

class AliHLTPHOSPhysicsAnalyzer;
class AliHLTPHOSRecPointHeaderStruct;
class AliHLTPHOSRecPointDataStruct;
class AliHLTPHOSCaloClusterHeaderStruct;
class AliHLTPHOSCaloClusterDataStruct;
class AliPHOSGeometry;

/** 
 * @class AliHLTPHOSClusterAnalyser
 * ClusterAnalyser for PHOS HLT. Algorithms for center of gravity
 * and moment calculations are based on classes from the PHOS
 * offline analysis directory
 *
 * @ingroup alihlt_phos
 */
class AliHLTPHOSClusterAnalyser : public AliHLTPHOSBase
{
public:

  /** Constructor */
  AliHLTPHOSClusterAnalyser();

  /** Destructor */
  virtual ~AliHLTPHOSClusterAnalyser();
  
  /** Copy constructor */
  AliHLTPHOSClusterAnalyser(const AliHLTPHOSClusterAnalyser &) : 
    AliHLTPHOSBase(),
    fLogWeight(0),
    fRecPointDataPtr(0),
    fNRecPoints(0),
    fCaloClusterDataPtr(0),
    fCaloClusterHeaderPtr(0),
    fPHOSGeometry(0),
    fAnalyzerPtr(0),
    fDoClusterFit(false),
    fHaveCPVInfo(false),
    fDoPID(false),
    fHaveDistanceToBadChannel(false)
    
  {
    //Copy constructor not implemented
  }
  
  /** Assignment */
  AliHLTPHOSClusterAnalyser & operator = (const AliHLTPHOSClusterAnalyser)
    {
      //Assignment
      return *this; 
    }
  
  /**
   * Set the rec point data buffer
   * @param recPointDataPtr is a pointer to the rec points
   */
  void SetRecPointDataPtr(AliHLTPHOSRecPointHeaderStruct *recPointDataPtr);


  /** 
   * Set the calo cluster output buffer
   * @param caloClusterDataPtr is a pointer to the calo cluster buffer
   */
  void SetCaloClusterDataPtr(AliHLTPHOSCaloClusterDataStruct *caloClusterDataPtr);

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
  Int_t CalculateClusterMoments(AliHLTPHOSRecPointDataStruct *recPointPtr, AliHLTPHOSCaloClusterDataStruct* clusterPtr);

  /** 
   * Deconvolute the clusters in an AliHLTPHOSRecPointContainerStruct
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
  Int_t FitCluster(AliHLTPHOSRecPointDataStruct* /*recPointPtr*/) { return 0; }

  /**
   * Get the distance to the nearest CPV rec point
   * param recPointPtr is the pointer to the emc rec point
   * @return the distance
   */
  Float_t GetCPVDistance(AliHLTPHOSRecPointDataStruct* /*recPointPtr*/) { return 0; };

  /**
   * Do partice identification
   * param clusterPtr is the pointer to the emc cluster
   * @return 
   */
  Int_t DoParticleIdentification(AliHLTPHOSCaloClusterDataStruct* /*clusterPtr*/) { return 0; }
  
  /**
   * Get the distance to the neares bad channel
   * param clusterPtr is a pointer to the calo cluster
   * @return the distance
   */
  Float_t GetDistanceToBadChannel(AliHLTPHOSCaloClusterDataStruct* /*clusterPtr*/) { return 0; }

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
  AliHLTPHOSRecPointDataStruct *fRecPointDataPtr;         //! transient

  /** Number of rec points */
  Int_t fNRecPoints;                                      //COMMENT

  /** Pointer to the cluster buffer */
  AliHLTPHOSCaloClusterDataStruct *fCaloClusterDataPtr;   //! transient

  /** Pointer to the cluster header */
  AliHLTPHOSCaloClusterHeaderStruct *fCaloClusterHeaderPtr;   //! transient

  /** Instance of the PHOS geometry */
  AliPHOSGeometry *fPHOSGeometry;                           //! transient

  //TODO: should not use PhysicsAnalyzer for global coord!
  /** */
  AliHLTPHOSPhysicsAnalyzer *fAnalyzerPtr;                  //! transient

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
