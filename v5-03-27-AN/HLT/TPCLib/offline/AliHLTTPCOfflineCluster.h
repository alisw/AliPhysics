// XEmacs -*-C++-*-
// $Id$

#ifndef AliHLTTPC_OFFLINECLUSTER_H
#define AliHLTTPC_OFFLINECLUSTER_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCOfflineCluster.h
    @author Kalliopi Kanaki
    @brief  Conversion of HLT cluster format to offline format and back
*/

#include "AliHLTLogging.h"

#include "AliTPCclusterMI.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCSpacePointData.h"


/**
 * @class AliHLTTPCOfflineCluster
 *
 * The conversion of HLT cluster format to offline one
 *
 * The conversion takes place cluster per cluster. The array
 * with the offline clusters is filled in the respective
 * conversion component (AliHLTTPCClusterConverter).
 * The same is valid for the transformation of the cluster
 * coordinates. This class only translates from the AliHLTTPCSpacePoint
 * data struct to the offline AliTPCclusterMI class.
 *
 */

class AliHLTTPCOfflineCluster : public AliTPCclusterMI {

 public:

  /** standard constructor */
  AliHLTTPCOfflineCluster();  
  /** overloading the constructor */
  AliHLTTPCOfflineCluster(const AliHLTTPCSpacePointData& hltCluster);  
  AliHLTTPCOfflineCluster(const AliTPCclusterMI& offCluster);
  /** assignment operator */
  AliHLTTPCOfflineCluster& operator=(const AliTPCclusterMI& offCluster);
  /** destructor */
  virtual ~AliHLTTPCOfflineCluster();
  
  AliTPCclusterMI*          ConvertHLTToOffline(AliHLTTPCSpacePointData  spacePoint);
  AliHLTTPCSpacePointData   ConvertOfflineToHLT(AliTPCclusterMI          *offCluster);
     
  ClassDef(AliHLTTPCOfflineCluster,0) 
};
#endif
