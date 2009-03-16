#ifndef ALIMUONVCLUSTERSERVER_H
#define ALIMUONVCLUSTERSERVER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONVClusterServer
/// \brief Interface of a cluster finder for combined tracking
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#include "TObject.h"
#endif

class AliMUONVDigitStore;
class AliMUONVClusterStore;
class AliMUONVTriggerTrackStore;
class AliMUONRecoParam;
class AliMpArea;
class TIter;

class AliMUONVClusterServer : public TObject
{
public:
  AliMUONVClusterServer();
  virtual ~AliMUONVClusterServer();
  
  /// Find and add clusters from a given region of a given chamber to the store.
  virtual Int_t Clusterize(Int_t chamberId, 
                           AliMUONVClusterStore& clusterStore,
                           const AliMpArea& area,
			   const AliMUONRecoParam* recoParam = 0x0) = 0;
  
  /// Specify an iterator to loop over the digits needed to perform our job.
  virtual void UseDigits(TIter& next, AliMUONVDigitStore* digitStore = 0x0) = 0;
  
  /// Use trigger tracks. Return kFALSE if not used.
  virtual Bool_t UseTriggerTrackStore(AliMUONVTriggerTrackStore* /*trackStore*/) { return kFALSE; }
  
  ClassDef(AliMUONVClusterServer,1) // Cluster server interface
};

#endif
