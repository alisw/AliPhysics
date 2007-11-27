#ifndef ALIMUONSIMPLECLUSTERSERVER_H
#define ALIMUONSIMPLECLUSTERSERVER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup reco
/// \class AliMUONSimpleClusterServer
/// \brief Implementation of AliMUONVClusterServer interface
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ALIMUONVCLUSTERSERVER_H
#  include "AliMUONVClusterServer.h"
#endif

class AliMUONVClusterFinder;
class AliMUONGeometryTransformer;
class TClonesArray;
class AliMpExMap;

class AliMUONSimpleClusterServer : public AliMUONVClusterServer
{
public:
  AliMUONSimpleClusterServer(AliMUONVClusterFinder& clusterFinder,
                             const AliMUONGeometryTransformer& transformer);
  
  virtual ~AliMUONSimpleClusterServer();
  
  Int_t Clusterize(Int_t chamberId,
                   AliMUONVClusterStore& clusterStore,
                   const AliMpArea& area);
  
  void UseDigitStore(const AliMUONVDigitStore& digitStore);
  
  void Print(Option_t* opt="") const;
  
private:
  /// Not implemented
  AliMUONSimpleClusterServer(const AliMUONSimpleClusterServer& rhs);
  /// Not implemented
  AliMUONSimpleClusterServer& operator=(const AliMUONSimpleClusterServer& rhs);
  
  Bool_t Overlap(Int_t detElemId, const AliMpArea& area, AliMpArea& deArea) const;
    
  void Global2Local(Int_t detElemId, const AliMpArea& globalArea, AliMpArea& localArea) const;

  TClonesArray* PadArray(Int_t detElemId, Int_t cathode) const;
  
private:
  AliMUONVClusterFinder& fClusterFinder; //!< the cluster finder
  const AliMUONGeometryTransformer& fTransformer; //!< the geometry transformer
  AliMpExMap* fPads[2]; ///< map of TClonesArray of AliMUONPads
  AliMpExMap* fDEAreas; ///< map of detection element areas in global coordinates
  
  ClassDef(AliMUONSimpleClusterServer,0) // Cluster server
};

#endif
