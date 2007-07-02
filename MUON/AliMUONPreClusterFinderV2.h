#ifndef ALIMUONPRECLUSTERFINDERV2_H
#define ALIMUONPRECLUSTERFINDERV2_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONPreClusterFinderV2
/// \brief A basic pre-cluster finder
/// 
// Author Laurent Aphecetche, Subatech

#ifndef AliMUONVCLUSTERFINDER_H
#  include "AliMUONVClusterFinder.h"
#endif

class TStopwatch;
class AliMUONPad;

class AliMUONPreClusterFinderV2 : public AliMUONVClusterFinder
{
public:
  AliMUONPreClusterFinderV2();
  virtual ~AliMUONPreClusterFinderV2();
  
  virtual Bool_t Prepare(const AliMpVSegmentation* segmentations[2],
                         const AliMUONVDigitStore& digitStore);
  
  virtual AliMUONCluster* NextCluster();

  virtual Bool_t UsePad(const AliMUONPad& pad);
  
private:
  /// Not implemented
  AliMUONPreClusterFinderV2(const AliMUONPreClusterFinderV2& rhs);
  /// Not implemented
  AliMUONPreClusterFinderV2& operator=(const AliMUONPreClusterFinderV2& rhs);

  void AddPad(AliMUONCluster& cluster, AliMUONPad* pad);
  
private:
  TClonesArray* fClusters; //!< the clusters we've found (owner)
  const AliMpVSegmentation** fSegmentations; //!< segmentations (not owner)
  TClonesArray* fPads[2]; //!< the pads corresponding to the digits (owner)
  Int_t fDetElemId; //!< which DE we're considering
  
  ClassDef(AliMUONPreClusterFinderV2,1) // A basic pre-cluster finder
};

#endif
