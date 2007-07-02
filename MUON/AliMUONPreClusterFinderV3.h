#ifndef ALIMUONPRECLUSTERFINDERV3_H
#define ALIMUONPRECLUSTERFINDERV3_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONPreClusterFinderV3
/// \brief A basic pre-cluster finder
/// 
// Author Laurent Aphecetche, Subatech

#ifndef AliMUONVCLUSTERFINDER_H
#  include "AliMUONVClusterFinder.h"
#endif

class TIterator;
class AliMUONPad;

class AliMUONPreClusterFinderV3 : public AliMUONVClusterFinder
{
public:
  AliMUONPreClusterFinderV3();
  virtual ~AliMUONPreClusterFinderV3();
  
  virtual Bool_t Prepare(const AliMpVSegmentation* segmentations[2],
                         const AliMUONVDigitStore& digitStore);
  
  virtual AliMUONCluster* NextCluster();

  virtual Bool_t UsePad(const AliMUONPad& pad);
  
private:
  /// Not implemented
  AliMUONPreClusterFinderV3(const AliMUONPreClusterFinderV3& rhs);
  /// Not implemented
  AliMUONPreClusterFinderV3& operator=(const AliMUONPreClusterFinderV3& rhs);

  void AddPad(AliMUONCluster& cluster, AliMUONPad* pad);
  void AddPreCluster(AliMUONCluster& cluster, AliMUONCluster* preCluster);
  void MakeCathodePreClusters(Int_t cathode);
  void MakeClusters();
  
  void DumpPreClusters() const;
  
private:
  TClonesArray* fClusters; //!< the clusters we've found (owner)
  const AliMpVSegmentation** fSegmentations; //!< segmentations (not owner)
  TClonesArray* fPads[2]; //!< the pads corresponding to the digits (owner)
  TClonesArray* fPreClusters[2]; //!< the preclusters per cathode (owner)
  Int_t fDetElemId; //!< which DE we're considering
  TIterator* fIterator; //!< iterator on fClusters
  
  ClassDef(AliMUONPreClusterFinderV3,1) // A basic pre-cluster finder
};

#endif
