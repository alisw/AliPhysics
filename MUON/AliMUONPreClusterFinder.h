#ifndef ALIMUONPRECLUSTERFINDER_H
#define ALIMUONPRECLUSTERFINDER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONPreClusterFinder
/// \brief A basic pre-cluster finder
/// 
// Author Laurent Aphecetche, Subatech

#ifndef AliMUONVCLUSTERFINDER_H
#  include "AliMUONVClusterFinder.h"
#endif
#ifndef ALI_MP_AREA_H
#  include "AliMpArea.h"
#endif

class TStopwatch;
class AliMUONPad;

class AliMUONPreClusterFinder : public AliMUONVClusterFinder
{
public:
  AliMUONPreClusterFinder();
  virtual ~AliMUONPreClusterFinder();
  
  using AliMUONVClusterFinder::Prepare;

  virtual Bool_t Prepare(Int_t detElemId,
                         TClonesArray* pads[2],
                         const AliMpArea& area);
  
  virtual AliMUONCluster* NextCluster();

  virtual Bool_t UsePad(const AliMUONPad& pad);
  
private:
  /// Not implemented
  AliMUONPreClusterFinder(const AliMUONPreClusterFinder& rhs);
  /// Not implemented
  AliMUONPreClusterFinder& operator=(const AliMUONPreClusterFinder& rhs);

  void AddPad(AliMUONCluster& cluster, AliMUONPad* pad);
  
  AliMUONPad* GetNextPad(Int_t cathode) const;

private:
  TClonesArray* fClusters; //!< the clusters we've found (owner)
  TClonesArray** fPads; //!< the pads corresponding to the digits (not owner)
  Int_t fDetElemId; //!< which DE we're considering
  AliMpArea fArea; //!< area into which to consider pads to *start* a cluster
  
  ClassDef(AliMUONPreClusterFinder,2) // A basic pre-cluster finder
};

#endif
