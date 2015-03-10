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
#ifndef ROOT_TClonesArray
#  include <TClonesArray.h>
#endif

class TStopwatch;
class AliMUONPad;
class TObjArray;

class AliMUONPreClusterFinder : public AliMUONVClusterFinder
{
public:
  AliMUONPreClusterFinder();
  virtual ~AliMUONPreClusterFinder();
  
  using AliMUONVClusterFinder::Prepare;

  virtual Bool_t Prepare(Int_t detElemId,
                         TObjArray* pads[2],
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

  /// Whether we should stop working...
  virtual Bool_t ShouldAbort() const { return fShouldAbort; }
  
  AliMUONCluster* NewCluster();
  void RemoveCluster(AliMUONCluster* cluster);
  
private:
  TClonesArray fClusters; //!<! the clusters we've found (owner)
  TObjArray** fPads; //!<! the pads corresponding to the digits (not owner)
  Int_t fDetElemId; //!<! which DE we're considering
  AliMpArea fArea; //!<! area into which to consider pads to *start* a cluster
  Bool_t fShouldAbort; //!<! to indicate clustering should stop right now
  
  ClassDef(AliMUONPreClusterFinder,4) // A basic pre-cluster finder
};

#endif
