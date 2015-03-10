#ifndef ALIMUONCLUSTERFINDERSIMPLEFIT_H
#define ALIMUONCLUSTERFINDERSIMPLEFIT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONClusterFinderSimpleFit
/// \brief Basic cluster finder 
/// 
// Author Laurent Aphecetche, Subatech

#ifndef AliMUONVCLUSTERFINDER_H
#  include "AliMUONVClusterFinder.h"
#endif

class AliMUONMathieson;

class AliMUONClusterFinderSimpleFit : public AliMUONVClusterFinder
{
public:
  AliMUONClusterFinderSimpleFit(AliMUONVClusterFinder* clusterFinder);
  virtual ~AliMUONClusterFinderSimpleFit();
  
  using AliMUONVClusterFinder::Prepare;

  virtual Bool_t Prepare(Int_t detElemId,
                         TObjArray* pads[2],
                         const AliMpArea& area);
  
  virtual AliMUONCluster* NextCluster();
  
  virtual void SetChargeHints(Double_t /*lowestPadCharge*/, Double_t lowestClusterCharge) { 
    fLowestClusterCharge=lowestClusterCharge; 
  }
  
private:
  /// Not implemented
  AliMUONClusterFinderSimpleFit(const AliMUONClusterFinderSimpleFit& rhs);
  /// Not implemented
  AliMUONClusterFinderSimpleFit& operator=(const AliMUONClusterFinderSimpleFit& rhs);

  void ComputePosition(AliMUONCluster& cluster);

private:
    AliMUONVClusterFinder* fClusterFinder; //!<! the preclustering we use
  AliMUONMathieson* fMathieson; //!<! Mathieson to compute the charge repartition
  Double_t fLowestClusterCharge; //!<! minimum cluster charge we allow
  
  ClassDef(AliMUONClusterFinderSimpleFit,2) // Basic cluster finder
};

#endif
