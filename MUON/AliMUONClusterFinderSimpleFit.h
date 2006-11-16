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
  AliMUONClusterFinderSimpleFit();
  virtual ~AliMUONClusterFinderSimpleFit();
  
  virtual Bool_t Prepare(const AliMpVSegmentation* segmentations[2],
                         TClonesArray* digits[2]);
  
  virtual AliMUONCluster* NextCluster();
  
private:
  AliMUONClusterFinderSimpleFit(const AliMUONClusterFinderSimpleFit& rhs);
  AliMUONClusterFinderSimpleFit& operator=(const AliMUONClusterFinderSimpleFit& rhs);
  void ComputePosition(AliMUONCluster& cluster);

private:
    AliMUONVClusterFinder* fClusterFinder; //!< the preclustering we use
  AliMUONMathieson* fMathieson; //!< Mathieson to compute the charge repartition
  
  ClassDef(AliMUONClusterFinderSimpleFit,1) // 
};

#endif
