#ifndef ALIMUONCLUSTERFINDERCOG_H
#define ALIMUONCLUSTERFINDERCOG_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONClusterFinderCOG
/// \brief A very basic (and mostly useless, probably) cluster finder
/// 
// Author Laurent Aphecetche, Subatech

#ifndef AliMUONVCLUSTERFINDER_H
#  include "AliMUONVClusterFinder.h"
#endif

class AliMUONMathieson;

class AliMUONClusterFinderCOG : public AliMUONVClusterFinder
{
public:
  AliMUONClusterFinderCOG();
  virtual ~AliMUONClusterFinderCOG();
  
  virtual Bool_t Prepare(const AliMpVSegmentation* segmentations[2],
                         TClonesArray* digits[2]);
  
  virtual AliMUONCluster* NextCluster();
  
private:
  AliMUONClusterFinderCOG(const AliMUONClusterFinderCOG& rhs);
  AliMUONClusterFinderCOG& operator=(const AliMUONClusterFinderCOG& rhs);
  void ComputePosition(AliMUONCluster& cluster);

private:
    AliMUONVClusterFinder* fPreClusterFinder; ///< the preclustering we use

  ClassDef(AliMUONClusterFinderCOG,1) // 
};

#endif
