#ifndef ALIMUONVCLUSTERFINDER_H
#define ALIMUONVCLUSTERFINDER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONClusterFinder
/// \brief Interface of a cluster finder.
/// 
//  Author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONCluster;
class AliMpVSegmentation;
class TClonesArray;
class AliMUONPad;

class AliMUONVClusterFinder : public TObject
{
public:
  AliMUONVClusterFinder();
  virtual ~AliMUONVClusterFinder();
  
  virtual Bool_t Prepare(const AliMpVSegmentation* segmentations[2],
                         TClonesArray* digits[2]) = 0;
  
  virtual AliMUONCluster* NextCluster() = 0;
  
  /** Add a pad to the list of pads to be considered for clustering.
    Typical usage is to "put-back-in-business" a pad that was part 
    of a previous cluster (returned by NextCluster) but was externally
    identified of not being part of that cluster, so it must be reuseable.
    Might not be implemented by all cluster finders...
    (in which case it must returns kFALSE)
    */
  virtual Bool_t UsePad(const AliMUONPad& pad);
  
  ClassDef(AliMUONVClusterFinder,0) // Interface of a MUON cluster finder.
};

#endif
