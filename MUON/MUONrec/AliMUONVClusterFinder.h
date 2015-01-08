#ifndef ALIMUONVCLUSTERFINDER_H
#define ALIMUONVCLUSTERFINDER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONVClusterFinder
/// \brief Interface of a cluster finder.
/// 
//  Author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONCluster;
class AliMUONRecoParam;
class AliMpVSegmentation;
class AliMUONPad;
class AliMpArea;

class AliMUONVClusterFinder : public TObject
{
public:
  AliMUONVClusterFinder();
  virtual ~AliMUONVClusterFinder();
  
  /// \todo add comment

  virtual Bool_t NeedSegmentation() const { return kFALSE; }
  
  virtual Bool_t Prepare(Int_t detElemId,
                         TObjArray* pads[2],
                         const AliMpArea& area);

  virtual Bool_t Prepare(Int_t detElemId,
                         TObjArray* pads[2],
                         const AliMpArea& area,
                         const AliMpVSegmentation* segmentations[2]);  
  
  /// \todo add comment
  virtual AliMUONCluster* NextCluster() = 0;
  
  /** Add a pad to the list of pads to be considered for clustering.
    Typical usage is to "put-back-in-business" a pad that was part 
    of a previous cluster (returned by NextCluster) but was externally
    identified of not being part of that cluster, so it must be reuseable.
    Might not be implemented by all cluster finders...
    (in which case it must returns kFALSE)
    */
  virtual Bool_t UsePad(const AliMUONPad& pad);
  
  /** Specify a couple of charge hints. We call them hints because some
   clustering need them and use them directly, other cook them before
   using them, and some others yet simply don't care about them.
   */
  virtual void SetChargeHints(Double_t /*lowestPadCharge*/, Double_t /*lowestClusterCharge*/) { }
  
  ClassDef(AliMUONVClusterFinder,0) // Interface of a MUON cluster finder.
};

#endif
