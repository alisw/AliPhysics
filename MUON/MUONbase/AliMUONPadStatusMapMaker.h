#ifndef ALIMUONPADSTATUSMAPMAKER_H
#define ALIMUONPADSTATUSMAPMAKER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup rec
/// \class AliMUONPadStatusMapMaker
/// \brief Convert a pad status container into a pad status *map* container
/// 
//  Author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONPadStatusMaker;
class AliMUONVCalibParam;
class AliMUONVStore;

class AliMUONPadStatusMapMaker : public TObject
{
public:
  AliMUONPadStatusMapMaker(const AliMUONPadStatusMaker& padStatusMaker,
                           Int_t mask,
                           Bool_t deferredInitialization=kTRUE);
  virtual ~AliMUONPadStatusMapMaker();
  
  /** Get access to internal status map store (for debug only, as it may not be complete,
    depending on whether you've already called StatusMap() for all possible de,manu,channel
    combinations or not...
    */
  AliMUONVStore* StatusMap() const { return fStatusMap; }
  
  Int_t StatusMap(Int_t detElemId, Int_t manuId, Int_t manuChannel) const;
  
  /// Return status bit map to tell a pad is bad
  static Int_t SelfDeadMask() { return fgkSelfDead; }
  
  void RefreshRejectProbabilities(); 
  
private:
  /// Not implemented
  AliMUONPadStatusMapMaker(const AliMUONPadStatusMapMaker&);
  /// Not implemented
  AliMUONPadStatusMapMaker& operator=(const AliMUONPadStatusMapMaker&);

private:
    
  AliMUONVCalibParam* ComputeStatusMap(Int_t detElemId, Int_t manuId) const;
  
private:
  
  static Int_t fgkSelfDead; //!<! status bit map to tell a pad is bad

  /// Bit numbers
  enum EBitNumbers
  {
    kLeftBottomBit = 6,
    kLeftBit = 7,
    kLeftTopBit = 8,
    kBottomBit = 11,
    kCenterBit = 12,
    kTopBit = 13,
    kRightBottomBit = 16,
    kRightBit = 17,
    kRightTopBit = 18
  };
  
  const AliMUONPadStatusMaker& fkStatusMaker; //!<! to access pad statuses
  Int_t fMask; //!<! mask to be tested
  mutable AliMUONVStore* fStatusMap; //!<! status map
  AliMUONVStore* fRejectProbabilities; //!<! reject probabilities (channel based, computed once per run)
  AliMUONVStore* fRejectList; //!<! reject list (which channels should be rejected, might change event-by-event for simulations)
  Bool_t fComputeOnDemand; //!<! whether we authorize to compute things on demand or not
  
  ClassDef(AliMUONPadStatusMapMaker,0) // Pad status map maker
};

#endif
