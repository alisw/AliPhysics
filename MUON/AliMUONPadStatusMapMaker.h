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

class AliMUONCalibrationData;
class AliMUONV2DStore;
class AliMUONVCalibParam;
class AliMpPad;
class AliMpVSegmentation;
class TObjArray;
class TVector2;

class AliMUONPadStatusMapMaker : public TObject
{
public:
  AliMUONPadStatusMapMaker(const AliMUONCalibrationData& calibData);
  virtual ~AliMUONPadStatusMapMaker();
  
  /// Return status bit map to tell a pad is bad
  static Int_t SelfDeadMask() { return fgkSelfDead; }
  
  AliMUONV2DStore* MakePadStatusMap(const AliMUONV2DStore& status,
                                    Int_t mask);
  
  static AliMUONV2DStore* MakeEmptyPadStatusMap();


private:
  /// Not implemented
  AliMUONPadStatusMapMaker(const AliMUONPadStatusMapMaker&);
  /// Not implemented
  AliMUONPadStatusMapMaker& operator=(const AliMUONPadStatusMapMaker&);

private:
  Int_t ComputeStatusMap(const AliMUONVCalibParam& neighbours, 
                        Int_t manuChannel, Int_t detElemId) const;
  
  Int_t GetPadStatus(Int_t detElemId, Int_t manuId, Int_t manuChannel) const;
  
private:
    static Int_t fgkSelfDead; //!< status bit map to tell a pad is bad
  const AliMUONV2DStore* fStatus; //!< status store
  Int_t fMask; //!< mask to be tested

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
  const AliMUONCalibrationData& fCalibrationData; //!< to access neighbourStore
  
  ClassDef(AliMUONPadStatusMapMaker,0) // Pad status map maker
};

#endif
