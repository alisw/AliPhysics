#ifndef ALIMUONPAINTERPADSTORE_H
#define ALIMUONPAINTERPADSTORE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup graphics
/// \class AliMUONPainterPadStore
/// \brief Container for pads
/// 
// Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMUONVCalibParam;
class AliMUONVStore;
class TArrayI;
class TVector2;

class AliMUONPainterPadStore : public TObject
{
public:
  AliMUONPainterPadStore();
  AliMUONPainterPadStore(TRootIOCtor* dummy);
  virtual ~AliMUONPainterPadStore();

  Int_t FindPadID(const TArrayI& pads, Double_t x, Double_t y) const;

  AliMUONVCalibParam* Get(Int_t detElemId, Int_t manuId) const;
  
  void GetBoundaries(const TArrayI& pads, Double_t& xmin, Double_t& ymin,
                     Double_t& xmax, Double_t& ymax) const;
    
  void GetPadGeometry(Int_t padID, TVector2& position, TVector2& dimensions) const;
  
  void PrintPads(const TArrayI& pads) const;

  Int_t GetSize() const;
  
private:
  /// not implemented
  AliMUONPainterPadStore(const AliMUONPainterPadStore& rhs);
  /// not implemented
  AliMUONPainterPadStore& operator=(const AliMUONPainterPadStore& rhs);
private:
  AliMUONVStore* fPadStore; ///< the pad container
  
  ClassDef(AliMUONPainterPadStore,1) // A pad container
};

#endif
