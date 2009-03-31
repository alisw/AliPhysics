#ifndef ALIMPFASTSEGMENTATION_H
#define ALIMPFASTSEGMENTATION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup management
/// \class AliMpFastSegmentation
/// \brief Fast version of AliMpVSegmentation
/// 
// author Laurent Aphecetche

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

#ifndef ALI_MP_V_SEGMENTATION_H
#  include "AliMpVSegmentation.h"
#endif

#ifndef ROOT_TObjArray
#  include "TObjArray.h"
#endif

#ifndef ROOT_TExMap
#  include "TExMap.h"
#endif

#ifndef ROOT_TVector2
#  include "TVector2.h"
#endif

class AliMpMotifPosition;

class AliMpFastSegmentation : public AliMpVSegmentation
{
public:
  AliMpFastSegmentation(AliMpVSegmentation* seg);
  virtual ~AliMpFastSegmentation();
  
  virtual AliMpVPadIterator* CreateIterator(const AliMpArea& area) const;
  virtual AliMpVPadIterator* CreateIterator() const;
  
  virtual Int_t GetNeighbours(const AliMpPad& pad, TObjArray& neighbours,
                              Bool_t includeSelf=kFALSE,
                              Bool_t includeVoid=kFALSE) const;
  
  virtual Bool_t HasPadByIndices(Int_t ix, Int_t iy) const;
  virtual Bool_t HasPadByLocation(Int_t manuId, Int_t manuChannel) const;
  
  virtual AliMpPad PadByLocation(Int_t manuId, Int_t manuChannel, Bool_t warning = true) const;
  virtual AliMpPad PadByIndices (Int_t ix, Int_t iy, Bool_t warning = true) const;
  virtual AliMpPad PadByPosition(const TVector2& position, Bool_t warning = true) const;
  
  virtual Int_t  MaxPadIndexX() const;
  virtual Int_t  MaxPadIndexY() const;
  virtual Int_t  NofPads() const;
  
  virtual void GetAllElectronicCardIDs(TArrayI& ecn) const;
  
  virtual Int_t GetNofElectronicCards() const;
  
  virtual AliMp::PlaneType PlaneType() const;
  
  virtual TVector2 Dimensions() const;
    
  virtual TVector2 Position() const;

  virtual AliMpMotifPosition* MotifPosition(Int_t manuId) const;

  virtual Bool_t HasMotifPosition(Int_t manuId) const;

  virtual void Print(Option_t* opt="") const;

  /// Return helper class 
  AliMpVSegmentation* GetHelper() const { return fHelper; }
  
  /// Return segmentation station type
  AliMp::StationType StationType() const { return fHelper->StationType(); }
  
private:
  AliMpFastSegmentation(const AliMpFastSegmentation& rhs);
  AliMpFastSegmentation& operator=(const AliMpFastSegmentation& rhs);

  virtual AliMpMotifPosition* InternalMotifPosition(Int_t index) const;

private:
  AliMpVSegmentation* fHelper; ///< helper class (owner)
  TObjArray fMotifPositions; ///< array of AliMpMotifPositions (not owner)
  mutable TExMap fIxIy; ///< map of (ix,iy) -> index in array above
  mutable TExMap fManuId; ///< map of (manuid) -> index in array above
  TVector2 fPosition; ///< to compute pad positions
  
  ClassDef(AliMpFastSegmentation,1) // Variant implementation for AliMpVSegmentation
};

#endif
