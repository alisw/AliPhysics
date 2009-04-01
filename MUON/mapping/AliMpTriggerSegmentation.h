/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$
// $MpId: AliMpTriggerSegmentation.h,v 1.8 2006/05/24 13:58:27 ivana Exp $

/// \ingroup mptrigger
/// \class AliMpTriggerSegmentation
/// \brief Implementation of AliMpVSegmentation for trigger slats.
///
//  Author: Laurent Aphecetche

#ifndef ALI_MP_TRIGGER_SEGMENTATION_H
#define ALI_MP_TRIGGER_SEGMENTATION_H

#ifndef ROOT_TString
#include "TString.h"
#endif

#ifndef ALI_MP_V_SEGMENTATION_H
#include "AliMpVSegmentation.h"
#endif

#ifndef ALI_MP_PAD_H
#include "AliMpPad.h"
#endif

class AliMpMotifPosition;
class AliMpPCB;
class AliMpTrigger;

class AliMpTriggerSegmentation : public AliMpVSegmentation
{
public:
  AliMpTriggerSegmentation();
  AliMpTriggerSegmentation(const AliMpTrigger* slat, Bool_t own = false);
  virtual ~AliMpTriggerSegmentation();
  
  virtual AliMpVPadIterator* CreateIterator(const AliMpArea& area) const;
  virtual AliMpVPadIterator* CreateIterator() const;
  virtual Int_t GetNeighbours(const AliMpPad& pad, TObjArray& neighbours,
                              Bool_t includeSelf=kFALSE,
                              Bool_t includeVoid=kFALSE) const;
  const char* GetName() const;
  
  Int_t MaxPadIndexX() const;
  Int_t MaxPadIndexY() const;
  Int_t NofPads() const { return fNofStrips; }
    
  virtual AliMpPad PadByLocation(Int_t manuId, Int_t manuChannel, 
                                 Bool_t warning) const;
  
  virtual AliMpPad PadByIndices(Int_t ix, Int_t iy,  
                                Bool_t warning) const;
  
  virtual AliMpPad PadByPosition(Double_t x, Double_t y,
                                 Bool_t warning) const;
  
  const AliMpTrigger* Slat() const;
   
  virtual void GetAllElectronicCardIDs(TArrayI& ecn) const;
  
  virtual AliMp::PlaneType PlaneType() const;
   
  virtual AliMp::StationType StationType() const;
 
  virtual Double_t  GetDimensionX() const;
  virtual Double_t  GetDimensionY() const;
  
  virtual Int_t GetNofElectronicCards() const;
  
  virtual Double_t  GetPositionX() const;
  virtual Double_t  GetPositionY() const;
  
  virtual Bool_t HasMotifPosition(Int_t manuId) const;
  
  virtual AliMpMotifPosition* MotifPosition(Int_t manuId) const;
  
private:
  /// Not implemented
  AliMpTriggerSegmentation(const AliMpTriggerSegmentation& right);
  /// Not implemented
  AliMpTriggerSegmentation&  operator = (const AliMpTriggerSegmentation& right);
  
  const AliMpTrigger* fkSlat;  ///< Slat
  Bool_t              fIsOwner;///< Trigger slat ownership     
  Int_t fNofStrips; ///< Number of strips in this slat

  ClassDef(AliMpTriggerSegmentation,3) // Segmentation for slat trigger stations
};

/// Return station type
inline AliMp::StationType AliMpTriggerSegmentation::StationType() const
{ return AliMp::kStationTrigger; }


#endif
