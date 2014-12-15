#ifndef ALIMPPADUID_H
#define ALIMPPADUID_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup management
/// \class AliMpPadUID
/// \brief Unique ID for pads
/// 
//  Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMpPadUID : public TObject
{
public:
  AliMpPadUID(UInt_t uid=0);
  AliMpPadUID(Int_t detElemId, Int_t manuId, Int_t manuChannel);
  virtual ~AliMpPadUID();
  
  /// Get detection element
  Int_t DetElemId() const { return AliMpPadUID::DetElemId(GetUniqueID()); }
  
  /// Get manuId
  Int_t ManuId() const { return AliMpPadUID::ManuId(GetUniqueID()); }
  
  /// Get manu channel
  Int_t ManuChannel() const { return AliMpPadUID::ManuChannel(GetUniqueID()); }
  
  static UInt_t BuildUniqueID(Int_t detElemId, Int_t manuId, 
                              Int_t manuChannel);
  
  static Int_t DetElemId(UInt_t uniqueID);
  
  static Int_t ManuChannel(UInt_t uniqueID);

  static Int_t ManuId(UInt_t uniqueID);
  
  ClassDef(AliMpPadUID,1) // Unique ID for MUON tracker pad
};

#endif
