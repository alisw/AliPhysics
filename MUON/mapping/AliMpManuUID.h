#ifndef ALIMPMANUUID_H
#define ALIMPMANUUID_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup management
/// \class AliMpManuUID
/// \brief Unique ID for manus
/// 
//  Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class AliMpManuUID : public TObject
{
public:
  AliMpManuUID();
  AliMpManuUID(Int_t detElemId, Int_t manuId);
  virtual ~AliMpManuUID();
  
  /// Get detection element
  Int_t DetElemId() const { return AliMpManuUID::DetElemId(GetUniqueID()); }

  /// Get manu identifier
  Int_t ManuId() const { return AliMpManuUID::ManuId(GetUniqueID()); }
  
  static UInt_t BuildUniqueID(Int_t detElemId, Int_t manuId);
  
  static Int_t DetElemId(UInt_t uniqueID);
  
  static Int_t ManuId(UInt_t uniqueID);
  
  ClassDef(AliMpManuUID,2) // Unique ID for MUON tracker manus
};

#endif
