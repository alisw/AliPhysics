#ifndef ALIMUONHVNAMER_H
#define ALIMUONHVNAMER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup shuttle
/// \class AliMUONHVNamer
/// \brief Collection of methods usefull to HV subprocessor for MUON TRK
/// 
/// \author Laurent Aphecetche, Subtech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class TObjArray;

class AliMUONHVNamer : public TObject
{
public:
  AliMUONHVNamer();
  virtual ~AliMUONHVNamer();
  
  const char* DCSHVChannelName(Int_t detElemId, Int_t sector=0) const;
  
  const char* DCSHVSwitchName(Int_t detElemId, Int_t pcbNumber) const;

  Int_t DCS2DE(Int_t chamberId, Int_t side, Int_t dcsNumber) const;
  
  Int_t DetElemId2DCS(Int_t detElemId, Int_t& side) const;
    
  Int_t DetElemIdFromDCSAlias(const char* dcsAlias) const;
  
  /// Returns the index of PCB (within a St345 slat) for a given manu number.
  Int_t ManuId2PCBIndex(Int_t detElemId, Int_t manuId) const;
  
  /// Return the HV-sector number (within a St12 quadrant) for a given manu number.
  Int_t ManuId2Sector(Int_t detElemId, Int_t manuId) const;
  
  Int_t NumberOfPCBs(Int_t detElemId) const;

  TObjArray* GenerateAliases() const;
  TObjArray* CompactAliases() const;
  
private:
  AliMUONHVNamer(const AliMUONHVNamer& right);
  AliMUONHVNamer&  operator = (const AliMUONHVNamer& right);
    
  static const char* fgHVChannelSt345Pattern[]; // HV Channel name template
  static const char* fgHVChannelSt12Pattern[]; // HV Channel name template
  static const char* fgHVSwitchSt345Pattern; // HV Switch name template
  
  ClassDef(AliMUONHVNamer,0) // Utility class for coding/decoding DCS HV aliases
};

#endif
