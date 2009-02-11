#ifndef ALIMPDCSNAMER_H
#define ALIMPDCSNAMER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup management
/// \class AliMpDCSNamer
/// \brief Collection of methods usefull to DCS handling for MUON TRK and TRG
/// 
//  Author Laurent Aphecetche and Diego Stocco, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif

class TObjArray;

class AliMpDCSNamer : public TObject
{
public:
  AliMpDCSNamer();
  AliMpDCSNamer(const char* detName);
  
  virtual ~AliMpDCSNamer();

  Bool_t SetDetector(const char* detName);
  
  const char* DCSChannelName(Int_t detElemId, Int_t sector=0, Int_t dcsMeasure=0) const;
  
  const char* DCSSwitchName(Int_t detElemId, Int_t pcbNumber) const;

  Int_t DCS2DE(Int_t chamberId, Int_t side, Int_t dcsNumber) const;
  
  Int_t DetElemId2DCS(Int_t detElemId, Int_t& side, Int_t& chId) const;
    
  Int_t DCSIndexFromDCSAlias(const char* dcsAlias) const;
  
  Int_t DetElemIdFromDCSAlias(const char* dcsAlias) const;

  Int_t DCSvariableFromDCSAlias(const char* dcsAlias) const;
  
  Int_t ManuId2Index(Int_t detElemId, Int_t manuId) const;
  
  /// Returns the index of PCB (within a St345 slat) for a given manu number.
  Int_t ManuId2PCBIndex(Int_t detElemId, Int_t manuId) const;
  
  /// Return the HV-sector number (within a St12 quadrant) for a given manu number.
  Int_t ManuId2Sector(Int_t detElemId, Int_t manuId) const;
  
  Int_t NumberOfPCBs(Int_t detElemId) const;

  TObjArray* GenerateAliases() const;
  TObjArray* CompactAliases() const;
  void AliasesAsLdif(const char* ldiffile) const;

  enum 
  {
    kDCSHV,    ///< High Voltage
    kDCSI,     ///< Currents
    kNDCSMeas  ///< Number of measured quantities
  };

  enum 
  {
    kTrackerDet, ///< Namer for tracker
    kTriggerDet  ///< Namer for trigger
  };
  
  
private:
  /// Not implemented
  AliMpDCSNamer(const AliMpDCSNamer& right);
  /// Not implemented
  AliMpDCSNamer&  operator = (const AliMpDCSNamer& right);

  Bool_t CheckConsistency(Int_t detElemId) const;
    
  static const char* fgkDCSChannelSt345Pattern[]; ///< DCS Tracker Channel name template
  static const char* fgkDCSChannelSt12Pattern[]; ///< DCS Tracker Channel name template
  static const char* fgkDCSSwitchSt345Pattern; ///< DCS Tracker Switch name template
  static const char* fgkDCSSideTrackerName[]; ///< DCS Tracker Name of the side written in DCS

  static const char* fgkDCSChannelTriggerPattern[]; ///< DCS Trigger Channel name template
  static const char* fgkDCSSideTriggerName[]; ///< DCS Trigger Name of the side written in DCS
  static const char* fgkDCSMeasureName[]; ///< DCS Trigger Name of the measure (HV or current) written in DCS

  static const char* fgkDetectorName[]; ///< Name of detector (Tracker or Trigger)

  Int_t fDetector; ///< Detector type (either tracker or trigger)
  
  ClassDef(AliMpDCSNamer,0) // Utility class for coding/decoding DCS aliases
};

#endif
