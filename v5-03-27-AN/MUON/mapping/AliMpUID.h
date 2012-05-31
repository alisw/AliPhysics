#ifndef ALIMPUID_H
#define ALIMPUID_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup management
/// \class AliMpUID
/// \brief Global (string-eable) ID of a tracker channel
/// 
//  Author Laurent Aphecetche, Subatech

#ifndef ROOT_TObject
#  include "TObject.h"
#endif
#ifndef ROOT_TString
#  include "TString.h"
#endif
#ifndef ALI_MP_CATHOD_TYPE_H
#  include "AliMpCathodType.h"
#endif

class AliMpUID : public TObject
{
public:
  AliMpUID();
  AliMpUID(AliMp::CathodType cathodeType, Int_t station, Int_t chamber=-1, Int_t de=-1, 
            Int_t bp=-1, Int_t manu=-1, Int_t pcb=-1);
  AliMpUID(AliMp::CathodType cathodeType, const AliMpUID& b);
  AliMpUID(AliMp::CathodType cathodeType, const char* pathname);
  AliMpUID(const char* pathname);
  
  /// dtor
  virtual ~AliMpUID() {}
  
  TString Name() const;
  TString PathName() const;
  TString BaseName() const;
  TString DirName() const;
  
  Bool_t IsStation() const;
  Bool_t IsChamber() const;
  Bool_t IsDetectionElement() const;
  Bool_t IsBusPatch() const;
  Bool_t IsManu() const;
  Bool_t IsPCB() const;
  Bool_t IsValid() const;
  
  AliMp::CathodType CathodeId() const;
  /// Return station Id
  Int_t StationId() const { return fStationId; }
  /// Return chamber Id
  Int_t ChamberId() const { return fChamberId; }
  /// Return detection element Id
  Int_t DetElemId() const { return fDetElemId; }
  /// Return bus patch Id
  Int_t BusPatchId() const { return fBusPatchId; }
  /// Return manu Id
  Int_t ManuId() const { return fManuId; }
  /// Return PCB Id
  Int_t PCBId() const { return fPCBId; }
  
  virtual void Print(Option_t* opt="") const;

  /// Return our type (e.g. PCB, Chamber, DE, MANU, etc...)
  TString Type() const;

private:
  
  Bool_t CheckTemplate(const char* name, const char* templateName, Int_t& value);
  TString StripCathode(const char* name) const;
  
private:
  Int_t fCathodeId; ///< Cathode number
  Int_t fStationId; ///< Station id
  Int_t fChamberId; ///< Chamber id
  Int_t fDetElemId; ///< Detection element id
  Int_t fBusPatchId;///< Bus patch id
  Int_t fManuId;    ///< Manu id
  Int_t fPCBId;     ///< PCB id
  
  ClassDef(AliMpUID,1) // UID of a tracker channel 
};

#endif
