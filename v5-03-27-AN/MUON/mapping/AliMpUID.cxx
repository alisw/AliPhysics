/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

// $Id$

#include "AliMpUID.h"

#include "AliLog.h"
#include "Riostream.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TSystem.h"

///
/// station/chamber/de/bp/manu
///
/// station/chamber/pcb/manu

ClassImp(AliMpUID)

namespace
{
  const char* nameTemplateMANU = "MANU %d";
  const char* nameTemplateDE = "DE %d";
  const char* nameTemplateBP = "BusPatch %d";
  const char* nameTemplateCHAMBER = "Chamber %d";
  const char* nameTemplateSTATION = "Station %d";
  const char* nameTemplatePCB = "PCB %d";
  
  const char* pathTemplateMANU = "Cathode%d/Station%d/Chamber%d/DE%04d/BUSPATCH%04d/MANU%04d";
  const char* pathTemplateBP = "Cathode%d/Station%d/Chamber%d/DE%04d/BUSPATCH%04d";
  const char* pathTemplateDE = "Cathode%d/Station%d/Chamber%d/DE%04d";
  const char* pathTemplateCHAMBER = "Cathode%d/Station%d/Chamber%d";
  const char* pathTemplateSTATION = "Cathode%d/Station%d";

  const char* pathTemplateMANUPCB = "Cathode%d/Station%d/Chamber%d/DE%04d/PCB%d/MANU%04d";
  const char* pathTemplatePCB = "Cathode%d/Station%d/Chamber%d/DE%04d/PCB%d";
}

//_____________________________________________________________________________
AliMpUID::AliMpUID()
: 
fCathodeId(-1),
fStationId(-1),
fChamberId(-1),
fDetElemId(-1),
fBusPatchId(-1),
fManuId(-1),
fPCBId(-1)
{
  /// empty ctor
}

//_____________________________________________________________________________
AliMpUID::AliMpUID(AliMp::CathodType cathodeType, Int_t station, Int_t chamber, Int_t de, Int_t bp, Int_t manu, Int_t pcb)
: 
fCathodeId(cathodeType),
fStationId(station),
fChamberId(chamber),
fDetElemId(de),
fBusPatchId(bp),
fManuId(manu),
fPCBId(pcb)
{
  /// default ctor
}

//_____________________________________________________________________________
AliMpUID::AliMpUID(AliMp::CathodType cathodeType, const AliMpUID& b)
: 
fCathodeId(cathodeType),
fStationId(b.StationId()),
fChamberId(b.ChamberId()),
fDetElemId(b.DetElemId()),
fBusPatchId(b.BusPatchId()),
fManuId(b.ManuId()),
fPCBId(b.PCBId())
{
  /// build the id from b, but using the given cathodeType
}

//_____________________________________________________________________________
AliMpUID::AliMpUID(AliMp::CathodType cathodeType, const char* pathname)
:
fCathodeId(cathodeType),
fStationId(-1),
fChamberId(-1),
fDetElemId(-1),
fBusPatchId(-1),
fManuId(-1),
fPCBId(-1)
{
  /// build id from path, but using the given cathodeType
  
  if ( CheckTemplate(pathname,pathTemplateMANUPCB,fPCBId) && fPCBId >= 0 ) return;
  if ( CheckTemplate(pathname,pathTemplatePCB,fPCBId) && fPCBId >= 0 ) return;
  
  if ( CheckTemplate(pathname,pathTemplateMANU,fBusPatchId) ) return;
  if ( CheckTemplate(pathname,pathTemplateBP,fBusPatchId) ) return;
  if ( CheckTemplate(pathname,pathTemplateDE,fBusPatchId) ) return;
  if ( CheckTemplate(pathname,pathTemplateCHAMBER,fBusPatchId) ) return;
  if ( CheckTemplate(pathname,pathTemplateSTATION,fBusPatchId) ) return;
}


//_____________________________________________________________________________
AliMpUID::AliMpUID(const char* pathname)
:
fCathodeId(2),
fStationId(-1),
fChamberId(-1),
fDetElemId(-1),
fBusPatchId(-1),
fManuId(-1),
fPCBId(-1)
{
  /// Build id from path
  
  if ( CheckTemplate(pathname,pathTemplateMANUPCB,fPCBId) && fPCBId >= 0 ) return;
  if ( CheckTemplate(pathname,pathTemplatePCB,fPCBId) && fPCBId >= 0 ) return;
  
  if ( CheckTemplate(pathname,pathTemplateMANU,fBusPatchId) ) return;
  if ( CheckTemplate(pathname,pathTemplateBP,fBusPatchId) ) return;
  if ( CheckTemplate(pathname,pathTemplateDE,fBusPatchId) ) return;
  if ( CheckTemplate(pathname,pathTemplateCHAMBER,fBusPatchId) ) return;
  if ( CheckTemplate(pathname,pathTemplateSTATION,fBusPatchId) ) return;
}

//_____________________________________________________________________________
TString
AliMpUID::BaseName() const
{
  /// Get the basename
  return gSystem->BaseName(PathName().Data());
}

//_____________________________________________________________________________
AliMp::CathodType 
AliMpUID::CathodeId() const
{
  /// return cathode id (not always valid)
  return AliMp::GetCathodType(fCathodeId);
}

//_____________________________________________________________________________
Bool_t
AliMpUID::CheckTemplate(const char* name, const char* pathTemplateName, Int_t& value)
{
  /// Check a name against a template
  
  if ( TString(name).Contains("Cathode") ) 
  {
    sscanf(name,pathTemplateName,&fCathodeId,&fStationId,&fChamberId,&fDetElemId,&value,&fManuId);
  }
  else
  {
    TString templ(pathTemplateName);
    Int_t i = templ.Index("/");
    templ = templ(i+1,templ.Length()-i-1);
    sscanf(name,templ.Data(),&fStationId,&fChamberId,&fDetElemId,&value,&fManuId);
  }
  return IsValid();
}

//_____________________________________________________________________________
TString
AliMpUID::DirName() const
{
  /// Get dirname
  return gSystem->DirName(PathName().Data());
}

//_____________________________________________________________________________
Bool_t 
AliMpUID::IsStation() const
{
  /// Whether we identify a station
  return fCathodeId >= 0 && fStationId >= 0 && fChamberId == -1 ;
}

//_____________________________________________________________________________
Bool_t 
AliMpUID::IsChamber() const
{
  /// Whether we identify a chamber

  return fCathodeId >= 0 && fStationId >= 0 && fChamberId >= 0 && fDetElemId == -1;
}

//_____________________________________________________________________________
Bool_t 
AliMpUID::IsDetectionElement() const
{
  /// whether we identify a detection element
  return fCathodeId >= 0 &&  fStationId >= 0 && fChamberId >= 0 && fDetElemId >= 0 && fBusPatchId==-1 && fPCBId == -1;
}

//_____________________________________________________________________________
Bool_t 
AliMpUID::IsBusPatch() const
{
  /// whether we identify a bus patch
  return fCathodeId >= 0 && fStationId >= 0 && fChamberId >= 0 && fDetElemId >= 0 && fBusPatchId>=0 && fManuId ==-1;
}

//_____________________________________________________________________________
Bool_t AliMpUID::IsManu() const
{
  /// whether we identify a manu
  return 
  fCathodeId >= 0 && 
  fStationId >= 0 && 
  fChamberId >= 0 && 
  fDetElemId >= 0 && 
  ( fBusPatchId>=0 || fPCBId >=0 ) && 
  fManuId >=0;
}

//_____________________________________________________________________________
Bool_t AliMpUID::IsPCB() const
{
  /// Whether we identify a PCB
  return fCathodeId >= 0 && fPCBId >= 0 && fManuId == -1;
}

//_____________________________________________________________________________
Bool_t AliMpUID::IsValid() const
{
  /// Whether we're a valid UID...
  return IsStation() || IsChamber() || IsDetectionElement() || IsBusPatch() || IsManu() || IsPCB();
}

//_____________________________________________________________________________
TString 
AliMpUID::Name() const
{
  /// Get our name
  if ( IsManu() ) 
  {
    return Form(nameTemplateMANU,ManuId());
  }
  
  if ( IsPCB() ) 
  {
    return Form(nameTemplatePCB,PCBId());
  }
  
  if ( IsBusPatch() ) 
  {
    return Form(nameTemplateBP,BusPatchId());
  }
  
  if ( IsDetectionElement() ) 
  {
    return Form(nameTemplateDE,DetElemId());
  }
  
  if ( IsChamber() ) 
  {
    return Form(nameTemplateCHAMBER,ChamberId());
  }
  
  if ( IsStation() ) 
  {
    return Form(nameTemplateSTATION,StationId());
  }
  
  return "INVALID NAME";
}

//_____________________________________________________________________________
TString 
AliMpUID::PathName() const
{
  /// Get our pathname
  if ( IsManu() ) 
  {
    if ( fPCBId >= 0 ) 
    {
      return StripCathode(Form(pathTemplateMANUPCB,CathodeId(),StationId(),ChamberId(),DetElemId(),PCBId(),ManuId()));
    }
    else
    {
      return StripCathode(Form(pathTemplateMANU,CathodeId(),StationId(),ChamberId(),DetElemId(),BusPatchId(),ManuId()));
    }
  }
  
  if ( IsPCB() ) 
  {
    return StripCathode(Form(pathTemplatePCB,CathodeId(),StationId(),ChamberId(),DetElemId(),PCBId()));
  }
  
  if ( IsBusPatch() ) 
  {
    return StripCathode(Form(pathTemplateBP,CathodeId(),StationId(),ChamberId(),DetElemId(),BusPatchId()));
  }
  
  if ( IsDetectionElement() ) 
  {
    return StripCathode(Form(pathTemplateDE,CathodeId(),StationId(),ChamberId(),DetElemId()));
  }
  
  if ( IsChamber() ) 
  {
    return StripCathode(Form(pathTemplateCHAMBER,CathodeId(),StationId(),ChamberId()));
  }
  
  if ( IsStation() ) 
  {
    return StripCathode(Form(pathTemplateSTATION,CathodeId(),StationId()));
  }
  
  return "INVALID PATHNAME";
}

//_____________________________________________________________________________
void 
AliMpUID::Print(Option_t*) const
{
  /// Printout
  cout << Name().Data() << " (" << PathName().Data() << ")" << endl;
}

//_____________________________________________________________________________
TString
AliMpUID::StripCathode(const char* name) const
{
  /// Remove cathode information if both cathodes are present
  
  TString rv(name);
  
  if ( fCathodeId == 2 ) 
  {
    rv.ReplaceAll("Cathode2/","");
  }
  
  return rv;
}

//_____________________________________________________________________________
TString
AliMpUID::Type() const
{
  /// Remove cathode information if both cathodes are present
  TString n(Name());
  TObjArray* s = n.Tokenize(" ");
  TString rv(static_cast<TObjString*>(s->At(0))->String());
  delete s;
  return rv;
}

