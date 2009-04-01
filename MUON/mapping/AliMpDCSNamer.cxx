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

#include "AliMpDCSNamer.h"

#include "AliCodeTimer.h"
#include "AliLog.h"
#include "AliMpArea.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpHelper.h"
#include "AliMpMotifMap.h"
#include "AliMpMotifPosition.h"
#include "AliMpSector.h"
#include "AliMpSegmentation.h"
#include "AliMpSlat.h"
#include "AliMpConstants.h"
#include <Riostream.h>
#include <TMap.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TString.h>
#include <TSystem.h>

//-----------------------------------------------------------------------------
/// \class AliMpDCSNamer
/// 
/// A utility class to manage DCS aliases names, in particular the
/// two conventions used to number the detection elements within a detector.
///
/// \author: Laurent Aphecetche and Diego Stocco, Subatech
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMpDCSNamer)
/// \endcond

const char* AliMpDCSNamer::fgkDCSChannelSt345Pattern[] = 
{ "MchHvLvLeft/Chamber%02dLeft/Slat%02d.actual.vMon",
  "MchHvLvRight/Chamber%02dRight/Slat%02d.actual.vMon" 
};

const char* AliMpDCSNamer::fgkDCSChannelSt12Pattern[] = 
{
  "MchHvLvLeft/Chamber%02dLeft/Quad%dSect%d.actual.vMon",
  "MchHvLvRight/Chamber%02dRight/Quad%dSect%d.actual.vMon",
};

const char* AliMpDCSNamer::fgkDCSSideTrackerName[] = { "Left", "Right" };


const char* AliMpDCSNamer::fgkDCSSwitchSt345Pattern = "MchDE%04dsw%d.inValue";

const char* AliMpDCSNamer::fgkDCSChannelTriggerPattern[] = {"MTR_%3sSIDE_MT%2i_RPC%i_HV.%s", "MTR_%2sSIDE_MT%2i_RPC%i_HV.%s"};
const char* AliMpDCSNamer::fgkDCSSideTriggerName[] = { "OUT", "IN" };
const char* AliMpDCSNamer::fgkDCSMeasureName[] = { "actual.iMon", "vEff" };

const char* AliMpDCSNamer::fgkDetectorName[] = { "TRACKER", "TRIGGER" };

//_____________________________________________________________________________
AliMpDCSNamer::AliMpDCSNamer():
fDetector(-1)
{
  SetDetector("TRACKER");
  /// default ctor 
}

//_____________________________________________________________________________
AliMpDCSNamer::AliMpDCSNamer(const char* detName):
fDetector(-1)
{
  /// ctor taking the detector name as argument (either trigger or tracker)
  SetDetector(detName);
}

//_____________________________________________________________________________
AliMpDCSNamer::~AliMpDCSNamer()
{
  /// dtor
}

//_____________________________________________________________________________
Bool_t AliMpDCSNamer::SetDetector(const char* detName)
{
  /// Set the detector type
  /// \param detName = tracker, trigger

  TString sDetName(detName);
  Bool_t isOk(kTRUE);
  sDetName.ToUpper();
  if(sDetName.Contains(fgkDetectorName[kTrackerDet]))
    fDetector = kTrackerDet;
  else if(sDetName.Contains(fgkDetectorName[kTriggerDet]))
    fDetector = kTriggerDet;
  else {
    AliWarning("Detector name must be either tracker or trigger. Default tracker selected");
    isOk = kFALSE;
  }
  return isOk;
}


//_____________________________________________________________________________
void 
AliMpDCSNamer::AliasesAsLdif(const char* ldiffile) const
{
/// Export the aliases in LDIF format

  ofstream out(ldiffile);
  
  TObjArray* a = CompactAliases();
  
  TIter next(a);
  TObjString* s;

  // Some header. host name and port probably not up to date.
  TString detName = (fDetector == kTriggerDet) ? "MTR" : "MCH";

  out << "#" << detName.Data() << " config" << endl
      << "dn: det=" << detName.Data() <<",o=alice,dc=cern,dc=ch" << endl
      << "objectClass: AliShuttleDetector" << endl
      << "det: " << detName.Data() << endl
      << "StrictRunOrder: 1" << endl
      << "responsible: aphecetc@in2p3.fr" << endl
      << "DCSHost: aldcs053.cern.ch" << endl
      << "DCSPort: 4242" <<endl;
  
  while ( ( s = (TObjString*)(next()) ) )
  {
    out << "DCSalias: " << s->String().Data() << endl;
  }
  
  out.close();
  
  delete a;
}

//_____________________________________________________________________________
TObjArray*
AliMpDCSNamer::CompactAliases() const
{
  /// Generate a compact list of aliases, for Shuttle test
  /// This one is completely hand-made, in contrast with GenerateAliases()
  /// method

  TObjArray* a = new TObjArray;
  a->SetOwner(kTRUE);

  switch(fDetector){
  case kTrackerDet:
    // St 12 (DCS Channels)
    a->Add(new TObjString("MchHvLvRight/Chamber[00..03]Right/Quad0Sect[0..2].actual.vMon"));
    a->Add(new TObjString("MchHvLvLeft/Chamber[00..03]Left/Quad1Sect[0..2].actual.vMon"));
    a->Add(new TObjString("MchHvLvLeft/Chamber[00..03]Left/Quad2Sect[0..2].actual.vMon"));
    a->Add(new TObjString("MchHvLvRight/Chamber[00..03]Right/Quad3Sect[0..2].actual.vMon"));
  
    // St345 (DCS Channels)
  
    a->Add(new TObjString("MchHvLvRight/Chamber[04..09]Right/Slat[00..08].actual.vMon"));
    a->Add(new TObjString("MchHvLvLeft/Chamber[04..09]Left/Slat[00..08].actual.vMon"));

    a->Add(new TObjString("MchHvLvRight/Chamber[06..09]Right/Slat[09..12].actual.vMon"));
    a->Add(new TObjString("MchHvLvLeft/Chamber[06..09]Left/Slat[09..12].actual.vMon"));
    break;
    
  case kTriggerDet:
    a->Add(new TObjString("MTR_OUTSIDE_MT[11..12]Right/RPC[1..9]_HV.imon"));
    a->Add(new TObjString("MTR_OUTSIDE_MT[21..22]Right/RPC[1..9]_HV.imon"));
    a->Add(new TObjString("MTR_INSIDE_MT[11..12]Right/RPC[1..9]_HV.imon"));
    a->Add(new TObjString("MTR_INSIDE_MT[21..22]Right/RPC[1..9]_HV.imon"));

    a->Add(new TObjString("MTR_OUTSIDE_MT[11..12]Right/RPC[1..9]_HV.vmon"));
    a->Add(new TObjString("MTR_OUTSIDE_MT[21..22]Right/RPC[1..9]_HV.vmon"));
    a->Add(new TObjString("MTR_INSIDE_MT[11..12]Right/RPC[1..9]_HV.vmon"));
    a->Add(new TObjString("MTR_INSIDE_MT[21..22]Right/RPC[1..9]_HV.vmon"));
  }
  

  if(fDetector == kTrackerDet){
    // St345 (DCS Switches)
    AliMpDEIterator it;
  
    it.First();
  
    while (!it.IsDone())
    {
      Int_t detElemId = it.CurrentDEId();
      if ( AliMpDEManager::GetStationType(detElemId) == AliMp::kStation345 )
      {
	a->Add(new TObjString(Form("MchDE%04dsw[0..%d].inValue",detElemId,NumberOfPCBs(detElemId)-1)));
      }
      it.Next();
    }
  }
  
  return a;
}

//_____________________________________________________________________________
Int_t 
AliMpDCSNamer::DCS2DE(Int_t chId, Int_t side, Int_t dcsNumber) const
{
  /// Convert DCS Tracker "slat number" (old convention) to DE (new) convention.
  ///
  /// \param chamberId : chamber number (starting at 0)
  /// \param side : 0 for Left, 1 for Right
  /// \param dcsNumber : slat number in DCS convention
  ///
  /// note that dcsNumber should be >=0 and < number of DEs/2 in chamber

  Int_t de(-1);
  Int_t chamberId = chId;

  if(fDetector == kTrackerDet){ // Tracker

    Int_t nofDE = AliMpDEManager::GetNofDEInChamber(chamberId);
  
    Int_t half = nofDE/2;
  
    dcsNumber = half - dcsNumber;
  
    Int_t quarter = nofDE/4;
    Int_t threeQuarter = half + quarter;
  
    if ( side == 0 ) // left
    {
      de = threeQuarter + 1 - dcsNumber;
    }
    else if ( side == 1 ) // right
    {
      if ( dcsNumber <= quarter )
      {
	de = dcsNumber + threeQuarter;
      }
      else
      {
	de = dcsNumber - quarter - 1;
      }
    }
  }
  else { // Trigger

    if ( chId < 19 ) chamberId = chId - 1;
    else chamberId = chId - 9;

    Int_t nofDE = AliMpDEManager::GetNofDEInChamber(chamberId);

    if ( side == 0 ) // left -> Outside
    {
      de = 14 - dcsNumber;
    }
    else if ( side == 1 ) // right -> Inside
    {
      de = (13 + dcsNumber) % nofDE;
    }
  }
  
  return (chamberId+1)*100 + de;
}


//_____________________________________________________________________________
Int_t
AliMpDCSNamer::DetElemId2DCS(Int_t detElemId, Int_t& side, Int_t &chId) const
{
  /// Convert DE to DCS "slat number"
  /// @see DCS2DE

  CheckConsistency(detElemId);
  
  Int_t chamberId = AliMpDEManager::GetChamberId(detElemId);
  if ( chamberId < 0 ) 
  {
    AliDebug(1,Form("DetElemId %d invalid",detElemId));
    return -1;
  }
  Int_t dcsNumber = (detElemId-(chamberId+1)*100);

  switch ( AliMpDEManager::GetStationType(detElemId) )
  {
    case AliMp::kStation12:
    {
      switch (dcsNumber)
      {
        case 0:
        case 3:
          side = 1; // right
          break;
        case 1:
        case 2:
          side = 0; // left
        default:
          break;
      }
    }
      break;
    case AliMp::kStation345:
    {
      Int_t nofDE = AliMpDEManager::GetNofDEInChamber(chamberId);
      
      Int_t quarter = nofDE/4;
      
      Int_t half = nofDE/2;
      
      Int_t threeQuarter = half + quarter;  
      
      side = -1;
      
      if ( dcsNumber <= quarter ) 
      {
        dcsNumber += quarter + 1 ;
        side = 1; // right
      }
      else if ( dcsNumber <= threeQuarter )
      {
        dcsNumber = ( threeQuarter - dcsNumber + 1 );
        side = 0; // left
      }
      else if ( dcsNumber > threeQuarter ) 
      {
        dcsNumber = dcsNumber - threeQuarter;
        side = 1; // right
      }
      else
      {
        AliFatal("oups");
      }  
      // dcs convention change : numbering from top, not from bottom
      dcsNumber = half-dcsNumber;
    }
      break;
    case AliMp::kStationTrigger:
    {
      if (chamberId < AliMpConstants::NofChambers()-2)
	chId = chamberId + 1;
      else chId = 23 + chamberId - AliMpConstants::NofChambers();

      Int_t nofDE = AliMpDEManager::GetNofDEInChamber(chamberId);

      if ( dcsNumber >=5 && dcsNumber <= 13 ) {
	side = 0;
	dcsNumber = 14 - dcsNumber;
      }
      else {
	side = 1;
	dcsNumber = (5 + dcsNumber) % nofDE;
      }
      AliDebug(10, Form("detElemId %i  -> MT%i_side%i_L%i", detElemId, chId, side, dcsNumber));
    }
      break;
    default:
      break;
  }

  return dcsNumber;
}


//_____________________________________________________________________________
const char* 
  AliMpDCSNamer::DCSChannelName(Int_t detElemId, Int_t sector, Int_t dcsMeasure) const
{
  /// Return the alias name of the DCS Channel for a given DCS area 
  /// \param detElemId 
  /// \param sector = 0,1 or 2 for St12, and is unused for st345 and trigger
  /// \param dcsMeasure = kDCSHV, kDCSI and is unused for tracker
  
  Int_t chamberId = AliMpDEManager::GetChamberId(detElemId);
  if ( chamberId < 0 ) return 0x0;

  Int_t side(-1), chId(-1);
  Int_t dcsNumber = DetElemId2DCS(detElemId,side,chId);

  switch (AliMpDEManager::GetStationType(detElemId))
  {
    case AliMp::kStation12:
      return Form(fgkDCSChannelSt12Pattern[side],chamberId,dcsNumber,sector);
      break;
    case AliMp::kStation345:
      return Form(fgkDCSChannelSt345Pattern[side],chamberId,dcsNumber);
      break;
    case AliMp::kStationTrigger:
      return Form(fgkDCSChannelTriggerPattern[side],fgkDCSSideTriggerName[side],chId,dcsNumber,fgkDCSMeasureName[dcsMeasure]);
      break;
    default:
      return 0x0;
      break;
  }
}

//_____________________________________________________________________________
const char* 
AliMpDCSNamer::DCSSwitchName(Int_t detElemId, Int_t pcbNumber) const
{
  /// Return the alias name of the DCS Switch for a given PCB 
  /// within a slat of St345
  
  if (AliMpDEManager::GetStationType(detElemId) == AliMp::kStation345)
  {
    return Form(fgkDCSSwitchSt345Pattern,detElemId,pcbNumber);
  }
  return 0x0;
}

//_____________________________________________________________________________
Int_t 
AliMpDCSNamer::DCSIndexFromDCSAlias(const char* dcsAlias) const
{
  /// Converts the dcs alias to a hv index 
  ///
  /// dcsAlias has one of the following 3 forms :
  ///
  /// MchHvLv[Left|Right]/Chamber##[Left|Right]/Chamber##[Left|Right]Slat##.actual.vMon
  ///
  /// MchHvLv[Left|Right]/Chamber##[Left|Right]/Chamber##[Left|Right]Quad#Sect#.actual.vMon
  ///
  /// MchDE####dsw#.inValue
  
  TString sDcsAlias(dcsAlias);
  Int_t de(-1);
  Int_t sw(-1);
  
  int side(-1);
  
  if ( sDcsAlias.Contains("Left") )
  {
    side = 0;
  }
  else if ( sDcsAlias.Contains("Right") )
  {
    side = 1;
  }
  else
  {
    /// it's a switch
    sscanf(sDcsAlias.Data(),fgkDCSSwitchSt345Pattern,&de,&sw);
    return sw;
  }
  
  int n1(-1);
  int n3(-1);
  int n4(-1);
  
  if ( sDcsAlias.Contains("Quad") )
  {
    sscanf(sDcsAlias.Data(),fgkDCSChannelSt12Pattern[side],&n1,&n3,&n4);    
    return n4;
  }
  
  return -2;
}

//_____________________________________________________________________________
Int_t 
AliMpDCSNamer::DetElemIdFromDCSAlias(const char* dcsAlias) const
{
  /// Converts the dcs alias to a detection element identifier
  ///
  /// dcsAlias has one of the following forms :
  ///
  /// MchHvLv[Left|Right]/Chamber##[Left|Right]/Chamber##[Left|Right]Slat##.actual.vMon
  ///
  /// MchHvLv[Left|Right]/Chamber##[Left|Right]/Chamber##[Left|Right]Quad#Sect#.actual.vMon
  ///
  /// MTR_Side[OUTSIDE|INSIDE]_MTChamber##_RPC#_HV.Type[actual.iMon|vEff]
  
  AliDebug(1,Form("dcsAlias=%s",dcsAlias));
  
  TString sDcsAlias(dcsAlias);
  
  int side(-1);

  const char** sideName = (fDetector == kTriggerDet) ? fgkDCSSideTriggerName : fgkDCSSideTrackerName;

  for(Int_t iside=0; iside<2; iside++){
    if ( sDcsAlias.Contains(sideName[iside]) ) {
      side = iside;
      break;
    }
  }
  if(side<0) return -2;
  
  int n1(-1);
  int n3(-1);
  int n4(-1);
  char type[15];
  char cside[4];
  int detElemId(-1);
  
  if ( sDcsAlias.Contains("Slat") )
  {
    sscanf(sDcsAlias.Data(),fgkDCSChannelSt345Pattern[side],&n1,&n3);
    detElemId = DCS2DE(n1,side,n3);
    AliDebug(1,Form("Slat side=%d n1=%d n3=%d de=%d",side,n1,n3,detElemId));
  }
  else if ( sDcsAlias.Contains("Quad") )
  {
    sscanf(sDcsAlias.Data(),fgkDCSChannelSt12Pattern[side],&n1,&n3,&n4);    
    detElemId = 100*(n1+1) + n3;
    AliDebug(1,Form("Quad side=%d n1=%d n3=%d n4=%d de=%d",side,n1,n3,n4,detElemId));
  }
  else if ( sDcsAlias.Contains("MT") )
  {
    sscanf(sDcsAlias.Data(),fgkDCSChannelTriggerPattern[side],cside,&n1,&n3,type);
    detElemId = DCS2DE(n1,side,n3);
    AliDebug(1,Form("Slat side=%d n1=%d n3=%d de=%d",side,n1,n3,detElemId));
  }
  else
  {
    return -3;
  }
  
  if ( !AliMpDEManager::IsValidDetElemId(detElemId)  )
  {
    AliError(Form("Invalid aliasName %s",dcsAlias));
    return -1;
  }
  
  return detElemId;
}

//_____________________________________________________________________________
Int_t AliMpDCSNamer::DCSvariableFromDCSAlias(const char* dcsAlias) const
{
  /// Get DCS variable from an alias (trigger)
  
  TString sDcsAlias(dcsAlias);

  Int_t dcsMeasurement = -1;

  for(Int_t iMeas=0; iMeas<kNDCSMeas; iMeas++){
    if ( sDcsAlias.Contains(fgkDCSMeasureName[iMeas]) ) {
      dcsMeasurement = iMeas;
      break;
    }
  }

  return dcsMeasurement;
}


//_____________________________________________________________________________
TObjArray*
AliMpDCSNamer::GenerateAliases() const
{
  /// Generate DCS alias names, for MUON Tracker High Voltage system.
  /// or for MUON Trigger HV and current system.
  ///
  /// We first generate aliases of DCS channels :
  ///
  /// St 1 ch  1 : 12 channels
  ///      ch  2 : 12 channels 
  /// St 2 ch  3 : 12 channels
  ///      ch  4 : 12 channels
  /// St 3 ch  5 : 18 channels
  ///      ch  6 : 18 channels
  /// St 4 ch  7 : 26 channels
  ///      ch  8 : 26 channels
  /// St 5 ch  9 : 26 channels
  ///      ch 10 : 26 channels
  ///
  /// then aliases of DCS switches (only for St345) : 1 switch per PCB.
  ///
  /// Returns a TObjArray of TObjString(=alias name)
  
  TObjArray* aliases = new TObjArray;
  aliases->SetOwner(kTRUE);

  Int_t nMeasures = (fDetector == kTriggerDet) ? kNDCSMeas : 1;
  
  for(Int_t iMeas=0; iMeas<nMeasures; iMeas++){

    AliMpDEIterator it;
  
    it.First();
  
    while (!it.IsDone())
    {
      Int_t detElemId = it.CurrentDEId();
      switch (fDetector){
      case kTrackerDet:
      {
	switch ( AliMpDEManager::GetStationType(detElemId) )
	{
	case AliMp::kStation12:
	  for ( int sector = 0; sector < 3; ++sector)
	  {
	    aliases->Add(new TObjString(DCSChannelName(detElemId,sector)));
	  }
	  break;
	case AliMp::kStation345:
	  aliases->Add(new TObjString(DCSChannelName(detElemId)));
	  for ( Int_t i = 0; i < NumberOfPCBs(detElemId); ++i )
	  {
	    aliases->Add(new TObjString(DCSSwitchName(detElemId,i)));
	  }
	  break;
	default:
	  break;
	}
      }
      break;
      case kTriggerDet:
      {
	switch ( AliMpDEManager::GetStationType(detElemId) )
	{
	case AliMp::kStationTrigger:
	  AliDebug(10,Form("Current DetElemId %i",detElemId));
	  aliases->Add(new TObjString(DCSChannelName(detElemId,0,iMeas)));
	  break;
	default:
	  break;
	}
      }
      break;
      }
      it.Next();
    } // loop on detElemId
  } // Loop on measurement type

  return aliases;
}

//_____________________________________________________________________________
Int_t 
AliMpDCSNamer::ManuId2Index(Int_t detElemId, Int_t manuId) const
{
  /// Convert (de,manu) to hv index, depending on the station
  
  AliMp::StationType stationType = AliMpDEManager::GetStationType(detElemId);
  if ( stationType == AliMp::kStation345 ) 
  {
    return ManuId2PCBIndex(detElemId,manuId);
  }
  else if ( stationType == AliMp::kStation12 ) 
  {
    return ManuId2Sector(detElemId,manuId);
  }
  return -1;
}

//_____________________________________________________________________________
Int_t 
AliMpDCSNamer::ManuId2PCBIndex(Int_t detElemId, Int_t manuId) const
{
  /// Returns the index of PCB (within a St345 slat) for a given manu number.
  /// Returns -1 if (detElemId,manuId) is incorrect
  
  AliCodeTimerAuto("")
  
  const AliMpSlat* slat 
    = AliMpSegmentation::Instance()->GetSlatByElectronics(detElemId, manuId);
  if ( ! slat ) return -1;
  
  return slat->FindPCBIndexByMotifPositionID(manuId);
}

//_____________________________________________________________________________
Int_t 
AliMpDCSNamer::ManuId2Sector(Int_t detElemId, Int_t manuId) const
{
  /// Return the DCS-sector number (within a St12 quadrant) for a given manu number.
  
  AliCodeTimerAuto("")
  
  const AliMpSector* sector 
    = AliMpSegmentation::Instance()->GetSectorByElectronics(detElemId, manuId);
  if ( ! sector ) return -1;
  
  const AliMpMotifMap* motifMap = sector->GetMotifMap();
  const AliMpMotifPosition* motifPos = motifMap->FindMotifPosition(manuId);

  Double_t lowerLeftX 
    = motifPos->GetPositionX()-motifPos->GetDimensionX();
  
  Double_t x = lowerLeftX*10.0; // cm -> mm
  Int_t isector(-1);

  AliMq::Station12Type stationType = AliMpDEManager::GetStation12Type(detElemId);
  
  if ( stationType == AliMq::kStation1 ) 
  {
    if ( x < -10 ) AliFatal("");
    
    if ( x < 291.65 ) isector = 0;
    else if ( x < 585.65 ) isector = 1;
    else if ( x < 879.65 ) isector = 2;
  }
  else
  {
    if ( x < -140 ) AliFatal("");
    
    if ( x < 283.75 ) isector = 0;
    else if ( x < 603.75 ) isector = 1;
    else if ( x < 1158.75 ) isector = 2;
  }
  
  return isector;
}

//_____________________________________________________________________________
Int_t 
AliMpDCSNamer::NumberOfPCBs(Int_t detElemId) const
{
  /// Returns the number of PCB in a given detection element
  /// Only works for St345
  
  AliMp::StationType stationType = AliMpDEManager::GetStationType(detElemId);
  if ( stationType != AliMp::kStation345 )
  {
    return 0;
  }
  else
  {
    const AliMpSlat* slat 
      = AliMpSegmentation::Instance()->GetSlat(detElemId, AliMp::kCath0);
    return slat->GetSize();
  }
}

//_____________________________________________________________________________
Bool_t AliMpDCSNamer::CheckConsistency(Int_t detElemId) const
{
  //
  /// Check that the required detElemId either belongs to tracker or trigger
  /// consistently with the initial definition of the namer
  //

  Bool_t isConsistent(kFALSE);
  TString requestInfo;
  switch(AliMpDEManager::GetStationType(detElemId))
  {
  case AliMp::kStation12:
  case AliMp::kStation345:
    if (fDetector == kTrackerDet) isConsistent = kTRUE;
    requestInfo = "TRACKER";
    break;
  case AliMp::kStationTrigger:
    if (fDetector == kTriggerDet) isConsistent = kTRUE;
    requestInfo = "TRIGGER";
    break;
  default:
    break;
  }

  if(!isConsistent) AliWarning(Form("Requesting information for %s station but class initialized for %s",requestInfo.Data(), fgkDetectorName[fDetector]));

  return isConsistent;
}
