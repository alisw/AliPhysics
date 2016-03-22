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
#include <cassert>

//-----------------------------------------------------------------------------
/// \class AliMpDCSNamer
///
/// A utility class to manage DCS aliases names, in particular the
/// two conventions used to number the detection elements within a detector.
///
/// \author: Laurent Aphecetche and Diego Stocco, Subatech
//-----------------------------------------------------------------------------

using std::cout;
using std::endl;
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
  "MchHvLvRight/Chamber%02dRight/Quad%dSect%d.actual.vMon"
};

const char* AliMpDCSNamer::fgkDCSQuadrantPattern[] =
{
  "MchHvLvLeft/Chamber%02dLeft/Quad%d",
  "MchHvLvRight/Chamber%02dRight/Quad%d"
};

const char* AliMpDCSNamer::fgkDCSChamberPattern[] =

{
  "MchHvLvLeft/Chamber%02dLeft",
  "MchHvLvRight/Chamber%02dRight"
};


const char* AliMpDCSNamer::fgkDCSMCHLVGroupPattern[] =
{
 "MchHvLvLeft/Chamber%02dLeft/Group%d%s",
 "MchHvLvRight/Chamber%02dRight/Group%d%s",
};

const char* AliMpDCSNamer::fgkDCSSideTrackerName[] = { "Left", "Right" };


const char* AliMpDCSNamer::fgkDCSSwitchSt345Pattern = "MchDE%04dsw%d.inValue";

const char* AliMpDCSNamer::fgkDCSChannelTriggerPatternRead[] = {"MTR_%3sSIDE_MT%2i_RPC%i_HV.%14s", "MTR_%2sSIDE_MT%2i_RPC%i_HV.%14s"};
const char* AliMpDCSNamer::fgkDCSChannelTriggerPattern[] = {"MTR_%3sSIDE_MT%2i_RPC%i_HV.%s", "MTR_%2sSIDE_MT%2i_RPC%i_HV.%s"};
const char* AliMpDCSNamer::fgkDCSSideTriggerName[] = { "OUT", "IN" };
const char* AliMpDCSNamer::fgkDCSMeasureName[] = { "vEff", "actual.iMon" };

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
    // St 12 (DCS HV Channels)
    a->Add(new TObjString("MchHvLvRight/Chamber[00..03]Right/Quad0Sect[0..2].actual.vMon"));
    a->Add(new TObjString("MchHvLvLeft/Chamber[00..03]Left/Quad1Sect[0..2].actual.vMon"));
    a->Add(new TObjString("MchHvLvLeft/Chamber[00..03]Left/Quad2Sect[0..2].actual.vMon"));
    a->Add(new TObjString("MchHvLvRight/Chamber[00..03]Right/Quad3Sect[0..2].actual.vMon"));

    // St345 (DCS HV Channels)

    a->Add(new TObjString("MchHvLvRight/Chamber[04..09]Right/Slat[00..08].actual.vMon"));
    a->Add(new TObjString("MchHvLvLeft/Chamber[04..09]Left/Slat[00..08].actual.vMon"));

    a->Add(new TObjString("MchHvLvRight/Chamber[06..09]Right/Slat[09..12].actual.vMon"));
    a->Add(new TObjString("MchHvLvLeft/Chamber[06..09]Left/Slat[09..12].actual.vMon"));

    // (LV groups)

    // St12 have 4 LV groups

    a->Add(new TObjString("MchHvLvLeft/Chamber[01..10]Left/Group[1..4]dig"));
    a->Add(new TObjString("MchHvLvLeft/Chamber[01..10]Left/Group[1..4]ann"));
    a->Add(new TObjString("MchHvLvLeft/Chamber[01..10]Left/Group[1..4]anp"));
    a->Add(new TObjString("MchHvLvRight/Chamber[01..10]Right/Group[1..4]dig"));
    a->Add(new TObjString("MchHvLvRight/Chamber[01..10]Right/Group[1..4]ann"));
    a->Add(new TObjString("MchHvLvRight/Chamber[01..10]Right/Group[1..4]anp"));

    // St3 has 5 LV groups

    a->Add(new TObjString("MchHvLvLeft/Chamber[05..06]Left/Group[5..5]dig"));
    a->Add(new TObjString("MchHvLvLeft/Chamber[05..06]Left/Group[5..5]ann"));
    a->Add(new TObjString("MchHvLvLeft/Chamber[05..06]Left/Group[5..5]anp"));
    a->Add(new TObjString("MchHvLvRight/Chamber[05..06]Right/Group[5..5]dig"));
    a->Add(new TObjString("MchHvLvRight/Chamber[05..06]Right/Group[5..5]ann"));
    a->Add(new TObjString("MchHvLvRight/Chamber[05..06]Right/Group[5..5]anp"));

    // St4-5 have 7 LV groups
    a->Add(new TObjString("MchHvLvLeft/Chamber[07..10]Left/Group[5..7]dig"));
    a->Add(new TObjString("MchHvLvLeft/Chamber[07..10]Left/Group[5..7]ann"));
    a->Add(new TObjString("MchHvLvLeft/Chamber[07..10]Left/Group[5..7]anp"));
    a->Add(new TObjString("MchHvLvRight/Chamber[07..10]Right/Group[5..7]dig"));
    a->Add(new TObjString("MchHvLvRight/Chamber[07..10]Right/Group[5..7]ann"));
    a->Add(new TObjString("MchHvLvRight/Chamber[07..10]Right/Group[5..7]anp"));

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
      if (nofDE>0)
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
      chId = chamberId;

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
      chId = chamberId;

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
	if (nofDE>0)
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
TString
AliMpDCSNamer::DCSNameFromAlias(const char* dcsAlias) const
{
  /// Convert a (possibly partial) aliasname to a name (only for MCH HV)

  TString salias(dcsAlias);

  if ( !salias.Contains("MchHvLv") ) return dcsAlias; // not MCH
  if ( salias.Contains("Group")) return dcsAlias; // not MCH HV (but LV)

  Int_t quadrantNumber(-1);
  Int_t chamberNumber(-1);
  Int_t side(-1);

  if ( salias.Contains("Left")) side = 0;
  if ( salias.Contains("Right")) side = 1;

  if ( side < 0 ) return "";

  TString channelName;

  if ( salias.Contains("Slat") )
  {
    Int_t slatNumber(-1);
    sscanf(salias.Data(),fgkDCSChannelSt345Pattern[side],&chamberNumber,&slatNumber);
    ++chamberNumber;
    ++slatNumber;
    channelName = TString::Format(fgkDCSChannelSt345Pattern[side],chamberNumber,slatNumber);
  }
  else if ( salias.Contains("Sect") )
  {
    Int_t sectorNumber(-1);
    sscanf(salias.Data(),fgkDCSChannelSt12Pattern[side],&chamberNumber,&quadrantNumber,&sectorNumber);
    ++chamberNumber;
    ++quadrantNumber;
    ++sectorNumber;
    channelName =  TString::Format(fgkDCSChannelSt12Pattern[side],chamberNumber,quadrantNumber,sectorNumber);
  }
  else if ( salias.Contains("Quad") )
  {
    sscanf(salias.Data(),fgkDCSQuadrantPattern[side],&chamberNumber,&quadrantNumber);
    ++chamberNumber;
    ++quadrantNumber;
    channelName =  TString::Format(fgkDCSQuadrantPattern[side],chamberNumber,quadrantNumber);
  }
  else if ( salias.Contains("Chamber") )
  {
    sscanf(salias.Data(),fgkDCSChamberPattern[side],&chamberNumber);
    ++chamberNumber;
    channelName =  TString::Format(fgkDCSChamberPattern[side],chamberNumber);
  }

  if ( TString(dcsAlias).Contains("iMon") )
  {
    channelName.ReplaceAll("vMon","iMon");
  }

  return channelName;
}

//_____________________________________________________________________________
TString
AliMpDCSNamer::DCSAliasFromName(const char* dcsName) const
{
  /// Convert a (possibly partial) dcsname to an alias (only for MCH HV)

  TString sname(dcsName);

  if ( !sname.Contains("MchHvLv") ) return dcsName;
  if ( sname.Contains("Group")) return dcsName; // not MCH HV (but LV)

  Int_t quadrantNumber(-1);
  Int_t chamberNumber(-1);
  Int_t side(-1);

  if ( sname.Contains("Left")) side = 0;
  if ( sname.Contains("Right")) side = 1;

  if ( side < 0 ) return "";

  TString channelName;

  if ( sname.Contains("Slat") )
  {
    Int_t slatNumber(-1);
    sscanf(sname.Data(),fgkDCSChannelSt345Pattern[side],&chamberNumber,&slatNumber);
    --chamberNumber;
    --slatNumber;
    channelName = TString::Format(fgkDCSChannelSt345Pattern[side],chamberNumber,slatNumber);
  }
  else if ( sname.Contains("Sect") )
  {
    Int_t sectorNumber(-1);
    sscanf(sname.Data(),fgkDCSChannelSt12Pattern[side],&chamberNumber,&quadrantNumber,&sectorNumber);
    --chamberNumber;
    --quadrantNumber;
    --sectorNumber;
    channelName =  TString::Format(fgkDCSChannelSt12Pattern[side],chamberNumber,quadrantNumber,sectorNumber);
  }
  else if ( sname.Contains("Quad") )
  {
    sscanf(sname.Data(),fgkDCSQuadrantPattern[side],&chamberNumber,&quadrantNumber);
    --chamberNumber;
    --quadrantNumber;
    channelName =  TString::Format(fgkDCSQuadrantPattern[side],chamberNumber,quadrantNumber);
  }
  else if ( sname.Contains("Chamber") )
  {
    sscanf(sname.Data(),fgkDCSChamberPattern[side],&chamberNumber);
    --chamberNumber;
    channelName =  TString::Format(fgkDCSChamberPattern[side],chamberNumber);
  }

  if ( TString(dcsName).Contains("iMon") )
  {
    channelName.ReplaceAll("vMon","iMon");
  }

  return channelName;
}

//_____________________________________________________________________________
Bool_t AliMpDCSNamer::DecodeDCSMCHLVAlias(const char* dcsAlias, Int_t*& detElemId, Int_t& numberOfDetectionElements, AliMp::PlaneType& planeType ) const
{
  /// Decode a MCH LV dcs alias in order to get :
  /// - the list of detection elements powered by this LV (between 1 for St 1-2 and max 4 DEs for St345 per LV group)
  /// - the plane type powered by this LV (only for St 1 and 2)

  TString salias(dcsAlias);

  detElemId = 0x0;
  planeType = AliMp::kBendingPlane;

  if (!salias.Contains("Group"))
  {
    // not a MCH LV alias
    return kFALSE;
  }

  Int_t side(-1);

  if ( salias.Contains("Left"))
  {
    side = 0;
  }
  else if ( salias.Contains("Right"))
  {
    side = 1;
  }
  else {
    AliError(Form("unexpected alias=%s",salias.Data()));
    return kFALSE;
  }

  Bool_t left = ( side == 0 );

  Int_t chamberNumber, groupNumber;
  TString voltageType;

  sscanf(salias.Data(),fgkDCSMCHLVGroupPattern[side],&chamberNumber,&groupNumber,&voltageType[0]);

  if ( chamberNumber >= 1 && chamberNumber <= 4 )
  {
    Int_t deOffset = 0;
    if ( groupNumber == 1 )
    {
      if ( left )
      {
        deOffset = 1;
        planeType = AliMp::kBendingPlane;
      }
      else
      {
        deOffset = 0;
        planeType = AliMp::kNonBendingPlane;
      }
    }
    else if ( groupNumber == 2 )
    {
      if ( left )
      {
        deOffset = 2;
        planeType = AliMp::kNonBendingPlane;
      }
      else
      {
        deOffset = 3;
        planeType = AliMp::kBendingPlane;
      }
    }
    else if ( groupNumber == 3 )
    {
      if ( left )
      {
        deOffset = 1;
        planeType = AliMp::kNonBendingPlane;
      }
      else
      {
        deOffset = 0;
        planeType = AliMp::kBendingPlane;
      }
    }
    else if ( groupNumber == 4 )
    {
      if ( left )
      {
        deOffset = 2;
        planeType = AliMp::kBendingPlane;
      }
      else
      {
        deOffset = 3;
        planeType = AliMp::kNonBendingPlane;
      }
    }
    else
    {
      AliError(Form("Got incorrect group number=%d from alias=%s",groupNumber,salias.Data()));
      return kFALSE;
    }

    numberOfDetectionElements=1;
    detElemId = new Int_t[numberOfDetectionElements];
    detElemId[0] = chamberNumber*100 + deOffset;
  }
  else if ( chamberNumber >= 5 && chamberNumber <= 10 )
  {
    Int_t* dcsSlatNumber(0x0);

    if ( chamberNumber >= 5 && chamberNumber <= 6 )
    {
      if ( groupNumber == 1 )
      {
        numberOfDetectionElements=3;
        dcsSlatNumber=new Int_t[numberOfDetectionElements];
        dcsSlatNumber[0]=1;
        dcsSlatNumber[1]=2;
        dcsSlatNumber[2]=3;
      }
      else if ( groupNumber == 5 )
      {
        numberOfDetectionElements=3;
        dcsSlatNumber=new Int_t[numberOfDetectionElements];
        dcsSlatNumber[0]=7;
        dcsSlatNumber[1]=8;
        dcsSlatNumber[2]=9;
      }
      else if ( groupNumber > 1 && groupNumber < 5 )
      {
        numberOfDetectionElements=1;
        dcsSlatNumber=new Int_t[numberOfDetectionElements];
        dcsSlatNumber[0]=groupNumber+2;
      }
      else
      {
        AliError(Form("Got incorrect group number=%d from alias=%s",groupNumber,salias.Data()));
        return kFALSE;
      }
    }
    else if ( chamberNumber >= 7 )
    {
      if ( groupNumber == 1 )
      {
        numberOfDetectionElements=4;
        dcsSlatNumber=new Int_t[numberOfDetectionElements];
        dcsSlatNumber[0]=1;
        dcsSlatNumber[1]=2;
        dcsSlatNumber[2]=3;
        dcsSlatNumber[3]=4;
      }
      else if ( groupNumber == 7 )
      {
        numberOfDetectionElements=4;
        dcsSlatNumber=new Int_t[numberOfDetectionElements];
        dcsSlatNumber[0]=10;
        dcsSlatNumber[1]=11;
        dcsSlatNumber[2]=12;
        dcsSlatNumber[3]=13;
      }
      else if ( groupNumber > 1 && groupNumber < 7 )
      {
        numberOfDetectionElements=1;
        dcsSlatNumber=new Int_t[numberOfDetectionElements];
        dcsSlatNumber[0]=groupNumber+3;
      }
      else
      {
        AliError(Form("Got incorrect group number=%d from alias=%s",groupNumber,salias.Data()));
        return kFALSE;
      }
    }

    detElemId = new Int_t[numberOfDetectionElements];
    for ( int i = 0; i < numberOfDetectionElements; ++i )
    {
      detElemId[i] = DCS2DE(chamberNumber-1,side,dcsSlatNumber[i]-1);
    }
  }
  else
  {
    AliError(Form("Got an invalid chamberNumber=%d from alias=%s",chamberNumber,salias.Data()));
    return kFALSE;
  }

  return kTRUE;
}

//_____________________________________________________________________________
TString
AliMpDCSNamer::DCSMCHLVAliasName(Int_t detElemId, Int_t voltageType, AliMp::PlaneType planeType) const
{
  TString voltage;

  if ( voltageType == -1 ) {
    voltage = "ann";
  } else if ( voltageType == 0 ) {
    voltage = "dig";
  } else if ( voltageType == 1 ) {
    voltage = "anp";
  }

  if ( !voltage.Length() )
  {
    AliError(Form("Incorrect voltageType=%d. Expected -1,0 or 1.",voltageType));
    return "";
  }

  Int_t chamberId = AliMpDEManager::GetChamberId(detElemId);
  if ( chamberId < 0 )
  {
    AliError(Form("Got an incorrect chamberId=%d from detElemId=%d",chamberId,detElemId));
    return "";
  }
  Int_t stationId = 1 + chamberId / 2;

  Int_t side(-1);
  Int_t cham;
  Int_t dcsNumber = DetElemId2DCS(detElemId, side, cham);

  Int_t group(0);

  switch (AliMpDEManager::GetStationType(detElemId))
  {
    case AliMp::kStation12:
      {
        // For Chamber 1 to 4 Left the relationship between DCS GUI names and groups is:
        // Quad2B    --> Group3 = DE x01 Non Bending
        // Quad2F    --> Group1 = DE x01 Bending
        // Quad3B    --> Group4 = DE x02 Bending
        // Quad3F    --> Group2 = DE x02 Non Bending
        // for Chamber 1 to 4 Right the relationship is:
        // Quad1B    --> Group3 = DE x00 Bending
        // Quad1F    --> Group1 = DE x00 Non  Bending
        // Quad4B    --> Group4 = DE x03 Non Bending
        // Quad4F    --> Group2 = DE x03 Bending
        // where x = 1,2,3,4
        // and Quad#B = Back = towards IP = cath1
        // while Quad#F = Front = towards muon trigger = cath0
        //

        Int_t remnant = detElemId % 100;
        switch (remnant) {
          case 0: // DE x00
            group = ( planeType == AliMp::kBendingPlane ) ? 3 : 1;
            break;
          case 1: // DE x01
            group = ( planeType == AliMp::kBendingPlane ) ? 1 : 3;
            break;
          case 2: // DE x02
            group = ( planeType == AliMp::kBendingPlane ) ? 4 : 2;
            break;
          case 3: // DE x03
            group = ( planeType == AliMp::kBendingPlane ) ? 2 : 4;
            break;
          default:
            AliFatal("");
            break;
        }
      }
      break;
    case AliMp::kStation345:
      {
        Int_t dcsSlatNumber = 1 + dcsNumber;
        if ( stationId == 3 ) {
          switch (dcsSlatNumber) {
            case 1:
            case 2:
            case 3:
              group = 1;
              break;
            case 4:
            case 5:
            case 6:
              group = dcsSlatNumber - 2;
              break;
            case 7:
            case 8:
            case 9:
              group = 5;
              break;
            default:
              break;
          }
        }
        else
        {
          switch (dcsSlatNumber) {
            case 1:
            case 2:
            case 3:
            case 4:
              group = 1;
              break;
            case 5:
            case 6:
            case 7:
            case 8:
            case 9:
              group = dcsSlatNumber - 3;
              break;
            case 10:
            case 11:
            case 12:
            case 13:
              group = 7;
              break;
            default:
              break;
          }
        }
      }
      break;
    default:
      break;
  }

  if ( group == 0 ) {
    AliError(Form("Could not get LV group id for detection element %d",detElemId));
    return "";
  }

  TString aliasName;

  aliasName.Form(fgkDCSMCHLVGroupPattern[side],chamberId+1,group,voltage.Data());

  return aliasName;
}

//_____________________________________________________________________________
TString
AliMpDCSNamer::DCSAliasName(Int_t detElemId, Int_t sector, Int_t dcsMeasure) const
{
  /// Return the alias name of the DCS Channel for a given DCS area
  /// \param detElemId
  /// \param sector = 0,1 or 2 for St12, and is unused for st345 and trigger
  /// \param dcsMeasure = kDCSHV, kDCSI

  Int_t chamberId = AliMpDEManager::GetChamberId(detElemId);
  if ( chamberId < 0 ) return "";

  Int_t side(-1), chId(-1);
  Int_t dcsNumber = DetElemId2DCS(detElemId,side,chId);

  TString aliasName;

  switch (AliMpDEManager::GetStationType(detElemId))
  {
    case AliMp::kStation12:
      aliasName.Form(fgkDCSChannelSt12Pattern[side],chamberId,dcsNumber,sector);
      break;
    case AliMp::kStation345:
      aliasName.Form(fgkDCSChannelSt345Pattern[side],chamberId,dcsNumber);
      break;
    case AliMp::kStationTrigger:
      return TString::Format(fgkDCSChannelTriggerPattern[side],fgkDCSSideTriggerName[side],chId,dcsNumber,fgkDCSMeasureName[dcsMeasure]);
      break;
    default:
      return "";
      break;
  }

  if ( dcsMeasure == AliMpDCSNamer::kDCSI )
  {
    aliasName.ReplaceAll("vMon","iMon");
  }

  return aliasName;

}

//_____________________________________________________________________________
TString
AliMpDCSNamer::DCSSwitchAliasName(Int_t detElemId, Int_t pcbNumber) const
{
  /// Return the alias name of the DCS Switch for a given PCB
  /// within a slat of St345

  if (AliMpDEManager::GetStationType(detElemId) == AliMp::kStation345)
  {
    return TString::Format(fgkDCSSwitchSt345Pattern,detElemId,pcbNumber);
  }
  return "";
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
    sscanf(sDcsAlias.Data(),fgkDCSChannelTriggerPatternRead[side],cside,&n1,&n3,type);
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
AliMpDCSNamer::GenerateAliases(const char* pattern) const
{
  /// Generate DCS alias names, for MUON Tracker High and Low Voltage systems
  /// or for MUON Trigger HV and current system.
  ///
  /// For MCH we first generate 188 aliases of HV DCS channels :
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
  /// then 600 aliases of DCS switches (only for St345) : 1 switch per PCB.
  ///
  /// and finally 324 LV groups (108 per voltage x 3 voltages)
  ///
  /// St 1 ch  1 left or right : 4 groups
  ///      ch  2 left or right : 4 groups
  /// St 2 ch  3 left or right : 4 groups
  ///      ch  4 left or right : 4 groups
  /// St 3 ch  5 left or right : 5 groups
  ///      ch  6 left or right : 5 groups
  /// St 4 ch  7 left or right : 7 groups
  ///      ch  8 left or right : 7 groups
  /// St 5 ch  9 left or right : 7 groups
  ///      ch 10 left or right : 7 groups
  ///
  /// Returns a TObjArray of TObjString(=alias name)

  TObjArray* aliases = new TObjArray;
  aliases->SetOwner(kTRUE);

  Int_t nMeasures = (fDetector == kTriggerDet) ? kNDCSMeas : 1;

  for(Int_t iMeas=0; iMeas<nMeasures; iMeas++){

    AliMpDEIterator it;

    it.First();

    Int_t voltageType[] = { -1,0,1 };

    while (!it.IsDone())
    {
      Int_t detElemId = it.CurrentDEId();

      switch (fDetector){
        case kTrackerDet:
        {
          switch ( AliMpDEManager::GetStationType(detElemId) )
          {
            case AliMp::kStation12:
            {
              for ( int sector = 0; sector < 3; ++sector)
              {
                // HV voltage
                aliases->Add(new TObjString(DCSAliasName(detElemId,sector)));
                // HV current
                aliases->Add(new TObjString(DCSAliasName(detElemId,sector,AliMpDCSNamer::kDCSI)));
              }

              AliMp::PlaneType planeType[] = { AliMp::kBendingPlane, AliMp::kNonBendingPlane };

              // LV groups, one per voltage (analog negative, digital, analog positive)
              // per plane (bending / non-bending)
              for ( int pt = 0; pt < 2; ++pt )
              {
                for ( int i = 0; i < 3; ++i )
                {
                  TString name = DCSMCHLVAliasName(detElemId,voltageType[i],planeType[pt]);
                  aliases->Add(new TObjString(name));
                }
              }
            }
            break;
            case AliMp::kStation345:
            {
              // HV voltage
              aliases->Add(new TObjString(DCSAliasName(detElemId)));
              // HV current
              aliases->Add(new TObjString(DCSAliasName(detElemId,0,AliMpDCSNamer::kDCSI)));
              // HV switches
              for ( Int_t i = 0; i < NumberOfPCBs(detElemId); ++i )
              {
                aliases->Add(new TObjString(DCSSwitchAliasName(detElemId,i)));
              }
              // LV groups, one per voltage (analog negative, digital, analog positive)
              for ( int i = 0; i < 3; ++i )
              {
                TString name = DCSMCHLVAliasName(detElemId,voltageType[i]);
                // for Station345 some detection elements share the same voltage group,
                // so we must insure we're not adding the same one several times
                if (!aliases->FindObject(name))
                {
                  aliases->Add(new TObjString(name));
                }
              }
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
              aliases->Add(new TObjString(DCSAliasName(detElemId,0,iMeas)));
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

  if (pattern && strlen(pattern)>0)
  {
    // remove aliases not containing the input pattern
    TObjArray* tmp = new TObjArray;
    tmp->SetOwner(kTRUE);
    for ( Int_t i = 0; i <= aliases->GetLast(); ++i )
    {
       TString name = static_cast<TObjString*>(aliases->At(i))->String();
       if (name.Contains(pattern))
       {
         tmp->Add(new TObjString(name.Data()));
       }
    }
    delete aliases;
    aliases = tmp;
  }
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

  AliCodeTimerAuto("",0)

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

  AliCodeTimerAuto("",0)

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
    else if ( x < 606.25 ) isector = 1;
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

//_____________________________________________________________________________
Bool_t AliMpDCSNamer::TestMCHLV() const
{
  /// Circular test of DCSMCHLVAliasName and DecodeDCSMCHLVAlias methods

  TObjArray* aliases = GenerateAliases("Group");

  const char* voltageTypeName[] = { "ann","dig","anp"};
  Int_t voltageType[] = { -1,0,1 };
  Int_t n(0);

  for ( Int_t iv = 0; iv < 3; ++iv )
  {
    for ( Int_t i = 0; i < aliases->GetEntries(); ++i )
    {
      TString dcsAlias = static_cast<TObjString*>(aliases->At(i))->String();
      if ( !dcsAlias.Contains(voltageTypeName[iv])) continue;

      Int_t* detElemId(0x0);
      Int_t numberOfDetectionElements(0);
      AliMp::PlaneType planeType;

      Bool_t ok = DecodeDCSMCHLVAlias(dcsAlias.Data(), detElemId, numberOfDetectionElements, planeType);

      if (!ok)
      {
        AliError(Form("Could not decode alias=%s",dcsAlias.Data()));
        delete[] detElemId;
        return kFALSE;
      }

      for ( Int_t id = 0; id < numberOfDetectionElements; ++id )
      {
        TString check =
          DCSMCHLVAliasName(detElemId[id], voltageType[iv], planeType);

        if (check!=dcsAlias)
        {
          AliError(Form("%s != %s",check.Data(),dcsAlias.Data()));
          return kFALSE;
        }
        else
        {
          ++n;
        }
      }

      delete[] detElemId;
    }
  }

  delete aliases;

  AliInfo(Form("%d aliases successfully tested",n));

  return kTRUE;
}
