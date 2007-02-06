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

#include "AliMUONHVNamer.h"

#include "AliMpArea.h"
#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpHelper.h"
#include "AliMpSegmentation.h"
#include "AliMpSlat.h"
#include "AliMpSlatSegmentation.h"
#include "AliMpVPadIterator.h"

#include "AliLog.h"
#include "Riostream.h"
#include "TMap.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TString.h"
#include "TSystem.h"

/// $Id$
///
/// \class AliMUONHVNamer
/// 
/// A utility class to manage HV DCS aliases names, in particular the
/// two conventions used to number the detection elements within a detector.
///
// Author: Laurent Aphecetche, Subatech

/// \cond CLASSIMP
ClassImp(AliMUONHVNamer)
/// \endcond

const char* AliMUONHVNamer::fgHVChannelSt345Pattern[] = 
{ "MchHvLvLeft/Chamber%02dLeft/Slat%02d.actual.vMon",
  "MchHvLvRight/Chamber%02dRight/Slat%02d.actual.vMon" 
};

const char* AliMUONHVNamer::fgHVChannelSt12Pattern[] = 
{
  "MchHvLvLeft/Chamber%02dLeft/Quad%dSect%d.actual.vMon",
  "MchHvLvRight/Chamber%02dRight/Quad%dSect%d.actual.vMon",
};
  
const char* AliMUONHVNamer::fgHVSwitchSt345Pattern = "MchDE%dsw%d.inValue";

//_____________________________________________________________________________
AliMUONHVNamer::AliMUONHVNamer()
{
 /// default ctor 
}

//_____________________________________________________________________________
AliMUONHVNamer::~AliMUONHVNamer()
{
  /// dtor
}

//_____________________________________________________________________________
TObjArray*
AliMUONHVNamer::CompactAliases() const
{
  /// Generate a compact list of aliases, for Shuttle test
  /// This one is completely hand-made, in contrast with GenerateAliases()
  /// method

  TObjArray* a = new TObjArray;
  a->SetOwner(kTRUE);
  
  // St 12 (HV Channels)
  a->Add(new TObjString("MchHvLvRight/Chamber[01..04]Right/Quad1Sect[1..3].actual.vMon"));
  a->Add(new TObjString("MchHvLvLeft/Chamber[01..04]Left/Quad2Sect[1..3].actual.vMon"));
  a->Add(new TObjString("MchHvLvLeft/Chamber[01..04]Left/Quad3Sect[1..3].actual.vMon"));
  a->Add(new TObjString("MchHvLvRight/Chamber[01..04]Right/Quad4Sect[1..3].actual.vMon"));
  
  // St345 (HV Channels)
  
  a->Add(new TObjString("MchHvLvRight/Chamber[05..10]Right/Slat[01..09].actual.vMon"));
  a->Add(new TObjString("MchHvLvLeft/Chamber[05..10]Left/Slat[01..09].actual.vMon"));

  a->Add(new TObjString("MchHvLvRight/Chamber[07..10]Right/Slat[10..13].actual.vMon"));
  a->Add(new TObjString("MchHvLvLeft/Chamber[07..10]Left/Slat[10..13].actual.vMon"));

  // St345 (HV Switches)
  AliMpDEIterator it;
  
  it.First();
  
  while (!it.IsDone())
  {
    Int_t detElemId = it.CurrentDEId();
    if ( AliMpDEManager::GetStationType(detElemId) == AliMp::kStation345 )
    {
          a->Add(new TObjString(Form("MchDE%dsw[1..%d].inValue",detElemId,NumberOfPCBs(detElemId))));
     }
    it.Next();
  }
  return a;
}

//_____________________________________________________________________________
Int_t 
AliMUONHVNamer::DCS2DE(Int_t chamberId, Int_t side, Int_t dcsNumber) const
{
  /// Convert DCS "slat number" (old convention) to DE (new) convention.
  ///
  /// @param chamberId : chamber number (starting at 1)
  /// @param side : 0 for Left, 1 for Right
  /// @param dcsNumber : slat number in DCS HV convention
  ///
  /// note that dcsNumber should be >=1 and <= number of DEs/2 in chamber

  Int_t nofDE = AliMpDEManager::GetNofDEInChamber(chamberId);
  
  Int_t half = nofDE/2;
  Int_t quarter = nofDE/4;
  Int_t threeQuarter = half + quarter;
  
  Int_t de(-1);
  
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
  
  return chamberId*100 + de;
}

//_____________________________________________________________________________
Int_t
AliMUONHVNamer::DetElemId2DCS(Int_t detElemId, Int_t& side) const
{
  /// Convert DE to DCS "slat number"
  /// @see DCS2DE
  
  Int_t chamberId = 1 + AliMpDEManager::GetChamberId(detElemId);
  if ( chamberId < 1 ) 
  {
    AliDebug(1,Form("DetElemId %d invalid",detElemId));
    return -1;
  }
  Int_t dcsNumber = (detElemId-chamberId*100);

  switch ( AliMpDEManager::GetStationType(detElemId) )
  {
    case AliMp::kStation1:
    case AliMp::kStation2:
    {
      switch (dcsNumber)
      {
        case 0:
        case 3:
          side = 1; // right
          ++dcsNumber;
          break;
        case 1:
        case 2:
          side = 0; // left
          ++dcsNumber;
        default:
          break;
      }
    }
      break;
    case AliMp::kStation345:
    {
      Int_t nofDE = AliMpDEManager::GetNofDEInChamber(chamberId-1);
      
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
    }
        break;
    default:
      break;
  }
  return dcsNumber;
}

//_____________________________________________________________________________
const char* 
AliMUONHVNamer::DCSHVChannelName(Int_t detElemId, Int_t sector) const
{
  /// Return the alias name of the HV Channel for a given HV area 
  /// @param sector=0,1 or 2 for St12, and is unused for st345
  
  Int_t chamberId = 1 + AliMpDEManager::GetChamberId(detElemId);
  if ( chamberId < 1 ) return 0x0;

  Int_t side(-1);
  Int_t dcsNumber = DetElemId2DCS(detElemId,side);
                                  
  switch (AliMpDEManager::GetStationType(detElemId))
  {
    case AliMp::kStation1:
    case AliMp::kStation2:
      return Form(fgHVChannelSt12Pattern[side],chamberId,dcsNumber,sector+1);
      break;
    case AliMp::kStation345:
      return Form(fgHVChannelSt345Pattern[side],chamberId,dcsNumber);
      break;
    default:
      return 0x0;
      break;
  }
}

//_____________________________________________________________________________
const char* 
AliMUONHVNamer::DCSHVSwitchName(Int_t detElemId, Int_t pcbNumber) const
{
  /// Return the alias name of the HV Switch for a given PCB 
  /// within a slat of St345
  
  if (AliMpDEManager::GetStationType(detElemId) == AliMp::kStation345)
  {
    return Form(fgHVSwitchSt345Pattern,detElemId,pcbNumber+1);
  }
  return 0x0;
}

//_____________________________________________________________________________
Int_t 
AliMUONHVNamer::DetElemIdFromDCSAlias(const char* dcsAlias) const
{
  /// Converts the dcs alias to a detection element identifier
  ///
  /// dcsAlias has one of the following 2 forms :
  ///
  /// MchHvLv[Left|Right]/Chamber##[Left|Right]/Chamber##[Left|Right]Slat##.actual.vMon
  ///
  /// MchHvLv[Left|Right]/Chamber##[Left|Right]/Chamber##[Left|Right]Quad#Sect#.actual.vMon
  
  TString sDcsAlias(dcsAlias);
  
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
    return -2;
  }
  
  int n1(-1);
  int n3(-1);
  int n4(-1);
  int detElemId(-1);
  
  if ( sDcsAlias.Contains("Slat") )
  {
    sscanf(sDcsAlias.Data(),fgHVChannelSt345Pattern[side],&n1,&n3);
    detElemId = DCS2DE(n1,side,n3);
  }
  else if ( sDcsAlias.Contains("Quad") )
  {
    sscanf(sDcsAlias.Data(),fgHVChannelSt12Pattern[side],&n1,&n3,&n4);    
    detElemId = n3-1;
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
TObjArray*
AliMUONHVNamer::GenerateAliases() const
{
  /// Generate DCS alias names, for MUON Tracker High Voltage system.
  ///
  /// We first generate aliases of HV channels :
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
  /// then aliases of HV switches (only for St345) : 1 switch per PCB.
  ///
  /// Returns a TObjArray of TObjString(=alias name)
  
  TObjArray* aliases = new TObjArray;
  aliases->SetOwner(kTRUE);
  
  AliMpDEIterator it;
  
  it.First();
  
  while (!it.IsDone())
  {
    Int_t detElemId = it.CurrentDEId();
    switch ( AliMpDEManager::GetStationType(detElemId) )
    {
      case AliMp::kStation1:
      case AliMp::kStation2:
        for ( int sector = 0; sector < 3; ++sector)
        {
          aliases->Add(new TObjString(DCSHVChannelName(detElemId,sector)));
        }
        break;
      case AliMp::kStation345:
        aliases->Add(new TObjString(DCSHVChannelName(detElemId)));
        for ( Int_t i = 0; i < NumberOfPCBs(detElemId); ++i )
        {
          aliases->Add(new TObjString(DCSHVSwitchName(detElemId,i)));
        }
        break;
      default:
        break;
    }
    it.Next();
  }
  
  return aliases;
}

//_____________________________________________________________________________
Int_t 
AliMUONHVNamer::ManuId2PCBIndex(Int_t detElemId, Int_t manuId) const
{
  /// Returns the index of PCB (within a St345 slat) for a given manu number.
  /// Returns -1 if (detElemId,manuId) is incorrect
  
  const AliMpSlatSegmentation* seg = static_cast<const AliMpSlatSegmentation*>
    (AliMpSegmentation::Instance()->GetMpSegmentationByElectronics(detElemId,manuId));
  const AliMpSlat* slat = seg->Slat();
  AliMpVPadIterator* it = seg->CreateIterator();
  it->First();
  AliMpPad pad = it->CurrentItem();
  
  while ( !it->IsDone() && pad.GetLocation().GetFirst() != manuId ) 
  {
    it->Next();
    pad = it->CurrentItem();
  }
  
  Int_t pcbIndex = slat->FindPCBIndex(pad.Position().X()+slat->Position().X(),
                                      pad.Position().Y()+slat->Position().Y());
//  AliDebug(1,Form("pcbIndex %d",pcbIndex));
//  StdoutToAliDebug(1,pad.Print());
  return pcbIndex;
}

//_____________________________________________________________________________
Int_t 
AliMUONHVNamer::ManuId2Sector(Int_t detElemId, Int_t manuId) const
{
  /// Return the HV-sector number (within a St12 quadrant) for a given manu number.
  
  //FIXME: write me !
  return 0;
}

//_____________________________________________________________________________
Int_t 
AliMUONHVNamer::NumberOfPCBs(Int_t detElemId) const
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
    const AliMpSlatSegmentation* seg = static_cast<const AliMpSlatSegmentation*>
    (AliMpSegmentation::Instance()->GetMpSegmentation(detElemId,AliMp::kCath0));
    const AliMpSlat* slat = seg->Slat();
    return slat->GetSize();
  }
}
