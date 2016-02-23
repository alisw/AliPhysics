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

//-----------------------------------------------------------------------------
/// \class AliMUONHVSubprocessor
///
/// A subprocessor to read HV values for one run
///
/// It simply creates a copy of the dcsAliasMap w/o information
/// from the MUON TRK, and dumps this copy into the CDB
///
/// \author Laurent Aphecetche, Subatech
//-----------------------------------------------------------------------------

#include "AliMUONHVSubprocessor.h"
#include "AliMUONPreprocessor.h"

#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpDCSNamer.h"

#include "AliCDBMetaData.h"
#include "AliLog.h"

#include "Riostream.h"
#include "TMap.h"
#include "TObjString.h"

#include "AliMUONCalibrationData.h"

/// \cond CLASSIMP
ClassImp(AliMUONHVSubprocessor)
/// \endcond

//_____________________________________________________________________________
AliMUONHVSubprocessor::AliMUONHVSubprocessor(AliMUONPreprocessor* master, Bool_t includeHVcurrents)
: AliMUONVSubprocessor(master,
                       "HV",
                       "Get MUON Tracker HV values from DCS"), fIncludeHVCurrents(includeHVcurrents)
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONHVSubprocessor::~AliMUONHVSubprocessor()
{
  /// dtor
}

//_____________________________________________________________________________
UInt_t
AliMUONHVSubprocessor::Process(TMap* dcsAliasMap)
{
  /// Make another alias map from dcsAliasMap, considering only MUON TRK aliases.

  TMap hv;
  hv.SetOwner(kTRUE);
  
  hv.SetUniqueID(AliMUONCalibrationData::PatchHVDCSAliasesSt1WasAppliedMask());
  
  AliMpDCSNamer hvNamer("TRACKER");

  AliMpDEIterator deIt;

  deIt.First();
  
  TObjArray aliases;
  aliases.SetOwner(kTRUE);
  
  // we first generate a list of expected MCH DCS aliases we'll then look for
  
  while ( !deIt.IsDone() )
  {
    Int_t detElemId = deIt.CurrentDEId();
    
    switch ( AliMpDEManager::GetStationType(detElemId) )
    {
      case AliMp::kStation12:
      {
        for ( int i = 0; i <3; ++i)
        {
          aliases.Add(new TObjString(hvNamer.DCSAliasName(detElemId,i)));
          if ( fIncludeHVCurrents )
          {
            aliases.Add(new TObjString(hvNamer.DCSAliasName(detElemId,i,AliMpDCSNamer::kDCSI)));
          }
        }
      }
      break;
      case AliMp::kStation345:
      {
        aliases.Add(new TObjString(hvNamer.DCSAliasName(detElemId)));
        if ( fIncludeHVCurrents )
        {
          aliases.Add(new TObjString(hvNamer.DCSAliasName(detElemId,0,AliMpDCSNamer::kDCSI)));
        }
        for ( int i = 0; i < hvNamer.NumberOfPCBs(detElemId); ++i)
        {
          aliases.Add(new TObjString(hvNamer.DCSSwitchAliasName(detElemId,i)));
        }
      }
      break;
      default:
        break;
    };

    deIt.Next();
  }

  TIter next(&aliases);
  TObjString* alias;
  Bool_t kNoAliases(kTRUE);
  Int_t aliasNotFound(0);
  Int_t valueNotFound(0);
  
  TObjArray temporaryOut; // due to a bug in PVSS some iMon aliases cause problem
  // we remove them for the moment.
  temporaryOut.SetOwner(kTRUE);
  
  temporaryOut.Add(new TObjString("MchHvLvLeft/Chamber06Left/Slat10.actual.iMon"));
  temporaryOut.Add(new TObjString("MchHvLvLeft/Chamber06Left/Slat11.actual.iMon"));
  temporaryOut.Add(new TObjString("MchHvLvLeft/Chamber06Left/Slat12.actual.iMon"));
  temporaryOut.Add(new TObjString("MchHvLvLeft/Chamber07Left/Slat10.actual.iMon"));
  temporaryOut.Add(new TObjString("MchHvLvLeft/Chamber07Left/Slat11.actual.iMon"));
  temporaryOut.Add(new TObjString("MchHvLvLeft/Chamber07Left/Slat12.actual.iMon"));
  temporaryOut.Add(new TObjString("MchHvLvLeft/Chamber08Left/Slat10.actual.iMon"));
  temporaryOut.Add(new TObjString("MchHvLvLeft/Chamber08Left/Slat11.actual.iMon"));
  temporaryOut.Add(new TObjString("MchHvLvLeft/Chamber08Left/Slat12.actual.iMon"));
  temporaryOut.Add(new TObjString("MchHvLvLeft/Chamber09Left/Slat10.actual.iMon"));
  temporaryOut.Add(new TObjString("MchHvLvLeft/Chamber09Left/Slat11.actual.iMon"));
  temporaryOut.Add(new TObjString("MchHvLvLeft/Chamber09Left/Slat12.actual.iMon"));
  temporaryOut.Add(new TObjString("MchHvLvRight/Chamber06Right/Slat10.actual.iMon"));
  temporaryOut.Add(new TObjString("MchHvLvRight/Chamber06Right/Slat11.actual.iMon"));
  temporaryOut.Add(new TObjString("MchHvLvRight/Chamber06Right/Slat12.actual.iMon"));
  temporaryOut.Add(new TObjString("MchHvLvRight/Chamber07Right/Slat10.actual.iMon"));
  temporaryOut.Add(new TObjString("MchHvLvRight/Chamber07Right/Slat11.actual.iMon"));
  temporaryOut.Add(new TObjString("MchHvLvRight/Chamber07Right/Slat12.actual.iMon"));
  temporaryOut.Add(new TObjString("MchHvLvRight/Chamber08Right/Slat10.actual.iMon"));
  temporaryOut.Add(new TObjString("MchHvLvRight/Chamber08Right/Slat11.actual.iMon"));
  temporaryOut.Add(new TObjString("MchHvLvRight/Chamber08Right/Slat12.actual.iMon"));
  temporaryOut.Add(new TObjString("MchHvLvRight/Chamber09Right/Slat10.actual.iMon"));
  temporaryOut.Add(new TObjString("MchHvLvRight/Chamber09Right/Slat11.actual.iMon"));
  temporaryOut.Add(new TObjString("MchHvLvRight/Chamber09Right/Slat12.actual.iMon"));
  
  while ( ( alias = static_cast<TObjString*>(next()) ) ) 
  {
    TString aliasName(alias->String());
    
    if ( temporaryOut.FindObject(alias->String().Data() ) )
    {
      // skip problematic aliases
      continue;
    }
    
    TPair* hvPair = static_cast<TPair*>(dcsAliasMap->FindObject(aliasName.Data()));
    if (!hvPair)
    {
      ++aliasNotFound;
    }
    else
    {
      kNoAliases = kFALSE;
      TObjArray* values = static_cast<TObjArray*>(hvPair->Value()->Clone());
      if (!values)
      {
        ++valueNotFound;
      }
      else
      {
        RemoveValuesOutsideRun(values);
        hv.Add(new TObjString(aliasName.Data()),values);
      }
    }
  }
  
  if ( kNoAliases ) 
  {
    Master()->Log("ERROR : no DCS values found");
    return 1;
  }
  
  if ( aliasNotFound ) 
  {
    Master()->Log(Form("WARNING %d aliases not found",aliasNotFound));
  }
  
  if ( valueNotFound )
  {
    Master()->Log(Form("WARNING %d values not found",valueNotFound));
  }
  
  Master()->Log("INFO Aliases successfully read in");
  
  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("MUON TRK");
  metaData.SetComment("Computed by AliMUONHVSubprocessor $Id$");
  
  Bool_t validToInfinity(kFALSE);
  
  Bool_t result = Master()->Store("Calib","HV",&hv,&metaData,0,validToInfinity);
  
  return ( result != kTRUE); // return 0 if everything is ok
}

