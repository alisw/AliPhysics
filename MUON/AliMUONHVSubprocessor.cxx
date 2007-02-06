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

/// $Id$
/// \class AliMUONHVSubprocessor
///
/// A subprocessor to read HV values for one run
///
/// It simply creates a copy of the dcsAliasMap w/o information
/// from the MUON TRK, and dumps this copy into the CDB
///
// Author: Laurent Aphecetche, Subatech

#include "AliMUONHVSubprocessor.h"
#include "AliMUONHVNamer.h"
#include "AliMUONPreprocessor.h"

#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"

#include "AliCDBMetaData.h"
#include "AliLog.h"

#include "Riostream.h"
#include "TMap.h"
#include "TObjString.h"

/// \cond CLASSIMP
ClassImp(AliMUONHVSubprocessor)
/// \endcond

//_____________________________________________________________________________
AliMUONHVSubprocessor::AliMUONHVSubprocessor(AliMUONPreprocessor* master)
: AliMUONVSubprocessor(master,
                       "HV",
                       "Get MUON Tracker HV values from DCS")
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONHVSubprocessor::~AliMUONHVSubprocessor()
{
  /// dtor
}

//_____________________________________________________________________________
void 
AliMUONHVSubprocessor::Initialize(Int_t run, 
                                  UInt_t startTime, UInt_t endTime)
{
  /// Initialisation for a given run
  AliDebug(1,Form("run %d startTime %ld endtime %ld",run,startTime,endTime));
}

//_____________________________________________________________________________
UInt_t
AliMUONHVSubprocessor::Process(TMap* dcsAliasMap)
{
  /// Make another alias map from dcsAliasMap, considering only MUON TRK aliases.

  TMap hv;
  hv.SetOwner(kTRUE);
  
  AliMUONHVNamer hvNamer;

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
      case AliMp::kStation1:
      case AliMp::kStation2:
      {
        for ( int i = 0; i <3; ++i)
        {
          aliases.Add(new TObjString(hvNamer.DCSHVChannelName(detElemId,i)));
        }
      }
      break;
      case AliMp::kStation345:
      {
        aliases.Add(new TObjString(hvNamer.DCSHVChannelName(detElemId)));
        for ( int i = 0; i < hvNamer.NumberOfPCBs(detElemId); ++i)
        {
          aliases.Add(new TObjString(hvNamer.DCSHVSwitchName(detElemId,i)));
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
  
  while ( ( alias = static_cast<TObjString*>(next()) ) ) 
  {
    TString aliasName(alias->String());
    TPair* hvPair = static_cast<TPair*>(dcsAliasMap->FindObject(aliasName.Data()));
    if (!hvPair)
    {
      AliError(Form("Did not find expected alias (%s)",aliasName.Data()));
    }
    else
    {
      TObjArray* values = static_cast<TObjArray*>(hvPair->Value());
      if (!values)
      {
        AliError(Form("Could not get values for alias (%s)",aliasName.Data()));
      }
      else
      {
        //FIXME : should insure here that values are only the ones within run
        //limits (startTime<timestamp<endTime)
        hv.Add(new TObjString(aliasName.Data()),values);
      }
    }
  }
  
  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("MUON TRK");
  metaData.SetComment("Computed by AliMUONHVSubprocessor $Id$");
  
  Bool_t validToInfinity(kFALSE);
  
  UInt_t result = Master()->Store("Calib","HV",&hv,&metaData,0,validToInfinity);
  
  return result;
}

