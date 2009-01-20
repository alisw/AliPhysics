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
/// \class AliMUONTriggerDCSSubprocessor
///
/// A subprocessor to read Trigger DCS values for one run
///
/// It simply creates a copy of the dcsAliasMap w/o information
/// from the MUON TRG, and dumps this copy into the CDB
///
/// \author Diego Stocco, Subatech
//-----------------------------------------------------------------------------

#include "AliMUONTriggerDCSSubprocessor.h"
#include "AliMUONPreprocessor.h"

#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpConstants.h"
#include "AliMpDCSNamer.h"

#include "AliCDBMetaData.h"
#include "AliLog.h"

#include "Riostream.h"
#include "TMap.h"
#include "TObjString.h"

/// \cond CLASSIMP
ClassImp(AliMUONTriggerDCSSubprocessor)
/// \endcond

//_____________________________________________________________________________
AliMUONTriggerDCSSubprocessor::AliMUONTriggerDCSSubprocessor(AliMUONPreprocessor* master)
: AliMUONVSubprocessor(master,
                       "TriggerDCS",
                       "Get MUON Trigger HV and Current values from DCS")
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONTriggerDCSSubprocessor::~AliMUONTriggerDCSSubprocessor()
{
  /// dtor
}

//_____________________________________________________________________________
UInt_t
AliMUONTriggerDCSSubprocessor::Process(TMap* dcsAliasMap)
{
  /// Make another alias map from dcsAliasMap, considering only MUON TRK aliases.

  TMap dcsMap;
  dcsMap.SetOwner(kTRUE);
  
  AliMpDCSNamer dcsMapNamer("TRIGGER");

  AliMpDEIterator deIt;

  deIt.First();
  
  TObjArray aliases;
  aliases.SetOwner(kTRUE);
  
  // we first generate a list of expected MTR DCS aliases we'll then look for
  
  while ( !deIt.IsDone() )
  {
    Int_t detElemId = deIt.CurrentDEId();
    
    if ( AliMpDEManager::GetStationType(detElemId) == AliMp::kStationTrigger) {

      for(Int_t iMeas=0; iMeas<AliMpDCSNamer::kNDCSMeas; iMeas++){
	aliases.Add(new TObjString(dcsMapNamer.DCSChannelName(detElemId, 0, iMeas)));
      }

    }

    deIt.Next();
  }

  TIter next(&aliases);
  TObjString* alias;
  Bool_t kNoAliases(kTRUE);
  Int_t aliasNotFound(0);
  Int_t valueNotFound(0);
  
  while ( ( alias = static_cast<TObjString*>(next()) ) ) 
  {
    TString aliasName(alias->String());
    TPair* dcsMapPair = static_cast<TPair*>(dcsAliasMap->FindObject(aliasName.Data()));
    if (!dcsMapPair)
    {
      ++aliasNotFound;
    }
    else
    {
      kNoAliases = kFALSE;
      if (!dcsMapPair->Value())
      {
        ++valueNotFound;
      }
      else
      {
	TObjArray* values = static_cast<TObjArray*>(dcsMapPair->Value()->Clone());
        //FIXME : should insure here that values are only the ones within run
        //limits (startTime<timestamp<endTime)
        dcsMap.Add(new TObjString(aliasName.Data()),values);
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
  metaData.SetResponsible("MUON TRG");
  metaData.SetComment("Computed by AliMUONTriggerDCSSubprocessor $Id$");
  
  Bool_t validToInfinity(kFALSE);
  
  Bool_t result = Master()->Store("Calib","TriggerDCS",&dcsMap,&metaData,0,validToInfinity);
  
  return ( result != kTRUE); // return 0 if everything is ok
}

