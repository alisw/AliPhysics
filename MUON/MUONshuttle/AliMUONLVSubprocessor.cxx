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
/// \class AliMUONLVSubprocessor
///
/// A subprocessor to read LV values for one run
///
/// It simply creates extract from the dcsAliasMap the information
/// from the MUON TRK Low Voltages, and dumps this into the CDB
///
/// \author Laurent Aphecetche, Subatech
//-----------------------------------------------------------------------------

#include "AliMUONLVSubprocessor.h"
#include "AliMUONPreprocessor.h"

#include "AliMpDCSNamer.h"

#include "AliCDBMetaData.h"
#include "AliLog.h"

#include "Riostream.h"
#include "TMap.h"
#include "TObjString.h"

#include "AliMUONCalibrationData.h"

/// \cond CLASSIMP
ClassImp(AliMUONLVSubprocessor)
/// \endcond

//_____________________________________________________________________________
AliMUONLVSubprocessor::AliMUONLVSubprocessor(AliMUONPreprocessor* master)
: AliMUONVSubprocessor(master,
                       "LV",
                       "Get MUON Tracker LV values from DCS")
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONLVSubprocessor::~AliMUONLVSubprocessor()
{
  /// dtor
}

//_____________________________________________________________________________
UInt_t
AliMUONLVSubprocessor::Process(TMap* dcsAliasMap)
{
  /// Make another alias map from dcsAliasMap, considering only MUON TRK LV aliases.

  TMap lv;
  lv.SetOwner(kTRUE);

  AliMpDCSNamer hvNamer("TRACKER");

  // we first generate a list of expected MCH DCS aliases we'll then look for
  TObjArray* aliases = hvNamer.GenerateAliases("Group");

  Master()->Log(Form("INFO : will look for %d LV aliases",aliases->GetEntries()));

  TIter next(aliases);
  TObjString* alias;
  Bool_t kNoAliases(kTRUE);
  Int_t valueNotFound(0);
  TList aliasesNotFound;
  aliasesNotFound.SetOwner(kTRUE);

  while ( ( alias = static_cast<TObjString*>(next()) ) )
  {
    TString aliasName(alias->String());

    TPair* lvPair = static_cast<TPair*>(dcsAliasMap->FindObject(aliasName.Data()));
    if (!lvPair)
    {
      aliasesNotFound.Add(new TObjString(aliasName));
    }
    else
    {
      kNoAliases = kFALSE;
      TObjArray* values = static_cast<TObjArray*>(lvPair->Value()->Clone());
      if (!values)
      {
        ++valueNotFound;
      }
      else
      {
        RemoveValuesOutsideRun(values);
        lv.Add(new TObjString(aliasName.Data()),values);
      }
    }
  }

  if ( kNoAliases )
  {
    Master()->Log("ERROR : no DCS values found");
    delete aliases;
    return 1;
  }

  if ( aliasesNotFound.GetEntries() )
  {
    Master()->Log(Form("WARNING %d aliases not found : ",aliasesNotFound.GetEntries()));
    TIter nextNotFound(&aliasesNotFound);
    TObjString* str;
    TString msg;
    while (( str = static_cast<TObjString*>(nextNotFound())))
    {
        msg += str->String();
        msg += "\n";
    }
    Master()->Log(msg.Data());
  }

  if ( valueNotFound )
  {
    Master()->Log(Form("WARNING %d values not found",valueNotFound));
  }

  Master()->Log(Form("INFO %d/%d aliases successfully read in.",aliases->GetEntries()-aliasesNotFound.GetEntries(),aliases->GetEntries()));

  AliCDBMetaData metaData;
  metaData.SetBeamPeriod(0);
  metaData.SetResponsible("MUON TRK");
  metaData.SetComment("Computed by AliMUONLVSubprocessor $Id$");

  Bool_t validToInfinity(kFALSE);

  Bool_t result = Master()->Store("Calib","LV",&lv,&metaData,0,validToInfinity);

  delete aliases;

  return ( result != kTRUE); // return 0 if everything is ok
}
