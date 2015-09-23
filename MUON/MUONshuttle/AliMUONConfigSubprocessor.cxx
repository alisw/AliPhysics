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

#include "AliMUONConfigSubprocessor.h"

#include "AliCDBMetaData.h"
#include "AliCDBEntry.h"
#include "AliDAQ.h"
#include "AliLog.h"
#include "AliMUON2DMap.h"
#include "AliMUON2DStoreValidator.h"
#include "AliMUONCalibParamNF.h"
#include "AliMUONPreprocessor.h"
#include "AliMUONTrackerIO.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include "TObjString.h"
#include <Riostream.h>
#include <TList.h>
#include <TObjString.h>
#include <TSystem.h>
#include <sstream>

//-----------------------------------------------------------------------------
/// \class AliMUONConfigSubprocessor
///
/// Implementation of AliMUONVSubprocessor class to deal with MUON TRK readout configuration.
///
/// Configs are read in from an ascii file, with the format :               \n
///---------------------------------------------------------------------------\n
/// BUS_PATCH MANU_ADDR
///---------------------------------------------------------------------------\n
///
/// \author L. Aphecetche
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONConfigSubprocessor)
/// \endcond

//_____________________________________________________________________________
AliMUONConfigSubprocessor::AliMUONConfigSubprocessor(AliMUONPreprocessor* master)
: AliMUONVSubprocessor(master,
                       "Config",
                       "Upload MUON Tracker readout configuration to OCDB"),
fConfig(0x0),
fConfigChanged(kFALSE)
{
  /// default ctor
}

//_____________________________________________________________________________
AliMUONConfigSubprocessor::~AliMUONConfigSubprocessor()
{
  /// dtor
  delete fConfig;
}

//_____________________________________________________________________________
Bool_t
AliMUONConfigSubprocessor::HasConfigChanged(const AliMUONVStore& newConfig) const
{
  /// Check whether the config changed. 
  /// Any error will return kTRUE to trig a config upload (safer way).
  
  AliCDBEntry* entry = Master()->GetFromOCDB("Calib","Config");
  if (!entry)
  {
    AliError("Could not get MUON/Calib/Config entry for current run ! That's not OK !");
    return kTRUE;
  }
  AliMUONVStore* oldConfig = dynamic_cast<AliMUONVStore*>(entry->GetObject());
  if (!oldConfig)
  {
    AliError("Could not get MUON/Calib/Config object for current run (wrong type ?) ! That's not OK !");
    return kTRUE;
  }
  
  if ( oldConfig->GetSize() != newConfig.GetSize() ) 
  {
    return kTRUE;
  }
  
  TIter next(oldConfig->CreateIterator());
  AliMUONVCalibParam* old;
  
  while ( ( old = static_cast<AliMUONVCalibParam*>(next()) ) )
  {
    Int_t detElemId = old->ID0();
    Int_t manuId = old->ID1();
    
    if ( ! newConfig.FindObject(detElemId,manuId) )
    {
      // not found in new. Configurations are different. Return right now.
      return kTRUE;
    }
  }
  
  // all tests OK. Configuration has not changed.
  return kFALSE;
}


//_____________________________________________________________________________
Bool_t 
AliMUONConfigSubprocessor::Initialize(Int_t run, UInt_t startTime, UInt_t endTime)
{
  /// When starting a new run, reads in the config ASCII files.
  
  const Int_t kSystem = AliMUONPreprocessor::kDAQ;
  const char* kIdConf = "CONFIG";

  delete fConfig;
  fConfig = new AliMUON2DMap(kTRUE);
  
  TList* sources = Master()->GetFileSources(kSystem,kIdConf);
  TIter nextConf(sources);
  Int_t nconf(0);
  Int_t nconfFiles(0);
  TObjString* o;
  
  while ( ( o = static_cast<TObjString*>(nextConf()) ) )
  {
    TString fileName(Master()->GetFile(kSystem,kIdConf,o->GetName()));
    Int_t ok = ReadConfigFile(fileName.Data());
    if (ok>0)
    {
      nconf += ok;
      ++nconfFiles;
    }
  }
  
  delete sources;
  
  if ( nconfFiles == 0 )
  {
    delete fConfig;
    fConfig = 0x0;
    fConfigChanged = kFALSE;
    AliInfo("No configuration files found for this run. That might be fine. Moving on...");
    return kTRUE;
  }

  fConfigChanged = HasConfigChanged(*fConfig);
  
  if (!fConfigChanged)
  {
    AliInfo("Will not upload the Configuration as it has not changed");
  }
  
  return kTRUE;
}

//_____________________________________________________________________________
UInt_t 
AliMUONConfigSubprocessor::Process(TMap* /*dcsAliasMap*/)
{
  /// Store the config into the CDB
  /// So far there is no way this sub processor can fail...
  
  if (!fConfig)
  {
    return 0;
  }
    
  if ( !fConfigChanged )
  {
    return 0;
  }
  
  Master()->Log("Storing readout configuration, as it has changed");
  
  AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("MUON TRK");
  TString comment("Computed by AliMUONConfigSubprocessor $Id$");
  comment.ReplaceAll("$","");
	metaData.SetComment(comment.Data());
  
  Bool_t validToInfinity = kFALSE;
	Bool_t result = Master()->Store("Calib", "Config", fConfig, &metaData, 0, validToInfinity);
  return ( result != kTRUE ); // return 0 if everything is ok.
}

//_____________________________________________________________________________
Int_t
AliMUONConfigSubprocessor::ReadConfigFile(const char* filename)
{
  /// Read the configuration from an ASCII file.                              
  /// Format of that file is one line per manu :                              
  /// BUS_PATCH MANU_ADDR
  /// Return kFALSE if reading was not successfull.                           
  ///
  
  TString sFilename(gSystem->ExpandPathName(filename));
  
  Master()->Log(Form("Reading %s",sFilename.Data()));
  
  Int_t n = AliMUONTrackerIO::ReadConfig(sFilename.Data(),*fConfig);
  
  switch (n)
  {
    case -1:
      Master()->Log(Form("Could not open %s",sFilename.Data()));
      break;
  }
  
  return n;
}


//_____________________________________________________________________________
void
AliMUONConfigSubprocessor::Print(Option_t* opt) const
{
  /// ouput to screen
  if (fConfig) fConfig->Print("",opt);
}

