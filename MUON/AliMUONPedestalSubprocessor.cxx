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

#include "AliMUONPedestalSubprocessor.h"

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
/// \class AliMUONPedestalSubprocessor
///
/// Implementation of AliMUONVSubprocessor class to deal with MUON TRK pedestals.
///
/// Pedestals are read in from an ascii file, with the format :               \n
///---------------------------------------------------------------------------\n
/// BUS_PATCH MANU_ADDR CHANNEL      MEAN       SIGMA                         \n
///---------------------------------------------------------------------------\n
///
/// \author L. Aphecetche
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONPedestalSubprocessor)
/// \endcond

//_____________________________________________________________________________
AliMUONPedestalSubprocessor::AliMUONPedestalSubprocessor(AliMUONPreprocessor* master)
: AliMUONVSubprocessor(master,
                       "Pedestals",
                       "Upload MUON Tracker pedestals to OCDB"),
fPedestals(0x0),
fConfig(0x0),
fConfigChanged(kFALSE),
fTooFewEvents(kFALSE)
{
  /// default ctor
}

//_____________________________________________________________________________
AliMUONPedestalSubprocessor::~AliMUONPedestalSubprocessor()
{
  /// dtor
  delete fPedestals;
  delete fConfig;
}

//_____________________________________________________________________________
Bool_t
AliMUONPedestalSubprocessor::HasConfigChanged(const AliMUONVStore& newConfig) const
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
AliMUONPedestalSubprocessor::Initialize(Int_t run, UInt_t startTime, UInt_t endTime)
{
  /// When starting a new run, reads in the pedestals ASCII files.
  
  const Int_t kSystem = AliMUONPreprocessor::kDAQ;
  const char* kId = "PEDESTALS";
  
  delete fPedestals;
  fPedestals = new AliMUON2DMap(kTRUE);
  
  delete fConfig;
  fConfig = new AliMUON2DMap(kTRUE);
  
  fTooFewEvents = kFALSE;
  
  Master()->Log(Form("Reading pedestal files for Run %d startTime %ld endTime %ld",
                     run,startTime,endTime));
  
  TList* sources = Master()->GetFileSources(kSystem,kId);
  TIter next(sources);
  TObjString* o(0x0);
  Int_t n(0);
  Int_t npedFiles(0);
  
  while ( ( o = static_cast<TObjString*>(next()) ) )
  {
    TString fileName(Master()->GetFile(kSystem,kId,o->GetName()));
    Int_t ok = ReadPedestalFile(fileName.Data());
    if (ok>0)
    {
      n += ok;
      ++npedFiles;
    }
  }

  delete sources;
  
  if (!n)
  {
    Master()->Log("Failed to read any pedestals");

    delete fPedestals;
    fPedestals = 0;
    delete fConfig;
    fConfig = 0;

    // OK, we did not get our pedestals. Check if the ped run itself
    // was bad, i.e. too few events
    TString nevents(Master()->GetRunParameter("totalEvents"));
    
    if ( nevents.Atoi() < 50 ) 
    {
      Master()->Log(Form("The run had only %d events, so the failure to read pedestals is normal",nevents.Atoi()));
      // too few events, failure is normal, returns OK.
      fTooFewEvents = kTRUE;
      return kTRUE;
    }
    
    // no ped, but run looks clean, that's an error
    return kFALSE;
  }
  
  const char* kIdConf = "CONFIG";

  sources = Master()->GetFileSources(kSystem,kIdConf);
  TIter nextConf(sources);
  Int_t nconf(0);
  Int_t nconfFiles(0);
  
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
  
  if ( npedFiles != nconfFiles )
  {
    Master()->Log(Form("ERROR : Number of config files (%d) different from number of pedestal files (%d)",nconfFiles,npedFiles));
    delete fPedestals;
    fPedestals = 0;
    delete fConfig;
    fConfig = 0;
    return kFALSE;
  }
  
  fConfigChanged = HasConfigChanged(*fConfig);
  
  return kTRUE;
}

//_____________________________________________________________________________
UInt_t 
AliMUONPedestalSubprocessor::Process(TMap* /*dcsAliasMap*/)
{
  /// Store the pedestals into the CDB
  
  if (!fPedestals || !fConfig) 
  {
    if ( fTooFewEvents ) 
    {
      // ped run was too short, no reason to complain about that, it's "normal" 
      // not to have pedestals in that case.
      return 0;
    }
    else
    {
      // this is the only reason to fail for the moment : getting no pedestal or no config
      // at all.
      return 1;
    }
  }
    
  AliMUON2DStoreValidator validator;

  Master()->Log("Validating");

  TObjArray* chambers = validator.Validate(*fPedestals,fConfig);
  
  if (chambers)
  {
    // we hereby report what's missing, but this is not a reason to fail ;-)
    // the only reason to fail would be if we got no pedestal at all
    TList lines;
    lines.SetOwner(kTRUE);
  
    validator.Report(lines,*chambers);
  
    TIter next(&lines);
    TObjString* line;
  
    while ( ( line = static_cast<TObjString*>(next()) ) )
    {
      Master()->Log(line->String().Data());
    }
  }
  
  Master()->Log("Storing pedestals...");
  if ( fConfigChanged ) 
  {
    Master()->Log("...and configuration, as it has changed");
  }
  
  AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("MUON TRK");
  TString comment("Computed by AliMUONPedestalSubprocessor $Id$");
  comment.ReplaceAll("$","");
	metaData.SetComment(comment.Data());
  
  Bool_t validToInfinity = kTRUE;
	Bool_t result = Master()->Store("Calib", "Pedestals", fPedestals, &metaData, 0, validToInfinity);
  if ( fConfigChanged ) 
  {
    result = result && Master()->Store("Calib", "Config", fConfig, &metaData, 0, validToInfinity);
  }
  return ( result != kTRUE ); // return 0 if everything is ok.  
}

//_____________________________________________________________________________
Int_t
AliMUONPedestalSubprocessor::ReadPedestalFile(const char* filename)
{
  /// Read the pedestals from an ASCII file.                                  \n
  /// Format of that file is one line per channel :                           \n
  ///-------------------------------------------------------------------------\n
  /// BUS_PATCH MANU_ADDR CHANNEL      MEAN       SIGMA                       \n
  ///-------------------------------------------------------------------------\n
  ///                                                                         \n
  /// Return kFALSE if reading was not successfull.                           \n
  ///
  
  TString sFilename(gSystem->ExpandPathName(filename));
  
  Master()->Log(Form("Reading %s",sFilename.Data()));
  
  Int_t n = AliMUONTrackerIO::ReadPedestals(sFilename.Data(),*fPedestals);
  
  switch (n)
  {
    case -1:
      Master()->Log(Form("Could not open %s",sFilename.Data()));
      break;
  }
  
  return n;
}

//_____________________________________________________________________________
Int_t
AliMUONPedestalSubprocessor::ReadConfigFile(const char* filename)
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
AliMUONPedestalSubprocessor::Print(Option_t* opt) const
{
  /// ouput to screen
  if (fPedestals) fPedestals->Print("",opt);
  if (fConfig) fConfig->Print("",opt);
}

