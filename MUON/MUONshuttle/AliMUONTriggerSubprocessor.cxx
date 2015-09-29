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

#include "AliMUONTriggerSubprocessor.h"

#include "AliCDBMetaData.h"
#include "AliLog.h"
#include "AliMUON1DArray.h"
#include "AliMUONGlobalCrateConfig.h"
#include "AliMUONRegionalTriggerConfig.h"
#include "AliMUONCalibParamNI.h"
#include "AliMUONPreprocessor.h"
#include "AliMUONTriggerIO.h"
#include "AliMUONTriggerLut.h"
#include "AliMpConstants.h"
#include <Riostream.h>
#include <TList.h>
#include <TObjString.h>
#include <TSystem.h>
#include <TClonesArray.h>

/// \class AliMUONTriggerSubprocessor
///
/// Implementation of AliMUONVSubprocessor for MUON trigger system
///
/// Reads masks and LUT online files to feed the OCDB
///
/// \author L. Aphecetche

using std::ifstream;
/// \cond CLASSIMP
ClassImp(AliMUONTriggerSubprocessor)
/// \endcond

//_____________________________________________________________________________
AliMUONTriggerSubprocessor::AliMUONTriggerSubprocessor(AliMUONPreprocessor* master)
: AliMUONVSubprocessor(master,
                       "Triggers",
                       "Upload MUON Trigger masks and LUT to OCDB"),
fRegionalConfig(0x0),
fLocalMasks(0x0),
fGlobalConfig(0x0),
fLUT(0x0),
fTrigScalers(0x0),
fRegionalConfigToOCDB(kFALSE),
fLocalMasksToOCDB(kFALSE),
fGlobalConfigToOCDB(kFALSE),
fLUTToOCDB(kFALSE),
fTrigScalersToOCDB(kFALSE)
{
  /// default ctor
}

//_____________________________________________________________________________
AliMUONTriggerSubprocessor::~AliMUONTriggerSubprocessor()
{
  /// dtor
  delete fRegionalConfig;
  delete fLocalMasks;
  delete fGlobalConfig;
  delete fLUT;
  delete [] fTrigScalers;
}

//_____________________________________________________________________________
TString
AliMUONTriggerSubprocessor::GetFileName(const char* fid) const
{
  /// Get filename for a given id
  
  const Int_t kSystem = AliMUONPreprocessor::kDAQ;
  
  TList* sources = Master()->GetFileSources(kSystem,fid);
  if ( sources && sources->GetSize() == 1 ) 
  {
    return Master()->GetFile(kSystem,fid,static_cast<TObjString*>(sources->First())->GetName());
  }
  return "";
}

//_____________________________________________________________________________
Bool_t 
AliMUONTriggerSubprocessor::Initialize(Int_t run, UInt_t startTime, UInt_t endTime)
{
  /// When starting a new run, reads in the trigger online files.
  
  // First thing to do (after cleanup, that is), is to check whether the DA
  // was alive or not. If it was, it should have put the "MtgCurrent.dat" file
  // on the FXS. 
  //
  
  Bool_t da(kTRUE);
  
  TString exportedFiles(gSystem->ExpandPathName(GetFileName("EXPORTED").Data()));
  
  if ( exportedFiles == "" ) 
  {
    Master()->Log("exported files not specified !");
    da = kFALSE;
  }
  
  if ( gSystem->AccessPathName(exportedFiles.Data(),kFileExists) ) // mind the strange return value convention of that method !
  {
    Master()->Log(Form("%s is not there !",exportedFiles.Data()));
    da = kFALSE;
  }
  
  if (!da)
  {
    Master()->Log("FATAL ERROR : DA does not seem to have been run !!!");
    Master()->Invalidate();
    return kFALSE;
  }
  
  // OK. We have an exportedFiles.dat file at hand.
  // Read it to know what other files should be there.
  
  Bool_t regionalFile(kFALSE);
  Bool_t globalFile(kFALSE);
  Bool_t localFile(kFALSE);
  Bool_t lutFile(kFALSE);
  Bool_t trigScalFile(kFALSE);
  
  WhichFilesToRead(GetFileName("EXPORTED").Data(),
                   globalFile,regionalFile,localFile,lutFile,trigScalFile);
  
  if ((globalFile+regionalFile+localFile+lutFile+trigScalFile) == 0) {
    Master()->Log("No file(s) to be processed for this run. Exiting.");
    return kTRUE;
  }
  
  delete fRegionalConfig; fRegionalConfig = 0x0;
  delete fLocalMasks; fLocalMasks = 0x0;
  delete fGlobalConfig; fGlobalConfig = 0x0;
  delete fLUT; fLUT = 0x0;
  delete fTrigScalers; fTrigScalers = 0x0;
  
  Master()->Log(Form("Reading trigger masks for Run %d startTime %u endTime %u",
                     run,startTime,endTime));
    
  Int_t check = 
    TestFile("REGIONAL",regionalFile) + 
    TestFile("LOCAL",localFile) + 
    TestFile("GLOBAL",globalFile) +
    TestFile("LUT",lutFile) +
    TestFile("TRIGSCAL",trigScalFile);

  if ( check ) 
  {
    Master()->Log("Could not read some input file(s). Exiting");
    Master()->Invalidate();
    return kFALSE;
  }
  
  if ( regionalFile ) globalFile = kTRUE;

  if ( regionalFile ) fRegionalConfig = new AliMUONRegionalTriggerConfig();
  if ( localFile ) fLocalMasks = new AliMUON1DArray(AliMpConstants::TotalNofLocalBoards()+1);
  if ( globalFile )   fGlobalConfig   = new AliMUONGlobalCrateConfig();

  AliMUONTriggerIO tio;
  
  Bool_t ok = tio.ReadConfig(localFile ? GetFileName("LOCAL").Data() : "",
                             regionalFile ? GetFileName("REGIONAL").Data() : "",
                             globalFile ? GetFileName("GLOBAL").Data() : "",
                             fLocalMasks,fRegionalConfig,fGlobalConfig);
  
  if (!ok)
  {
    Master()->Log("ERROR : ReadConfig failed");
    delete fLocalMasks;
    delete fRegionalConfig;
    delete fGlobalConfig;
    fLocalMasks = 0x0;
    fRegionalConfig = 0x0;
    fGlobalConfig = 0x0;
  }

  if ( lutFile ) 
  {
    fLUT = new AliMUONTriggerLut;
    
    Master()->Log(Form("Reading trigger LUT for Run %d startTime %u endTime %u",
                       run,startTime,endTime));
  
    ok = tio.ReadLUT(GetFileName("LUT").Data(),*fLUT);

    if (!ok)
    {
      Master()->Log("ERROR : ReadLUT failed");
      delete fLUT;
      fLUT = 0x0;
    }
  }

  if ( trigScalFile ) {
    fTrigScalers = new TClonesArray("AliMUONTriggerScalers",10);
    ok = tio.ReadTrigScalers(GetFileName("TRIGSCAL").Data(),*fTrigScalers);
    if (!ok)
    {
      Master()->Log("ERROR : ReadTrigScalers failed");
      delete fTrigScalers;
      fTrigScalers = 0x0;
    }
  }

  return kTRUE;

}

//_____________________________________________________________________________
UInt_t 
AliMUONTriggerSubprocessor::Process(TMap* /*dcsAliasMap*/)
{
  /// Store the trigger masks into the CDB
  
  if ( !fGlobalConfig && !fRegionalConfig && !fLocalMasks && !fLUT )
  {
    // nothing to do
    return 0;
  }
  
  Master()->Log(Form("N global = %d N regional = %d N local %d N lut %d",                     
                     (fGlobalConfig ? 1 : 0 ),
                     (fRegionalConfig ? fRegionalConfig->GetNofTriggerCrates() : 0 ),
                     (fLocalMasks ? fLocalMasks->GetSize() : 0 ),
                     (fLUT ? 1 : 0)));
  
  AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("MUON TRG");
	metaData.SetComment("Computed by AliMUONTriggerSubprocessor $Id$");
  
  Bool_t validToInfinity = kFALSE;

  Bool_t result1(kTRUE);
  Bool_t result2(kTRUE);
  Bool_t result3(kTRUE);
  Bool_t result4(kTRUE);
  Bool_t result5(kTRUE);  

  if ( fGlobalConfig && fGlobalConfigToOCDB ) 
  {
    result1 = Master()->Store("Calib", "GlobalTriggerCrateConfig", fGlobalConfig, 
                              &metaData, 0, validToInfinity);
  }
  
  if ( fRegionalConfig && fRegionalConfig->GetNofTriggerCrates() > 0 && fRegionalConfigToOCDB )
  {
    result2 = Master()->Store("Calib", "RegionalTriggerConfig", fRegionalConfig, 
                              &metaData, 0, validToInfinity);
  }
  
  if ( fLocalMasks && fLocalMasks->GetSize() > 0 && fLocalMasksToOCDB ) 
  {
    result3 = Master()->Store("Calib", "LocalTriggerBoardMasks", fLocalMasks, 
                              &metaData, 0, validToInfinity);
  }

  if ( fLUT && fLUTToOCDB )
  {
    result4 = Master()->Store("Calib", "TriggerLut", fLUT, 
                              &metaData, 0, validToInfinity);
  }

  if ( fTrigScalers && fTrigScalersToOCDB )
  {
    result5 = Master()->Store("Calib","TriggerScalers", fTrigScalers,
			      &metaData, 0, validToInfinity);
  }

  return ( result1 != kTRUE || result2 != kTRUE || result3 != kTRUE || result4 != kTRUE || result5 != kTRUE); // return 0 if everything is ok.  
}

//_____________________________________________________________________________
Int_t
AliMUONTriggerSubprocessor::TestFile(const char* baseName, Bool_t shouldBeThere) const
{
  /// Check if required file can be accessed 
  if (!shouldBeThere) return 0; // all is ok
  
  TString fileName(gSystem->ExpandPathName(GetFileName(baseName).Data()));
  
  if ( gSystem->AccessPathName(fileName.Data(),kFileExists) ) 
  {
    // mind the strange return value convention of that method !
    Master()->Log(Form("File %s does not exist",fileName.Data()));    
    return 1; // this is an error
  }
  return 0; // all is ok.
}

//_____________________________________________________________________________
void
AliMUONTriggerSubprocessor::WhichFilesToRead(const char* exportedFiles,
                                             Bool_t& globalFile,
                                             Bool_t& regionalFile,
                                             Bool_t& localFile,
                                             Bool_t& lutFile,
					     Bool_t& trigScalFile)
{
  /// From the exportedFiles file, determine which other files will need
  /// to be read in
  ifstream in(gSystem->ExpandPathName(exportedFiles));
  
  globalFile = kFALSE;
  regionalFile = kFALSE;
  localFile = kFALSE;
  lutFile = kFALSE;
  trigScalFile = kFALSE;
  
  char line[1024];
  TObjArray *oa;
  TObjString *os;
  
  while ( in.getline(line,1024,'\n') )
  {
    TString sline(line);
    sline.ToUpper();
    
    if ( sline.Contains("REGIONAL") ) {
      regionalFile = kTRUE;
      oa = sline.Tokenize(" ");
      if (oa->GetLast() >= 1) {
	// new format, use the integer flag after the file name
	os = (TObjString*)oa->At(1);
	fRegionalConfigToOCDB |= (os->GetString()).Atoi();
      } else {
	// old format, transfer anyway
	fRegionalConfigToOCDB |= 1;
      }
    }
    if ( sline.Contains("LUT") ) {
      lutFile = kTRUE;
      oa = sline.Tokenize(" ");
      if (oa->GetLast() >= 1) {
	// new format, use the integer flag after the file name
	os = (TObjString*)oa->At(1);
	fLUTToOCDB |= (os->GetString()).Atoi();
      } else {
	// old format, transfer anyway
	fLUTToOCDB |= 1;
      }
    }
    if ( sline.Contains("LOCAL") && sline.Contains("MASK") ) {
      localFile = kTRUE;
      oa = sline.Tokenize(" ");
      if (oa->GetLast() >= 1) {
	// new format, use the integer flag after the file name
	os = (TObjString*)oa->At(1);
	fLocalMasksToOCDB |= (os->GetString()).Atoi();
      } else {
	// old format, transfer anyway
	fLocalMasksToOCDB |= 1;
      }
    }
    if ( sline.Contains("GLOBAL") ) {
      globalFile = kTRUE;
      oa = sline.Tokenize(" ");
      if (oa->GetLast() >= 1) {
	// new format, use the integer flag after the file name
	os = (TObjString*)oa->At(1);
	fGlobalConfigToOCDB |= (os->GetString()).Atoi();
      } else {
	// old format, transfer anyway
	fGlobalConfigToOCDB |= 1;
      }
    }
    if ( sline.Contains("TRIGSCAL") ) {
      trigScalFile = kTRUE;
      oa = sline.Tokenize(" ");
      if (oa->GetLast() >= 1) {
	// new format, use the integer flag after the file name
	os = (TObjString*)oa->At(1);
	fTrigScalersToOCDB |= (os->GetString()).Atoi();
      } else {
	// old format, transfer anyway
	fTrigScalersToOCDB |= 1;
      }
    }

  }
    
  if ( regionalFile || localFile || globalFile || lutFile || trigScalFile) 
  {
    Master()->Log(Form("Will have to read the following file types:"));
    if ( globalFile) Master()->Log(Form("GLOBAL (%d)",fGlobalConfigToOCDB));
    if ( regionalFile ) Master()->Log(Form("REGIONAL (%d)",fRegionalConfigToOCDB));
    if ( localFile) Master()->Log(Form("LOCAL (%d)",fLocalMasksToOCDB));
    if ( lutFile ) Master()->Log(Form("LUT (%d)",fLUTToOCDB));
    if ( trigScalFile ) Master()->Log(Form("TRIGSCAL (%d)",fTrigScalersToOCDB));
  }
}

