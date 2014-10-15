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

#include "AliMUONGainSubprocessor.h"

#include "AliCDBMetaData.h"
#include "AliLog.h"
#include "AliMUON2DMap.h"
#include "AliMUON2DStoreValidator.h"
#include "AliMUONCalibParamNF.h"
#include "AliMUONConstants.h"
#include "AliMUONPreprocessor.h"
#include "AliMUONTrackerIO.h"
#include "AliMpConstants.h"
#include "AliMpDDLStore.h"
#include <Riostream.h>
#include <TList.h>
#include <TObjString.h>
#include <TObjString.h>
#include <TSystem.h>
#include <sstream>

//-----------------------------------------------------------------------------
/// \class AliMUONGainSubprocessor
///
/// Implementation of AliMUONVSubprocessor class to deal with MUON TRK Gains.
///
/// Gains are read in from an ascii file, with the format :                               
///
///---------------------------------------------------------------------------
///
///BUS_PATCH   MANU   CHANNEL    a0   a1   thres Qual
///
///---------------------------------------------------------------------------
///
/// \author L. Aphecetche
///
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONGainSubprocessor)
/// \endcond

//_____________________________________________________________________________
AliMUONGainSubprocessor::AliMUONGainSubprocessor(AliMUONPreprocessor* master)
: AliMUONVSubprocessor(master,
                       "Gains",
                       "Upload MUON Tracker Gains to OCDB"),
fGains(0x0),
fSkip(kFALSE),
fComment("")
{
  /// default ctor
}

//_____________________________________________________________________________
AliMUONGainSubprocessor::~AliMUONGainSubprocessor()
{
  /// dtor
  delete fGains;
}

//_____________________________________________________________________________
Bool_t 
AliMUONGainSubprocessor::Initialize(Int_t run, UInt_t startTime, UInt_t endTime)
{
  /// When starting a new run, reads in the Gains ASCII files.
  
  const Int_t kSystem = AliMUONPreprocessor::kDAQ;
  const char* kId = "GAINS";
  
  fComment = "";
  fSkip = kFALSE;
  
  delete fGains;
  fGains = new AliMUON2DMap(kTRUE);
  
  Master()->Log(Form("Reading Gain files for Run %d startTime %u endTime %u",
                     run,startTime,endTime));
  
  TList* sources = Master()->GetFileSources(kSystem,kId);
  TIter next(sources);
  TObjString* o(0x0);
  Int_t n(0);
  
  while ( ( o = static_cast<TObjString*>(next()) ) )
  {
    TString fileName(Master()->GetFile(kSystem,kId,o->GetName()));
    Int_t ok = ReadFile(fileName.Data());
    if (ok>0)
    {
      n += ok;
    }
    else if ( ok == AliMUONTrackerIO::kDummyFile )
    {
      // not an interesting file.
      fSkip = kTRUE;
      break;
    }
  }
  
  delete sources;

  if ( fSkip ) 
  {
    delete fGains;
    fGains = 0x0;
  }
  
  if (!n && !fSkip)
  {
    Master()->Log("Failed to read any Gains");
    delete fGains;
    fGains = 0x0;
    return kFALSE;
  }
  
  return kTRUE;
}

//_____________________________________________________________________________
UInt_t 
AliMUONGainSubprocessor::Process(TMap* /*dcsAliasMap*/)
{
  /// Store the Gains into the CDB
  
  if (!fGains) 
  {
    // this is the only reason to fail for the moment : getting no Gain
    // at all, except if the input file was a dummy one
    if ( fSkip ) 
    {
      return 0;
    }
    else
    {
      return 1;
    }
  }
    
  AliMUON2DStoreValidator validator;

  Master()->Log("Validating");

  TObjArray* chambers = validator.Validate(*fGains);
  
  if (chambers)
  {
    // we hereby report what's missing, but this is not a reason to fail ;-)
    // the only reason to fail would be if we got no Gain at all
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
  
  Master()->Log("Storing Gains");
  
  AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("MUON TRK");
	metaData.SetComment(Form("Computed by AliMUONGainSubprocessor "
                           "$Id$ ; %s",fComment.Data()));
  
  Bool_t validToInfinity = kTRUE;
	Bool_t result = Master()->Store("Calib", "Gains", fGains, &metaData, 0, validToInfinity);
  
  return ( result != kTRUE ); // return 0 if everything is ok  
}

//_____________________________________________________________________________
Int_t
AliMUONGainSubprocessor::ReadFile(const char* filename)
{
  /// Read the Gains from an ASCII file.                                  
  /// Format of that file is one line per channel :                       
  ///-------------------------------------------------------------------------
  ///BUS_PATCH   MANU   CHANNEL  a0  a1 thres Qual
  ///-------------------------------------------------------------------------
  ///                                                                         
  /// Return kFALSE if reading was not successfull.                           
  ///

  TString sFilename(gSystem->ExpandPathName(filename));
  
  Master()->Log(Form("Reading %s",sFilename.Data()));
    
  Int_t n = AliMUONTrackerIO::ReadGains(sFilename.Data(),*fGains,fComment);

  switch (n)
  {
    case AliMUONTrackerIO::kCannotOpenFile:
      Master()->Log(Form("Could not open %s",sFilename.Data()));
      break;
    case AliMUONTrackerIO::kFormatError:
      Master()->Log(Form("File %s is not of the expected format",sFilename.Data()));
      break;
    case AliMUONTrackerIO::kDummyFile:
      Master()->Log(Form("File %s is a dummy one. That's fine. I won't do anything then ;-)",sFilename.Data()));
      break;
    default:
      break;
  }
  
  return n;
}
