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
#include "AliMUONCalibParamNI.h"
#include "AliMUONPreprocessor.h"
#include "AliMUONTriggerIO.h"
#include "AliMUONTriggerLut.h"
#include "AliMpConstants.h"
#include <Riostream.h>
#include <TList.h>
#include <TObjString.h>
#include <TSystem.h>

/// \class AliMUONTriggerSubprocessor
///
/// Implementation of AliMUONVSubprocessor for MUON trigger system
///
/// Reads masks and LUT online files to feed the OCDB
///
/// \author L. Aphecetche

/// \cond CLASSIMP
ClassImp(AliMUONTriggerSubprocessor)
/// \endcond

//_____________________________________________________________________________
AliMUONTriggerSubprocessor::AliMUONTriggerSubprocessor(AliMUONPreprocessor* master)
: AliMUONVSubprocessor(master,
                       "Triggers",
                       "Upload MUON Trigger masks and LUT to OCDB"),
fRegionalMasks(0x0),
fLocalMasks(0x0),
fGlobalMasks(0x0),
fLUT(0x0)
{
  /// default ctor
}

//_____________________________________________________________________________
AliMUONTriggerSubprocessor::~AliMUONTriggerSubprocessor()
{
  /// dtor
  delete fRegionalMasks;
  delete fLocalMasks;
  delete fGlobalMasks;
  delete fLUT;
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
void 
AliMUONTriggerSubprocessor::Initialize(Int_t run, UInt_t startTime, UInt_t endTime)
{
  /// When starting a new run, reads in the trigger online files.
  
  Master()->Log(Form("Reading trigger masks for Run %d startTime %ld endTime %ld",
                     run,startTime,endTime));
    
  delete fRegionalMasks;
  delete fLocalMasks;
  delete fGlobalMasks;
  
  fRegionalMasks = new AliMUON1DArray(16);
  fLocalMasks = new AliMUON1DArray(AliMpConstants::NofLocalBoards()+1);
  fGlobalMasks = 0x0; // new AliMUONCalibParamNI(1,16,1,0,0);

  AliMUONTriggerIO tio;
  
  Bool_t ok = tio.ReadMasks(GetFileName("LOCAL").Data(),
                            GetFileName("REGIONAL").Data(),
                            GetFileName("GLOBAL").Data(),
                            fLocalMasks,fRegionalMasks,fGlobalMasks);
  
  if (!ok)
  {
    Master()->Log("ERROR : ReadMasks failed");
    delete fLocalMasks;
    delete fRegionalMasks;
    delete fGlobalMasks;
    fLocalMasks = 0x0;
    fRegionalMasks = 0x0;
    fGlobalMasks = 0x0;
  }

  delete fLUT;
  fLUT = new AliMUONTriggerLut;
    
  Master()->Log(Form("Reading trigger LUT for Run %d startTime %ld endTime %ld",
                     run,startTime,endTime));
  
  tio.ReadLUT(GetFileName("LUT").Data(),*fLUT);

  if (!ok)
  {
    Master()->Log("ERROR : ReadLUT failed");
    delete fLUT;
    fLUT = 0x0;
  }
}

//_____________________________________________________________________________
UInt_t 
AliMUONTriggerSubprocessor::Process(TMap* /*dcsAliasMap*/)
{
  /// Store the trigger masks into the CDB
  
  if ( !fGlobalMasks && !fRegionalMasks && !fLocalMasks && !fLUT )
  {
    // nothing to do
    return 0;
  }
  
  Master()->Log(Form("N global = %d N regional = %d N local %d",                     
                     (fGlobalMasks ? 1 : 0 ),
                     (fRegionalMasks ? fRegionalMasks->GetSize() : 0 ),
                     (fLocalMasks ? fLocalMasks->GetSize() : 0 )));
  
  AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("MUON TRG");
	metaData.SetComment("Computed by AliMUONTriggerSubprocessor $Id$");
  
  Bool_t validToInfinity = kTRUE;

	Bool_t result1(kTRUE);
  Bool_t result2(kTRUE);
  Bool_t result3(kTRUE);
  Bool_t result4(kTRUE);
  
  if ( fGlobalMasks ) 
  {
    result1 = Master()->Store("Calib", "GlobalTriggerBoardMasks", fGlobalMasks, 
                              &metaData, 0, validToInfinity);
  }
  
  if ( fRegionalMasks && fRegionalMasks->GetSize() > 0 )
  {
    result2 = Master()->Store("Calib", "RegionalTriggerBoardMasks", fRegionalMasks, 
                              &metaData, 0, validToInfinity);
  }
  
  if ( fLocalMasks && fLocalMasks->GetSize() > 0 ) 
  {
    result3 = Master()->Store("Calib", "LocalTriggerBoardMasks", fLocalMasks, 
                              &metaData, 0, validToInfinity);
  }

  if ( fLUT )
  {
    result4 = Master()->Store("Calib", "TriggerLut", fLUT, 
                              &metaData, 0, validToInfinity);
  }
  
  return ( result1 != kTRUE && result2 != kTRUE && result3 != kTRUE && result4 != kTRUE ); // return 0 if everything is ok.  
}
