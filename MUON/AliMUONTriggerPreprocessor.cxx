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

#include "AliMUONTriggerPreprocessor.h"

#include "AliLog.h"
#include "AliMUONTriggerSubprocessor.h"
#include "AliShuttleInterface.h"
#include "Riostream.h"

//-----------------------------------------------------------------------------
/// \class AliMUONTriggerPreprocessor
///
/// Shuttle preprocessor for MUON trigger. The real worker
/// class is AliMUONTriggerSubprocessor
/// 
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONTriggerPreprocessor)
/// \endcond

//_____________________________________________________________________________
AliMUONTriggerPreprocessor::AliMUONTriggerPreprocessor(AliShuttleInterface* shuttle)
: AliMUONPreprocessor("MTR",shuttle),
fTriggerSubprocessor(new AliMUONTriggerSubprocessor(this))
{
  /// ctor. 
}

//_____________________________________________________________________________
AliMUONTriggerPreprocessor::~AliMUONTriggerPreprocessor()
{
  /// dtor
  delete fTriggerSubprocessor;
}

//_____________________________________________________________________________
void
AliMUONTriggerPreprocessor::Initialize(Int_t run, UInt_t startTime, UInt_t endTime)
{
  /// Re-register the subprocessor(s) depending on the actual runType

  ClearSubprocessors();
  
  fIsValid = kTRUE;
  
  TString runType = GetRunType();
  
  if ( runType == "PHYSICS" ||
       runType == "ELECTRONICS_CALIBRATION_RUN" ||
       runType == "DETECTOR_CALIBRATION_RUN" ) 
  {
    Add(fTriggerSubprocessor);
  }
  else
  {
    Log(Form("ERROR-Unknown RunType=%",runType.Data()));
    fIsValid = kFALSE;
  }
  
  AliMUONPreprocessor::Initialize(run,startTime,endTime);  
}
