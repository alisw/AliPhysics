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
#include "AliMUONTriggerDCSSubprocessor.h"
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
  fTriggerSubprocessor(new AliMUONTriggerSubprocessor(this)),
  fTriggerDCSSubprocessor(new AliMUONTriggerDCSSubprocessor(this))
{
  /// ctor. 
  AddRunType("PHYSICS");
  AddRunType("CALIBRATION");
}

//_____________________________________________________________________________
AliMUONTriggerPreprocessor::~AliMUONTriggerPreprocessor()
{
  /// dtor
  delete fTriggerSubprocessor;
  delete fTriggerDCSSubprocessor;
}

//_____________________________________________________________________________
void
AliMUONTriggerPreprocessor::Initialize(Int_t run, UInt_t startTime, UInt_t endTime)
{
  /// Re-register the subprocessor(s) depending on the actual runType

  ClearSubprocessors();
  
  fIsValid = kTRUE;
  fIsApplicable = kTRUE;
  
  TString runType = GetRunType();
  
  if ( runType == "PHYSICS" || runType == "CALIBRATION" )
  {
    Add(fTriggerSubprocessor);
    Add(fTriggerDCSSubprocessor,kTRUE); // uses DCS
  }
  else
  {
    fIsApplicable = kFALSE;
  }
  
  AliMUONPreprocessor::Initialize(run,startTime,endTime);  
}
