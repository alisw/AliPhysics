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

#include "AliMUONTrackerPreprocessor.h"

#include "AliMUONPedestalSubprocessor.h"
#include "AliMUONHVSubprocessor.h"
#include "AliMUONGMSSubprocessor.h"

#include "AliLog.h"
#include "AliShuttleInterface.h"
#include "Riostream.h"
#include "TObjArray.h"

/// \class AliMUONTrackerPreprocessor
///
/// Shuttle preprocessor for MUON tracker
/// 
/// It's simply a manager class that deals with a list of sub-tasks 
/// (of type AliMUONVSubprocessor).
///
/// \author Laurent Aphecetche

/// \cond CLASSIMP
ClassImp(AliMUONTrackerPreprocessor)
/// \endcond

//_____________________________________________________________________________
AliMUONTrackerPreprocessor::AliMUONTrackerPreprocessor(AliShuttleInterface* shuttle)
: AliMUONPreprocessor("MCH",shuttle)
{
  /// ctor. 
}

//_____________________________________________________________________________
AliMUONTrackerPreprocessor::~AliMUONTrackerPreprocessor()
{
  /// dtor
}

//_____________________________________________________________________________
void
AliMUONTrackerPreprocessor::Initialize(Int_t run, UInt_t startTime, UInt_t endTime)
{
  DeleteSubprocessors();
  
  TString runType = GetRunType();
  if ( runType == "PEDESTAL_RUN" ) // FIXME : check the name
  {
    Add(new AliMUONPedestalSubprocessor(this)); // to be called only for pedestal runs
    Log("INFO-Will run Pedestal subprocessor");
  }
  else if ( runType == "ELECTRONICS_CALIBRATION_RUN" ) // FIXME : check the name
  {
    Log("WARNING-Subprocessor for gains not yet implemented");
    //fSubprocessors->Add(new AliMUONGainSubprocessor(this)); // to be called only for gain runs
  }
  else if ( runType == "GMS" ) // FIXME : check the name
  {
    Add(new AliMUONGMSSubprocessor(this));
    Log("INFO-Will run GMS subprocessor");
  }
  else if ( runType == "PHYSICS" ) // FIXME : check the name
  {
    Add(new AliMUONHVSubprocessor(this)); // to be called only for physics runs
    Log("INFO-Will run HV subprocessor");
  }
  else
  {
    Log(Form("ERROR-Unknown RunType=%",runType.Data()));
  }
  
  AliMUONPreprocessor::Initialize(run,startTime,endTime);
  
}
