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

#include "AliMUONPreprocessor.h"

#include "AliMUONPedestalSubprocessor.h"
#include "AliMUONHVSubprocessor.h"
#include "AliMUONGMSSubprocessor.h"

#include "AliLog.h"
#include "AliShuttleInterface.h"
#include "Riostream.h"
#include "TObjArray.h"

/// \class AliMUONPreprocessor
///
/// Shuttle preprocessor for MUON subsystems (TRK and TRG)
/// 
/// It's simply a manager class that deals with a list of sub-tasks 
/// (of type AliMUONVSubprocessor).
///
/// \author Laurent Aphecetche

/// \cond CLASSIMP
ClassImp(AliMUONPreprocessor)
/// \endcond

const TString  AliMUONPreprocessor::fgkTrackerDetName = "MCH";
const TString  AliMUONPreprocessor::fgkTriggerDetName = "MTR";

//_____________________________________________________________________________
AliMUONPreprocessor::AliMUONPreprocessor(const TString& detName, 
                                         AliShuttleInterface* shuttle) 
: AliPreprocessor(detName.Data(),shuttle), 
//  fSubprocessors(new TObjArray[kLast])
  fSubprocessors(new TObjArray())
{
  /// ctor. Builds the list of subtasks
  /// Tests detector wrt to tracker or trigger to 
  /// instantiate the correct list of subtasks, which should be : 
  /// Tracker : 
  /// - pedestals
  /// - gains
  /// - deadchannels
  /// - gms
  ///
  /// Trigger : 
  /// - masks
  /// - lut
  
  if ( detName == fgkTrackerDetName ) 
  { 
    TString runType = shuttle->GetRunParameter("RunType");
    if ( runType == "PEDESTAL_RUN" ) // FIXME : check the name 
    {
      fSubprocessors->Add(new AliMUONPedestalSubprocessor(this)); // to be called only for pedestal runs
      Log("INFO-Will run Pedestal subprocessor");
    }
    else if ( runType == "ELECTRONICS_CALIBRATION_RUN" ) // FIXME : check the name 
    {  
      Log("WARNING-Subprocessor for gains not yet implemented");
      //fSubprocessors->Add(new AliMUONGainSubprocessor(this)); // to be called only for gain runs
    }
    else if ( runType == "GMS" ) // FIXME : check the name 
    {
      fSubprocessors->Add(new AliMUONGMSSubprocessor(this));
      Log("INFO-Will run GMS subprocessor");
    }
    else if ( runType == "PHYSICS" ) // FIXME : check the name
    {
      fSubprocessors->Add(new AliMUONHVSubprocessor(this)); // to be called only for physics runs
      Log("INFO-Will run HV subprocessor");
    }
    else
    {
      Log(Form("ERROR-Unknown RunType=%",runType.Data()));
    }
  }
  else if ( detName == fgkTriggerDetName ) 
  {
    Log("WARNING-Trigger subprocessors not yet implemented.");
  }  
  else 
  { 
    // Wrong detector name
    Log("ERROR-Wrong detector name.");
  }  
}

//_____________________________________________________________________________
AliMUONPreprocessor::~AliMUONPreprocessor()
{
  /// dtor
  delete fSubprocessors;
}

//_____________________________________________________________________________
void
AliMUONPreprocessor::Initialize(Int_t run, UInt_t startTime, UInt_t endTime)
{
  /// loop over subtasks and initialize them
  for ( Int_t i = 0; i <= fSubprocessors->GetLast(); ++i )
  {
    Subprocessor(i)->Initialize(run,startTime,endTime);
  }
}

//_____________________________________________________________________________
UInt_t
AliMUONPreprocessor::Process(TMap* dcsAliasMap)
{
  /// loop over subtasks to make them work
  UInt_t rv(0);
  
  for ( Int_t i = 0; i <= fSubprocessors->GetLast(); ++i )
  {
    rv += Subprocessor(i)->Process(dcsAliasMap);
  }
  return rv;
}

//_____________________________________________________________________________
void
AliMUONPreprocessor::Print(Option_t* opt) const
{
  /// output to screen
  cout << "<AliMUONPreprocessor> subprocessors :" << endl;
  for ( Int_t i=0; i <= fSubprocessors->GetLast(); ++i )
  {
    Subprocessor(i)->Print(opt);
  }
}

//_____________________________________________________________________________
AliMUONVSubprocessor*
AliMUONPreprocessor::Subprocessor(Int_t i) const
{
  /// return i-th subprocessor
  if ( i >= 0 && i <= fSubprocessors->GetLast() )
  {
    return static_cast<AliMUONVSubprocessor*>(fSubprocessors->At(i));
  }
  return 0x0;
}
