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
#include "AliMUONLVSubprocessor.h"
#include "AliMUONGMSSubprocessor.h"
#include "AliMUONOccupancySubprocessor.h"
#include "AliLog.h"
#include "AliMUONBusPatchEvolutionSubprocessor.h"
#include "AliShuttleInterface.h"
#include "Riostream.h"
#include "TObjArray.h"
#include "AliMUONConfigSubprocessor.h"

//-----------------------------------------------------------------------------
/// \class AliMUONTrackerPreprocessor
///
/// Shuttle preprocessor for MUON tracker
///
/// It's simply a manager class that deals with a list of sub-tasks
/// (of type AliMUONVSubprocessor).
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONTrackerPreprocessor)
/// \endcond

//_____________________________________________________________________________
AliMUONTrackerPreprocessor::AliMUONTrackerPreprocessor(AliShuttleInterface* shuttle)
: AliMUONPreprocessor("MCH",shuttle),
fPedestalSubprocessor(new AliMUONPedestalSubprocessor(this)),
fGMSSubprocessor(new AliMUONGMSSubprocessor(this)),
fHVSubprocessor(new AliMUONHVSubprocessor(this,kTRUE)),
fOccupancySubprocessor(new AliMUONOccupancySubprocessor(this)),
fBusPatchEvolutionSubprocessor(new AliMUONBusPatchEvolutionSubprocessor(this)),
fConfigSubprocessor(new AliMUONConfigSubprocessor(this)),
fLVSubprocessor(new AliMUONLVSubprocessor(this))
{
  /// ctor.

    AddRunType("PEDESTAL");
    AddRunType("CALIBRATION");
    AddRunType("GMS");
    AddRunType("PHYSICS");
}

//_____________________________________________________________________________
AliMUONTrackerPreprocessor::~AliMUONTrackerPreprocessor()
{
  /// dtor

  delete fPedestalSubprocessor;
  delete fGMSSubprocessor;
  delete fHVSubprocessor;
  delete fOccupancySubprocessor;
  delete fBusPatchEvolutionSubprocessor;
  delete fConfigSubprocessor;
  delete fLVSubprocessor;
}

//_____________________________________________________________________________
void
AliMUONTrackerPreprocessor::Initialize(Int_t run, UInt_t startTime, UInt_t endTime)
{
  /// Re-register the subprocessor(s) depending on the actual runType

  ClearSubprocessors();

  TString runType = GetRunType();

  fIsValid = kTRUE;
  fIsApplicable = kTRUE;

  if ( runType == "PEDESTAL" )
  {
    Add(fPedestalSubprocessor); // to be called only for pedestal runs
    Log("INFO-Will run Pedestal subprocessor");
  }
  else if ( runType == "GMS" )
  {
    Add(fGMSSubprocessor);
    Log("INFO-Will run GMS subprocessor");
  }
  else if ( runType == "PHYSICS" )
  {
    Bool_t useDCS(kTRUE);

    Add(fHVSubprocessor,useDCS); // to be called only for physics runs
    Add(fLVSubprocessor,useDCS);
    Add(fOccupancySubprocessor);
    Add(fBusPatchEvolutionSubprocessor);
    Add(fConfigSubprocessor);

    Log("INFO-Will run LV subprocessor");
    Log("INFO-Will run HV subprocessor");
    if ( static_cast<AliMUONHVSubprocessor*>(fHVSubprocessor)->IncludeHVCurrent() )
    {
      Log("INFO-HV subprocessor will store HV currents in addition to the voltages");
    }
    Log("INFO-Will run Occupancy subprocessor");
    Log("INFO-Will run Bus Patch Evolution subprocessor");
    Log("INFO-Will run Config subprocessor");

  }
  else
  {
    fIsApplicable = kFALSE;
  }

  AliMUONPreprocessor::Initialize(run,startTime,endTime);

}
