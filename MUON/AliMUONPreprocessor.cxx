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

#include "AliMpSegmentation.h"
#include "AliMpDDLStore.h"

#include "AliLog.h"
#include "AliShuttleInterface.h"
#include "Riostream.h"
#include "TObjArray.h"

//-----------------------------------------------------------------------------
/// \class AliMUONPreprocessor
///
/// Shuttle preprocessor for MUON subsystems (TRK and TRG)
/// 
/// It's simply a manager class that deals with a list of sub-tasks 
/// (of type AliMUONVSubprocessor).
///
/// \author Laurent Aphecetche
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONPreprocessor)
/// \endcond

//_____________________________________________________________________________
AliMUONPreprocessor::AliMUONPreprocessor(const char* detName, AliShuttleInterface* shuttle)
: AliPreprocessor(detName, shuttle),
  fSubprocessors(new TObjArray()),
  fProcessDCS(kFALSE),
  fIsValid(kFALSE)
{
  /// ctor
}

//_____________________________________________________________________________
AliMUONPreprocessor::~AliMUONPreprocessor()
{
  /// dtor
  delete fSubprocessors;
}

//_____________________________________________________________________________
void
AliMUONPreprocessor::ClearSubprocessors()
{
  /// Empty our subprocessor list
  fSubprocessors->Clear();
  fProcessDCS = kFALSE;
  fIsValid = kFALSE;
}

//_____________________________________________________________________________
void
AliMUONPreprocessor::Add(AliMUONVSubprocessor* sub, Bool_t processDCS)
{
  /// Add a subprocessor to our list of workers
  fSubprocessors->Add(sub);
  if ( processDCS == kTRUE ) fProcessDCS = processDCS;
}

//_____________________________________________________________________________
void
AliMUONPreprocessor::Initialize(Int_t run, UInt_t startTime, UInt_t endTime)
{
  /// Load mapping and initialize subtasks  

  // Delete previous mapping
  delete AliMpSegmentation::Instance(false);
  delete AliMpDDLStore::Instance(false);
  
  // Load mapping from CDB for this run
  
  AliCDBEntry* cdbEntry = GetFromOCDB("Calib", "Mapping");
  if (!cdbEntry)
  {
    Log("Could not get Mapping from OCDB !");
    fIsValid = kFALSE;
  }
  
  cdbEntry = GetFromOCDB("Calib", "DDLStore");
  if (!cdbEntry)
  {
    Log("Could not get DDLStore from OCDB");
    fIsValid = kFALSE;
  }

  if (IsValid())
  {
    // loop over subtasks and initialize them
    for ( Int_t i = 0; i <= fSubprocessors->GetLast(); ++i )
    {
      Subprocessor(i)->Initialize(run,startTime,endTime);
    }
  }

  Log(Form("Initialize was %s",( IsValid() ? "fine" : "NOT OK")));
}

//_____________________________________________________________________________
UInt_t
AliMUONPreprocessor::Process(TMap* dcsAliasMap)
{
  /// loop over subtasks to make them work
  
  if (!IsValid())
  {
    Log("Will not run as not properly initialized");
    return 99;
  }
  
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
