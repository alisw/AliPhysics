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
fPedestals(0x0)
{
  /// default ctor
}

//_____________________________________________________________________________
AliMUONPedestalSubprocessor::~AliMUONPedestalSubprocessor()
{
  /// dtor
  delete fPedestals;
}

//_____________________________________________________________________________
void 
AliMUONPedestalSubprocessor::Initialize(Int_t run, UInt_t startTime, UInt_t endTime)
{
  /// When starting a new run, reads in the pedestals ASCII files.
  
  const Int_t kSystem = AliMUONPreprocessor::kDAQ;
  const char* kId = "PEDESTALS";
  
  delete fPedestals;
  fPedestals = new AliMUON2DMap(kTRUE);
  
  Master()->Log(Form("Reading pedestal files for Run %d startTime %ld endTime %ld",
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
  }
  
  if (!n)
  {
    Master()->Log("Failed to read any pedestals");
    delete fPedestals;
    fPedestals = 0;
  }
  delete sources;
}

//_____________________________________________________________________________
UInt_t 
AliMUONPedestalSubprocessor::Process(TMap* /*dcsAliasMap*/)
{
  /// Store the pedestals into the CDB
  
  if (!fPedestals) 
  {
    // this is the only reason to fail for the moment : getting no pedestal
    // at all.
    return 1;
  }
    
  AliMUON2DStoreValidator validator;

  Master()->Log("Validating");

  TObjArray* chambers = validator.Validate(*fPedestals);
  
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
  
  Master()->Log("Storing pedestals");
  
  AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("MUON TRK");
  TString comment("Computed by AliMUONPedestalSubprocessor $Id$");
  comment.ReplaceAll("$","");
	metaData.SetComment(comment.Data());
  
  Bool_t validToInfinity = kTRUE;
	Bool_t result = Master()->Store("Calib", "Pedestals", fPedestals, &metaData, 0, validToInfinity);
  
  return ( result != kTRUE ); // return 0 if everything is ok.  
}

//_____________________________________________________________________________
Int_t
AliMUONPedestalSubprocessor::ReadFile(const char* filename)
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
void
AliMUONPedestalSubprocessor::Print(Option_t* opt) const
{
  /// ouput to screen
  if (fPedestals) fPedestals->Print("",opt);
}
