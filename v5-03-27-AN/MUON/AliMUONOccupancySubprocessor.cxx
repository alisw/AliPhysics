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

#include "AliMUONOccupancySubprocessor.h"

#include "AliCDBMetaData.h"
#include "AliMUON2DMap.h"
#include "AliMUONPreprocessor.h"
#include "AliMUONTrackerIO.h"
#include "TObjString.h"
#include "TSystem.h"

//-----------------------------------------------------------------------------
/// \class AliMUONOccupancySubprocessor
///
/// Implementation of AliMUONVSubprocessor class to deal with MUON TRK occupancy.
///
/// Values to compute the occupancy are read in from an ascii file,
/// with the format :               \n
///---------------------------------------------------------------------------\n
/// BUS_PATCH MANU_ADDR SUM_N NEVENTS
///---------------------------------------------------------------------------\n
///
/// \author L. Aphecetche
//-----------------------------------------------------------------------------

/// \cond CLASSIMP
ClassImp(AliMUONOccupancySubprocessor)
/// \endcond

//_____________________________________________________________________________
AliMUONOccupancySubprocessor::AliMUONOccupancySubprocessor(AliMUONPreprocessor* master)
: AliMUONVSubprocessor(master,"Occupancy","Upload MUON Tracker occupancy to OCDB"),
fOccupancyMap(0x0)
{
  /// Default ctor
}

//_____________________________________________________________________________
AliMUONOccupancySubprocessor::~AliMUONOccupancySubprocessor()
{
  /// dtor
  delete fOccupancyMap;
}

//_____________________________________________________________________________
Bool_t 
AliMUONOccupancySubprocessor::Initialize(Int_t run, UInt_t startTime, UInt_t endTime)
{
  /// When starting a new run, reads in the occupancy ASCII files.
  
  const Int_t kSystem = AliMUONPreprocessor::kDAQ;
  const char* kId = "OCCUPANCY";
  
  delete fOccupancyMap;
  fOccupancyMap = new AliMUON2DMap(kTRUE);
  
  Master()->Log(Form("Reading occupancy file for Run %d startTime %u endTime %u",
                     run,startTime,endTime));
  
  TList* sources = Master()->GetFileSources(kSystem,kId);
  TIter next(sources);
  TObjString* o(0x0);
  Int_t n(0);
  
  while ( ( o = static_cast<TObjString*>(next()) ) )
  {
    TString fileName(Master()->GetFile(kSystem,kId,o->GetName()));
    Int_t ok = ReadFile(fileName.Data());
    if ( ok>0 || ok == AliMUONTrackerIO::kNoInfoFile )
    {
      n += ok;
    }
  }
  
  delete sources;

  if (!n)
  {
    Master()->Log("Failed to read any occupancy");
    delete fOccupancyMap;
    fOccupancyMap = 0;
    return kFALSE;
  }
  return kTRUE;
}

//_____________________________________________________________________________
UInt_t 
AliMUONOccupancySubprocessor::Process(TMap* /*dcsAliasMap*/)
{
  /// Store the occupancy map into the CDB
  
  if (!fOccupancyMap) 
  {
    // this is the only reason to fail for the moment : getting no occupancy
    // at all.
    return 1;
  }
  
  if ( fOccupancyMap->GetSize() )
  {
    Master()->Log("Storing occupancy map");
  
    AliCDBMetaData metaData;
    metaData.SetBeamPeriod(0);
    metaData.SetResponsible("MUON TRK");
    TString comment("Computed by AliMUONOccupancySubprocessor $Id$");
    comment.ReplaceAll("$","");
    metaData.SetComment(comment.Data());
    
    Bool_t validToInfinity = kFALSE;
    Bool_t result = Master()->Store("Calib", "OccupancyMap", fOccupancyMap, &metaData, 0, validToInfinity);
  
    return ( result != kTRUE ); // return 0 if everything is ok.  
  }
  else
  {
    Master()->Log("No occupancy map to store");
    return 0;
  }
}

//_____________________________________________________________________________
Int_t
AliMUONOccupancySubprocessor::ReadFile(const char* filename)
{
  /// Read the occupancy from an ASCII file.                                  \n
  /// Return kFALSE if reading was not successfull.                           \n
  ///
  
  TString sFilename(gSystem->ExpandPathName(filename));
  
  Master()->Log(Form("Reading %s",sFilename.Data()));
  
  Int_t n = AliMUONTrackerIO::ReadOccupancy(sFilename.Data(),*fOccupancyMap);
  
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
AliMUONOccupancySubprocessor::Print(Option_t* opt) const
{
  /// ouput to screen
  if (fOccupancyMap) fOccupancyMap->Print("",opt);
}

