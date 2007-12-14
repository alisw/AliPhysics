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
fGains(0x0)
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
void 
AliMUONGainSubprocessor::Initialize(Int_t run, UInt_t startTime, UInt_t endTime)
{
  /// When starting a new run, reads in the Gains ASCII files.
  
  const Int_t kSystem = AliMUONPreprocessor::kDAQ;
  const char* kId = "GAINS";
  
  delete fGains;
  fGains = new AliMUON2DMap(kTRUE);
  
  Master()->Log(Form("Reading Gain files for Run %d startTime %ld endTime %ld",
                     run,startTime,endTime));
  
  TList* sources = Master()->GetFileSources(kSystem,kId);
  TIter next(sources);
  TObjString* o(0x0);
  Int_t n(0);
  
  while ( ( o = static_cast<TObjString*>(next()) ) )
  {
    TString fileName(Master()->GetFile(kSystem,kId,o->GetName()));
    Int_t ok = ReadFile(fileName.Data());
    if (!ok)
    {
      Master()->Log(Form("Could not read file %s",fileName.Data()));
    }
    else
    {
      n += ok;
    }
  }
  
  if (!n)
  {
    Master()->Log("Failed to read any Gains");
    delete fGains;
    fGains = 0;
  }
  delete sources;
}

//_____________________________________________________________________________
UInt_t 
AliMUONGainSubprocessor::Process(TMap* /*dcsAliasMap*/)
{
  /// Store the Gains into the CDB
  
  if (!fGains) 
  {
    // this is the only reason to fail for the moment : getting no Gain
    // at all.
    return 1;
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
	metaData.SetComment("Computed by AliMUONGainSubprocessor $Id$");
  
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
  
  std::ifstream in(sFilename.Data());
  if (!in.good()) 
  {
    return 0;
  }
  char line[1024];
  Int_t busPatchID, manuID, manuChannel;
  Float_t a0, a1;
  Int_t thres;
  UInt_t qual;
  const Int_t kSaturation(3000); // FIXME: how to get this number ?
  
  static const Int_t kNchannels(AliMpConstants::ManuNofChannels());
  Int_t n(0);
  
  while ( in.getline(line,1024) )
  {
    if ( strlen(line) < 10 ) continue;
    if ( line[0] == '/' && line[1] == '/' ) continue;

    sscanf(line,"%d %d %d %f %f %d %x",&busPatchID,&manuID,&manuChannel,
           &a0,&a1,&thres,&qual); 
    AliDebug(3,Form("line=%s",line));
    Int_t detElemID = AliMpDDLStore::Instance()->GetDEfromBus(busPatchID);
    AliDebug(3,Form("BUSPATCH %3d DETELEMID %4d MANU %3d CH %3d A0 %7.2f "
                    "A1 %e THRES %5d QUAL %x",
                    busPatchID,detElemID,manuID,manuChannel,a0,a1,thres,qual));
    if ( qual == 0 ) continue;
    
    AliMUONVCalibParam* gain = 
      static_cast<AliMUONVCalibParam*>(fGains->FindObject(detElemID,manuID));
    
    if (!gain) 
    {
      gain = new AliMUONCalibParamNF(5,kNchannels,detElemID,manuID,0);
      fGains->Add(gain);
    }
    gain->SetValueAsFloat(manuChannel,0,a0);
    gain->SetValueAsFloat(manuChannel,1,a1);
    gain->SetValueAsInt(manuChannel,2,thres);
    gain->SetValueAsInt(manuChannel,3,qual);
    gain->SetValueAsInt(manuChannel,4,kSaturation);
    ++n;
  }
  in.close();
  return n;
}
