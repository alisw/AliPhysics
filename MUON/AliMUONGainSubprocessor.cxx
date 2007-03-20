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

#include "AliMUON2DMap.h"
#include "AliMUONCalibParamNF.h"
#include "AliMUONConstants.h"
#include "AliMUONObjectPair.h"
#include "AliMUONGainSubprocessor.h"
#include "AliMUONPreprocessor.h"
#include "AliMUONVDataIterator.h"
#include "AliMpDDLStore.h"
#include "AliMUON2DStoreValidator.h"

#include "AliCDBMetaData.h"
#include "AliLog.h"

#include <TObjString.h>
#include <Riostream.h>
#include <TList.h>
#include <TObjString.h>
#include <TSystem.h>

#include <sstream>

///
/// \class AliMUONGainSubprocessor
///
/// Implementation of AliMUONVSubprocessor class to deal with MUON TRK Gains.
///
/// Gains are read in from an ascii file, with the format :                   \n            
///---------------------------------------------------------------------------\n
///BUS_PATCH   MANU   CHANNEL    Ped.     a0        a1         a2         xlim        P(chi2)    P(chi2)_2  \n 
///---------------------------------------------------------------------------\n
///
/// \author L. Aphecetche
///

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
    return 0;
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
	UInt_t result = Master()->Store("Calib", "Gains", fGains, &metaData, 0, validToInfinity);
  
  return result;  
}

//_____________________________________________________________________________
Int_t
AliMUONGainSubprocessor::ReadFile(const char* filename)
{
  /// Read the Gains from an ASCII file.                                  \n
  /// Format of that file is one line per channel :                           \n
  ///-------------------------------------------------------------------------\n
  ///BUS_PATCH   MANU   CHANNEL    Ped.     a0        a1         a2         xlim        P(chi2)    P(chi2)_2  \n 
  ///-------------------------------------------------------------------------\n
  ///                                                                         \n
  /// Return kFALSE if reading was not successfull.                           \n
  ///
  
  TString sFilename(gSystem->ExpandPathName(filename));
  
  Master()->Log(Form("Reading %s",sFilename.Data()));
  
  std::ifstream in(sFilename.Data());
  if (!in.good()) 
  {
    return 0;
  }
  char line[256];
  Int_t busPatchID, manuID, manuChannel;
  Float_t ped;
  Float_t a0, a1, a2;
  Float_t xlim;
  Float_t chi2, chi22;
  
  static const Int_t kNchannels(64);
  static Bool_t replace(kFALSE);
  Int_t n(0);
  
  while ( in.getline(line,256) )
  {
    if ( strlen(line) < 10 ) continue;
    if ( line[0] == '/' && line[1] == '/' ) continue;
    std::istringstream sin(line);
    AliDebug(3,Form("line=%s",line));
    sin >> busPatchID >> manuID >> manuChannel >> ped >> a0 >> a1 >> a2 >> xlim >> chi2 >> chi22;
    Int_t detElemID = AliMpDDLStore::Instance()->GetDEfromBus(busPatchID);
    AliDebug(3,Form("BUSPATCH %3d DETELEMID %4d MANU %3d CH %3d PED %7.2f A0 %7.2f A1 %7.2f A2 %7.2f"
                    " XLIM %7.2f CHI2 %7.2f CHI22 %7.2f",
                    busPatchID,detElemID,manuID,manuChannel,ped,a0,a1,a2,xlim,chi2,chi22));
    if ( a0==a1 && a1==a2 && a0==-2) continue;
    
    AliMUONVCalibParam* gain = 
      static_cast<AliMUONVCalibParam*>(fGains->Get(detElemID,manuID));
    
    if (!gain) 
    {
      gain = new AliMUONCalibParamNF(6,kNchannels,0);//AliMUONVCalibParam::InvalidFloatValue());  
      fGains->Set(detElemID,manuID,gain,replace);
    }
    gain->SetValueAsFloat(manuChannel,0,a0);
    gain->SetValueAsFloat(manuChannel,1,a1);
    gain->SetValueAsFloat(manuChannel,2,a2);
    gain->SetValueAsFloat(manuChannel,3,xlim);
    gain->SetValueAsFloat(manuChannel,4,chi2);
    gain->SetValueAsFloat(manuChannel,5,chi22);
    ++n;
  }
  in.close();
  return n;
}


//_____________________________________________________________________________
void
AliMUONGainSubprocessor::Print(Option_t* opt) const
{
  /// ouput to screen
  AliMUONVDataIterator* it = fGains->Iterator();
  AliMUONObjectPair* p;

  while ( ( p = static_cast<AliMUONObjectPair*>(it->Next() ) ) )
  {
    AliMUONVCalibParam* value = static_cast<AliMUONVCalibParam*>(p->Value());
    value->Print(opt);
  }
}
