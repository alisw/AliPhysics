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
#include "AliMUONCalibParam2F.h"
#include "Riostream.h"
#include "TObjString.h"
#include "TSystem.h"
#include <sstream>
#include "AliMUONVDataIterator.h"
#include "AliMUONConstants.h"
#include "AliMUONObjectPair.h"
#include "AliMpDDLStore.h"
#include "AliMUONPreprocessor.h"

///
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
///

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
  // default ctor
}

//_____________________________________________________________________________
AliMUONPedestalSubprocessor::~AliMUONPedestalSubprocessor()
{
  // dtor
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
  
  AliInfo(Form("Reading pedestal files for Run %d startTime %ld endTime %ld",
               run,startTime,endTime));
  
  TList* sources = Master()->GetFileSources(kSystem,kId);
  TIter next(sources);
  TObjString* o(0x0);
  while ( ( o = static_cast<TObjString*>(next()) ) )
  {
    TString fileName(Master()->GetFile(kSystem,kId,o->GetName()));
    Bool_t ok = ReadFile(fileName.Data());
    if (!ok)
    {
      AliError(Form("Could not read file %s",fileName.Data()));
    }
  }
  delete sources;
}

//_____________________________________________________________________________
UInt_t 
AliMUONPedestalSubprocessor::Process(TMap* /*dcsAliasMap*/)
{
  // Store the pedestals into the CDB
  
  if (!fPedestals) return 0;
  
  AliInfo("Validating pedestals");
  AliMUON2DStoreValidator validator;
  TObjArray* missing =
    validator.Validate(*fPedestals,AliMUONCalibParam2F::InvalidFloatValue());  
  
  if (missing)  
  {
    validator.Report(*missing);
//    AliError("Will not write into CDB as some pieces are missing...");
//    return 0;
  }
  
  AliInfo("Storing pedestals");
  
  AliCDBMetaData metaData;
	metaData.SetBeamPeriod(0);
	metaData.SetResponsible("MUON TRK");
	metaData.SetComment("Computed by AliMUONPedestalSubprocessor $Id$");
  
	UInt_t result = Master()->Store("Calib", "Pedestals", fPedestals, &metaData, 0, 0);
  
  return result;  
}

//_____________________________________________________________________________
Bool_t
AliMUONPedestalSubprocessor::ReadFile(const char* filename)
{
  // Read the pedestals from an ASCII file.
  // Format of that file is one line per channel :
  //---------------------------------------------------------------------------
  // BUS_PATCH MANU_ADDR CHANNEL      MEAN       SIGMA
  //---------------------------------------------------------------------------
  //
  // Return kFALSE if reading was not successfull.
  //
  
  AliInfo(Form("Reading %s",filename));
  
  std::ifstream in(gSystem->ExpandPathName(filename));
  if (!in.good()) return kFALSE;
  
  char line[80];
  Int_t busPatchID, manuID, manuChannel;
  Float_t pedMean, pedSigma;
  static const Int_t kNchannels(64);
  static Bool_t replace(kFALSE);
  
  while ( in.getline(line,80) )
  {
    if ( line[0] == '/' && line[1] == '/' ) continue;
    std::istringstream sin(line);
    sin >> busPatchID >> manuID >> manuChannel >> pedMean >> pedSigma;
    Int_t detElemID = AliMpDDLStore::Instance()->GetDEfromBus(busPatchID);
    AliDebug(3,Form("BUSPATCH %3d DETELEMID %4d MANU %3d CH %3d MEAN %7.2f SIGMA %7.2f",
             busPatchID,detElemID,manuID,manuChannel,pedMean,pedSigma));
    
    AliMUONVCalibParam* ped = 
      static_cast<AliMUONVCalibParam*>(fPedestals->Get(detElemID,manuID));
    
    if (!ped) 
    {
      ped = new AliMUONCalibParam2F(kNchannels,AliMUONCalibParam2F::InvalidFloatValue());  
      fPedestals->Set(detElemID,manuID,ped,replace);
    }
    ped->SetValueAsFloat(manuChannel,0,pedMean);
    ped->SetValueAsFloat(manuChannel,1,pedSigma);
  }
  in.close();
  return kTRUE;
}


//_____________________________________________________________________________
void
AliMUONPedestalSubprocessor::Print(Option_t* opt) const
{
  /// ouput to screen
  AliMUONVDataIterator* it = fPedestals->Iterator();
  AliMUONObjectPair* p;

  while ( ( p = static_cast<AliMUONObjectPair*>(it->Next() ) ) )
  {
    AliMUONVCalibParam* value = static_cast<AliMUONVCalibParam*>(p->Value());
    value->Print(opt);
  }
}
