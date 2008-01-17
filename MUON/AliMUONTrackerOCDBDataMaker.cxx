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

#include "AliMUONTrackerOCDBDataMaker.h"

#include "AliMUONTrackerData.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliMUONCalibrationData.h"
#include "AliMUONVStore.h"
#include "AliLog.h"
#include <TString.h>

///\class AliMUONTrackerOCDBDataMaker
///
/// Producer of AliMUONVTrackerData from OCDB data
///
///\author Laurent Aphecetche, Subatech

///\cond CLASSIMP
ClassImp(AliMUONTrackerOCDBDataMaker)
///\endcond

//_____________________________________________________________________________
AliMUONTrackerOCDBDataMaker::AliMUONTrackerOCDBDataMaker(const char* ocdbPath,
                                                           Int_t runNumber,
                                                           const char* type)
: AliMUONVTrackerDataMaker(),
  fIsValid(kTRUE),
  fIsOwner(kTRUE),
  fData(0x0),
  fSource(Form("%s-%010d-%s",ocdbPath,runNumber,type))
{
    /// Ctor
    AliCDBStorage* storage = AliCDBManager::Instance()->GetDefaultStorage();
    
    AliCDBManager::Instance()->SetDefaultStorage(ocdbPath);
    
    AliMUONVStore* store(0x0);
    
    TString stype(type);
    stype.ToUpper();
    
    if ( stype == "PEDESTALS" )
    {
      fData = new AliMUONTrackerData(Form("PED%d",runNumber),"Pedestals",2,kFALSE);
      fData->SetDimensionName(0,"Mean");
      fData->SetDimensionName(1,"Sigma");
      store = AliMUONCalibrationData::CreatePedestals(runNumber);
    }
    else if ( stype == "GAINS" ) 
    {
      fData = new AliMUONTrackerData(Form("GAIN%d",runNumber),"Gains",5,kFALSE);
      fData->SetDimensionName(0,"a0");
      fData->SetDimensionName(1,"a1");
      fData->SetDimensionName(2,"thres");
      fData->SetDimensionName(3,"qual");
      fData->SetDimensionName(4,"sat");
      store = AliMUONCalibrationData::CreateGains(runNumber);
    }
    else if ( stype == "CAPACITANCES" )
    {
      fData = new AliMUONTrackerData(Form("CAPA%d",runNumber),"Capacitances",2,kFALSE);
      fData->SetDimensionName(0,"Capa");
      fData->SetDimensionName(1,"Injection gain");
      store = AliMUONCalibrationData::CreateCapacitances(runNumber);
    }
    
    AliCDBManager::Instance()->SetDefaultStorage(storage);

    if (!store)
    {
      fIsValid = kFALSE;
      delete fData;
      fData = 0x0;
      AliError("Could not create store");
      return;
    }
    
    fData->Add(*store);
}

//_____________________________________________________________________________
AliMUONTrackerOCDBDataMaker::~AliMUONTrackerOCDBDataMaker()
{
  /// dtor
  if ( fIsOwner) delete fData;
}
