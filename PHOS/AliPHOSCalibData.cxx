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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// class for PHOS calibration                                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "TRandom.h"
#include "AliPHOSCalibData.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliPHOSEmcCalibData.h"
#include "AliPHOSCpvCalibData.h"
#include "AliCDBMetaData.h"

ClassImp(AliPHOSCalibData)

//________________________________________________________________
  AliPHOSCalibData::AliPHOSCalibData(): 
    TNamed(), 
    fCalibDataEmc(0x0), 
    fCalibDataCpv(0x0),
    fEmcDataPath("PHOS/Calib/EmcGainPedestals"),
    fCpvDataPath("PHOS/Calib/CpvGainPedestals")
{
  // Default constructor
  
  AliCDBEntry* entryEmc = AliCDBManager::Instance()->Get(fEmcDataPath.Data());
  if(entryEmc)
    fCalibDataEmc = (AliPHOSEmcCalibData*)entryEmc->GetObject();
  
  AliCDBEntry* entryCpv = AliCDBManager::Instance()->Get(fCpvDataPath.Data());
  if(entryCpv)
    fCalibDataCpv = (AliPHOSCpvCalibData*)entryCpv->GetObject();

}

//________________________________________________________________
AliPHOSCalibData::AliPHOSCalibData(Int_t runNumber) :
  TNamed("phosCalib","PHOS Calibration Data Manager"),
  fCalibDataEmc(0x0), fCalibDataCpv(0x0),
  fEmcDataPath("PHOS/Calib/EmcGainPedestals"),
  fCpvDataPath("PHOS/Calib/CpvGainPedestals")
{
  // Constructor
  AliCDBEntry* entryEmc = AliCDBManager::Instance()->Get(fEmcDataPath.Data(),runNumber);
  if(entryEmc)
    fCalibDataEmc = (AliPHOSEmcCalibData*)entryEmc->GetObject();
  
  AliCDBEntry* entryCpv = AliCDBManager::Instance()->Get(fCpvDataPath.Data(),runNumber);
  if(entryCpv)
    fCalibDataCpv = (AliPHOSCpvCalibData*)entryCpv->GetObject();

}

//________________________________________________________________
AliPHOSCalibData::AliPHOSCalibData(AliPHOSCalibData & phosCDB) :
  TNamed(phosCDB),
  fCalibDataEmc(phosCDB.fCalibDataEmc),
  fCalibDataCpv(phosCDB.fCalibDataCpv),
  fEmcDataPath(phosCDB.fEmcDataPath),
  fCpvDataPath(phosCDB.fCpvDataPath)
{
  // Copy constructor
}
//________________________________________________________________
AliPHOSCalibData::~AliPHOSCalibData()
{
  // Destructor
 
}

AliPHOSCalibData & AliPHOSCalibData::operator = (const AliPHOSCalibData & rhs)
{
  //Copy-assignment. Does not delete anything (see destructor)
  //compiler generated is ok, but ... because -Weffc++ and pointer
  //members we have to define it explicitly.
  TNamed::operator=(rhs);
  fCalibDataEmc = rhs.fCalibDataEmc;
  fCalibDataCpv = rhs.fCalibDataCpv;
  fEmcDataPath  = rhs.fEmcDataPath;
  fCpvDataPath  = rhs.fCpvDataPath;
  
  return *this;
}

//________________________________________________________________
void AliPHOSCalibData::Reset()
{
  // Set all pedestals to 0 and all ADC channels to 1

  fCalibDataEmc->Reset();
  fCalibDataCpv->Reset();
}

//________________________________________________________________
void  AliPHOSCalibData::Print(Option_t *option) const
{
  // Print EMC and CPV calibration containers
  // Input: option="ped"  to print pedestals
  //        option="gain" to print calibration coefficients
  if (fCalibDataEmc) fCalibDataEmc->Print(option);
  if (fCalibDataCpv) fCalibDataCpv->Print(option);
}

//________________________________________________________________
void AliPHOSCalibData::CreateNew()
{
  // Create new EMC and CPV calibration containers with ideal coefficients

  if(fCalibDataEmc) delete fCalibDataEmc;
  fCalibDataEmc = new AliPHOSEmcCalibData("PHOS-EMC");

  if(fCalibDataCpv) delete fCalibDataCpv;
  fCalibDataCpv = new AliPHOSCpvCalibData("PHOS-CPV");

}

//________________________________________________________________
Bool_t AliPHOSCalibData::WriteEmc(Int_t firstRun, Int_t lastRun, AliCDBMetaData *md)
{
  // Write EMC calibration container to CDB

  if(!fCalibDataEmc) return kFALSE;

  AliCDBStorage* storage = AliCDBManager::Instance()->GetSpecificStorage("PHOS");
  if(!storage)
    storage = AliCDBManager::Instance()->GetDefaultStorage();

  if(storage) { 
    AliCDBId id(fEmcDataPath.Data(),firstRun,lastRun);
    storage->Put(fCalibDataEmc,id, md);
    return kTRUE;
  }
  else
    return kFALSE;

}

//________________________________________________________________
Bool_t AliPHOSCalibData::WriteCpv(Int_t firstRun, Int_t lastRun, AliCDBMetaData *md)
{
  // Write CPV calibration container to CDB

  if(!fCalibDataCpv) return kFALSE;
  
  AliCDBStorage* storage = AliCDBManager::Instance()->GetSpecificStorage("PHOS");
  if(!storage)
    storage = AliCDBManager::Instance()->GetDefaultStorage();

  if(storage) { 
    AliCDBId id(fCpvDataPath.Data(),firstRun,lastRun);
    storage->Put(fCalibDataCpv,id, md);
    return kTRUE;
  }
  else
    return kFALSE;

}

//________________________________________________________________
Float_t AliPHOSCalibData::GetADCchannelEmc(Int_t module, Int_t column, Int_t row) const
{
  // Return EMC calibration coefficient
  // for channel defined by (module,column,row)
  // module, column,raw should follow the internal PHOS convention:
  // module 1:5, column 1:56, row 1:64
  // if CBD instance exists, the value is taken from CDB.
  // Otherwise it is an ideal one

  if(fCalibDataEmc) 
    return fCalibDataEmc->GetADCchannelEmc(module,column,row);
  else
    return 1.0; // default width of one EMC ADC channel in GeV
}

//________________________________________________________________
Float_t AliPHOSCalibData::GetADCpedestalEmc(Int_t module, Int_t column, Int_t row) const
{
  // Return EMC pedestal for channel defined by (module,column,row)
  // module, column,raw should follow the internal PHOS convention:
  // module 1:5, column 1:56, row 1:64
  // if CBD instance exists, the value is taken from CDB.
  // Otherwise it is an ideal one

  if(fCalibDataEmc) 
    return fCalibDataEmc->GetADCpedestalEmc(module,column,row);
  else
    return 0.0; // default EMC ADC pedestal
}

//________________________________________________________________
void AliPHOSCalibData::SetADCchannelEmc(Int_t module, Int_t column, Int_t row, Float_t value)
{
  // Set EMC calibration coefficient for (module,column,row)

  if(!fCalibDataEmc)
    fCalibDataEmc = new AliPHOSEmcCalibData("PHOS-EMC");

  fCalibDataEmc->SetADCchannelEmc(module,column,row,value);
}

//________________________________________________________________
void AliPHOSCalibData::SetADCpedestalEmc(Int_t module, Int_t column, Int_t row, Float_t value)
{
  // Set EMC pedestal for (module,column,row)

  if(!fCalibDataEmc)
    fCalibDataEmc = new AliPHOSEmcCalibData("PHOS-EMC");

  fCalibDataEmc->SetADCpedestalEmc(module,column,row,value);
}

//________________________________________________________________
Float_t AliPHOSCalibData::GetADCchannelCpv(Int_t module, Int_t column, Int_t row) const
{
  // Return CPV calibration coefficient
  // for channel defined by (module,column,row)
  // module, column,raw should follow the internal CPV convention:
  // module 1:5, column 1:56, row 1:128
  // if CBD instance exists, the value is taken from CDB.
  // Otherwise it is an ideal one

  if(fCalibDataCpv) 
    return fCalibDataCpv->GetADCchannelCpv(module,column,row);
  else
    return 0.0012; // default width of one ADC channel in CPV arbitrary units
}

//________________________________________________________________
Float_t AliPHOSCalibData::GetADCpedestalCpv(Int_t module, Int_t column, Int_t row) const
{
  // Return CPV pedestal
  // for channel defined by (module,column,row)
  // module, column,raw should follow the internal CPV convention:
  // module 1:5, column 1:56, row 1:128
  // if CBD instance exists, the value is taken from CDB.
  // Otherwise it is an ideal one

  if(fCalibDataCpv) 
    return fCalibDataCpv->GetADCpedestalCpv(module,column,row);
  else
    return 0.012; // default CPV ADC pedestal
}

//________________________________________________________________
void AliPHOSCalibData::SetADCchannelCpv(Int_t module, Int_t column, Int_t row, Float_t value)
{
  // Set CPV calibration coefficient for (module,column,row)

  if(!fCalibDataCpv)
    fCalibDataCpv = new AliPHOSCpvCalibData("PHOS-CPV");

  fCalibDataCpv->SetADCchannelCpv(module,column,row,value);
}

//________________________________________________________________
void AliPHOSCalibData::SetADCpedestalCpv(Int_t module, Int_t column, Int_t row, Float_t value)
{
  // Set CPV pedestal for (module,column,row)

  if(!fCalibDataCpv)
    fCalibDataCpv = new AliPHOSCpvCalibData("PHOS-CPV");

  fCalibDataCpv->SetADCpedestalCpv(module,column,row,value);
}

//________________________________________________________________
void AliPHOSCalibData::RandomEmc(Float_t ccMin, Float_t ccMax)
{
  // Create decalibrated EMC with calibration coefficients and pedestals
  // randomly distributed within hard-coded limits
  // Default spread of calibration parameters is Cmax/Cmin = 4, (Cmax-Cmin)/2 = 1

  if(fCalibDataEmc) delete fCalibDataEmc;
  fCalibDataEmc = new AliPHOSEmcCalibData("PHOS-EMC");

  TRandom rn;
  rn.SetSeed(0); //the seed is set to the current  machine clock
  
  Float_t adcChannelEmc,adcPedestalEmc;

  for(Int_t module=1; module<6; module++) {
    for(Int_t column=1; column<57; column++) {
      for(Int_t row=1; row<65; row++) {
        adcChannelEmc =rn.Uniform(ccMin,ccMax);
        adcPedestalEmc=rn.Uniform(0.0,0.0); // 0 spread of pedestals
        fCalibDataEmc->SetADCchannelEmc(module,column,row,adcChannelEmc);
        fCalibDataEmc->SetADCpedestalEmc(module,column,row,adcPedestalEmc);
      }
    }
  }

}

//________________________________________________________________
void AliPHOSCalibData::RandomCpv(Float_t ccMin, Float_t ccMax)
{
  // Create decalibrated CPV with calibration coefficients and pedestals
  // randomly distributed within hard-coded limits
  // Default spread of calibration parameters is  0.0012 +- 25%

  if(fCalibDataCpv) delete fCalibDataCpv;
  fCalibDataCpv = new AliPHOSCpvCalibData("PHOS-CPV");

  TRandom rn;
  rn.SetSeed(0); //the seed is set to the current  machine clock
  
  Float_t adcChannelCpv,adcPedestalCpv;

  for(Int_t module=1; module<6; module++) {
    for(Int_t column=1; column<57; column++) {
      for(Int_t row=1; row<129; row++) {
	adcChannelCpv =rn.Uniform(ccMin,ccMax);
        adcPedestalCpv=rn.Uniform(0.0048,0.0192); // Ped[max]/Ped[min] = 4, <Ped> = 0.012
        fCalibDataCpv->SetADCchannelCpv(module,column,row,adcChannelCpv);
        fCalibDataCpv->SetADCpedestalCpv(module,column,row,adcPedestalCpv);
      }
    }
  }
}
