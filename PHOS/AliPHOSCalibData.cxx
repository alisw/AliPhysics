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
#include "AliLog.h"
#include "AliPHOSCalibData.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliPHOSEmcCalibData.h"
#include "AliPHOSCpvCalibData.h"
#include "AliPHOSEmcBadChannelsMap.h"
#include "AliCDBMetaData.h"

ClassImp(AliPHOSCalibData)

//________________________________________________________________
  AliPHOSCalibData::AliPHOSCalibData(): 
    TNamed(), 
    fCalibDataEmc(0x0), 
    fCalibDataCpv(0x0),
    fEmcBadChannelsMap(0x0),
    fEmcDataPath("PHOS/Calib/EmcGainPedestals"),
    fCpvDataPath("PHOS/Calib/CpvGainPedestals"),
    fEmcBadChannelsMapPath("PHOS/Calib/EmcBadChannels")
{
  // Default constructor.
  // Open CDB entry, get EMC and CPV calibration data and bad channel map.
  // If EMC or CPV calibration data does not exist, stop the run
 
}

//________________________________________________________________
AliPHOSCalibData::AliPHOSCalibData(Int_t runNumber) :
  TNamed("phosCalib","PHOS Calibration Data Manager"),
  fCalibDataEmc(0x0), fCalibDataCpv(0x0), fEmcBadChannelsMap(0x0),
  fEmcDataPath("PHOS/Calib/EmcGainPedestals"),
  fCpvDataPath("PHOS/Calib/CpvGainPedestals"),
  fEmcBadChannelsMapPath("PHOS/Calib/EmcBadChannels")
{
  // Constructor
  // Open CDB entry, get EMC and CPV calibration data and bad channel map.
  // If EMC or CPV calibration data does not exist, stop the run

  AliCDBEntry* entryEmc = AliCDBManager::Instance()->Get(fEmcDataPath.Data(),runNumber);
  if(entryEmc)
    fCalibDataEmc = (AliPHOSEmcCalibData*)entryEmc->GetObject();

  if(!fCalibDataEmc)
    AliFatal("Calibration parameters for PHOS EMC not found. Stop reconstruction!\n");
  
  AliCDBEntry* entryCpv = AliCDBManager::Instance()->Get(fCpvDataPath.Data(),runNumber);
  if(entryCpv)
    fCalibDataCpv = (AliPHOSCpvCalibData*)entryCpv->GetObject();

  if(!fCalibDataCpv)
    AliFatal("Calibration parameters for PHOS CPV not found. Stop reconstruction!\n");
  
  AliCDBEntry* entryEmcBadMap = AliCDBManager::Instance()->
    Get(fEmcBadChannelsMapPath.Data(),runNumber);
  if(entryEmcBadMap)
    fEmcBadChannelsMap = (AliPHOSEmcBadChannelsMap*)entryEmcBadMap->GetObject(); 

}

//________________________________________________________________
AliPHOSCalibData::AliPHOSCalibData(AliPHOSCalibData & phosCDB) :
  TNamed(phosCDB),
  fCalibDataEmc(phosCDB.fCalibDataEmc),
  fCalibDataCpv(phosCDB.fCalibDataCpv),
  fEmcBadChannelsMap(phosCDB.fEmcBadChannelsMap),
  fEmcDataPath(phosCDB.fEmcDataPath),
  fCpvDataPath(phosCDB.fCpvDataPath),
  fEmcBadChannelsMapPath(phosCDB.fEmcBadChannelsMapPath)
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
  fEmcBadChannelsMap = rhs.fEmcBadChannelsMap;
  fEmcDataPath  = rhs.fEmcDataPath;
  fCpvDataPath  = rhs.fCpvDataPath;
  fEmcBadChannelsMapPath = rhs.fEmcBadChannelsMapPath;
  
  return *this;
}

//________________________________________________________________
void AliPHOSCalibData::Reset()
{
  // Set all pedestals to 0 and all ADC channels to 1,
  // and all channels are good (alive)

  fCalibDataEmc     ->Reset();
  fCalibDataCpv     ->Reset();
  fEmcBadChannelsMap->Reset();
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

  if(fEmcBadChannelsMap) delete fEmcBadChannelsMap;
  fEmcBadChannelsMap = new AliPHOSEmcBadChannelsMap();

}

//________________________________________________________________
Bool_t AliPHOSCalibData::WriteEmc(Int_t firstRun, Int_t lastRun, AliCDBMetaData *md)
{
  // Write EMC calibration container to CDB

  if(!fCalibDataEmc) return kFALSE;

  AliCDBStorage* storage = AliCDBManager::Instance()->GetSpecificStorage("PHOS/*");
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
  
  AliCDBStorage* storage = AliCDBManager::Instance()->GetSpecificStorage("PHOS/*");
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
Bool_t AliPHOSCalibData::WriteEmcBadChannelsMap(Int_t firstRun,Int_t lastRun,AliCDBMetaData *md)
{
  //Write EMC bad channels map into CDB.

  if(!fEmcBadChannelsMap) return kFALSE;
  
  AliCDBStorage* storage = AliCDBManager::Instance()->GetSpecificStorage("PHOS/*");
  if(!storage)
    storage = AliCDBManager::Instance()->GetDefaultStorage();

  if(storage) { 
    AliCDBId id(fEmcBadChannelsMapPath.Data(),firstRun,lastRun);
    storage->Put(fEmcBadChannelsMap,id, md);
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
void AliPHOSCalibData::SetADCchannelEmc(Int_t module, Int_t column, Int_t row, Float_t value)
{
  // Set EMC calibration coefficient for (module,column,row)

  if(!fCalibDataEmc)
    fCalibDataEmc = new AliPHOSEmcCalibData("PHOS-EMC");

  fCalibDataEmc->SetADCchannelEmc(module,column,row,value);
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
void AliPHOSCalibData::SetADCpedestalEmc(Int_t module, Int_t column, Int_t row, Float_t value)
{
  // Set EMC pedestal for (module,column,row)

  if(!fCalibDataEmc)
    fCalibDataEmc = new AliPHOSEmcCalibData("PHOS-EMC");

  fCalibDataEmc->SetADCpedestalEmc(module,column,row,value);
}

//________________________________________________________________
Float_t AliPHOSCalibData::GetHighLowRatioEmc(Int_t module, Int_t column, Int_t row) const
{
  // Return EMC calibration coefficient
  // for channel defined by (module,column,row)
  // module, column,raw should follow the internal PHOS convention:
  // module 1:5, column 1:56, row 1:64
  // if CBD instance exists, the value is taken from CDB.
  // Otherwise it is an ideal one
 
  if(fCalibDataEmc)
    return fCalibDataEmc->GetHighLowRatioEmc(module,column,row);
  else
    return 1.0; // default width of one EMC ADC channel in GeV
}
 
//________________________________________________________________
void AliPHOSCalibData::SetHighLowRatioEmc(Int_t module, Int_t column, Int_t row, Float_t value)
{
  // Set EMC calibration coefficient for (module,column,row)
 
  if(!fCalibDataEmc)
    fCalibDataEmc = new AliPHOSEmcCalibData("PHOS-EMC");
 
  fCalibDataEmc->SetHighLowRatioEmc(module,column,row,value);
}
 
//________________________________________________________________
Float_t AliPHOSCalibData::GetTimeShiftEmc(Int_t module, Int_t column, Int_t row) const
{ 
  // Return EMC calibration coefficient 
  // for channel defined by (module,column,row)                                
  // module, column,raw should follow the internal PHOS convention:            
  // module 1:5, column 1:56, row 1:64 
  // if CBD instance exists, the value is taken from CDB. 
  // Otherwise it is an ideal one  
  
  if(fCalibDataEmc)
    return fCalibDataEmc->GetTimeShiftEmc(module,column,row);
  else
    return 1.0; // default width of one EMC ADC channel in GeV
}
 
//________________________________________________________________
void AliPHOSCalibData::SetTimeShiftEmc(Int_t module, Int_t column, Int_t row, Float_t value)
{
  // Set EMC calibration coefficient for (module,column,row)
 
  if(!fCalibDataEmc)
    fCalibDataEmc = new AliPHOSEmcCalibData("PHOS-EMC");
 
  fCalibDataEmc->SetTimeShiftEmc(module,column,row,value);
}
//________________________________________________________________
Float_t AliPHOSCalibData::GetSampleTimeStep() const 
{
  //Get conversion coeff. from sample time step to seconds.
  //Negative value means that it is not used in reconstruction
  //but only in simulation of raw.
  if(fCalibDataEmc)
    return fCalibDataEmc->GetSampleTimeStep();
  else
    return 0.0; // default width of one EMC ADC channel in GeV
}
//________________________________________________________________
void   AliPHOSCalibData::SetSampleTimeStep(Float_t step)
{
  //Set conversion coeff. from sample time step to seconds.
  //Negative value means that it is not used in reconstruction
  //but only in simulation of raw.
  if(!fCalibDataEmc)
    fCalibDataEmc = new AliPHOSEmcCalibData("PHOS-EMC");

  fCalibDataEmc->SetSampleTimeStep(step) ;
}
//________________________________________________________________
Int_t AliPHOSCalibData::GetAltroOffsetEmc(Int_t module, Int_t column, Int_t row) const
{
  // Return ALTRO pedestal coefficient
  // for channel defined by (module,column,row)
  // module, column,raw should follow the internal PHOS convention:
  // module 1:5, column 1:56, row 1:64
  // if CBD instance exists, the value is taken from CDB.
  // Otherwise it is an ideal one
 
  if(fCalibDataEmc)
    return fCalibDataEmc->GetAltroOffsetEmc(module,column,row);
  else
    return 0; // default width of one EMC ADC channel in GeV
}
 
//________________________________________________________________
void AliPHOSCalibData::SetAltroOffsetEmc(Int_t module, Int_t column, Int_t row, Int_t value)
{
  // Set altro offset for (module,column,row)
 
  if(!fCalibDataEmc)
    fCalibDataEmc = new AliPHOSEmcCalibData("PHOS-EMC");
 
  fCalibDataEmc->SetAltroOffsetEmc(module,column,row,value);
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
//________________________________________________________________
Bool_t AliPHOSCalibData::IsBadChannelEmc(Int_t module, Int_t col, Int_t row) const
{
  //If no bad channels map found, channel considered good

  if(fEmcBadChannelsMap) 
    return fEmcBadChannelsMap->IsBadChannel(module,col,row);
  else
    return kFALSE;
}

//________________________________________________________________
Int_t AliPHOSCalibData::GetNumOfEmcBadChannels() const
{
  if(fEmcBadChannelsMap)
    return fEmcBadChannelsMap->GetNumOfBadChannels();
  else
    return 0;
}
//________________________________________________________________
void AliPHOSCalibData::EmcBadChannelIds(Int_t *badIds)
{
  //Fill array badIds by the Ids of EMC bad channels.
  //Array badIds of length GetNumOfBadChannels() should be prepared in advance. 

  if(fEmcBadChannelsMap)              
    fEmcBadChannelsMap->BadChannelIds(badIds);
}

//________________________________________________________________
Float_t AliPHOSCalibData::GetADCchannelEmcDecalib(Int_t module, Int_t column, Int_t row) const
{
  // Return random EMC (de)calibration factor O(1) for channel defined by (module,column,row). 
  // Used in simulation.
  
  // module, column,raw should follow the internal PHOS convention:
  // module 1:5, column 1:56, row 1:64
  // if CBD instance exists, the value is taken from CDB.
  // Otherwise it is an ideal one (no decalibration).
  
  if(fCalibDataEmc) 
    return fCalibDataEmc->GetADCchannelEmcDecalib(module,column,row);
  else
    return 1.0; // no decalibration by default
}

//________________________________________________________________
void AliPHOSCalibData::SetADCchannelEmcDecalib(Int_t module, Int_t column, Int_t row, Float_t value)
{
  // Set EMC (de)calibration factor for (module,column,row).
  // Used in simulation.
  
  if(!fCalibDataEmc)
    fCalibDataEmc = new AliPHOSEmcCalibData("PHOS-EMC");
  
  fCalibDataEmc->SetADCchannelEmcDecalib(module,column,row,value);
}
