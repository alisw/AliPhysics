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

#include "TMath.h"
#include "TRandom.h"
#include "AliPHOSCalibData.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"

ClassImp(AliPHOSCalibData)

//________________________________________________________________
  AliPHOSCalibData::AliPHOSCalibData(): 
    TNamed(), fCalibDataEmc(0x0), fCalibDataCpv(0x0)
{
  // Default constructor  

  fEmcDataPath="PHOS/Calib/EmcGainPedestals";
  fCpvDataPath="PHOS/Calib/CpvGainPedestals";

}

//________________________________________________________________
AliPHOSCalibData::AliPHOSCalibData(Int_t runNumber) :
  TNamed("phosCalib","PHOS Calibration Data Manager"),
  fCalibDataEmc(0x0), fCalibDataCpv(0x0)
{
  // Constructor
  
  fEmcDataPath="PHOS/Calib/EmcGainPedestals";
  fCpvDataPath="PHOS/Calib/CpvGainPedestals";

  AliCDBEntry* entryEmc = AliCDBManager::Instance()->Get(fEmcDataPath.Data(),runNumber);
  if(entryEmc)
    fCalibDataEmc = (AliPHOSEmcCalibData*)entryEmc->GetObject();
  
  AliCDBEntry* entryCpv = AliCDBManager::Instance()->Get(fCpvDataPath.Data(),runNumber);
  if(entryCpv)
    fCalibDataCpv = (AliPHOSCpvCalibData*)entryCpv->GetObject();

}

//________________________________________________________________
AliPHOSCalibData::~AliPHOSCalibData()
{
  // Destructor
 
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
  if (fCalibDataEmc) fCalibDataEmc->Print(option);
  if (fCalibDataCpv) fCalibDataCpv->Print(option);
}

//________________________________________________________________
void AliPHOSCalibData::CreateNew()
{
  if(fCalibDataEmc) delete fCalibDataEmc;
  fCalibDataEmc = new AliPHOSEmcCalibData("PHOS-EMC");

  if(fCalibDataCpv) delete fCalibDataCpv;
  fCalibDataCpv = new AliPHOSCpvCalibData("PHOS-CPV");

}

//________________________________________________________________
Bool_t AliPHOSCalibData::WriteEmc(Int_t firstRun, Int_t lastRun, AliCDBMetaData *md)
{

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
  //module, column,raw should follow the internal PHOS convention:
  //module 1:5, column 1:56, row 1:64

  if(fCalibDataEmc) 
    return fCalibDataEmc->GetADCchannelEmc(module,column,row);
  else
    return 0.0015; // default width of one EMC ADC channel in GeV
}

Float_t AliPHOSCalibData::GetADCpedestalEmc(Int_t module, Int_t column, Int_t row) const
{
  if(fCalibDataEmc) 
    return fCalibDataEmc->GetADCpedestalEmc(module,column,row);
  else
    return 0.005; // default EMC ADC pedestal
}

void AliPHOSCalibData::SetADCchannelEmc(Int_t module, Int_t column, Int_t row, Float_t value)
{
  if(!fCalibDataEmc)
    fCalibDataEmc = new AliPHOSEmcCalibData("PHOS-EMC");

  fCalibDataEmc->SetADCchannelEmc(module,column,row,value);
}

void AliPHOSCalibData::SetADCpedestalEmc(Int_t module, Int_t column, Int_t row, Float_t value)
{
  if(!fCalibDataEmc)
    fCalibDataEmc = new AliPHOSEmcCalibData("PHOS-EMC");

  fCalibDataEmc->SetADCpedestalEmc(module,column,row,value);
}

//________________________________________________________________
Float_t AliPHOSCalibData::GetADCchannelCpv(Int_t module, Int_t column, Int_t row) const
{
  //module, column,raw should follow the internal CPV convention:
  //module 1:5, column 1:64, row 1:128

  if(fCalibDataCpv) 
    return fCalibDataCpv->GetADCchannelCpv(module,column,row);
  else
    return 0.0012; // default width of one ADC channel in CPV 'popugais'
}

Float_t AliPHOSCalibData::GetADCpedestalCpv(Int_t module, Int_t column, Int_t row) const
{
  if(fCalibDataCpv) 
    return fCalibDataCpv->GetADCpedestalCpv(module,column,row);
  else
    return 0.012; // default CPV ADC pedestal
}

void AliPHOSCalibData::SetADCchannelCpv(Int_t module, Int_t column, Int_t row, Float_t value)
{
  if(!fCalibDataCpv)
    fCalibDataCpv = new AliPHOSCpvCalibData("PHOS-CPV");

  fCalibDataCpv->SetADCchannelCpv(module,column,row,value);
}

void AliPHOSCalibData::SetADCpedestalCpv(Int_t module, Int_t column, Int_t row, Float_t value)
{
  if(!fCalibDataCpv)
    fCalibDataCpv = new AliPHOSCpvCalibData("PHOS-CPV");

  fCalibDataCpv->SetADCpedestalCpv(module,column,row,value);
}

//________________________________________________________________
void AliPHOSCalibData::RandomEmc()
{

  if(fCalibDataEmc) delete fCalibDataEmc;
  fCalibDataEmc = new AliPHOSEmcCalibData("PHOS-EMC");

  TRandom rn;
  rn.SetSeed(0); //the seed is set to the current  machine clock
  
  Float_t ADCchanelEmc,ADCpedestalEmc;

  for(Int_t module=1; module<6; module++) {
    for(Int_t column=1; column<57; column++) {
      for(Int_t row=1; row<65; row++) {
        ADCchanelEmc=rn.Uniform(0.00075,0.00375); // Cmax/Cmin = 5, (Cmax-Cmin)/2 = 0.0015
        ADCpedestalEmc=rn.Uniform(0.0045,0.0055); //+-10% spread of pedestals from 0.005
        fCalibDataEmc->SetADCchannelEmc(module,column,row,ADCchanelEmc);
        fCalibDataEmc->SetADCpedestalEmc(module,column,row,ADCpedestalEmc);
      }
    }
  }

}

//________________________________________________________________
void AliPHOSCalibData::RandomCpv()
{

  if(fCalibDataCpv) delete fCalibDataCpv;
  fCalibDataCpv = new AliPHOSCpvCalibData("PHOS-CPV");

  TRandom rn;
  rn.SetSeed(0); //the seed is set to the current  machine clock
  
  Float_t ADCchanelCpv,ADCpedestalCpv;

  for(Int_t module=1; module<6; module++) {
    for(Int_t column=1; column<65; column++) {
      for(Int_t row=1; row<129; row++) {
	ADCchanelCpv=TMath::Abs(rn.Uniform(0.0009,0.0015)); // 0.0012 +- 25%
        ADCpedestalCpv=rn.Uniform(0.0048,0.0192); // Ped[max]/Ped[min] = 4, <Ped> = 0.012
        fCalibDataCpv->SetADCchannelCpv(module,column,row,ADCchanelCpv);
        fCalibDataCpv->SetADCpedestalCpv(module,column,row,ADCpedestalCpv);
      }
    }
  }
}
