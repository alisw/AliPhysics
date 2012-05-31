// Script to create simulation parameters and store them into CDB
// Author: Gustavo Conesa Balbastre

/* $Id: AliEMCALSetSimParamCDB.C 37597 2009-12-04 13:48:42Z gconesab $ */

#if !defined(__CINT__)
#include "TControlBar.h"
#include "TString.h"

#include "AliEMCALSimParam.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#endif


void AliEMCALSetSimParamCDB()
{
  
  // Create an object AliEMCALRecParam and store it to OCDB
  
  //Activate CDB storage
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  
  // Create simulation parameter object and set parameter values
  AliEMCALSimParam *simParam = new AliEMCALSimParam();
  //Digits 
//   simParam->SetDigitThreshold(3) ;
//   simParam->SetPinNoise(0.012) ;
//  simParam->SetTimeDelay(600e-9) ;      
//   simParam->SetTimeResolution(0.6e-9) ; 
//   simParam->SetNADCED( (Int_t) TMath::Power(2,16)) ;     
//   simParam->SetMeanPhotonElectron(4400);
  //SDigits
//   simParam->SetA(0) ;          
//   simParam->SetB(1e6) ;             
//   simParam->SetECPrimaryThreshold(0.05);

  // Store calibration data into database  
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("G. Conesa");
  md->SetComment("Simulation Parameters: EMCAL");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);
  
  AliCDBId id("EMCAL/Calib/SimParam",0,AliCDBRunRange::Infinity());
  cdb->GetDefaultStorage()->Put(simParam, id, md);
  
  return;
}



