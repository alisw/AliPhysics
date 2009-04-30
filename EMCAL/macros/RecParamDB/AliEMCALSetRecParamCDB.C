// Script to create reconstruction parameters and store them into CDB
// Author: Yuri Kharlov

/* $Id$ */

#if !defined(__CINT__)
#include "TControlBar.h"
#include "TString.h"

#include "AliEMCALRecParam.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#endif


void AliEMCALSetRecParamCDB(AliRecoParam::EventSpecie_t default = AliRecoParam::kDefault)
{

  // Create an object AliEMCALRecParam and store it to OCDB

  //Activate CDB storage
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");


  // Create reconstruction parameter object and set parameter values
  TObjArray* recParamArray = new TObjArray();

  {
    //default
    AliEMCALRecParam *recParamDB = AliEMCALRecParam::GetDefaultParameters();
    
    //Clusterization
    recParamDB->SetClusteringThreshold(0.5);
    recParamDB->SetW0(4.5);
    recParamDB->SetMinECut(0.45);
    recParamDB->SetUnfold(kFALSE);
    recParamDB->SetLocMaxCut(0.03);
    
    //Track matching
    recParamDB->SetTrkCutX(6.0);
    recParamDB->SetTrkCutY(6.0);
    recParamDB->SetTrkCutZ(6.0);
    recParamDB->SetTrkCutR(10.0);
    recParamDB->SetTrkCutAlphaMin(-50.0);
    recParamDB->SetTrkCutAlphaMax( 50.0);
    recParamDB->SetTrkCutNITS(3.0);
    recParamDB->SetTrkCutNTPC(20.0);
    recParamDB->SetTrkCutAngle(10000.0);      // i.e. exclude this for the moment
    
    //PID
    
    // as a first step, all array elements are initialized to 0.0
    Int_t i, j;
    for (i = 0; i < 6; i++) {
      for (j = 0; j < 6; j++) {
	
	recParamDB->SetGamma(i,j,0.);
	recParamDB->SetHadron(i,j,0.);
	recParamDB->SetPiZero5to10(i,j, 0.);
	recParamDB->SetPiZero10to60(i,j,0.);
      }
    } 
    
    recParamDB->SetGamma(0,0,  0.038022);
    recParamDB->SetGamma(0,1, -0.0001883);
    recParamDB->SetGamma(0,2,  5.449e-06);
    
    recParamDB->SetGamma(1,0,  0.207313);
    recParamDB->SetGamma(1,1, -0.000978);
    recParamDB->SetGamma(1,2,  0.00001634);
    
    recParamDB->SetGamma(2,0,  0.043364);
    recParamDB->SetGamma(2,1, -0.0002048);
    recParamDB->SetGamma(2,2,  8.661e-06);
    recParamDB->SetGamma(2,3, -1.353e-07);
    
    recParamDB->SetGamma(3,0,  0.265004);
    recParamDB->SetGamma(3,1,  0.061298);
    recParamDB->SetGamma(3,2, -0.003203);
    recParamDB->SetGamma(3,3,  4.73e-05);
    
    recParamDB->SetGamma(4,0,  0.243579);
    recParamDB->SetGamma(4,1, -1.614e-05);
    
    recParamDB->SetGamma(5,0,  0.002942);
    recParamDB->SetGamma(5,1, -3.976e-05);
    
    recParamDB->SetHadron(0,0,  0.011945 / 3.);
    recParamDB->SetHadron(0,1,  0.000386 / 3.);
    recParamDB->SetHadron(0,2, -0.000014 / 3.);
    recParamDB->SetHadron(0,3,  1.336e-07 / 3.);
    
    recParamDB->SetHadron(1,0,  0.496544);
    recParamDB->SetHadron(1,1, -0.003226);
    recParamDB->SetHadron(1,2,  0.00001678);
    
    recParamDB->SetHadron(2,0,  0.144838);
    recParamDB->SetHadron(2,1, -0.002954);
    recParamDB->SetHadron(2,2,  0.00008754);
    recParamDB->SetHadron(2,3, -7.587e-07);
    
    recParamDB->SetHadron(3,0,  1.264461 / 7.);
    recParamDB->SetHadron(3,1,  0.002097 / 7.);
    
    recParamDB->SetHadron(4,0,  0.261950);
    recParamDB->SetHadron(4,1, -0.001078);
    recParamDB->SetHadron(4,2,  0.00003237);
    recParamDB->SetHadron(4,3, -3.241e-07);
    recParamDB->SetHadron(4,4,  0.);
    recParamDB->SetHadron(4,5,  0.);
    recParamDB->SetHadron(5,0,  0.010317);
    recParamDB->SetHadron(5,1,  0.);
    recParamDB->SetHadron(5,2,  0.);
    recParamDB->SetHadron(5,3,  0.);
    recParamDB->SetHadron(5,4,  0.);
    recParamDB->SetHadron(5,5,  0.);
    
    recParamDB->SetPiZero5to10(0,0, 0.009138);
    recParamDB->SetPiZero5to10(0,1, 0.0006377);
    
    recParamDB->SetPiZero5to10(1,0, 0.08);
    
    recParamDB->SetPiZero5to10(2,0, -0.061119);
    recParamDB->SetPiZero5to10(2,1,  0.019013);
    
    recParamDB->SetPiZero5to10(3,0,  0.2);
    
    recParamDB->SetPiZero5to10(4,0,  0.252044);
    recParamDB->SetPiZero5to10(4,1, -0.002315);
    
    recParamDB->SetPiZero5to10(5,0,  0.002942);
    recParamDB->SetPiZero5to10(5,1, -3.976e-05);
    
    recParamDB->SetPiZero10to60(0,0,  0.009138);
    recParamDB->SetPiZero10to60(0,1,  0.0006377);
    
    recParamDB->SetPiZero10to60(1,0,  1.272837);
    recParamDB->SetPiZero10to60(1,1, -0.069708);
    recParamDB->SetPiZero10to60(1,2,  0.001568);
    recParamDB->SetPiZero10to60(1,3, -1.162e-05);
    
    recParamDB->SetPiZero10to60(2,0,  0.139703);
    recParamDB->SetPiZero10to60(2,1,  0.003687);
    recParamDB->SetPiZero10to60(2,2, -0.000568);
    recParamDB->SetPiZero10to60(2,3,  1.498e-05);
    recParamDB->SetPiZero10to60(2,4, -1.174e-07);
    
    recParamDB->SetPiZero10to60(3,0, -0.826367);
    recParamDB->SetPiZero10to60(3,1,  0.096951);
    recParamDB->SetPiZero10to60(3,2, -0.002215);
    recParamDB->SetPiZero10to60(3,3,  2.523e-05);
    
    recParamDB->SetPiZero10to60(4,0,  0.249890);
    recParamDB->SetPiZero10to60(4,1, -0.000063);
    
    recParamDB->SetPiZero10to60(5,0,  0.002942);
    recParamDB->SetPiZero10to60(5,1, -3.976e-05);
    
    // raw signal fitting
    recParamDB->SetHighLowGainFactor(16.);
    recParamDB->SetOrderParameter(2);
    recParamDB->SetTau(2.35);
    recParamDB->SetNoiseThreshold(3);
    recParamDB->SetNPedSamples(5);

    //Add to the recParamArray
    recParamDB->SetEventSpecie(AliRecoParam::kDefault);
    recParamArray->AddLast(recParamDB);
  }

  //Add other options here, if desired, for
  //Cosmic, LowMult and HighMult type events
  //and add them to the array


  {
    //For now, default is Pb+Pb, but let's add it again as
    //the "high mult" version too...
    AliEMCALRecParam *recParamDB = AliEMCALRecParam::GetDefaultParameters();

    recParamDB->SetEventSpecie(AliRecoParam::kHighMult);
    recParamArray->AddLast(recParamDB);
  }


  {
    //Low multiplicity parameter modifications:
    AliEMCALRecParam *recParamDB = AliEMCALRecParam::GetDefaultParameters();

    recParamDB->SetClusteringThreshold(0.2); // 200 MeV
    recParamDB->SetMinECut(0.01);  //10 MeV
    recParamDB->SetEventSpecie(AliRecoParam::kLowMult);
    recParamArray->AddLast(recParamDB);
    
  }


  //Set the default version in the array
  Bool_t defaultIsSet = kFALSE;
  for(Int_t i = 0; i < recParamArray->GetEntriesFast(); i++) {
    AliDetectorRecoParam* param = (AliDetectorRecoParam*)recParamArray->UncheckedAt(i);
    if(!param) continue;
    if(default & param->GetEventSpecie()) {
      param->SetAsDefault();
      defaultIsSet = kTRUE;
    }
  }

  if(!defaultIsSet) {
    AliError("The default reconstruction parameters are not set!  Exiting...");
    return;
  }

  // Store calibration data into database  
  AliCDBMetaData *md = new AliCDBMetaData();
  md->SetResponsible("J. Klay");
  md->SetComment("Reconstruction Parameters: EMCAL");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);
  
  AliCDBId id("EMCAL/Calib/RecoParam",0,AliCDBRunRange::Infinity());
  cdb->GetDefaultStorage()->Put(recParamArray, id, md);

  return;
}
