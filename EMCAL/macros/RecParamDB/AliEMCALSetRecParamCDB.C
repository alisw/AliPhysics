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


void AliEMCALSetRecParamCDB(AliRecoParam::EventSpecie_t default = AliRecoParam::kLowMult)
{
  
  // Create an object AliEMCALRecParam and store it to OCDB
  
  //Activate CDB storage
  AliCDBManager* cdb = AliCDBManager::Instance();
  if(!cdb->IsDefaultStorageSet()) cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  
  // Create reconstruction parameter object and set parameter values
  TObjArray* recParamArray = new TObjArray();
  
  {
    //default
    //AliEMCALRecParam *recParamDB = AliEMCALRecParam::GetDefaultParameters();
    AliEMCALRecParam *recParamDB = GetLowMultiplicityParameters();
    recParamDB->SetName("Default - p+p");
    recParamDB->SetTitle("Default - p+p");
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
    //AliEMCALRecParam *recParamDB = AliEMCALRecParam::GetHighFluxParam();
    AliEMCALRecParam *recParamDB = GetHighMultiplicityParameters();
    recParamDB->SetName("High Flux - Pb+Pb");
    recParamDB->SetTitle("High Flux - Pb+Pb");
    recParamDB->SetEventSpecie(AliRecoParam::kHighMult);
    recParamArray->AddLast(recParamDB);
  }
  
  {
    //Low multiplicity parameter modifications:
    //AliEMCALRecParam *recParamDB = AliEMCALRecParam::GetLowFluxParam();
    AliEMCALRecParam *recParamDB = GetLowMultiplicityParameters();
    recParamDB->SetName("Low Flux - p+p");
    recParamDB->SetTitle("Low Flux - p+p");
    recParamDB->SetEventSpecie(AliRecoParam::kLowMult);
    recParamArray->AddLast(recParamDB);
    
  }
  
  {
    //Cosmic parameter modifications (same as low multiplicity):
    //AliEMCALRecParam *recParamDB = AliEMCALRecParam::GetLowFluxParam();
    AliEMCALRecParam *recParamDB = GetLowMultiplicityParameters();
    recParamDB->SetName("Cosmic");
    recParamDB->SetTitle("Cosmic");
    recParamDB->SetEventSpecie(AliRecoParam::kCosmic);
    recParamArray->AddLast(recParamDB);
    
  }
  
  {
    //Calib parameter modifications (same as low multiplicity):
    //AliEMCALRecParam *recParamDB = AliEMCALRecParam::GetLowFluxParam();
    AliEMCALRecParam *recParamDB = GetLowMultiplicityParameters();
    recParamDB->SetName("Calibration - LED");
    recParamDB->SetTitle("Calibration - LED");
    recParamDB->SetEventSpecie(AliRecoParam::kCalib);
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

//-----------------------------------------------------------------------------
AliEMCALRecParam* GetHighMultiplicityParameters()
{
  //Set here the high flux/multiplicity ("Pb+Pb") parameters for the reconstruction
  //Right now it should be the same settings as with
  //AliEMCALRecParam *recParamDB = AliEMCALRecParam::GetHighFluxParam();
  //or
  //AliEMCALRecParam *recParamDB = AliEMCALRecParam::GetDefaultParameters();
  
  AliEMCALRecParam* params =  AliEMCALRecParam::GetDefaultParameters();
  //Clusterization

  // params->SetClusterizerFlag(AliEMCALRecParam::kClusterizerv1);
  // params->SetClusteringThreshold(0.5);
  // params->SetMinECut(0.45);

  params->SetClusterizerFlag(AliEMCALRecParam::kClusterizerv2);
  params->SetClusteringThreshold(0.1); // 100 MeV                                             
  params->SetMinECut(0.05);  //50 MeV    
   
  params->SetUnfold(kFALSE);

  params->SetW0(4.5);

  params->SetTimeCut(250e-9);//250 ns
  params->SetTimeMin(425e-9);//425 ns
  params->SetTimeMax(825e-9);//825 ns

  //Track matching

  params->SetTrkCutNITS(1.0);
  params->SetTrkCutNTPC(20.0);
  params->SetExtrapolateStep(20.);
  
  //PID
	
  // as a first step, all array elements are initialized to 0.0
  Int_t i, j;
  for (i = 0; i < 6; i++) {
    params->SetGammaEnergyProb(i,0.);
    params->SetHadronEnergyProb(i,0.); 
    params->SetPiZeroEnergyProb(i,0.);
    
    for (j = 0; j < 6; j++) {
      
      params->SetGamma(i,j,0.);
      params->SetHadron(i,j,0.);
      params->SetPiZero(i,j, 0.);
      params->SetGamma1to10(i,j,0.);
      params->SetHadron1to10(i,j,0.);
    }
  } 
  
  params->SetGamma(0,0, -7.656908e-01); 
  params->SetGamma(0,1,  2.352536e-01); 
  params->SetGamma(0,2,  1.555996e-02);
  params->SetGamma(0,3,  2.243525e-04);
  params->SetGamma(0,4, -2.560087e-06);
  
  params->SetGamma(1,0,  6.500216e+00);
  params->SetGamma(1,1, -2.564958e-01);
  params->SetGamma(1,2,  1.967894e-01);
  params->SetGamma(1,3, -3.982273e-04);
  params->SetGamma(1,4,  2.797737e-06);
  
  params->SetGamma(2,0,  2.416489e+00);
  params->SetGamma(2,1, -1.601258e-01);
  params->SetGamma(2,2,  3.126839e-02);
  params->SetGamma(2,3,  3.387532e-04);
  params->SetGamma(2,4, -4.089145e-06);
  
  params->SetGamma(3,0,  0.);
  params->SetGamma(3,1, -2.696008e+00);
  params->SetGamma(3,2,  6.920305e-01);
  params->SetGamma(3,3, -2.281122e-03);
  params->SetGamma(3,4,  0.);
  
  params->SetGamma(4,0,  2.281564e-01);
  params->SetGamma(4,1, -7.575040e-02);
  params->SetGamma(4,2,  3.813423e-01);
  params->SetGamma(4,3, -1.243854e-04);
  params->SetGamma(4,4,  1.232045e-06);
  
  params->SetGamma(5,0, -3.290107e-01);
  params->SetGamma(5,1,  3.707545e-02);
  params->SetGamma(5,2,  2.917397e-03);
  params->SetGamma(5,3,  4.695306e-05);
  params->SetGamma(5,4, -3.572981e-07);
  
  params->SetHadron(0,0,   1.519112e-01);
  params->SetHadron(0,1, -8.267603e-02);
  params->SetHadron(0,2,  1.914574e-02);
  params->SetHadron(0,3, -2.677921e-04);
  params->SetHadron(0,4,  5.447939e-06);
  
  params->SetHadron(1,0,  0.);
  params->SetHadron(1,1, -7.549870e-02); 
  params->SetHadron(1,2,  3.930087e-01);
  params->SetHadron(1,3, -2.368500e-03); 
  params->SetHadron(1,4,  0.);
  
  params->SetHadron(2,0,  0.);
  params->SetHadron(2,1, -2.463152e-02);
  params->SetHadron(2,2,  1.349257e-01);
  params->SetHadron(2,3, -1.089440e-03);
  params->SetHadron(2,4,  0.);
  
  params->SetHadron(3,0, 0.);
  params->SetHadron(3,1, 5.101560e-01);
  params->SetHadron(3,2, 1.458679e-01);
  params->SetHadron(3,3, 4.903068e-04);
  params->SetHadron(3,4, 0.);
  
  params->SetHadron(4,0, 0.);
  params->SetHadron(4,1, -6.693943e-03); 
  params->SetHadron(4,2,  2.444753e-01);
  params->SetHadron(4,3, -5.553749e-05);
  params->SetHadron(4,4, 0.);
  
  params->SetHadron(5,0, -4.414030e-01);
  params->SetHadron(5,1, 2.292277e-01);
  params->SetHadron(5,2, -2.433737e-02);
  params->SetHadron(5,3,  1.758422e-03);
  params->SetHadron(5,4, -3.001493e-05);
  
  params->SetPiZero(0,0,  5.072157e-01);
  params->SetPiZero(0,1, -5.352747e-01);
  params->SetPiZero(0,2,  8.499259e-02);
  params->SetPiZero(0,3, -3.687401e-03);
  params->SetPiZero(0,4,  5.482280e-05);
  
  params->SetPiZero(1,0,  4.590137e+02); 
  params->SetPiZero(1,1, -7.079341e+01);
  params->SetPiZero(1,2,  4.990735e+00);
  params->SetPiZero(1,3, -1.241302e-01);
  params->SetPiZero(1,4,  1.065772e-03);
  
  params->SetPiZero(2,0,  1.376415e+02);
  params->SetPiZero(2,1, -3.031577e+01);
  params->SetPiZero(2,2,  2.474338e+00);
  params->SetPiZero(2,3, -6.903410e-02);
  params->SetPiZero(2,4,  6.244089e-04);
  
  params->SetPiZero(3,0, 0.);
  params->SetPiZero(3,1,  1.145983e+00);
  params->SetPiZero(3,2, -2.476052e-01);
  params->SetPiZero(3,3,  1.367373e-02);
  params->SetPiZero(3,4, 0.);
  
  params->SetPiZero(4,0, -2.097586e+02);
  params->SetPiZero(4,1,  6.300800e+01);
  params->SetPiZero(4,2, -4.038906e+00);
  params->SetPiZero(4,3,  1.088543e-01);
  params->SetPiZero(4,4, -9.362485e-04);
  
  params->SetPiZero(5,0, -1.671477e+01); 
  params->SetPiZero(5,1,  2.995415e+00);
  params->SetPiZero(5,2, -6.040360e-02);
  params->SetPiZero(5,3, -6.137459e-04);
  params->SetPiZero(5,4,  1.847328e-05);
  
  // High flux ones pp 
  
  params->SetHadronEnergyProb(0, 0.);
  params->SetHadronEnergyProb(1, 0.);
  params->SetHadronEnergyProb(2,  6.188452e-02);
  params->SetHadronEnergyProb(3,  2.030230e+00);
  params->SetHadronEnergyProb(4, -6.402242e-02);
  
  // raw signal fitting
  params->SetHighLowGainFactor(16.);
  params->SetOrderParameter(2);
  params->SetTau(2.35);
  params->SetNoiseThreshold(3);
  params->SetNPedSamples(5);
  params->SetRemoveBadChannels(kFALSE);
  params->SetFittingAlgorithm(0);//(AliCaloConstants::kStandard);  
  params->SetFALTROUsage(kTRUE); 
  params->SetLEDFit(kFALSE);   

  return params ;
}	

//-----------------------------------------------------------------------------
AliEMCALRecParam* GetLowMultiplicityParameters()
{
  // Set here the low flux/multiplicity ("p+p") parameters for the reconstruction
  //Right now it should be the same settings as with
  //AliEMCALRecParam *recParamDB = AliEMCALRecParam::GetLowFluxParam();
  
  AliEMCALRecParam* params =  AliEMCALRecParam::GetDefaultParameters();
  //params->SetClusterizerFlag(AliEMCALRecParam::kClusterizerNxN);
  params->SetClusterizerFlag(AliEMCALRecParam::kClusterizerv1);
  params->SetClusteringThreshold(0.1); // 100 MeV                                             
  params->SetMinECut(0.05);  //50 MeV       	

  params->SetUnfold(kFALSE);

  params->SetW0(4.5);

  params->SetTimeCut(250e-9);//250 ns
  params->SetTimeMin(425e-9);//425 ns
  params->SetTimeMax(825e-9);//825 ns

  // Track Matching
  
  params->SetTrkCutNITS(1.0);
  params->SetTrkCutNTPC(20.0);
  params->SetExtrapolateStep(20.);  
  
  //PID parameters for pp  implemented 
  // as a first step, all array elements are initialized to 0.0
  Int_t i, j;
  for (i = 0; i < 6; i++) {
    for (j = 0; j < 6; j++) {
      params->SetGamma(i,j,0.);
      params->SetGamma1to10(i,j,0.);
      params->SetHadron(i,j,0.);
      params->SetHadron1to10(i,j,0.);
      params->SetPiZero(i,j,0.);
      
    }
    params->SetGammaEnergyProb(i,0.); // not yet implemented
    params->SetHadronEnergyProb(i,0.);
    params->SetPiZeroEnergyProb(i,0.); // not yet implemented
  }
  
  params->SetGamma(0,0, -7.656908e-01);
  params->SetGamma(0,1,  2.352536e-01);
  params->SetGamma(0,2,  1.555996e-02);
  params->SetGamma(0,3,  2.243525e-04);
  params->SetGamma(0,4, -2.560087e-06);
  
  params->SetGamma(1,0,  6.500216e+00);
  params->SetGamma(1,1, -2.564958e-01);
  params->SetGamma(1,2,  1.967894e-01);
  params->SetGamma(1,3, -3.982273e-04);
  params->SetGamma(1,4,  2.797737e-06);
  
  params->SetGamma(2,0,  2.416489e+00);
  params->SetGamma(2,1, -1.601258e-01);
  params->SetGamma(2,2,  3.126839e-02);
  params->SetGamma(2,3,  3.387532e-04);
  params->SetGamma(2,4, -4.089145e-06);
  
  params->SetGamma(3,0,0.);
  params->SetGamma(3,1,-2.696008e+00); 
  params->SetGamma(3,2, 6.920305e-01);
  params->SetGamma(3,3,-2.281122e-03);
  params->SetGamma(3,4,0.);
  
  params->SetGamma(4,0,  2.281564e-01); 
  params->SetGamma(4,1, -7.575040e-02);
  params->SetGamma(4,2,  3.813423e-01);
  params->SetGamma(4,3, -1.243854e-04);
  params->SetGamma(4,4,  1.232045e-06);
  
  params->SetGamma(5,0, -3.290107e-01);
  params->SetGamma(5,1,  3.707545e-02);
  params->SetGamma(5,2,  2.917397e-03);
  params->SetGamma(5,3,  4.695306e-05);
  params->SetGamma(5,4, -3.572981e-07);
  
  params->SetHadron(0,0,  9.482243e-01); 
  params->SetHadron(0,1, -2.780896e-01);
  params->SetHadron(0,2,  2.223507e-02);
  params->SetHadron(0,3,  7.294263e-04);
  params->SetHadron(0,4, -5.665872e-06); 
  
  params->SetHadron(1,0,  0.);
  params->SetHadron(1,1,  0.);
  params->SetHadron(1,2,  2.483298e-01);
  params->SetHadron(1,3,  0.);
  params->SetHadron(1,4,  0.);
  
  params->SetHadron(2,0, -5.601199e+00);
  params->SetHadron(2,1,  2.097382e+00);
  params->SetHadron(2,2, -2.307965e-01);
  params->SetHadron(2,3,  9.206871e-03);
  params->SetHadron(2,4, -8.887548e-05);
  
  params->SetHadron(3,0,  6.543101e+00);
  params->SetHadron(3,1, -2.305203e+00);
  params->SetHadron(3,2,  2.761673e-01);
  params->SetHadron(3,3, -5.465855e-03);
  params->SetHadron(3,4,  2.784329e-05);
  
  params->SetHadron(4,0, -2.443530e+01);
  params->SetHadron(4,1,  8.902578e+00);
  params->SetHadron(4,2, -5.265901e-01);
  params->SetHadron(4,3,  2.549111e-02);
  params->SetHadron(4,4, -2.196801e-04);
  
  params->SetHadron(5,0,  2.102007e-01);
  params->SetHadron(5,1, -3.844418e-02);
  params->SetHadron(5,2,  1.234682e-01);
  params->SetHadron(5,3, -3.866733e-03);
  params->SetHadron(5,4,  3.362719e-05);
  
  params->SetPiZero(0,0,  5.07215e-01);
  params->SetPiZero(0,1, -5.35274e-01);
  params->SetPiZero(0,2,  8.49925e-02);
  params->SetPiZero(0,3, -3.68740e-03);
  params->SetPiZero(0,4,  5.48228e-05);
  
  params->SetPiZero(1,0,  4.590137e+02);
  params->SetPiZero(1,1, -7.079341e+01);
  params->SetPiZero(1,2,  4.990735e+00);
  params->SetPiZero(1,3, -1.241302e-01);
  params->SetPiZero(1,4,  1.065772e-03);
  
  params->SetPiZero(2,0,  1.376415e+02); 
  params->SetPiZero(2,1, -3.031577e+01);
  params->SetPiZero(2,2,  2.474338e+00);
  params->SetPiZero(2,3, -6.903410e-02); 
  params->SetPiZero(2,4,  6.244089e-04);
  
  params->SetPiZero(3,0,  0.);
  params->SetPiZero(3,1,  1.145983e+00);
  params->SetPiZero(3,2, -2.476052e-01);
  params->SetPiZero(3,3,  1.367373e-02);
  params->SetPiZero(3,4,  0.);
  
  params->SetPiZero(4,0, -2.09758e+02);
  params->SetPiZero(4,1,  6.30080e+01);
  params->SetPiZero(4,2, -4.03890e+00);
  params->SetPiZero(4,3,  1.08854e-01);
  params->SetPiZero(4,4, -9.36248e-04);
  
  params->SetPiZero(5,0, -1.671477e+01);
  params->SetPiZero(5,1,  2.995415e+00);
  params->SetPiZero(5,2, -6.040360e-02);
  params->SetPiZero(5,3,  -6.137459e-04);
  params->SetPiZero(5,4, 1.847328e-05);
  
  //     params->SetHadronEnergyProb(0,0.);
  //     params->SetHadronEnergyProb(1,0.);
  //     params->SetHadronEnergyProb(2,1.);
  //     params->SetHadronEnergyProb(3,0.);
  //     params->SetHadronEnergyProb(4,0.);
  
  params->SetHadronEnergyProb(0,  4.767543e-02);
  params->SetHadronEnergyProb(1, -1.537523e+00);
  params->SetHadronEnergyProb(2,  2.956727e-01);
  params->SetHadronEnergyProb(3, -3.051022e+01);
  params->SetHadronEnergyProb(4, -6.036931e-02);

  // raw signal fitting
  params->SetHighLowGainFactor(16.);
  params->SetOrderParameter(2);
  params->SetTau(2.35);
  params->SetNoiseThreshold(3);
  params->SetNPedSamples(5);
  params->SetRemoveBadChannels(kFALSE);
  params->SetFittingAlgorithm(0);//(AliCaloConstants::kStandard);  
  params->SetFALTROUsage(kTRUE); 
  params->SetLEDFit(kFALSE);   

  return params;
  
}


