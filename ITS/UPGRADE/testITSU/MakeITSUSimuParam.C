void MakeITSUSimuParam(const char* cdbURI="local://") {
//========================================================================
//
// Steering macro for ITS simulation parameters
//
// Author: L.Molnar
// Contact: levente.molnar@cern.ch
//
//========================================================================

  const char* macroname = "MakeITSUSimuParam.C";
  //
  gSystem->Load("libITSUpgradeBase.so");
  gSystem->Load("libITSUpgradeSim.so");
  //
  // Activate CDB storage and load geometry from CDB
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(cdbURI);

  AliITSUSimuParam* itsSimuParam = new AliITSUSimuParam();
  //
  // Add spread function parameterization data
  AliParamList* parData = 0;
  int offs = 0;
  //
  //------------------------ parameterization data for segmentation 0 ----------------------
  parData = new AliParamList(AliITSUSimulationPix::kNG2Par); // 2 common + 9 params for double gaussian
  parData->SetUniqueID(0); // this is a function for detId=0
  parData->SetID(AliITSUSimulationPix::kSpreadFunDoubleGauss2D); // and uses double gaussian
  parData->SetNameTitle("Monopix_seg0","double gaussian for segmentation 0");
  //
  // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
  // injected one to consider
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,3,"nPixX");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,3,"nPixZ"); 
  //
  // now set the parameters according selected function
  parData->SetParameter(AliITSUSimulationPix::kG2MeanX0  , -0.1e-4  , "G1 Mean_x");
  parData->SetParameter(AliITSUSimulationPix::kG2SigX0   ,  8.125e-4, "G1 Sigma_x");
  parData->SetParameter(AliITSUSimulationPix::kG2MeanZ0  ,  2.011e-4, "G1 Mean_z"); 
  parData->SetParameter(AliITSUSimulationPix::kG2SigZ0   ,  8.125e-4, "G1 Sigma_z"); 
  parData->SetParameter(AliITSUSimulationPix::kG2MeanX1  , -0.069e-4, "G2 Mean_x");
  parData->SetParameter(AliITSUSimulationPix::kG2SigX1   , 15.050e-4, "G2 Sigma_x");
  parData->SetParameter(AliITSUSimulationPix::kG2MeanZ1  , -8.425e-4, "G2 Mean_z");
  parData->SetParameter(AliITSUSimulationPix::kG2SigZ1   , 15.050e-4, "G2 Sigma_z"); 
  // scaling of 2nd gaussian amplitude wrt 1st one
  parData->SetParameter(AliITSUSimulationPix::kG2ScaleG2 , 3.904037/59.468672, "G2 A2/A1");  
  // 
  itsSimuParam->AddRespFunParam(parData);
  //
  //
  //------------------------ parameterization data for segmentation 1 ----------------------
  parData = new AliParamList(AliITSUSimulationPix::kNG2Par); // 2 common + 9 params for double gaussian
  parData->SetUniqueID(1); // this is a function for detId=1
  parData->SetID(AliITSUSimulationPix::kSpreadFunDoubleGauss2D); // and uses double gaussian
  parData->SetNameTitle("Monopix_seg1","double gaussian for segmentation 1");
  //
  // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
  // injected one to consider
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,3,"nPixX");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,3,"nPixZ"); 
  // 
  // now set the parameters according selected function
  parData->SetParameter(AliITSUSimulationPix::kG2MeanX0  , -0.1e-4  , "G1 Mean_x");
  parData->SetParameter(AliITSUSimulationPix::kG2SigX0   ,  8.125e-4, "G1 Sigma_x");
  parData->SetParameter(AliITSUSimulationPix::kG2MeanZ0  ,  2.011e-4, "G1 Mean_z"); 
  parData->SetParameter(AliITSUSimulationPix::kG2SigZ0   ,  8.125e-4, "G1 Sigma_z"); 
  parData->SetParameter(AliITSUSimulationPix::kG2MeanX1  , -0.069e-4, "G2 Mean_x");
  parData->SetParameter(AliITSUSimulationPix::kG2SigX1   , 15.050e-4, "G2 Sigma_x");
  parData->SetParameter(AliITSUSimulationPix::kG2MeanZ1  , -8.425e-4, "G2 Mean_z");
  parData->SetParameter(AliITSUSimulationPix::kG2SigZ1   , 15.050e-4, "G2 Sigma_z"); 
  // scaling of 2nd gaussian amplitude wrt 1st one
  parData->SetParameter(AliITSUSimulationPix::kG2ScaleG2 , 3.904037/59.468672, "G2 A2/A1");  
  // 
  itsSimuParam->AddRespFunParam(parData);
  //
  //
  //------------------------ parameterization data for segmentation 2 ----------------------
  parData = new AliParamList(AliITSUSimulationPix::kNG1Par); // 2 common + 3 params for double gaussian
  parData->SetUniqueID(2); // this is a function for detId=2
  parData->SetID(AliITSUSimulationPix::kSpreadFunGauss2D); // and uses double gaussian
  parData->SetNameTitle("Monopix_seg2","single gaussian for segmentation 1");
  //
  // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
  // injected one to consider
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,3,"nPixX");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,3,"nPixZ"); 
  // 
  // now set the parameters according selected function
  parData->SetParameter(AliITSUSimulationPix::kG1MeanX  , -0.1e-4  , "Mean_x");
  parData->SetParameter(AliITSUSimulationPix::kG1SigX   ,  8.125e-4, "Sigma_x");
  parData->SetParameter(AliITSUSimulationPix::kG1MeanZ  ,  2.011e-4, "Mean_z"); 
  parData->SetParameter(AliITSUSimulationPix::kG1SigZ   ,  8.125e-4, "Sigma_z"); 
  // 
  itsSimuParam->AddRespFunParam(parData);
  //
  //
  // save in CDB storage
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("ITS Upgrade Project");
  md->SetComment("Simulation parameters for ITS Upgrade.");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);
  AliCDBId id("ITS/Calib/SimuParam",0,AliCDBRunRange::Infinity());
  cdb->GetDefaultStorage()->Put(itsSimuParam,id, md);
  //
  return;
}

