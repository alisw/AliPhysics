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
  gSystem->Load("libITSUpgradeRec.so");
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
  parData = new AliParamList(AliITSUSimulationPix::kParamStart+11); // 2 common + 9 params for double gaussian
  parData->SetUniqueID(0); // this is a function for detId=0
  parData->SetID(AliITSUSimulationPix::kSpreadDoubleGauss); // and uses double gaussian
  parData->SetNameTitle("Monopix_seg0","double gaussian for segmentation 0");
  //
  // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
  // injected one to consider
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,3,"nPixX");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,3,"nPixZ"); 
  // 
  // now set the parameters according selected function
  offs = AliITSUSimulationPix::kParamStart;
  parData->SetParameter(offs++,-0.1e-4  , "G1 Mean_x");
  parData->SetParameter(offs++, 8.125e-4, "G1 Sigma_x");
  parData->SetParameter(offs++, 2.011e-4, "G1 Mean_z"); 
  parData->SetParameter(offs++, 8.125e-4, "G1 Sigma_z"); 
  parData->SetParameter(offs++,-0.069e-4, "G2 Mean_x");
  parData->SetParameter(offs++,15.050e-4, "G2 Sigma_x");
  parData->SetParameter(offs++,-8.425e-4, "G2 Mean_z");
  parData->SetParameter(offs++,15.050e-4, "G2 Sigma_z"); 
  parData->SetParameter(offs++,3.904037/59.468672, "G2 A2/A1");  // scaling of 2nd gaussian amplitude wrt 1st one
  // 
  itsSimuParam->AddRespFunParam(parData);
  //
  //
  //------------------------ parameterization data for segmentation 1 ----------------------
  parData = new AliParamList(AliITSUSimulationPix::kParamStart+9); // 2 common + 9 params for double gaussian
  parData->SetUniqueID(1); // this is a function for detId=1
  parData->SetID(AliITSUSimulationPix::kSpreadDoubleGauss); // and uses double gaussian
  parData->SetNameTitle("Monopix_seg1","double gaussian for segmentation 1");
  //
  // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
  // injected one to consider
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,3,"nPixX");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,3,"nPixZ"); 
  // 
  // now set the parameters according selected function
  offs = AliITSUSimulationPix::kParamStart;
  parData->SetParameter(offs++,-0.1e-4  , "G1 Mean_x");
  parData->SetParameter(offs++, 8.125e-4, "G1 Sigma_x");
  parData->SetParameter(offs++, 2.011e-4, "G1 Mean_z"); 
  parData->SetParameter(offs++, 8.125e-4, "G1 Sigma_z"); 
  parData->SetParameter(offs++,-0.069e-4, "G2 Mean_x");
  parData->SetParameter(offs++,15.050e-4, "G2 Sigma_x");
  parData->SetParameter(offs++,-8.425e-4, "G2 Mean_z");
  parData->SetParameter(offs++,15.050e-4, "G2 Sigma_z"); 
  parData->SetParameter(offs++,3.904037/59.468672, "G2 A2/A1");  // scaling of 2nd gaussian amplitude wrt 1st one
  // 
  itsSimuParam->AddRespFunParam(parData);
  //
  //
  //------------------------ parameterization data for segmentation 2 ----------------------
  parData = new AliParamList(AliITSUSimulationPix::kParamStart+9); // 2 common + 9 params for double gaussian
  parData->SetUniqueID(2); // this is a function for detId=2
  parData->SetID(AliITSUSimulationPix::kSpreadDoubleGauss); // and uses double gaussian
  parData->SetNameTitle("Monopix_seg1","double gaussian for segmentation 1");
  //
  // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
  // injected one to consider
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,3,"nPixX");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,3,"nPixZ"); 
  // 
  // now set the parameters according selected function
  offs = AliITSUSimulationPix::kParamStart;
  parData->SetParameter(offs++,-0.1e-4  , "G1 Mean_x");
  parData->SetParameter(offs++, 8.125e-4, "G1 Sigma_x");
  parData->SetParameter(offs++, 2.011e-4, "G1 Mean_z"); 
  parData->SetParameter(offs++, 8.125e-4, "G1 Sigma_z"); 
  parData->SetParameter(offs++,-0.069e-4, "G2 Mean_x");
  parData->SetParameter(offs++,15.050e-4, "G2 Sigma_x");
  parData->SetParameter(offs++,-8.425e-4, "G2 Mean_z");
  parData->SetParameter(offs++,15.050e-4, "G2 Sigma_z"); 
  parData->SetParameter(offs++,3.904037/59.468672, "G2 A2/A1");  // scaling of 2nd gaussian amplitude wrt 1st one
  // 
  itsSimuParam->AddRespFunParam(parData);
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

