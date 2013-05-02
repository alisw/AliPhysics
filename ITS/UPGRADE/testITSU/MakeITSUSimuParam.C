//RS: before running MakeITSUSimuParam call ".x LoadLibs.C"

// Reference sensor parameterizations
AliITSUSimuParam*  MakeITSUSimuParam_M32P26Func(Double_t qscale, Double_t sncut); // don't forget to turn off diode shift for P26 in Config.C
AliITSUSimuParam*  MakeITSUSimuParam_M32P26Map(Double_t qscale, Double_t sncut);// don't forget to turn off diode shift for P26 in Config.C
AliITSUSimuParam*  MakeITSUSimuParam_M32terP31Map(Double_t qscale, Double_t sncut);
// Irradiated sensor parameterizations P31
AliITSUSimuParam*  MakeITSUSimuParam_M32terP31Map1MRad(Double_t qscale, Double_t sncut);
AliITSUSimuParam*  MakeITSUSimuParam_M32terP31Map300kRad3e12(Double_t qscale, Double_t sncut);
AliITSUSimuParam*  MakeITSUSimuParam_M32terP31Map1MRad1e13(Double_t qscale, Double_t sncut);
// Irradiated sensor parameterizations P26
AliITSUSimuParam*  MakeITSUSimuParam_M32P26Map300kRad3e12(Double_t qscale, Double_t sncut);
AliITSUSimuParam*  MakeITSUSimuParam_M32P26Map1MRad1e13(Double_t qscale, Double_t sncut);
AliITSUSimuParam*  MakeITSUSimuParam_M32P26Map300kRad3e12_20deg(Double_t qscale, Double_t sncut);
AliITSUSimuParam*  MakeITSUSimuParam_M32P26Map1MRad1e13_20deg(Double_t qscale, Double_t sncut);
void SetPSFParams(TString pixType, AliITSUParamList* parData);

// To turn noise generation ON set these values to 1
const int kAddNoise = -1;
const int kAddNoiseInAllMod = -1;

const char* inpPSFName = "$ALICE_ROOT/ITS/UPGRADE/misc/ITSU_pixel_response_PSFs.root";
const int kNLayers = 7; 

/*
  these are readout phases settings: 
  the rule is: 
  1) abs(kROShifts)<1: the modules of the layer are synchronized, the layer 
     phase is set to kROShifts*ROCycleLength of the modules of this layer
  2) abs(kROShifts)>1: each module within the layer will have random phase within its ROCycleLength
*/
const float kROShifts[kNLayers] = {0.5,0.5,0.5, -0.5,-0.5, 0.5,0.5};

void MakeITSUSimuParam(Double_t qscale=1.5, Double_t sncut=5.0,const char* cdbURI="local://") {
  //========================================================================
  //
  // Steering macro for ITS simulation parameters
  //
  // Author: L.Molnar
  // Contact: levente.molnar@cern.ch
  //
  //========================================================================
  AliITSUSimuParam *param = 0;
  //
  // Select only one parameterziation... and don't forget to set 18 um thickness in Config.C !!!
  
  //____ MIMOSA32 P26 Response parameterzied by fit functions
  //param = MakeITSUSimuParam_M32P26Func(qscale,sncut);
  
  //____ MIMOSA32 P26 Response parameterzied by map
  param = MakeITSUSimuParam_M32P26Map(qscale,sncut);

  //____ MIMOSA32ter P31 Response parameterzied by map //suggested!!! 
  //param = MakeITSUSimuParam_M32terP31Map(qscale,sncut);
  
  //____ MIMOSA32ter P31 Response parameterzied by map 1MRad irradiation
  //param = MakeITSUSimuParam_M32terP31Map1MRad(qscale,sncut);
  
  //____ MIMOSA32ter P31 Response parameterzied by map , 300kRad + 3e12 neq/cm2 irradiation
  //param = MakeITSUSimuParam_M32terP31Map300kRad3e12(qscale,sncut);
    
  //____ MIMOSA32ter P31 Response parameterzied by map , 1MRad+ 1e13 neq/cm2 irradiation
  //param = MakeITSUSimuParam_M32terP31Map1MRad1e13(qscale,sncut);
  
  
  //___ MIMOSA32 P26 , 300kRad + 3e12 neq/cm2 irradiation 30 deg
  // param = MakeITSUSimuParam_M32P26Map300kRad3e12(qscale,sncut);
  
  //____ MIMOSA32 P26 , 300kRad + 3e12 neq/cm2 irradiation 20 deg
  //param = MakeITSUSimuParam_M32P26Map300kRad3e12_20deg(qscale,sncut);
    
  //___ MIMOSA32 P26 ,  1MRad+ 1e13 neq/cm2 irradiation 30 deg
  //param = MakeITSUSimuParam_M32P26Map1MRad1e13(qscale,sncut);
    
  //___ MIMOSA32 P26 ,  1MRad+ 1e13 neq/cm2 irradiation 30 deg
  //param = MakeITSUSimuParam_M32P26Map1MRad1e13_20deg(qscale,sncut);
    
  
  
  param->Print();
  //
  // ----------------------------------------------------------
  // save in CDB storage
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetDefaultStorage(cdbURI);
  //
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("ITS Upgrade Project");
  md->SetComment("Simulation parameters for ITS Upgrade.");
  md->SetAliRootVersion(gSystem->Getenv("ARVERSION"));
  md->SetBeamPeriod(0);
  AliCDBId id("ITS/Calib/SimuParam",0,AliCDBRunRange::Infinity());
  cdb->GetDefaultStorage()->Put(param,id, md);
  //  
}

//________________________________________________________________________________________________________//
//   ||                                                                                                   //
//   || Paremeterization by fit functions of the MIMOSA32 P26 chip, 30 C, no irradiation                  //
//   ||                                                                                                   //
//  \||/                                                                                                  //
//___\/___________________________________________________________________________________________________//
AliITSUSimuParam* MakeITSUSimuParam_M32P26Func(Double_t qscale, Double_t sncut) 
{
  //const char* macroname = "MakeITSUSimuParam.C";
  //
  AliITSUSimuParam* itsSimuParam = new AliITSUSimuParam();
  //
  itsSimuParam->SetNLayers(kNLayers);
  for (int ilr=kNLayers;ilr--;) itsSimuParam->SetLrROCycleShift(kROShifts[ilr],ilr);
  //
  // Add spread function parameterization data
  AliITSUParamList* parData = 0;
  //  
  //------------------------ parameterization data for segmentation 0 ----------------------
  parData = new AliITSUParamList(AliITSUSimulationPix::kNG2Par); // 4 common + 9 params for double gaussian
  parData->SetUniqueID(0); // this is a function for detId=0
  //
  // and uses double gaussian for charge spread parameterization
  parData->SetNameTitle("Monopix_seg0","double gaussian for segmentation 0");
  parData->SetParameter(AliITSUSimulationPix::kChargeSpreadType,AliITSUSimulationPix::kSpreadFunDoubleGauss2D,"ChargeSpread"); 
  //
  // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
  // injected one to consider
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
  parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");  
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");

  parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
  parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
  // Noise
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,17.53,"pixNoiseMPV");
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,2.93,"pixNoiseSigma");  
 // and readout timing scheme
  parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
  parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
  //
  // now set the parameters according selected function
  parData->SetParameter(AliITSUSimulationPix::kG2MeanX0  , -5.63484e-01 * 1e-4, "G1 Mean_x");
  parData->SetParameter(AliITSUSimulationPix::kG2SigX0   ,  2.49464e+01 * 1e-4, "G1 Sigma_x");
  parData->SetParameter(AliITSUSimulationPix::kG2MeanZ0  , -2.76353e+00 * 1e-4, "G1 Mean_z"); 
  parData->SetParameter(AliITSUSimulationPix::kG2SigZ0   ,  2.59449e+01 * 1e-4, "G1 Sigma_z"); 
  parData->SetParameter(AliITSUSimulationPix::kG2MeanX1  ,  5.43664e-01 * 1e-4, "G2 Mean_x");
  parData->SetParameter(AliITSUSimulationPix::kG2SigX1   ,  7.97169e+00 * 1e-4, "G2 Sigma_x");
  parData->SetParameter(AliITSUSimulationPix::kG2MeanZ1  ,  1.76857e+00 * 1e-4, "G2 Mean_z");
  parData->SetParameter(AliITSUSimulationPix::kG2SigZ1   ,  1.01543e+01 * 1e-4, "G2 Sigma_z"); 
  // scaling of 2nd gaussian amplitude wrt 1st one
  parData->SetParameter(AliITSUSimulationPix::kG2ScaleG2 , 3.904037/59.468672, "G2 A2/A1");  
  // 
  itsSimuParam->AddRespFunParam(parData);
  //
  //
  //------------------------ parameterization data for segmentation 1 ----------------------
  parData = new AliITSUParamList(AliITSUSimulationPix::kNG2Par); // 4 common + 9 params for double gaussian
  parData->SetUniqueID(1); // this is a function for detId=1

  // and uses double gaussian for charge spread parameterization
  parData->SetNameTitle("Monopix_seg1","double gaussian for segmentation 1");
  parData->SetParameter(AliITSUSimulationPix::kChargeSpreadType,AliITSUSimulationPix::kSpreadFunDoubleGauss2D,"ChargeSpread"); 
  //
  // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
  // injected one to consider
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
  parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");  
  parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
  parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
  // Noise
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");

  parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,17.53,"pixNoiseMPV");
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,2.93,"pixNoiseSigma");  
 // and readout timing scheme
  parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
  parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
  //
  // now set the parameters according selected function
  parData->SetParameter(AliITSUSimulationPix::kG2MeanX0  , -5.63484e-01 * 1e-4, "G1 Mean_x");
  parData->SetParameter(AliITSUSimulationPix::kG2SigX0   ,  2.49464e+01 * 1e-4, "G1 Sigma_x");
  parData->SetParameter(AliITSUSimulationPix::kG2MeanZ0  , -2.76353e+00 * 1e-4, "G1 Mean_z"); 
  parData->SetParameter(AliITSUSimulationPix::kG2SigZ0   ,  2.59449e+01 * 1e-4, "G1 Sigma_z"); 
  parData->SetParameter(AliITSUSimulationPix::kG2MeanX1  ,  5.43664e-01 * 1e-4, "G2 Mean_x");
  parData->SetParameter(AliITSUSimulationPix::kG2SigX1   ,  7.97169e+00 * 1e-4, "G2 Sigma_x");
  parData->SetParameter(AliITSUSimulationPix::kG2MeanZ1  ,  1.76857e+00 * 1e-4, "G2 Mean_z");
  parData->SetParameter(AliITSUSimulationPix::kG2SigZ1   ,  1.01543e+01 * 1e-4, "G2 Sigma_z"); 
  // scaling of 2nd gaussian amplitude wrt 1st one
  parData->SetParameter(AliITSUSimulationPix::kG2ScaleG2 , 3.904037/59.468672, "G2 A2/A1");  
  // 
  itsSimuParam->AddRespFunParam(parData);
  //
  //
  //------------------------ parameterization data for segmentation 2 ----------------------
  parData = new AliITSUParamList(AliITSUSimulationPix::kNG2Par); // 4 common + 9 params for double gaussian
  parData->SetUniqueID(2); // this is a function for detId=2
  //
  parData->SetNameTitle("Monopix_seg2","double gaussian for segmentation 2");
  // and uses double gaussian for charge spread parameterization
  parData->SetParameter(AliITSUSimulationPix::kChargeSpreadType,AliITSUSimulationPix::kSpreadFunDoubleGauss2D,"ChargeSpread"); 
  //
  // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
  // injected one to consider
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
  parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");  
  parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
  parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
  // Noise
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");

  parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,17.53,"pixNoiseMPV");
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,2.93,"pixNoiseSigma");  
 // and readout timing scheme
  parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
  parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
  //
  // now set the parameters according selected function
  parData->SetParameter(AliITSUSimulationPix::kG2MeanX0  , -5.63484e-01 * 1e-4, "G1 Mean_x");
  parData->SetParameter(AliITSUSimulationPix::kG2SigX0   ,  2.49464e+01 * 1e-4, "G1 Sigma_x");
  parData->SetParameter(AliITSUSimulationPix::kG2MeanZ0  , -2.76353e+00 * 1e-4, "G1 Mean_z"); 
  parData->SetParameter(AliITSUSimulationPix::kG2SigZ0   ,  2.59449e+01 * 1e-4, "G1 Sigma_z"); 
  parData->SetParameter(AliITSUSimulationPix::kG2MeanX1  ,  5.43664e-01 * 1e-4, "G2 Mean_x");
  parData->SetParameter(AliITSUSimulationPix::kG2SigX1   ,  7.97169e+00 * 1e-4, "G2 Sigma_x");
  parData->SetParameter(AliITSUSimulationPix::kG2MeanZ1  ,  1.76857e+00 * 1e-4, "G2 Mean_z");
  parData->SetParameter(AliITSUSimulationPix::kG2SigZ1   ,  1.01543e+01 * 1e-4, "G2 Sigma_z"); 
  // scaling of 2nd gaussian amplitude wrt 1st one
  parData->SetParameter(AliITSUSimulationPix::kG2ScaleG2 , 3.904037/59.468672, "G2 A2/A1");  
  // 
  itsSimuParam->AddRespFunParam(parData);
  //
  //
  return itsSimuParam;
  
}
      
//________________________________________________________________________________________________________//
//   ||                                                                                                   //
//   || Paremeterization by charge spread map of the MIMOSA32 P26 chip, 30 C, no irradiation              //
//   ||                                                                                                   //
//  \||/                                                                                                  //
//___\/___________________________________________________________________________________________________//
AliITSUSimuParam* MakeITSUSimuParam_M32P26Map(Double_t qscale, Double_t sncut) 
{
  //const char* macroname = "MakeITSUSimuParam.C";
  //
  AliITSUSimuParam* itsSimuParam = new AliITSUSimuParam();
  //
  itsSimuParam->SetNLayers(kNLayers);
  for (int ilr=kNLayers;ilr--;) itsSimuParam->SetLrROCycleShift(kROShifts[ilr],ilr);
  //
  // Add spread function parameterization data
  AliITSUParamList* parData = 0;
  //
  //------------------------ parameterization data for segmentation 0 ----------------------
  parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
  parData->SetUniqueID(0);                                              // this is a function for detId=0
  parData->SetNameTitle("Monopix_seg1","PSF map for M32P26");
  SetPSFParams("hProfPSD_M32P26",parData);
  //
  // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
  // injected one to consider
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ"); 
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
  parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");  
  parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
  parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
  // Noise
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");

  parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,17.53,"pixNoiseMPV");
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,2.93,"pixNoiseSigma");  
 // and readout timing scheme
  parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
  parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
  // 
  itsSimuParam->AddRespFunParam(parData);
  //
  //
   //------------------------ parameterization data for segmentation 1 ----------------------
  parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
  parData->SetUniqueID(1);                                              // this is a function for detId=1
  parData->SetNameTitle("Monopix_seg1","PSF map for M32P26");
  SetPSFParams("hProfPSD_M32P26",parData);
  //
  // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
  // injected one to consider
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ"); 
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
  parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");  
  parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
  parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
  // Noise
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");

  parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,17.53,"pixNoiseMPV");
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,2.93,"pixNoiseSigma");  
 // and readout timing scheme
  parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
  parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
  //
  itsSimuParam->AddRespFunParam(parData);
  //
  //     
  //------------------------ parameterization data for segmentation 2 ----------------------
  parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
  parData->SetUniqueID(2);                                              // this is a function for detId=2
  parData->SetNameTitle("Monopix_seg2","PSF map for M32P26");
  SetPSFParams("hProfPSD_M32P26",parData);
  //
  // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
  // injected one to consider
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ"); 
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
  parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate"); 
  parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
  parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
  // Noise
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");

  parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,17.53,"pixNoiseMPV");
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,2.93,"pixNoiseSigma");  
  // and readout timing scheme
  parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
  parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
  //
  itsSimuParam->AddRespFunParam(parData);
  //
  return itsSimuParam;
}
    
//________________________________________________________________________________________________________//
//   ||                                                                                                   //
//   || Paremeterization by charge spread map of the MIMOSA32ter P31 chip, 30 C, no irradiation           //
//   ||                                                                                                   //
//  \||/                                                                                                  //
//___\/___________________________________________________________________________________________________//

AliITSUSimuParam* MakeITSUSimuParam_M32terP31Map(Double_t qscale, Double_t sncut) 
{  
  //const char* macroname = "MakeITSUSimuParam.C";
  //
  AliITSUSimuParam* itsSimuParam = new AliITSUSimuParam();
  //
  itsSimuParam->SetNLayers(kNLayers);
  for (int ilr=kNLayers;ilr--;) itsSimuParam->SetLrROCycleShift(kROShifts[ilr],ilr);
  //
  // Add spread function parameterization data
  AliITSUParamList* parData = 0;
  // 
  //------------------------ parameterization data for segmentation 0 ----------------------
  parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
  parData->SetUniqueID(0);                                              // this is a function for detId=0
  parData->SetNameTitle("Monopix_seg2","PSF map for M32terP31");
  SetPSFParams("hProfPSD_M32terP31",parData);
  //
  // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
  // injected one to consider
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ"); 
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
  parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");  
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");

  parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
  parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
  // Noise
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,20.62,"pixNoiseMPV");
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,3.55,"pixNoiseSigma");  
  // and readout timing scheme
  parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
  parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
  // 
  itsSimuParam->AddRespFunParam(parData);
  //
  //
  // 
  //------------------------ parameterization data for segmentation 1 ----------------------
  parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
  parData->SetUniqueID(1);                                              // this is a function for detId=1
  parData->SetNameTitle("Monopix_seg1","PSF map for M32terP31");
  SetPSFParams("hProfPSD_M32terP31",parData);

  //
  // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
  // injected one to consider
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ"); 
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
  parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");

  parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
  parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
   // Noise
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,20.62,"pixNoiseMPV");
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,3.55,"pixNoiseSigma");  
  // and readout timing scheme
  parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
  parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
  // 
  // now set the parameters according selected function
  itsSimuParam->AddRespFunParam(parData);
  //
  //
  // 
  //------------------------ parameterization data for segmentation 2 ----------------------
  parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
  parData->SetUniqueID(2);                                              // this is a function for detId=2
  parData->SetNameTitle("Monopix_seg2","PSF map for M32terP31");
  SetPSFParams("hProfPSD_M32terP31", parData );
  //
  // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
  // injected one to consider
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ"); 
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
  parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale"); //980./1080.
  parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");
  parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
  parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
   // Noise
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,20.62,"pixNoiseMPV");
  parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,3.55,"pixNoiseSigma");  
// and readout timing scheme
  parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
  parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)"); // cycle lenght is 10-30 us --> 10-30e-6 s
  //
  itsSimuParam->AddRespFunParam(parData);
  //
  return itsSimuParam;
}

//________________________________________________________________________________________________________//
//   ||                                                                                                   //
//   || Paremeterization by charge spread map of the MIMOSA32ter P31 chip, 30 C, 1 MRad irradiation       //
//   ||                                                                                                   //
//  \||/                                                                                                  //
//___\/___________________________________________________________________________________________________//
AliITSUSimuParam* MakeITSUSimuParam_M32terP31Map1MRad(Double_t qscale, Double_t sncut)
{
    //const char* macroname = "MakeITSUSimuParam.C";
    //
    AliITSUSimuParam* itsSimuParam = new AliITSUSimuParam();
    //
    itsSimuParam->SetNLayers(kNLayers);
    for (int ilr=kNLayers;ilr--;) itsSimuParam->SetLrROCycleShift(kROShifts[ilr],ilr);
    //
    // Add spread function parameterization data
    AliITSUParamList* parData = 0;
    //
    //------------------------ parameterization data for segmentation 0 ----------------------
    parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(0);                                              // this is a function for detId=0
    parData->SetNameTitle("Monopix_seg2","PSF map for M32terP31");
    SetPSFParams("hProfPSD_M32terP31_1MRad",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
    parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,21.843706,"pixNoiseMPV");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,3.417494,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //
    //------------------------ parameterization data for segmentation 1 ----------------------
    parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(1);                                              // this is a function for detId=1
    parData->SetNameTitle("Monopix_seg1","PSF map for M32terP31");
    SetPSFParams("hProfPSD_M32terP31_1MRad",parData);
    
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
    parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,21.843706,"pixNoiseMPV");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,3.417494,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    // now set the parameters according selected function
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //
    //------------------------ parameterization data for segmentation 2 ----------------------
    parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(2);                                              // this is a function for detId=2
    parData->SetNameTitle("Monopix_seg2","PSF map for M32terP31");
    SetPSFParams("hProfPSD_M32terP31_1MRad", parData );
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale"); //980./1080.
    parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,21.843706,"pixNoiseMPV");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,3.417494,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)"); // cycle lenght is 10-30 us --> 10-30e-6 s
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    return itsSimuParam;
}
//________________________________________________________________________________________________________//
//   ||                                                                                                   //
//   || Paremeterization by charge spread map of the MIMOSA32ter P31 chip, 30 C,                          //
//   ||          irradiation 300 kRad + 3e12 neq /cm2                                                     //  \||/                                                                                                  //
//___\/___________________________________________________________________________________________________//


AliITSUSimuParam* MakeITSUSimuParam_M32terP31Map300kRad3e12(Double_t qscale, Double_t sncut)
{
    //const char* macroname = "MakeITSUSimuParam.C";
    //
    AliITSUSimuParam* itsSimuParam = new AliITSUSimuParam();
    //
    itsSimuParam->SetNLayers(kNLayers);
    for (int ilr=kNLayers;ilr--;) itsSimuParam->SetLrROCycleShift(kROShifts[ilr],ilr);
    //
    // Add spread function parameterization data
    AliITSUParamList* parData = 0;
    //
    //------------------------ parameterization data for segmentation 0 ----------------------
    parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(0);                                              // this is a function for detId=0
    parData->SetNameTitle("Monopix_seg2","PSF map for M32terP31");
    SetPSFParams("hProfPSD_M32terP31_300kRad3e12neq",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
    parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,27.274548,"pixNoiseMPV");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,2.723949,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //
    //------------------------ parameterization data for segmentation 1 ----------------------
    parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(1);                                              // this is a function for detId=1
    parData->SetNameTitle("Monopix_seg1","PSF map for M32terP31");
    SetPSFParams("hProfPSD_M32terP31_300kRad3e12neq",parData);
    
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
    parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,27.274548,"pixNoiseMPV");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,2.723949,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    // now set the parameters according selected function
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //
    //------------------------ parameterization data for segmentation 2 ----------------------
    parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(2);                                              // this is a function for detId=2
    parData->SetNameTitle("Monopix_seg2","PSF map for M32terP31");
    SetPSFParams("hProfPSD_M32terP31_300kRad3e12neq", parData );
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale"); //980./1080.
    parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,27.274548,"pixNoiseMPV");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,2.723949,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)"); // cycle lenght is 10-30 us --> 10-30e-6 s
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    return itsSimuParam;
}

//________________________________________________________________________________________________________//
//   ||                                                                                                   //
//   || Paremeterization by charge spread map of the MIMOSA32ter P31 chip, 30 C,                          //
//   ||          irradiation 1MRad + 1e13 neq /cm2                                                        //  \||/                                                                                                  //
//___\/___________________________________________________________________________________________________//


AliITSUSimuParam* MakeITSUSimuParam_M32terP31Map1MRad1e13(Double_t qscale, Double_t sncut)
{
    //const char* macroname = "MakeITSUSimuParam.C";
    //
    AliITSUSimuParam* itsSimuParam = new AliITSUSimuParam();
    //
    itsSimuParam->SetNLayers(kNLayers);
    for (int ilr=kNLayers;ilr--;) itsSimuParam->SetLrROCycleShift(kROShifts[ilr],ilr);
    //
    // Add spread function parameterization data
    AliITSUParamList* parData = 0;
    //
    //------------------------ parameterization data for segmentation 0 ----------------------
    parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(0);                                              // this is a function for detId=0
    parData->SetNameTitle("Monopix_seg2","PSF map for M32terP31");
    SetPSFParams("hProfPSD_M32terP31_1MRad1e13neq",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
    parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,43.840798,"pixNoiseMPV");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,3.187048,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //
    //------------------------ parameterization data for segmentation 1 ----------------------
    parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(1);                                              // this is a function for detId=1
    parData->SetNameTitle("Monopix_seg1","PSF map for M32terP31");
    SetPSFParams("hProfPSD_M32terP31_1MRad1e13neq",parData);
    
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
    parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,43.840798,"pixNoiseMPV");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,3.187048,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    // now set the parameters according selected function
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //
    //------------------------ parameterization data for segmentation 2 ----------------------
    parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(2);                                              // this is a function for detId=2
    parData->SetNameTitle("Monopix_seg2","PSF map for M32terP31");
    SetPSFParams("hProfPSD_M32terP31_1MRad1e13neq", parData );
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale"); //980./1080.
    parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,43.840798,"pixNoiseMPV");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,3.187048,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)"); // cycle lenght is 10-30 us --> 10-30e-6 s
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    return itsSimuParam;
}



//________________________________________________________________________________________________________//
//   ||                                                                                                   //
//   || Paremeterization by charge spread map of the MIMOSA32 P26 chip, 30 C,                             //
//   ||    300kRad + 3e12 new/cm2                                                                         //
//  \||/                                                                                                  //
//___\/___________________________________________________________________________________________________//

AliITSUSimuParam* MakeITSUSimuParam_M32P26Map300kRad3e12(Double_t qscale, Double_t sncut)
{
    //const char* macroname = "MakeITSUSimuParam.C";
    //
    AliITSUSimuParam* itsSimuParam = new AliITSUSimuParam();
    //
    itsSimuParam->SetNLayers(kNLayers);
    for (int ilr=kNLayers;ilr--;) itsSimuParam->SetLrROCycleShift(kROShifts[ilr],ilr);
    //
    // Add spread function parameterization data
    AliITSUParamList* parData = 0;
    //
    //------------------------ parameterization data for segmentation 0 ----------------------
    parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(0);                                              // this is a function for detId=0
    parData->SetNameTitle("Monopix_seg1","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26_300kRad3e12neq",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
    parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,18.566901,"pixNoiseMPV");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,2.665233,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //------------------------ parameterization data for segmentation 1 ----------------------
    parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(1);                                              // this is a function for detId=1
    parData->SetNameTitle("Monopix_seg1","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26_300kRad3e12neq",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
    parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,18.566901,"pixNoiseMPV");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,2.665233,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //------------------------ parameterization data for segmentation 2 ----------------------
    parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(2);                                              // this is a function for detId=2
    parData->SetNameTitle("Monopix_seg2","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26_300kRad3e12neq",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
    parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,18.566901,"pixNoiseMPV");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,2.665233,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    return itsSimuParam;
}

//________________________________________________________________________________________________________//
//   ||                                                                                                   //
//   || Paremeterization by charge spread map of the MIMOSA32 P26 chip, 30 C,                             //
//   ||    1 MRad + 1e13 new/cm2                                                                         //
//  \||/                                                                                                  //
//___\/___________________________________________________________________________________________________//

AliITSUSimuParam* MakeITSUSimuParam_M32P26Map1MRad1e13(Double_t qscale, Double_t sncut)
{
    //const char* macroname = "MakeITSUSimuParam.C";
    //
    AliITSUSimuParam* itsSimuParam = new AliITSUSimuParam();
    //
    itsSimuParam->SetNLayers(kNLayers);
    for (int ilr=kNLayers;ilr--;) itsSimuParam->SetLrROCycleShift(kROShifts[ilr],ilr);
    //
    // Add spread function parameterization data
    AliITSUParamList* parData = 0;
    //
    //------------------------ parameterization data for segmentation 0 ----------------------
    parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(0);                                              // this is a function for detId=0
    parData->SetNameTitle("Monopix_seg1","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26_1MRad1e13neq",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
    parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,18.982178,"pixNoiseMPV");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,2.067505,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //------------------------ parameterization data for segmentation 1 ----------------------
    parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(1);                                              // this is a function for detId=1
    parData->SetNameTitle("Monopix_seg1","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26_1MRad1e13neq",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
    parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,18.982178,"pixNoiseMPV");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,2.067505,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //------------------------ parameterization data for segmentation 2 ----------------------
    parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(2);                                              // this is a function for detId=2
    parData->SetNameTitle("Monopix_seg2","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26_1MRad1e13neq",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
    parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,18.982178,"pixNoiseMPV");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,2.067505,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    return itsSimuParam;
}


//________________________________________________________________________________________________________//
//   ||                                                                                                   //
//   || Paremeterization by charge spread map of the MIMOSA32 P26 chip, 20 C,                             //
//   ||    1 MRad + 1e13 new/cm2                                                                         //
//  \||/                                                                                                  //
//___\/___________________________________________________________________________________________________//

AliITSUSimuParam* MakeITSUSimuParam_M32P26Map1MRad1e13_20deg(Double_t qscale, Double_t sncut)
{
    //const char* macroname = "MakeITSUSimuParam.C";
    //
    AliITSUSimuParam* itsSimuParam = new AliITSUSimuParam();
    //
    itsSimuParam->SetNLayers(kNLayers);
    for (int ilr=kNLayers;ilr--;) itsSimuParam->SetLrROCycleShift(kROShifts[ilr],ilr);
    //
    // Add spread function parameterization data
    AliITSUParamList* parData = 0;
    //
    //------------------------ parameterization data for segmentation 0 ----------------------
    parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(0);                                              // this is a function for detId=0
    parData->SetNameTitle("Monopix_seg1","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26_1MRad1e13neq_20deg",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
    parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,17.980760,"pixNoiseMPV");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,2.091568,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //------------------------ parameterization data for segmentation 1 ----------------------
    parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(1);                                              // this is a function for detId=1
    parData->SetNameTitle("Monopix_seg1","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26_1MRad1e13neq_20deg",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
    parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,17.980760,"pixNoiseMPV");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,2.091568,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //------------------------ parameterization data for segmentation 2 ----------------------
    parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(2);                                              // this is a function for detId=2
    parData->SetNameTitle("Monopix_seg2","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26_1MRad1e13neq_20deg",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
    parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,17.980760,"pixNoiseMPV");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,2.091568,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    return itsSimuParam;
}

//________________________________________________________________________________________________________//
//   ||                                                                                                   //
//   || Paremeterization by charge spread map of the MIMOSA32 P26 chip, 20 C,                             //
//   ||    300kRad + 3e12 new/cm2                                                                         //
//  \||/                                                                                                  //
//___\/___________________________________________________________________________________________________//

AliITSUSimuParam* MakeITSUSimuParam_M32P26Map300kRad3e12_20deg(Double_t qscale, Double_t sncut)
{
    //const char* macroname = "MakeITSUSimuParam.C";
    //
    AliITSUSimuParam* itsSimuParam = new AliITSUSimuParam();
    //
    itsSimuParam->SetNLayers(kNLayers);
    for (int ilr=kNLayers;ilr--;) itsSimuParam->SetLrROCycleShift(kROShifts[ilr],ilr);
    //
    // Add spread function parameterization data
    AliITSUParamList* parData = 0;
    //
    //------------------------ parameterization data for segmentation 0 ----------------------
    parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(0);                                              // this is a function for detId=0
    parData->SetNameTitle("Monopix_seg1","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26_300kRad3e12neq_20deg",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
    parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,18.218402,"pixNoiseMPV");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,2.616870,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //------------------------ parameterization data for segmentation 1 ----------------------
    parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(1);                                              // this is a function for detId=1
    parData->SetNameTitle("Monopix_seg1","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26_300kRad3e12neq_20deg",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
    parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,18.218402,"pixNoiseMPV");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,2.616870,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //------------------------ parameterization data for segmentation 2 ----------------------
    parData = new AliITSUParamList(AliITSUSimulationPix::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(2);                                              // this is a function for detId=2
    parData->SetNameTitle("Monopix_seg2","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26_300kRad3e12neq_20deg",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSUSimulationPix::kSpreadFunGlobalQScale,qscale,"globQscale");
    parData->SetParameter(AliITSUSimulationPix::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSUSimulationPix::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSUSimulationPix::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseMPV,18.218402,"pixNoiseMPV");
    parData->SetParameter(AliITSUSimulationPix::kPixNoiseSigma,2.616870,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSUSimulationPix::kReadOutSchemeType,AliITSUSimulationPix::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSUSimulationPix::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    return itsSimuParam;
}










//_______________________________________________________________
void SetPSFParams(TString pixType, AliITSUParamList* parData )
{
  //
  // Reads the PSF map and sets the parameters
  //
  printf("Loading charge spread histo %s from file %s\n",pixType.Data(),inpPSFName);
  TFile* fin = new TFile(inpPSFName);
  if (!fin) { 
    printf("NO parameters are set! Input file %s doesn't exist\n",inpPSFName); 
    exit(1);
  }
  //
  TH2* hProfWrk = 0;
  hProfWrk =  dynamic_cast<TH2*> fin->Get(pixType.Data());
  if(!hProfWrk) {
    printf("PSF map %s doesn't exist!!!\n",pixType.Data()); 
    exit(1);
  }
  hProfWrk = (TH2*) hProfWrk->Clone();
  hProfWrk->SetDirectory(0);
  fin->Close();  
  //
  parData->AddParamObject(hProfWrk);
  parData->SetParameter(AliITSUSimulationPix::kChargeSpreadType,AliITSUSimulationPix::kSpreadFunHisto,pixType.Data());
}
