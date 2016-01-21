//RS: before running MakeITSUSimuParam call ".x LoadLibs.C"

// Reference sensor parameterizations
AliITSMFTSimuParam*  MakeITSUSimuParam_M32P26Map(Double_t sncut);// don't forget to turn off diode shift for P26 in Config.C
AliITSMFTSimuParam*  MakeITSUSimuParam_M32terP31Map( Double_t sncut);
// Irradiated sensor parameterizations P31
AliITSMFTSimuParam*  MakeITSUSimuParam_M32terP31Map1MRad( Double_t sncut);
AliITSMFTSimuParam*  MakeITSUSimuParam_M32terP31Map300kRad3e12( Double_t sncut);
AliITSMFTSimuParam*  MakeITSUSimuParam_M32terP31Map1MRad1e13( Double_t sncut);
// Irradiated sensor parameterizations P26
AliITSMFTSimuParam*  MakeITSUSimuParam_M32P26Map300kRad3e12( Double_t sncut);
AliITSMFTSimuParam*  MakeITSUSimuParam_M32P26Map1MRad1e13( Double_t sncut);
AliITSMFTSimuParam*  MakeITSUSimuParam_M32P26Map300kRad3e12_20deg( Double_t sncut);
AliITSMFTSimuParam*  MakeITSUSimuParam_M32P26Map1MRad1e13_20deg( Double_t sncut);

// Reference sensor parameterizations with 2 2D Gaussian
AliITSMFTSimuParam*  MakeITSUSimuParam_M32P26MapRecenteredBroadened(Double_t broadening);//


// Digital sensor parameter response test
AliITSMFTSimuParam*  MakeITSUSimuParam_1stDigital();//



void SetPSFParams(TString pixType, AliITSMFTParamList* parData);

// To turn noise generation ON set these values to 1
const int kAddNoise = -1;
const int kAddNoiseInAllMod = -1;

const char* inpPSFName = "$ALICE_ROOT/ITSMFT/COMMONITSMFT/misc/ITSU_pixel_response_PSFs.root";
const int kNLayers = 7;


const int kSNcut = 5;   // Threshold/Noise cut. CAN BE CHANGED IN THE RANGE [5,10] no other values are allowed.
const int knSNcut = 6 ; // Number of tuned SNcut fixed to 6.

const double kBroadeningFactor = 1.0; // For the 2 2D  Gaussian parameterization, allowed range [0.5,2.0]


/*
 these are readout phases settings:
 the rule is:
 1) abs(kROShifts)<1: the modules of the layer are synchronized, the layer
 phase is set to kROShifts*ROCycleLength of the modules of this layer
 2) abs(kROShifts)>1: each module within the layer will have random phase within its ROCycleLength
 */
const float kROShifts[kNLayers] = {0.5,0.5,0.5, -0.5,-0.5, 0.5,0.5};

void MakeITSUSimuParam(Int_t simuType = 9, const char* cdbURI="local://") {
    //========================================================================
    //
    // Steering macro for ITS simulation parameters
    //
    // Author: L.Molnar
    // Contact: levente.molnar@cern.ch
    //
    //
    //****** DO NOT FORGET TO SET THE RIGHT GEOMETRY IN THE Config.C *****
    //
    //****** P26 chips Config.C parameters *****
    // col pitch 20 um
    // row pitch 20 um
    // sensor thickness 18 um
    // SET diode staggering to:   kDiodShiftM32terP31X[][] = {0.0, 0.0};
    // SET diode staggering to:   kDiodShiftM32terP31Z[][] = {0.0, 0.0};
    //
    //
    //****** P31 chips Config.C parameters  *****
    // col pitch 20 um
    // row pitch 33 um
    // sensor thickness 18 um
    // SET diode staggering to:   kDiodShiftM32terP31X[][] = {0.30,-0.19};
    // SET diode staggering to:   kDiodShiftM32terP31Z[][] = {0.0, 0.0};
    //
    //
    //========================================================================
    AliITSMFTSimuParam *param = 0;
    //
    // Select only one parameterziation... and don't forget to set 18 um thickness in Config.C !!!
    
    switch(simuType) {
            
        case 0:
            //____ MIMOSA32 P26 Response parameterzied by map ---> selected for TDR
            param = MakeITSUSimuParam_M32P26Map(kSNcut);
            break;
            //
        case 1:
            //____ MIMOSA32ter P31 Response parameterzied by map
            param = MakeITSUSimuParam_M32terP31Map(kSNcut);
            break;
            //
        case 2:
            //____ MIMOSA32ter P31 Response parameterzied by map 1MRad irradiation
            param = MakeITSUSimuParam_M32terP31Map1MRad(kSNcut);
            break;
            //
        case 3:
            //____ MIMOSA32ter P31 Response parameterzied by map , 300kRad + 3e12 neq/cm2 irradiation
            param = MakeITSUSimuParam_M32terP31Map300kRad3e12(kSNcut);
            break;
            //
        case 4:
            //____ MIMOSA32ter P31 Response parameterzied by map , 1MRad+ 1e13 neq/cm2 irradiation
            param = MakeITSUSimuParam_M32terP31Map1MRad1e13(kSNcut);
            break;
            //
        case 5:
            //___ MIMOSA32 P26 , 300kRad + 3e12 neq/cm2 irradiation 30 deg
            param = MakeITSUSimuParam_M32P26Map300kRad3e12(kSNcut);
            break;
            //
        case 6:
            //____ MIMOSA32 P26 , 300kRad + 3e12 neq/cm2 irradiation 20 deg
            param = MakeITSUSimuParam_M32P26Map300kRad3e12_20deg(kSNcut);
            break;
            //
        case 7:
            //___ MIMOSA32 P26 ,  1MRad+ 1e13 neq/cm2 irradiation 30 deg
            param = MakeITSUSimuParam_M32P26Map1MRad1e13(kSNcut);
            break;
            //
        case 8:
            //___ MIMOSA32 P26 ,  1MRad+ 1e13 neq/cm2 irradiation 20 deg
            param = MakeITSUSimuParam_M32P26Map1MRad1e13_20deg(kSNcut);
            break;
            //
        case 9:
            //____ MIMOSA32 P26 Response parameterzied 2 2D Gaussian and recentered to 0,0 and the sigmas are broadened for TDR6 vs TDR7 geometry study
            //____ Map only available for SNcut = 5 !!!
            //____ Introduced on the 31/03/2014
            param = MakeITSUSimuParam_M32P26MapRecenteredBroadened(kBroadeningFactor);
            break;
            //
        case 10:
            //____ 1st digital chip response test
            param = MakeITSUSimuParam_1stDigital();
            break;
            //
        default:
            break;
    }
    
    
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
//   || Paremeterization by charge spread map of the MIMOSA32 P26 chip, 30 C, no irradiation              //
//   ||                                                                                                   //
//  \||/                                                                                                  //
//___\/___________________________________________________________________________________________________//
AliITSMFTSimuParam* MakeITSUSimuParam_M32P26Map( Int_t sncut)
{
    //const char* macroname = "MakeITSUSimuParam.C";
    //
    AliITSMFTSimuParam* itsSimuParam = new AliITSMFTSimuParam();
    //
    itsSimuParam->SetNLayers(kNLayers);
    for (int ilr=kNLayers;ilr--;) itsSimuParam->SetLrROCycleShift(kROShifts[ilr],ilr);
    //
    Double_t qscale[knSNcut]={1.036868,	1.055369,	1.083679,	1.098877,	1.126203,	1.145552};
    if(sncut < 5 || sncut > 10 ) {printf("---> ERROR ERROR ERROR requested SNcut: %d is not valid... Check the macro header! \n",sncut); return;}
    //
    // Add spread function parameterization data
    AliITSMFTParamList* parData = 0;
    //
    //------------------------ parameterization data for segmentation 0 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(0);                                              // this is a function for detId=0
    parData->SetNameTitle("Monopix_seg1","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,17.53,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,2.93,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //------------------------ parameterization data for segmentation 1 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(1);                                              // this is a function for detId=1
    parData->SetNameTitle("Monopix_seg1","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,17.53,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,2.93,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //------------------------ parameterization data for segmentation 2 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(2);                                              // this is a function for detId=2
    parData->SetNameTitle("Monopix_seg2","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,17.53,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,2.93,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
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

AliITSMFTSimuParam* MakeITSUSimuParam_M32terP31Map( Double_t sncut)
{
    //const char* macroname = "MakeITSUSimuParam.C";
    //
    AliITSMFTSimuParam* itsSimuParam = new AliITSMFTSimuParam();
    //
    itsSimuParam->SetNLayers(kNLayers);
    for (int ilr=kNLayers;ilr--;) itsSimuParam->SetLrROCycleShift(kROShifts[ilr],ilr);
    //
    //
    Double_t qscale[knSNcut]={1.396168,	1.439231,	1.484984,	1.534129,	1.570807,	1.600674};
    if(sncut < 5 || sncut > 10 ) {printf("---> ERROR ERROR ERROR requested SNcut: %d is not valid... Check the macro header! \n",sncut); return;}
    //
    // Add spread function parameterization data
    AliITSMFTParamList* parData = 0;
    //
    //------------------------ parameterization data for segmentation 0 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(0);                                              // this is a function for detId=0
    parData->SetNameTitle("Monopix_seg2","PSF map for M32terP31");
    SetPSFParams("hProfPSD_M32terP31",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,20.62,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,3.55,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //
    //------------------------ parameterization data for segmentation 1 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(1);                                              // this is a function for detId=1
    parData->SetNameTitle("Monopix_seg1","PSF map for M32terP31");
    SetPSFParams("hProfPSD_M32terP31",parData);
    
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,20.62,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,3.55,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    // now set the parameters according selected function
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //
    //------------------------ parameterization data for segmentation 2 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(2);                                              // this is a function for detId=2
    parData->SetNameTitle("Monopix_seg2","PSF map for M32terP31");
    SetPSFParams("hProfPSD_M32terP31", parData );
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale"); //980./1080.
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,20.62,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,3.55,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)"); // cycle lenght is 10-30 us --> 10-30e-6 s
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
AliITSMFTSimuParam* MakeITSUSimuParam_M32terP31Map1MRad( Double_t sncut)
{
    //const char* macroname = "MakeITSUSimuParam.C";
    //
    AliITSMFTSimuParam* itsSimuParam = new AliITSMFTSimuParam();
    //
    itsSimuParam->SetNLayers(kNLayers);
    for (int ilr=kNLayers;ilr--;) itsSimuParam->SetLrROCycleShift(kROShifts[ilr],ilr);
    //
    //
    printf(" ---> WARNING WARNING WARNING --- Parameterization is not final, it is set to the reference sensor!"); return;
    Double_t qscale[knSNcut]={1.396168,	1.439231,	1.484984,	1.534129,	1.570807,	1.600674};
    if(sncut < 5 || sncut > 10 ) {printf("---> ERROR ERROR ERROR requested SNcut: %d is not valid... Check the macro header! \n",sncut); return;}
    //
    // Add spread function parameterization data
    AliITSMFTParamList* parData = 0;
    //
    //------------------------ parameterization data for segmentation 0 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(0);                                              // this is a function for detId=0
    parData->SetNameTitle("Monopix_seg2","PSF map for M32terP31");
    SetPSFParams("hProfPSD_M32terP31_1MRad",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,21.843706,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,3.417494,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //
    //------------------------ parameterization data for segmentation 1 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(1);                                              // this is a function for detId=1
    parData->SetNameTitle("Monopix_seg1","PSF map for M32terP31");
    SetPSFParams("hProfPSD_M32terP31_1MRad",parData);
    
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,21.843706,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,3.417494,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    // now set the parameters according selected function
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //
    //------------------------ parameterization data for segmentation 2 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(2);                                              // this is a function for detId=2
    parData->SetNameTitle("Monopix_seg2","PSF map for M32terP31");
    SetPSFParams("hProfPSD_M32terP31_1MRad", parData );
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale"); //980./1080.
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,21.843706,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,3.417494,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)"); // cycle lenght is 10-30 us --> 10-30e-6 s
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


AliITSMFTSimuParam* MakeITSUSimuParam_M32terP31Map300kRad3e12( Double_t sncut)
{
    //const char* macroname = "MakeITSUSimuParam.C";
    //
    AliITSMFTSimuParam* itsSimuParam = new AliITSMFTSimuParam();
    //
    itsSimuParam->SetNLayers(kNLayers);
    for (int ilr=kNLayers;ilr--;) itsSimuParam->SetLrROCycleShift(kROShifts[ilr],ilr);
    //
    //
    Double_t qscale[knSNcut]={1.532517,	1.617223,	1.641962,	1.714945,	1.73809,	1.749223};
    if(sncut < 5 || sncut > 10 ) {printf("---> ERROR ERROR ERROR requested SNcut: %d is not valid... Check the macro header! \n",sncut); return;}
    //
    // Add spread function parameterization data
    AliITSMFTParamList* parData = 0;
    //
    //------------------------ parameterization data for segmentation 0 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(0);                                              // this is a function for detId=0
    parData->SetNameTitle("Monopix_seg2","PSF map for M32terP31");
    SetPSFParams("hProfPSD_M32terP31_300kRad3e12neq",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,27.274548,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,2.723949,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //
    //------------------------ parameterization data for segmentation 1 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(1);                                              // this is a function for detId=1
    parData->SetNameTitle("Monopix_seg1","PSF map for M32terP31");
    SetPSFParams("hProfPSD_M32terP31_300kRad3e12neq",parData);
    
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,27.274548,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,2.723949,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    // now set the parameters according selected function
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //
    //------------------------ parameterization data for segmentation 2 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(2);                                              // this is a function for detId=2
    parData->SetNameTitle("Monopix_seg2","PSF map for M32terP31");
    SetPSFParams("hProfPSD_M32terP31_300kRad3e12neq", parData );
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale"); //980./1080.
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,27.274548,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,2.723949,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)"); // cycle lenght is 10-30 us --> 10-30e-6 s
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


AliITSMFTSimuParam* MakeITSUSimuParam_M32terP31Map1MRad1e13( Double_t sncut)
{
    //const char* macroname = "MakeITSUSimuParam.C";
    //
    AliITSMFTSimuParam* itsSimuParam = new AliITSMFTSimuParam();
    //
    itsSimuParam->SetNLayers(kNLayers);
    for (int ilr=kNLayers;ilr--;) itsSimuParam->SetLrROCycleShift(kROShifts[ilr],ilr);
    //
    //
    Double_t qscale[knSNcut]={1.530764,	1.542968,	1.602497,	1.621188,	1.68353	,1.639263};
    if(sncut < 5 || sncut > 10 ) {printf("---> ERROR ERROR ERROR requested SNcut: %d is not valid... Check the macro header! \n",sncut); return;}
    //
    // Add spread function parameterization data
    AliITSMFTParamList* parData = 0;
    //
    //------------------------ parameterization data for segmentation 0 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(0);                                              // this is a function for detId=0
    parData->SetNameTitle("Monopix_seg2","PSF map for M32terP31");
    SetPSFParams("hProfPSD_M32terP31_1MRad1e13neq",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,43.840798,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,3.187048,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //
    //------------------------ parameterization data for segmentation 1 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(1);                                              // this is a function for detId=1
    parData->SetNameTitle("Monopix_seg1","PSF map for M32terP31");
    SetPSFParams("hProfPSD_M32terP31_1MRad1e13neq",parData);
    
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,43.840798,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,3.187048,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    // now set the parameters according selected function
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //
    //------------------------ parameterization data for segmentation 2 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(2);                                              // this is a function for detId=2
    parData->SetNameTitle("Monopix_seg2","PSF map for M32terP31");
    SetPSFParams("hProfPSD_M32terP31_1MRad1e13neq", parData );
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale"); //980./1080.
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,43.840798,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,3.187048,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)"); // cycle lenght is 10-30 us --> 10-30e-6 s
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

AliITSMFTSimuParam* MakeITSUSimuParam_M32P26Map300kRad3e12( Double_t sncut)
{
    //const char* macroname = "MakeITSUSimuParam.C";
    //
    AliITSMFTSimuParam* itsSimuParam = new AliITSMFTSimuParam();
    //
    itsSimuParam->SetNLayers(kNLayers);
    for (int ilr=kNLayers;ilr--;) itsSimuParam->SetLrROCycleShift(kROShifts[ilr],ilr);
    //
    //
    Double_t qscale[knSNcut]={0.446298,	0.45217,	0.405819,	0.457178,	0.488428,	0.508903};
    if(sncut < 5 || sncut > 10 ) {printf("---> ERROR ERROR ERROR requested SNcut: %d is not valid... Check the macro header! \n",sncut); return;}
    //
    // Add spread function parameterization data
    AliITSMFTParamList* parData = 0;
    //
    //------------------------ parameterization data for segmentation 0 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(0);                                              // this is a function for detId=0
    parData->SetNameTitle("Monopix_seg1","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26_300kRad3e12neq",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,18.566901,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,2.665233,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //------------------------ parameterization data for segmentation 1 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(1);                                              // this is a function for detId=1
    parData->SetNameTitle("Monopix_seg1","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26_300kRad3e12neq",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,18.566901,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,2.665233,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //------------------------ parameterization data for segmentation 2 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(2);                                              // this is a function for detId=2
    parData->SetNameTitle("Monopix_seg2","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26_300kRad3e12neq",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,18.566901,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,2.665233,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
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

AliITSMFTSimuParam* MakeITSUSimuParam_M32P26Map1MRad1e13( Double_t sncut)
{
    //const char* macroname = "MakeITSUSimuParam.C";
    //
    AliITSMFTSimuParam* itsSimuParam = new AliITSMFTSimuParam();
    //
    itsSimuParam->SetNLayers(kNLayers);
    for (int ilr=kNLayers;ilr--;) itsSimuParam->SetLrROCycleShift(kROShifts[ilr],ilr);
    //
    //
    Double_t qscale[knSNcut]={0.406783,	0.418178,	0.410906,  0.406477,	0.355895,	0.479006};
    if(sncut < 5 || sncut > 10 ) {printf("---> ERROR ERROR ERROR requested SNcut: %d is not valid... Check the macro header! \n",sncut); return;}
    //
    // Add spread function parameterization data
    AliITSMFTParamList* parData = 0;
    //
    //------------------------ parameterization data for segmentation 0 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(0);                                              // this is a function for detId=0
    parData->SetNameTitle("Monopix_seg1","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26_1MRad1e13neq",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,18.982178,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,2.067505,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //------------------------ parameterization data for segmentation 1 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(1);                                              // this is a function for detId=1
    parData->SetNameTitle("Monopix_seg1","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26_1MRad1e13neq",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,18.982178,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,2.067505,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //------------------------ parameterization data for segmentation 2 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(2);                                              // this is a function for detId=2
    parData->SetNameTitle("Monopix_seg2","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26_1MRad1e13neq",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,18.982178,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,2.067505,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
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

AliITSMFTSimuParam* MakeITSUSimuParam_M32P26Map1MRad1e13_20deg( Double_t sncut)
{
    //const char* macroname = "MakeITSUSimuParam.C";
    //
    AliITSMFTSimuParam* itsSimuParam = new AliITSMFTSimuParam();
    //
    itsSimuParam->SetNLayers(kNLayers);
    for (int ilr=kNLayers;ilr--;) itsSimuParam->SetLrROCycleShift(kROShifts[ilr],ilr);
    //
    //
    printf(" ---> WARNING WARNING WARNING --- Parameterization is not final, it is set to the 30 deg irradiated sensor!"); return;
    Double_t qscale[knSNcut]={0.446298,	0.45217,	0.405819,	0.457178,	0.488428,	0.508903};
    if(sncut < 5 || sncut > 10 ) {printf("---> ERROR ERROR ERROR requested SNcut: %d is not valid... Check the macro header! \n",sncut); return;}
    //
    // Add spread function parameterization data
    AliITSMFTParamList* parData = 0;
    //
    //------------------------ parameterization data for segmentation 0 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(0);                                              // this is a function for detId=0
    parData->SetNameTitle("Monopix_seg1","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26_1MRad1e13neq_20deg",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,17.980760,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,2.091568,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //------------------------ parameterization data for segmentation 1 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(1);                                              // this is a function for detId=1
    parData->SetNameTitle("Monopix_seg1","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26_1MRad1e13neq_20deg",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,17.980760,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,2.091568,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //------------------------ parameterization data for segmentation 2 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(2);                                              // this is a function for detId=2
    parData->SetNameTitle("Monopix_seg2","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26_1MRad1e13neq_20deg",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,17.980760,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,2.091568,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
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

AliITSMFTSimuParam* MakeITSUSimuParam_M32P26Map300kRad3e12_20deg( Double_t sncut)
{
    //const char* macroname = "MakeITSUSimuParam.C";
    //
    AliITSMFTSimuParam* itsSimuParam = new AliITSMFTSimuParam();
    //
    itsSimuParam->SetNLayers(kNLayers);
    for (int ilr=kNLayers;ilr--;) itsSimuParam->SetLrROCycleShift(kROShifts[ilr],ilr);
    //
    //
    printf(" ---> WARNING WARNING WARNING --- Parameterization is not final, it is set to the 30 deg irradiated sensor!"); return;
    Double_t qscale[knSNcut]={0.446298,	0.45217	,0.405819,	0.457178,	0.488428,	0.508903};
    if(sncut < 5 || sncut > 10 ) {printf("---> ERROR ERROR ERROR requested SNcut: %d is not valid... Check the macro header! \n",sncut); return;}
    //
    // Add spread function parameterization data
    AliITSMFTParamList* parData = 0;
    //
    //------------------------ parameterization data for segmentation 0 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(0);                                              // this is a function for detId=0
    parData->SetNameTitle("Monopix_seg1","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26_300kRad3e12neq_20deg",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,18.218402,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,2.616870,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //------------------------ parameterization data for segmentation 1 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(1);                                              // this is a function for detId=1
    parData->SetNameTitle("Monopix_seg1","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26_300kRad3e12neq_20deg",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,18.218402,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,2.616870,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //------------------------ parameterization data for segmentation 2 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(2);                                              // this is a function for detId=2
    parData->SetNameTitle("Monopix_seg2","PSF map for M32P26");
    SetPSFParams("hProfPSD_M32P26_300kRad3e12neq_20deg",parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-4,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,18.218402,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,2.616870,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    return itsSimuParam;
}


//
//________________________________________________________________________________________________________//
//   ||                                                                                                   //
//   || Paremeterization by charge spread map of the MIMOSA32 P26 chip, 30 C, no irradiation              //
//   ||   - Map is fitted by two 2D Gaussian and recetnered to the x,z = [0,0]                            //
//  \||/  - The Gaussian sigma can be broadened with the factor "broadening" to mimic larger pixels       //
//   \/   - then 20 x0 20 um^2                                                                            //
//        - The fake rate is set to 10^-5                                                                 //
//________________________________________________________________________________________________________//

AliITSMFTSimuParam* MakeITSUSimuParam_M32P26MapRecenteredBroadened(Double_t broadening)
{
    //const char* macroname = "MakeITSUSimuParam.C";
    //
    const Int_t sncut = 5;
    TString histoName = Form("hProfPSD_M32P26_Cent_Broad_%d",TMath::Nint(broadening*100));
    
    
    //
    AliITSMFTSimuParam* itsSimuParam = new AliITSMFTSimuParam();
    //
    itsSimuParam->SetNLayers(kNLayers);
    for (int ilr=kNLayers;ilr--;) itsSimuParam->SetLrROCycleShift(kROShifts[ilr],ilr);
    //
    
    Double_t qscale[1]={1.036868};
    if(sncut != 5 ) {printf("---> ERROR ERROR ERROR requested SNcut: %d is not valid... Check the macro! \n",sncut); return;}
    //
    // Add spread function parameterization data
    AliITSMFTParamList* parData = 0;
    //
    //------------------------ parameterization data for segmentation 0 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(0);                                              // this is a function for detId=0
    parData->SetNameTitle("Monopix_seg1","PSF map for M32P26");
    SetPSFParams(histoName,parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-5,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,17.53,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,2.93,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //------------------------ parameterization data for segmentation 1 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(1);                                              // this is a function for detId=1
    parData->SetNameTitle("Monopix_seg1","PSF map for M32P26");
    SetPSFParams(histoName,parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-5,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,17.53,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,2.93,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //------------------------ parameterization data for segmentation 2 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(2);                                              // this is a function for detId=2
    parData->SetNameTitle("Monopix_seg2","PSF map for M32P26");
    SetPSFParams(histoName,parData);
    //
    // obligatory params for all AliITSUSimulationPix functions: number of pixels in X,Z around
    // injected one to consider
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,0,"DigitalSim");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNXoffs,2,"nPixX");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunParamNZoffs,2,"nPixZ");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunMinSteps,10,"nChargeSteps");
    parData->SetParameter(AliITSMFTSimuParam::kSpreadFunGlobalQScale,qscale[sncut-5],"globQscale");
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-5,"pixFakeRate");
    parData->SetParameter(AliITSMFTSimuParam::kPixSNDisrcCut,sncut,"pixSNDisrcCut");
    parData->SetParameter(AliITSMFTSimuParam::kPixMinElToAdd,1,"pixMinElToAdd");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,17.53,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,2.93,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    return itsSimuParam;
}

//
//________________________________________________________________________________________________________//
//   ||                                                                                                   //
//   ||     First digital chip response parameterization                                                  //
//   ||   - The fake rate is set to 10^-5                                                                 //
//  \||/                                                                                                  //
//   \/                                                                                                   //
//                                                                                                        //
//________________________________________________________________________________________________________//

AliITSMFTSimuParam* MakeITSUSimuParam_1stDigital()
{
    //
    AliITSMFTSimuParam* itsSimuParam = new AliITSMFTSimuParam();
    //
    itsSimuParam->SetNLayers(kNLayers);
    for (int ilr=kNLayers;ilr--;) itsSimuParam->SetLrROCycleShift(kROShifts[ilr],ilr);
    //
    // Add spread function parameterization data
    AliITSMFTParamList* parData = 0;
    //
    //------------------------ parameterization data for segmentation 0 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(0);                                              // this is a function for detId=0
    parData->SetNameTitle("Monopix_seg1","First digital param");
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,1,"DigitalSim");
    // obligatory params for all AliITSUSimulationPix functions: fake rate
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-5,"pixFakeRate");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    //
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,17.53,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,2.93,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //------------------------ parameterization data for segmentation 1 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(1);                                              // this is a function for detId=0
    parData->SetNameTitle("Monopix_seg1","First digital param");
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,1,"DigitalSim");
    //
    // obligatory params for all AliITSUSimulationPix functions: fake rate
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-5,"pixFakeRate");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    //
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,17.53,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,2.93,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);
    //
    //
    //------------------------ parameterization data for segmentation 2 ----------------------
    parData = new AliITSMFTParamList(AliITSMFTSimuParam::kNReservedParams);   // no custom params are needed
    parData->SetUniqueID(2);                                              // this is a function for detId=0
    parData->SetNameTitle("Monopix_seg1","First digital param");
    //
    parData->SetParameter(AliITSMFTSimuParam::kDigitalSim,1,"DigitalSim");
    //
    // obligatory params for all AliITSUSimulationPix functions: fake rate
    parData->SetParameter(AliITSMFTSimuParam::kPixFakeRate,1e-5,"pixFakeRate");
    // Noise
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseIsOn,kAddNoise,"pixNoiseIsOn");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseInAllMod,kAddNoiseInAllMod,"pixNoiseIsOnInAllMod");
    //
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseMPV,17.53,"pixNoiseMPV");
    parData->SetParameter(AliITSMFTSimuParam::kPixNoiseSigma,2.93,"pixNoiseSigma");
    // and readout timing scheme
    parData->SetParameter(AliITSMFTSimuParam::kReadOutSchemeType,AliITSMFTSimuParam::kReadOutRollingShutter,"ROType");
    parData->SetParameter(AliITSMFTSimuParam::kReadOutCycleLength,25e-6,"ROCycle(s)");
    //
    itsSimuParam->AddRespFunParam(parData);

    //
    return itsSimuParam;

}





//_______________________________________________________________
void SetPSFParams(TString pixType, AliITSMFTParamList* parData )
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
    parData->SetParameter(AliITSMFTSimuParam::kChargeSpreadType,AliITSMFTSimulationPix::kSpreadFunHisto,pixType.Data());
}
