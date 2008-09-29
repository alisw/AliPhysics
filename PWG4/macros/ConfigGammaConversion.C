
/** ------------------------------ Monte Carlo flag -----------------------------------------*/
Bool_t doMCTruth = kTRUE;
/** ---------------------------- end Monte Carlo flag ---------------------------------------*/

/** ------------------------- Choose KFParticle OR ESDTrack  --------------------------------*/
Bool_t useKFParticle = kTRUE;
Bool_t useESDTrack   = kFALSE;
/** ----------------------- end Choose KFParticle OR ESDTrack  -----------------------------*/


Bool_t calculateBackground = kTRUE;

/** ---------------------------------- define cuts here ------------------------------------*/

Int_t pidOfNegativeTrack=11;
Int_t pidOfPositiveTrack=-11;

Double_t maxRCut   = 200.;
Double_t etaCut    = 1.2;
Double_t ptCut     = 0.1;
Double_t chi2Cut   = 20.;

Double_t xVertexCut = 0.;
Double_t yVertexCut = 0.;
Double_t zVertexCut = 0.;

Double_t sigmaCutGammaMass=0.0001;

Bool_t useImprovedVertex = kTRUE;

// define masses of different particles, this will be used by the KF particle
// together with the width to set mass constraints. Units in GeV.
Double_t electronMass = 0.00051099892;
Double_t gammaMass    = 0.;
Double_t pi0Mass      = 0.1349766;
Double_t etaMass      = 0.54751;

// define the width constraint used by KF particle.
Double_t gammaWidth = 0.01;
Double_t pi0Width   = 0.01;
Double_t etaWidth   = 0.01;

// define the probability of track being an electron
Double_t probElectron = 0.5;

/** ----------------------------------end define cuts here----------------------------------*/

/** -------------------------------- Phi/R Mapping ---------------------------------------*/
Int_t nPhiIndex = 18;
Int_t nRIndex   = 40;

Double_t minRadius   = 0.;
Double_t maxRadius   = 200.;
Double_t minPhi      = -TMath::Pi();
Double_t maxPhi      = TMath::Pi();
/** ------------------------------- end Phi/R Mapping ------------------------------------*/

/** ------------------- define which histograms to plot here --------------------------------*/
/**   NB: to change the bin numbers, see below the histogram flags                           */
Bool_t plotMC_EP_R                             = kTRUE;
Bool_t plotMC_EP_Z_R                           = kTRUE;
Bool_t plotMC_EP_X_Y                           = kTRUE;
Bool_t plotMC_EP_OpeningAngle                  = kTRUE;

Bool_t plotMC_E_Energy                         = kTRUE;
Bool_t plotMC_E_Pt                             = kTRUE;
Bool_t plotMC_E_Eta                            = kTRUE;
Bool_t plotMC_E_Phi                            = kTRUE;

Bool_t plotMC_P_Energy                         = kTRUE;
Bool_t plotMC_P_Pt                             = kTRUE;
Bool_t plotMC_P_Eta                            = kTRUE;
Bool_t plotMC_P_Phi                            = kTRUE;

Bool_t plotMC_Gamma_Energy                     = kTRUE;
Bool_t plotMC_Gamma_Pt                         = kTRUE;
Bool_t plotMC_Gamma_Eta                        = kTRUE;
Bool_t plotMC_Gamma_Phi                        = kTRUE;

Bool_t plotMC_DirectGamma_Energy               = kTRUE;
Bool_t plotMC_DirectGamma_Pt                   = kTRUE;
Bool_t plotMC_DirectGamma_Eta                  = kTRUE;
Bool_t plotMC_DirectGamma_Phi                  = kTRUE;

Bool_t plotMC_Match_Gamma_Eta                  = kTRUE;
Bool_t plotMC_Match_Gamma_Phi                  = kTRUE;
Bool_t plotMC_Match_Gamma_Pt                   = kTRUE;
Bool_t plotMC_Match_Gamma_Energy               = kTRUE;
Bool_t plotMC_Match_Gamma_Mass                 = kTRUE;
Bool_t plotMC_Match_Gamma_OpeningAngle         = kTRUE;
Bool_t plotMC_Match_Gamma_R                    = kTRUE;
Bool_t plotMC_Match_Gamma_Z_R                  = kTRUE;
Bool_t plotMC_Match_Gamma_X_Y                  = kTRUE;

Bool_t plotMC_Pi0_Eta                          = kTRUE;
Bool_t plotMC_Pi0_Phi                          = kTRUE;
Bool_t plotMC_Pi0_Pt                           = kTRUE;
Bool_t plotMC_Pi0_Energy                       = kTRUE;
Bool_t plotMC_Pi0_Mass                         = kTRUE;
Bool_t plotMC_Pi0_OpeningAngle                 = kTRUE;
Bool_t plotMC_Pi0_R                            = kTRUE;
Bool_t plotMC_Pi0_Z_R                          = kTRUE;
Bool_t plotMC_Pi0_X_Y                          = kTRUE;

Bool_t plotMC_Eta_Eta                          = kTRUE;
Bool_t plotMC_Eta_Phi                          = kTRUE;
Bool_t plotMC_Eta_Pt                           = kTRUE;
Bool_t plotMC_Eta_Energy                       = kTRUE;
Bool_t plotMC_Eta_Mass                         = kTRUE;
Bool_t plotMC_Eta_OpeningAngleGamma            = kTRUE;
Bool_t plotMC_Eta_R                            = kTRUE;
Bool_t plotMC_Eta_Z_R                          = kTRUE;
Bool_t plotMC_Eta_X_Y                          = kTRUE;
    
// Histograms from esd tracks
Bool_t plotESD_EP_R                            = kTRUE;
Bool_t plotESD_EP_Z_R                          = kTRUE;
Bool_t plotESD_EP_X_Y                          = kTRUE;
Bool_t plotESD_EP_OpeningAngle                 = kTRUE;

Bool_t plotESD_E_Energy                        = kTRUE;
Bool_t plotESD_E_Pt                            = kTRUE;
Bool_t plotESD_E_Eta                           = kTRUE;
Bool_t plotESD_E_Phi                           = kTRUE;

Bool_t plotESD_P_Energy                        = kTRUE;
Bool_t plotESD_P_Pt                            = kTRUE;
Bool_t plotESD_P_Eta                           = kTRUE;
Bool_t plotESD_P_Phi                           = kTRUE;


Bool_t plotESD_Gamma_Energy                    = kTRUE;
Bool_t plotESD_Gamma_Pt                        = kTRUE;
Bool_t plotESD_Gamma_Eta                       = kTRUE;
Bool_t plotESD_Gamma_Phi                       = kTRUE;

Bool_t plotESD_Match_Gamma_OpeningAngle        = kTRUE;
Bool_t plotESD_Match_Gamma_Energy              = kTRUE;
Bool_t plotESD_Match_Gamma_Pt                  = kTRUE;
Bool_t plotESD_Match_Gamma_Eta                 = kTRUE;
Bool_t plotESD_Match_Gamma_Phi                 = kTRUE;
Bool_t plotESD_Match_Gamma_Mass                = kTRUE;
Bool_t plotESD_Match_Gamma_Width               = kTRUE;
Bool_t plotESD_Match_Gamma_Chi2                = kTRUE;
Bool_t plotESD_Match_Gamma_NDF                 = kTRUE;
Bool_t plotESD_Match_Gamma_R                   = kTRUE;
Bool_t plotESD_Match_Gamma_Z_R                 = kTRUE;
Bool_t plotESD_Match_Gamma_X_Y                 = kTRUE;

Bool_t plotESD_Pi0_OpeningAngleGamma           = kTRUE;
Bool_t plotESD_Pi0_Energy                      = kTRUE;
Bool_t plotESD_Pi0_Pt                          = kTRUE;
Bool_t plotESD_Pi0_Eta                         = kTRUE;
Bool_t plotESD_Pi0_Phi                         = kTRUE;
Bool_t plotESD_Pi0_Mass                        = kTRUE;
Bool_t plotESD_Pi0_R                           = kTRUE;
Bool_t plotESD_Pi0_Z_R                         = kTRUE;
Bool_t plotESD_Pi0_X_Y                         = kTRUE;

Bool_t plotESD_Eta_OpeningAngleGamma           = kTRUE;
Bool_t plotESD_Eta_Energy                      = kTRUE;
Bool_t plotESD_Eta_Pt                          = kTRUE;
Bool_t plotESD_Eta_Eta                         = kTRUE;
Bool_t plotESD_Eta_Phi                         = kTRUE;
Bool_t plotESD_Eta_Mass                        = kTRUE;
Bool_t plotESD_Eta_R                           = kTRUE;
Bool_t plotESD_Eta_Z_R                         = kTRUE;
Bool_t plotESD_Eta_X_Y                         = kTRUE;

Bool_t plotESD_Background_OpeningAngleGamma    = kTRUE;
Bool_t plotESD_Background_Energy               = kTRUE;
Bool_t plotESD_Background_Pt                   = kTRUE;
Bool_t plotESD_Background_Eta                  = kTRUE;
Bool_t plotESD_Background_Phi                  = kTRUE;
Bool_t plotESD_Background_Mass                 = kTRUE;
Bool_t plotESD_Background_R                    = kTRUE;
Bool_t plotESD_Background_Z_R                  = kTRUE;
Bool_t plotESD_Background_X_Y                  = kTRUE;

Bool_t plotMapping                             = kTRUE;       

Bool_t plotResolution_dPt                      = kTRUE;
Bool_t plotResolution_dR                       = kTRUE;
Bool_t plotResolution_dZ                       = kTRUE;
  
Bool_t plotResolution_dR_dPt                   = kTRUE;

Bool_t plotResolution_MC_Pt                    = kTRUE;
Bool_t plotResolution_MC_R                     = kTRUE;
Bool_t plotResolution_MC_Z                     = kTRUE;

Bool_t plotResolution_ESD_Pt                   = kTRUE;
Bool_t plotResolution_ESD_R                    = kTRUE;
Bool_t plotResolution_ESD_Z                    = kTRUE;

Bool_t plotNumberOfV0s                         = kTRUE;
Bool_t plotNumberOfSurvivingV0s                = kTRUE;

  //  debug histograms
Bool_t plotV0MassDebugCut1                     = kTRUE;
Bool_t plotV0MassDebugCut2                     = kTRUE;
Bool_t plotV0MassDebugCut3                     = kTRUE;
Bool_t plotV0MassDebugCut4                     = kTRUE;
Bool_t plotV0MassDebugCut5                     = kTRUE;
Bool_t plotV0MassDebugCut6                     = kTRUE;
Bool_t plotV0MassDebugCut7                     = kTRUE;
Bool_t plotV0MassDebugCut8                     = kTRUE;


/** ----------------- end define which histograms to plot here -------------------------------*/



/** ----------- Define the binning for the different plot types here -------------------------*/
//R-plots
Int_t nXBinsR = 1000;
Double_t firstXBinR = 0.;
Double_t lastXBinR = 250.;

//Z_R-plots
Int_t nXBinsZ_R = 2000;
Double_t firstXBinZ_R = -10.;
Double_t lastXBinZ_R = 10.;
Int_t nYBinsZ_R = 1000;
Double_t firstYBinZ_R = 0.;
Double_t lastYBinZ_R = 250.;

//X_Y-plots
Int_t nXBinsX_Y = 2000;
Double_t firstXBinX_Y = -250.;
Double_t lastXBinX_Y = 250.;
Int_t nYBinsX_Y = 2000;
Double_t firstYBinX_Y = -250.;
Double_t lastYBinX_Y = 250.;

//OpenAngle-plots
Int_t nXBinsOpeningAngle = 200;
Double_t firstXBinOpeningAngle = 0.;
Double_t lastXBinOpeningAngle = TMath::Pi()/2;

//Energy-plots
Int_t nXBinsEnergy = 500;
Double_t firstXBinEnergy = 0.;
Double_t lastXBinEnergy = 5.;

//Pt-plots
Int_t nXBinsPt = 500;
Double_t firstXBinPt = 0.;
Double_t lastXBinPt = 5.;

//Eta-plots
Int_t nXBinsEta = 400;
Double_t firstXBinEta = -2.;
Double_t lastXBinEta = 2.;

//Phi-plots
Int_t nXBinsPhi = 720;
Double_t firstXBinPhi = -TMath::Pi();
Double_t lastXBinPhi = TMath::Pi();

//Mapping-plots
Int_t nXBinsMapping = 40;
Double_t firstXBinMapping = -20.;
Double_t lastXBinMapping = 20.;
Int_t nYBinsMapping = 30;
Double_t firstYBinMapping = -1.5;
Double_t lastYBinMapping = 1.5;

//ResolutionPlots
//RES_dPt
Int_t nXBinsResdPt=500;
Int_t firstXBinResdPt= 0;
Int_t lastXBinResdPt=5;
Int_t nYBinsResdPt=1000;
Int_t firstYBinResdPt= -5;
Int_t lastYBinResdPt=5;

//RES_dR
Int_t nXBinsResdR=1000;
Int_t firstXBinResdR= 0;
Int_t lastXBinResdR=250;
Int_t nYBinsResdR=1000;
Int_t firstYBinResdR= -25;
Int_t lastYBinResdR=25;

//RES_dZ
Int_t nXBinsResdZ=2000;
Int_t firstXBinResdZ= -20;
Int_t lastXBinResdZ=20;
Int_t nYBinsResdZ=1000;
Int_t firstYBinResdZ= -20;
Int_t lastYBinResdZ=20;

//RES_dR_dPt
Int_t nXBinsResdR_dPt=1000;
Int_t firstXBinResdR_dPt= -22;
Int_t lastXBinResdR_dPt=22;
Int_t nYBinsResdR_dPt=1000;
Int_t firstYBinResdR_dPt= -5;
Int_t lastYBinResdR_dPt=5;


//RES_MC_Pt
Int_t nXBinsResPt=500;
Int_t firstXBinResPt= 0;
Int_t lastXBinResPt=5;

//RES_MC_R
Int_t nXBinsResR=1000;
Int_t firstXBinResR= 0;
Int_t lastXBinResR=250;

//RES_MC_Z
Int_t nXBinsResZ=1000;
Int_t firstXBinResZ= 0;
Int_t lastXBinResZ=250;

//GammaMass-plots
Int_t nXBinsGamma_Mass = 100;
Double_t firstXBinGamma_Mass = 0.;
Double_t lastXBinGamma_Mass = 1.;

//Pi0Mass-plots
Int_t nXBinsPi0_Mass = 100;
Double_t firstXBinPi0_Mass = 0.;
Double_t lastXBinPi0_Mass = 1.;

//EtaMass-plots
Int_t nXBinsEta_Mass = 100;
Double_t firstXBinEta_Mass = 0.;
Double_t lastXBinEta_Mass = 1.;

//Gamma_Width-plots
Int_t nXBinsGamma_Width = 100;
Double_t firstXBinGamma_Width = 0.;
Double_t lastXBinGamma_Width = 1.;

//Gamma_Chi2-plots
Int_t nXBinsGamma_Chi2 = 100;
Double_t firstXBinGamma_Chi2 = 0;
Double_t lastXBinGamma_Chi2 = 100.;

//Gamma_NDF-plots
Int_t nXBinsGamma_NDF = 10;
Double_t firstXBinGamma_NDF = 0.;
Double_t lastXBinGamma_NDF = 10.;
/** ---------- end Define the binning for the different plot types here ----------------------*/



/************************************************************************************************
 *                                                                                              *
 *                                                                                              *
 *                     EVERYTHING BELOW IS FOR DEVELOPERS ONLY                                  *
 *                                                                                              *
 *                                                                                              *
 ************************************************************************************************/

void ConfigGammaConversion(const char *chainName, const char *sample, int limit = 0){
  
  build();//build (if necessary) and load the libraries needed

  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C"); // load the CreateChain macro

  AliLog::SetGlobalLogLevel(AliLog::kError);

  //-------------------------------- Creating the histograms -------------------------------
  AliGammaConversionHistograms * histograms = new AliGammaConversionHistograms();

  //Change here
  // All the plots need to be given good names for the x-axis and y-axis so that people don't have to ask, or get confused.
  // I put in cm and cm in the first one as an example where to put in the axis titles.
  
  if(plotMC_EP_R != NULL){ histograms->Initialize_MC_EP_R(nXBinsR, firstXBinR, lastXBinR,"cm","cm");}
  if(plotMC_EP_Z_R != NULL){ histograms->Initialize_MC_EP_Z_R(nXBinsZ_R, firstXBinZ_R, lastXBinZ_R, nYBinsZ_R, firstYBinZ_R, lastYBinZ_R, "", "");}
  if(plotMC_EP_X_Y != NULL){ histograms->Initialize_MC_EP_X_Y(nXBinsX_Y, firstXBinX_Y, lastXBinX_Y, nYBinsX_Y, firstYBinX_Y, lastYBinX_Y, "", "");}
  if(plotMC_EP_OpeningAngle != NULL){ histograms->Initialize_MC_EP_OpeningAngle(nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}

  if(plotMC_E_Energy != NULL){ histograms->Initialize_MC_E_Energy(nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotMC_E_Pt != NULL){ histograms->Initialize_MC_E_Pt(nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotMC_E_Eta != NULL){ histograms->Initialize_MC_E_Eta(nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMC_E_Phi != NULL){ histograms->Initialize_MC_E_Phi(nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}

  if(plotMC_P_Energy != NULL){ histograms->Initialize_MC_P_Energy(nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotMC_P_Pt != NULL){ histograms->Initialize_MC_P_Pt(nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotMC_P_Eta != NULL){ histograms->Initialize_MC_P_Eta(nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMC_P_Phi != NULL){ histograms->Initialize_MC_P_Phi(nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}

  if(plotMC_Gamma_Energy != NULL){ histograms->Initialize_MC_Gamma_Energy(nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotMC_Gamma_Pt != NULL){ histograms->Initialize_MC_Gamma_Pt(nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotMC_Gamma_Eta != NULL){ histograms->Initialize_MC_Gamma_Eta(nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMC_Gamma_Phi != NULL){ histograms->Initialize_MC_Gamma_Phi(nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}

  if(plotMC_DirectGamma_Energy != NULL){ histograms->Initialize_MC_DirectGamma_Energy(nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotMC_DirectGamma_Pt != NULL){ histograms->Initialize_MC_DirectGamma_Pt(nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotMC_DirectGamma_Eta != NULL){ histograms->Initialize_MC_DirectGamma_Eta(nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMC_DirectGamma_Phi != NULL){ histograms->Initialize_MC_DirectGamma_Phi(nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}

  if(plotMC_Match_Gamma_Eta != NULL){ histograms->Initialize_MC_Match_Gamma_Eta(nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMC_Match_Gamma_Phi != NULL){ histograms->Initialize_MC_Match_Gamma_Phi(nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotMC_Match_Gamma_Pt != NULL){ histograms->Initialize_MC_Match_Gamma_Pt(nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotMC_Match_Gamma_Energy != NULL){ histograms->Initialize_MC_Match_Gamma_Energy(nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotMC_Match_Gamma_Mass != NULL){ histograms->Initialize_MC_Match_Gamma_Mass(nXBinsGamma_Mass, firstXBinGamma_Mass, lastXBinGamma_Mass, "", "");}
  if(plotMC_Match_Gamma_OpeningAngle != NULL){ histograms->Initialize_MC_Match_Gamma_OpeningAngle(nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}
  if(plotMC_Match_Gamma_R != NULL){ histograms->Initialize_MC_Match_Gamma_R(nXBinsR, firstXBinR, lastXBinR, "", "");}
  if(plotMC_Match_Gamma_Z_R != NULL){ histograms->Initialize_MC_Match_Gamma_Z_R(nXBinsZ_R, firstXBinZ_R, lastXBinZ_R, nYBinsZ_R, firstYBinZ_R, lastYBinZ_R, "", "");}
  if(plotMC_Match_Gamma_X_Y != NULL){ histograms->Initialize_MC_Match_Gamma_X_Y(nXBinsX_Y, firstXBinX_Y, lastXBinX_Y, nYBinsX_Y, firstYBinX_Y, lastYBinX_Y, "", "");}

  if(plotMC_Pi0_Eta != NULL){ histograms->Initialize_MC_Pi0_Eta(nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMC_Pi0_Phi != NULL){ histograms->Initialize_MC_Pi0_Phi(nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotMC_Pi0_Pt != NULL){ histograms->Initialize_MC_Pi0_Pt(nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotMC_Pi0_Energy != NULL){ histograms->Initialize_MC_Pi0_Energy(nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotMC_Pi0_Mass != NULL){ histograms->Initialize_MC_Pi0_Mass(nXBinsPi0_Mass, firstXBinPi0_Mass, lastXBinPi0_Mass, "", "");}
  if(plotMC_Pi0_OpeningAngle != NULL){ histograms->Initialize_MC_Pi0_OpeningAngleGamma(nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}
  if(plotMC_Pi0_R != NULL){ histograms->Initialize_MC_Pi0_R(nXBinsR, firstXBinR, lastXBinR, "", "");}
  if(plotMC_Pi0_Z_R != NULL){ histograms->Initialize_MC_Pi0_Z_R(nXBinsZ_R, firstXBinZ_R, lastXBinZ_R, nYBinsZ_R, firstYBinZ_R, lastYBinZ_R, "", "");}
  if(plotMC_Pi0_X_Y != NULL){ histograms->Initialize_MC_Pi0_X_Y(nXBinsX_Y, firstXBinX_Y, lastXBinX_Y, nYBinsX_Y, firstYBinX_Y, lastYBinX_Y, "", "");}

  if(plotMC_Eta_Eta != NULL){ histograms->Initialize_MC_Eta_Eta(nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMC_Eta_Phi != NULL){ histograms->Initialize_MC_Eta_Phi(nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotMC_Eta_Pt != NULL){ histograms->Initialize_MC_Eta_Pt(nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotMC_Eta_Energy != NULL){ histograms->Initialize_MC_Eta_Energy(nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotMC_Eta_Mass != NULL){ histograms->Initialize_MC_Eta_Mass(nXBinsEta_Mass, firstXBinEta_Mass, lastXBinEta_Mass, "", "");}
  if(plotMC_Eta_OpeningAngleGamma != NULL){ histograms->Initialize_MC_Eta_OpeningAngleGamma(nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}
  if(plotMC_Eta_R != NULL){ histograms->Initialize_MC_Eta_R(nXBinsR, firstXBinR, lastXBinR, "", "");}
  if(plotMC_Eta_Z_R != NULL){ histograms->Initialize_MC_Eta_Z_R(nXBinsZ_R, firstXBinZ_R, lastXBinZ_R, nYBinsZ_R, firstYBinZ_R, lastYBinZ_R, "", "");}
  if(plotMC_Eta_X_Y != NULL){ histograms->Initialize_MC_Eta_X_Y(nXBinsX_Y, firstXBinX_Y, lastXBinX_Y, nYBinsX_Y, firstYBinX_Y, lastYBinX_Y, "", "");}
    
  // Histograms from esd tracks
  if(plotESD_EP_R != NULL){ histograms->Initialize_ESD_EP_R(nXBinsR, firstXBinR, lastXBinR, "", "");}
  if(plotESD_EP_Z_R != NULL){ histograms->Initialize_ESD_EP_Z_R(nXBinsZ_R, firstXBinZ_R, lastXBinZ_R, nYBinsZ_R, firstYBinZ_R, lastYBinZ_R, "", "");}
  if(plotESD_EP_X_Y != NULL){ histograms->Initialize_ESD_EP_X_Y(nXBinsX_Y, firstXBinX_Y, lastXBinX_Y, nYBinsX_Y, firstYBinX_Y, lastYBinX_Y, "", "");}
  if(plotESD_EP_OpeningAngle != NULL){ histograms->Initialize_ESD_EP_OpeningAngle(nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}

  if(plotESD_E_Energy != NULL){ histograms->Initialize_ESD_E_Energy(nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotESD_E_Pt != NULL){ histograms->Initialize_ESD_E_Pt(nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotESD_E_Eta != NULL){ histograms->Initialize_ESD_E_Eta(nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotESD_E_Phi != NULL){ histograms->Initialize_ESD_E_Phi(nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}

  if(plotESD_P_Energy != NULL){ histograms->Initialize_ESD_P_Energy(nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotESD_P_Pt != NULL){ histograms->Initialize_ESD_P_Pt(nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotESD_P_Eta != NULL){ histograms->Initialize_ESD_P_Eta(nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotESD_P_Phi != NULL){ histograms->Initialize_ESD_P_Phi(nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}

  if(plotESD_Gamma_Energy != NULL){ histograms->Initialize_ESD_Gamma_Energy(nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotESD_Gamma_Pt != NULL){ histograms->Initialize_ESD_Gamma_Pt(nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotESD_Gamma_Eta != NULL){ histograms->Initialize_ESD_Gamma_Eta(nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotESD_Gamma_Phi != NULL){ histograms->Initialize_ESD_Gamma_Phi(nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}

  if(plotESD_Match_Gamma_OpeningAngle != NULL){ histograms->Initialize_ESD_Match_Gamma_OpeningAngle(nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}
  if(plotESD_Match_Gamma_Energy != NULL){ histograms->Initialize_ESD_Match_Gamma_Energy(nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotESD_Match_Gamma_Pt != NULL){ histograms->Initialize_ESD_Match_Gamma_Pt(nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotESD_Match_Gamma_Eta != NULL){ histograms->Initialize_ESD_Match_Gamma_Eta(nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotESD_Match_Gamma_Phi != NULL){ histograms->Initialize_ESD_Match_Gamma_Phi(nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotESD_Match_Gamma_Mass != NULL){ histograms->Initialize_ESD_Match_Gamma_Mass(nXBinsGamma_Mass, firstXBinGamma_Mass, lastXBinGamma_Mass, "", "");}
  if(plotESD_Match_Gamma_Width != NULL){ histograms->Initialize_ESD_Match_Gamma_Width(nXBinsGamma_Width, firstXBinGamma_Width, lastXBinGamma_Width, "", "");}
  if(plotESD_Match_Gamma_Chi2 != NULL){ histograms->Initialize_ESD_Match_Gamma_Chi2(nXBinsGamma_Chi2, firstXBinGamma_Chi2, lastXBinGamma_Chi2, "", "");}
  if(plotESD_Match_Gamma_NDF != NULL){ histograms->Initialize_ESD_Match_Gamma_NDF(nXBinsGamma_NDF, firstXBinGamma_NDF, lastXBinGamma_NDF, "", "");}
  if(plotESD_Match_Gamma_R != NULL){ histograms->Initialize_ESD_Match_Gamma_R(nXBinsR, firstXBinR, lastXBinR, "", "");}
  if(plotESD_Match_Gamma_Z_R != NULL){ histograms->Initialize_ESD_Match_Gamma_Z_R(nXBinsZ_R, firstXBinZ_R, lastXBinZ_R, nYBinsZ_R, firstYBinZ_R, lastYBinZ_R, "", "");}
  if(plotESD_Match_Gamma_X_Y != NULL){ histograms->Initialize_ESD_Match_Gamma_X_Y(nXBinsX_Y, firstXBinX_Y, lastXBinX_Y, nYBinsX_Y, firstYBinX_Y, lastYBinX_Y, "", "");}

  if(plotESD_Pi0_OpeningAngleGamma != NULL){ histograms->Initialize_ESD_Pi0_OpeningAngleGamma(nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}
  if(plotESD_Pi0_Energy != NULL){ histograms->Initialize_ESD_Pi0_Energy(nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotESD_Pi0_Pt != NULL){ histograms->Initialize_ESD_Pi0_Pt(nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotESD_Pi0_Eta != NULL){ histograms->Initialize_ESD_Pi0_Eta(nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotESD_Pi0_Phi != NULL){ histograms->Initialize_ESD_Pi0_Phi(nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotESD_Pi0_Mass != NULL){ histograms->Initialize_ESD_Pi0_Mass(nXBinsPi0_Mass, firstXBinPi0_Mass, lastXBinPi0_Mass, "", "");}
  if(plotESD_Pi0_R != NULL){ histograms->Initialize_ESD_Pi0_R(nXBinsR, firstXBinR, lastXBinR, "", "");}
  if(plotESD_Pi0_Z_R != NULL){ histograms->Initialize_ESD_Pi0_Z_R(nXBinsZ_R, firstXBinZ_R, lastXBinZ_R, nYBinsZ_R, firstYBinZ_R, lastYBinZ_R, "", "");}
  if(plotESD_Pi0_X_Y != NULL){ histograms->Initialize_ESD_Pi0_X_Y(nXBinsX_Y, firstXBinX_Y, lastXBinX_Y, nYBinsX_Y, firstYBinX_Y, lastYBinX_Y, "", "");}

  if(plotESD_Eta_OpeningAngleGamma != NULL){ histograms->Initialize_ESD_Eta_OpeningAngleGamma(nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}
  if(plotESD_Eta_Energy != NULL){ histograms->Initialize_ESD_Eta_Energy(nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotESD_Eta_Pt != NULL){ histograms->Initialize_ESD_Eta_Pt(nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotESD_Eta_Eta != NULL){ histograms->Initialize_ESD_Eta_Eta(nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotESD_Eta_Phi != NULL){ histograms->Initialize_ESD_Eta_Phi(nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotESD_Eta_Mass != NULL){ histograms->Initialize_ESD_Eta_Mass(nXBinsEta_Mass, firstXBinEta_Mass, lastXBinEta_Mass, "", "");}
  if(plotESD_Eta_R != NULL){ histograms->Initialize_ESD_Eta_R(nXBinsR, firstXBinR, lastXBinR, "", "");}
  if(plotESD_Eta_Z_R != NULL){ histograms->Initialize_ESD_Eta_Z_R(nXBinsZ_R, firstXBinZ_R, lastXBinZ_R, nYBinsZ_R, firstYBinZ_R, lastYBinZ_R, "", "");}
  if(plotESD_Eta_X_Y != NULL){ histograms->Initialize_ESD_Eta_X_Y(nXBinsX_Y, firstXBinX_Y, lastXBinX_Y, nYBinsX_Y, firstYBinX_Y, lastYBinX_Y, "", "");}

  if(plotESD_Background_OpeningAngleGamma != NULL){ histograms->Initialize_ESD_Background_OpeningAngleGamma(nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}
  if(plotESD_Background_Energy != NULL){ histograms->Initialize_ESD_Background_Energy(nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotESD_Background_Pt != NULL){ histograms->Initialize_ESD_Background_Pt(nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotESD_Background_Eta != NULL){ histograms->Initialize_ESD_Background_Eta(nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotESD_Background_Phi != NULL){ histograms->Initialize_ESD_Background_Phi(nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotESD_Background_Mass != NULL){ histograms->Initialize_ESD_Background_Mass(nXBinsEta_Mass, firstXBinEta_Mass, lastXBinEta_Mass, "", "");}
  if(plotESD_Background_R != NULL){ histograms->Initialize_ESD_Background_R(nXBinsR, firstXBinR, lastXBinR, "", "");}
  if(plotESD_Background_Z_R != NULL){ histograms->Initialize_ESD_Background_Z_R(nXBinsZ_R, firstXBinZ_R, lastXBinZ_R, nYBinsZ_R, firstYBinZ_R, lastYBinZ_R, "", "");}
  if(plotESD_Background_X_Y != NULL){ histograms->Initialize_ESD_Background_X_Y(nXBinsX_Y, firstXBinX_Y, lastXBinX_Y, nYBinsX_Y, firstYBinX_Y, lastYBinX_Y, "", "");}

  if(plotMapping != NULL){histograms->Initialize_MappingValues(nPhiIndex,nRIndex,nXBinsMapping,minRadius,maxRadius,nYBinsMapping,minPhi,maxPhi);}
  //create the mapping histograms
  if(plotMapping != NULL){histograms->Initialize_MappingHistograms(nPhiIndex,nRIndex,nXBinsMapping,firstXBinMapping,lastXBinMapping,nYBinsMapping,firstYBinMapping,lastYBinMapping,"", "");}




  
  if(plotResolution_dPt !=NULL){histograms->Initialize_Resolution_dPt(nXBinsResdPt, firstXBinResdPt, lastXBinResdPt, nYBinsResdPt, firstYBinResdPt, lastYBinResdPt, "", "");}
  if(plotResolution_dR !=NULL){histograms->Initialize_Resolution_dR(nXBinsResdR, firstXBinResdR, lastXBinResdR, nYBinsResdR, firstYBinResdR, lastYBinResdR, "", "");}
  if(plotResolution_dZ !=NULL){histograms->Initialize_Resolution_dZ(nXBinsResdZ, firstXBinResdZ, lastXBinResdZ, nYBinsResdZ, firstYBinResdZ, lastYBinResdZ, "", "");}
  
  if(plotResolution_dR_dPt !=NULL){histograms->Initialize_Resolution_dR_dPt(nXBinsResdR_dPt, firstXBinResdR_dPt, lastXBinResdR_dPt, nYBinsResdR_dPt, firstYBinResdR_dPt, lastYBinResdR_dPt, "", "");}
  
  if(plotResolution_MC_Pt !=NULL){histograms->Initialize_Resolution_MC_Pt(nXBinsResPt, firstXBinResPt, lastXBinResPt,"","");}
  if(plotResolution_MC_R !=NULL){histograms->Initialize_Resolution_MC_R(nXBinsResR, firstXBinResR, lastXBinResR,"","");}
  if(plotResolution_MC_Z !=NULL){histograms->Initialize_Resolution_MC_Z(nXBinsResZ, firstXBinResZ, lastXBinResZ,"","");}
  
  if(plotResolution_ESD_Pt !=NULL){histograms->Initialize_Resolution_ESD_Pt(nXBinsResPt, firstXBinResPt, lastXBinResPt,"","");}
  if(plotResolution_ESD_R !=NULL){histograms->Initialize_Resolution_ESD_R(nXBinsResR, firstXBinResR, lastXBinResR,"","");}
  if(plotResolution_ESD_Z !=NULL){histograms->Initialize_Resolution_ESD_Z(nXBinsResZ, firstXBinResZ, lastXBinResZ,"","");}
  
  if(plotNumberOfV0s != NULL){histograms->Initialize_NumberOfV0s(100, 0, 100,"","");}
  if(plotNumberOfSurvivingV0s != NULL){histograms->Initialize_NumberOfSurvivingV0s(100, 0, 100,"","");}

  //  debug histograms
  if(plotV0MassDebugCut1 != NULL){histograms->Initialize_V0MassDebugCut1(nXBinsGamma_Mass, firstXBinGamma_Mass, lastXBinGamma_Mass,"","");}
  if(plotV0MassDebugCut2 != NULL){histograms->Initialize_V0MassDebugCut2(nXBinsGamma_Mass, firstXBinGamma_Mass, lastXBinGamma_Mass,"","");}
  if(plotV0MassDebugCut3 != NULL){histograms->Initialize_V0MassDebugCut3(nXBinsGamma_Mass, firstXBinGamma_Mass, lastXBinGamma_Mass,"","");}
  if(plotV0MassDebugCut4 != NULL){histograms->Initialize_V0MassDebugCut4(nXBinsGamma_Mass, firstXBinGamma_Mass, lastXBinGamma_Mass,"","");}
  if(plotV0MassDebugCut5 != NULL){histograms->Initialize_V0MassDebugCut5(nXBinsGamma_Mass, firstXBinGamma_Mass, lastXBinGamma_Mass,"","");}
  if(plotV0MassDebugCut6 != NULL){histograms->Initialize_V0MassDebugCut6(nXBinsGamma_Mass, firstXBinGamma_Mass, lastXBinGamma_Mass,"","");}
  if(plotV0MassDebugCut7 != NULL){histograms->Initialize_V0MassDebugCut7(nXBinsGamma_Mass, firstXBinGamma_Mass, lastXBinGamma_Mass,"","");}
  if(plotV0MassDebugCut8 != NULL){histograms->Initialize_V0MassDebugCut8(nXBinsGamma_Mass, firstXBinGamma_Mass, lastXBinGamma_Mass,"","");}



  //------------------------------ end Creating the histograms -----------------------------

  // Create the Analysis manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager", "My Analysis");

  // Define Input Event Handler 
  AliESDInputHandler* inpHandler = new AliESDInputHandler();

  // Define Output Event Handler
  AliAODHandler* aodHandler = new AliAODHandler();
  aodHandler->SetOutputFileName("aodAliGammaConversion.root");
  
  // Define MC Truth Event Handler
  AliMCEventHandler* mcHandler = new AliMCEventHandler();
  
  // Add Handlers to the Task Manager
  mgr->SetInputEventHandler  (inpHandler);
  mgr->SetOutputEventHandler (aodHandler);
  mgr->SetMCtruthEventHandler(mcHandler);

  // Be sure you are told what you are doing
  mgr->SetDebugLevel(10);

  // Declare Common Input Tchain
  AliAnalysisDataContainer *cinput1 = mgr->CreateContainer("Chain",TChain::Class(),AliAnalysisManager::kInputContainer);

  // Common Output Tree in common ‘default’ output file
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("tree", TTree::Class(),AliAnalysisManager::kOutputContainer, "default");

  // Private output objects
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("histogramsAliGammaConversion", TList::Class(),AliAnalysisManager::kOutputContainer, "histogramsAliGammaConversion.root");


  //------------------------ END: Define input/output handlers ---------------------------------------------------


  //check for errors in the specified data
  if(useKFParticle == kTRUE && useESDTrack == kTRUE){
    //Print warning, cannot use both
  }
  if(useKFParticle == kFALSE && useESDTrack == kFALSE){
    //Print warning, one have to be specified
  }


  //Create the V0Reader
  AliV0Reader * v0Reader = new AliV0Reader();
  if(useKFParticle){
    v0Reader->UseKFParticle();
  }
  else if(useESDTrack){
    v0Reader->UseESDTrack();
  }
  v0Reader->SetNegativeTrackPID(pidOfNegativeTrack);
  v0Reader->SetPositiveTrackPID(pidOfPositiveTrack);
  v0Reader->SetMaxRCut(maxRCut);
  v0Reader->SetEtaCut(etaCut);
  v0Reader->SetPtCut(ptCut);
  v0Reader->SetChi2Cut(chi2Cut);
  v0Reader->SetPIDProbability(probElectron);
  v0Reader->SetXVertexCut(xVertexCut);
  v0Reader->SetYVertexCut(yVertexCut);
  v0Reader->SetZVertexCut(zVertexCut);
  v0Reader->SetSigmaMass(sigmaCutGammaMass);
  v0Reader->SetUseImprovedVertex(useImprovedVertex);

  // Create the GammaConversionTask
  AliAnalysisTaskGammaConversion *gammaconversion = new AliAnalysisTaskGammaConversion("GammaConversionTask");
  gammaconversion->SetDebugLevel(10);
  
  gammaconversion->SetV0Reader(v0Reader);
  gammaconversion->SetCalculateBackground(calculateBackground);
  gammaconversion->Init();
  
  gammaconversion->SetElectronMass(electronMass);
  gammaconversion->SetGammaMass(gammaMass);
  gammaconversion->SetPi0Mass(pi0Mass);
  gammaconversion->SetEtaMass(etaMass);

  gammaconversion->SetGammaWidth(gammaWidth);
  gammaconversion->SetPi0Width(pi0Width);
  gammaconversion->SetEtaWidth(etaWidth);

  // define the width constraint used by KF particle.
  Double_t gammaWidth = 0.01;
  Double_t pi0Width   = 0.01;
  Double_t etaWidth   = 0.01;

  gammaconversion->SetHistograms(histograms);
  v0Reader->SetHistograms(histograms);// also give the pointer to the v0reader, for debugging cuts
  
  gammaconversion->SetDoMCTruth(doMCTruth);


  // Add task to the manager 
  mgr->AddTask(gammaconversion);

  // Connect I/O to the task
  mgr->ConnectInput (gammaconversion, 0, cinput1);
  mgr->ConnectOutput(gammaconversion, 0, coutput1);
  mgr->ConnectOutput(gammaconversion, 1, coutput2);

  TChain* chain= CreateESDChain(sample);
  
  mgr->InitAnalysis();
  
  mgr->PrintStatus();
  
  mgr->StartAnalysis("local",chain);

}




void build() {
  TStopwatch timer;
  timer.Start();
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom");
  //  gSystem->Load("libANALYSISalice");

  //____________________________________________________//
  //____________________________________________________//
  //_____________Setting up STEERBase.par_______________//
  //____________________________________________________//
  setupPar("STEERBase");
  gSystem->Load("libSTEERBase.so");

  //____________________________________________________//
  //_____________Setting up ESD.par_____________________//
  //____________________________________________________//
  setupPar("ESD");
  gSystem->Load("libVMC.so");
  gSystem->Load("libESD.so");

  //____________________________________________________//
  //_____________Setting up AOD.par_____________________//
  //____________________________________________________//
  setupPar("AOD");
  gSystem->Load("libAOD.so");
                                                                
  //_____________________________________________________________//
  //_____________Setting up ANALYSIS.par_________________________//
  //_____________________________________________________________//
  setupPar("ANALYSIS");
  gSystem->Load("libANALYSIS.so");

  //_____________________________________________________________//
  //_____________Setting up ANALYSISalice.par_________________________//
  //_____________________________________________________________//
  setupPar("ANALYSISalice");
  gSystem->Load("libANALYSISalice.so");
                                                                                                                                  
  //_____________________________________________________________//
  //_____________Setting up PWG4Gamma.par_____________________//
  //_____________________________________________________________//
  //  setupPar("PWG4Gamma");
  //  gSystem->Load("libPWG4Gamma.so");
  setupPar("PWG4PartCorr");
  gSystem->Load("libPWG4PartCorr.so");
  //if head:: use PWG4PartCorr
             
  //gROOT->LoadMacro("AliAnalysisTaskPi0.cxx+");
  // gROOT->LoadMacro("AliAnalysisTaskPtMC.cxx+");
  //  gROOT->LoadMacro("AliAnalysisTaskPi0MC.cxx+");
  //  gROOT->LoadMacro("AliAnaScale.cxx+");


  //gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  //TChain* chain = CreateESDChain("files1.txt");



  //____________________________________________//


}

Int_t setupPar(const char* pararchivename) {
  ///////////////////
  // Setup PAR File//
  ///////////////////
  if (pararchivename) {
    char processline[1024];
    sprintf(processline,".! tar xvzf %s.par",pararchivename);
    gROOT->ProcessLine(processline);
    const char* ocwd = gSystem->WorkingDirectory();
    gSystem->ChangeDirectory(pararchivename);
                                                                                                                                               
    // check for BUILD.sh and execute
    if (!gSystem->AccessPathName("PROOF-INF/BUILD.sh")) {
      printf("*******************************\n");
      printf("*** Building PAR archive    ***\n");
      printf("*******************************\n");
                                                                                                                                               
      if (gSystem->Exec("PROOF-INF/BUILD.sh")) {
        Error("runAnalysis","Cannot Build the PAR Archive! - Abort!");
        return -1;
      }
    }
    // check for SETUP.C and execute
    if (!gSystem->AccessPathName("PROOF-INF/SETUP.C")) {
      printf("*******************************\n");
      printf("*** Setup PAR archive       ***\n");
      printf("*******************************\n");
      gROOT->Macro("PROOF-INF/SETUP.C");
    }
                                                                                                                                               
    gSystem->ChangeDirectory("../");
  }                                                                                                                                               
  return 1;
}
