/** VERSION NUMBER 0 */

Bool_t usePWG4PartCorr = kTRUE;


/** ------------------------------ Monte Carlo flag -----------------------------------------*/
Bool_t doMCTruth = kTRUE;
/** ---------------------------- end Monte Carlo flag ---------------------------------------*/

/** ------------------------- Choose KFParticle OR ESDTrack  --------------------------------*/
Bool_t useKFParticle = kTRUE;
Bool_t useESDTrack   = kFALSE;
/** ----------------------- end Choose KFParticle OR ESDTrack  -----------------------------*/


Bool_t calculateBackground = kTRUE;

Int_t numberOfFilesToAnalyze=0;

/** ---------------------------------- define cuts here ------------------------------------*/

Int_t pidOfNegativeTrack=11;
Int_t pidOfPositiveTrack=-11;

Double_t maxRCut   = 200.;
Double_t etaCut    = 1.2;
Double_t ptCut     = 0.1;
Double_t chi2CutConversion   = 20.;
Double_t chi2CutMeson   = 20.;

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
Bool_t plotMCEPR                                           = kTRUE;
Bool_t plotMCEPZR                                          = kTRUE;
Bool_t plotMCEPXY                                          = kTRUE;
Bool_t plotMCEPOpeningAngle                                = kTRUE;

Bool_t plotMCEEnergy                                       = kTRUE;
Bool_t plotMCEPt                                           = kTRUE;
Bool_t plotMCEEta                                          = kTRUE;
Bool_t plotMCEPhi                                          = kTRUE;

Bool_t plotMCPEnergy                                       = kTRUE;
Bool_t plotMCPPt                                           = kTRUE;
Bool_t plotMCPEta                                          = kTRUE;
Bool_t plotMCPPhi                                          = kTRUE;

Bool_t plotMCGammaEnergy                                   = kTRUE;
Bool_t plotMCGammaPt                                       = kTRUE;
Bool_t plotMCGammaEta                                      = kTRUE;
Bool_t plotMCGammaPhi                                      = kTRUE;

Bool_t plotMCDirectGammaEnergy                             = kTRUE;
Bool_t plotMCDirectGammaPt                                 = kTRUE;
Bool_t plotMCDirectGammaEta                                = kTRUE;
Bool_t plotMCDirectGammaPhi                                = kTRUE;

Bool_t plotMCMatchGammaEta                                 = kTRUE;
Bool_t plotMCMatchGammaPhi                                 = kTRUE;
Bool_t plotMCMatchGammaPt                                  = kTRUE;
Bool_t plotMCMatchGammaEnergy                              = kTRUE;
Bool_t plotMCMatchGammaMass                                = kTRUE;
Bool_t plotMCMatchGammaOpeningAngle                        = kTRUE;
Bool_t plotMCMatchGammaR                                   = kTRUE;
Bool_t plotMCMatchGammaZR                                  = kTRUE;
Bool_t plotMCMatchGammaXY                                  = kTRUE;

Bool_t plotMCPi0Eta                                        = kTRUE;
Bool_t plotMCPi0Phi                                        = kTRUE;
Bool_t plotMCPi0Pt                                         = kTRUE;
Bool_t plotMCPi0Energy                                     = kTRUE;
Bool_t plotMCPi0Mass                                       = kTRUE;
Bool_t plotMCPi0OpeningAngle                               = kTRUE;
Bool_t plotMCPi0R                                          = kTRUE;
Bool_t plotMCPi0ZR                                         = kTRUE;
Bool_t plotMCPi0XY                                         = kTRUE;

Bool_t plotMCEtaEta                                        = kTRUE;
Bool_t plotMCEtaPhi                                        = kTRUE;
Bool_t plotMCEtaPt                                         = kTRUE;
Bool_t plotMCEtaEnergy                                     = kTRUE;
Bool_t plotMCEtaMass                                       = kTRUE;
Bool_t plotMCEtaOpeningAngleGamma                          = kTRUE;
Bool_t plotMCEtaR                                          = kTRUE;
Bool_t plotMCEtaZR                                         = kTRUE;
Bool_t plotMCEtaXY                                         = kTRUE;
    
// Histograms from esd tracks
Bool_t plotESDEPR                                          = kTRUE;
Bool_t plotESDEPZR                                         = kTRUE;
Bool_t plotESDEPXY                                         = kTRUE;
Bool_t plotESDEPOpeningAngle                               = kTRUE;

Bool_t plotESDEEnergy                                      = kTRUE;
Bool_t plotESDEPt                                          = kTRUE;
Bool_t plotESDEEta                                         = kTRUE;
Bool_t plotESDEPhi                                         = kTRUE;

Bool_t plotESDPEnergy                                      = kTRUE;
Bool_t plotESDPPt                                          = kTRUE;
Bool_t plotESDPEta                                         = kTRUE;
Bool_t plotESDPPhi                                         = kTRUE;

Bool_t plotESDGammaEnergy                                  = kTRUE;
Bool_t plotESDGammaPt                                      = kTRUE;
Bool_t plotESDGammaEta                                     = kTRUE;
Bool_t plotESDGammaPhi                                     = kTRUE;

Bool_t plotESDMatchGammaOpeningAngle                       = kTRUE;
Bool_t plotESDMatchGammaEnergy                             = kTRUE;
Bool_t plotESDMatchGammaPt                                 = kTRUE;
Bool_t plotESDMatchGammaEta                                = kTRUE;
Bool_t plotESDMatchGammaPhi                                = kTRUE;
Bool_t plotESDMatchGammaMass                               = kTRUE;
Bool_t plotESDMatchGammaWidth                              = kTRUE;
Bool_t plotESDMatchGammaChi2                               = kTRUE;
Bool_t plotESDMatchGammaNDF                                = kTRUE;
Bool_t plotESDMatchGammaR                                  = kTRUE;
Bool_t plotESDMatchGammaZR                                 = kTRUE;
Bool_t plotESDMatchGammaXY                                 = kTRUE;

Bool_t plotESDTwoGammaCombinationOpeningAngleGamma         = kTRUE;
Bool_t plotESDTwoGammaCombinationEnergy                    = kTRUE;
Bool_t plotESDTwoGammaCombinationPt                        = kTRUE;
Bool_t plotESDTwoGammaCombinationEta                       = kTRUE;
Bool_t plotESDTwoGammaCombinationPhi                       = kTRUE;
Bool_t plotESDTwoGammaCombinationMass                      = kTRUE;
Bool_t plotESDTwoGammaCombinationR                         = kTRUE;
Bool_t plotESDTwoGammaCombinationZR                        = kTRUE;
Bool_t plotESDTwoGammaCombinationXY                        = kTRUE;

Bool_t plotESDBackgroundOpeningAngleGamma                  = kTRUE;
Bool_t plotESDBackgroundEnergy                             = kTRUE;
Bool_t plotESDBackgroundPt                                 = kTRUE;
Bool_t plotESDBackgroundEta                                = kTRUE;
Bool_t plotESDBackgroundPhi                                = kTRUE;
Bool_t plotESDBackgroundMass                               = kTRUE;
Bool_t plotESDBackgroundR                                  = kTRUE;
Bool_t plotESDBackgroundZR                                 = kTRUE;
Bool_t plotESDBackgroundXY                                 = kTRUE;

Bool_t plotMapping                                         = kTRUE;       

Bool_t plotResolutiondPt                                   = kTRUE;
Bool_t plotResolutiondR                                    = kTRUE;
Bool_t plotResolutiondZ                                    = kTRUE;
  
Bool_t plotResolutiondRdPt                                 = kTRUE;

Bool_t plotResolutionMCPt                                  = kTRUE;
Bool_t plotResolutionMCR                                   = kTRUE;
Bool_t plotResolutionMCZ                                   = kTRUE;

Bool_t plotResolutionESDPt                                 = kTRUE;
Bool_t plotResolutionESDR                                  = kTRUE;
Bool_t plotResolutionESDZ                                  = kTRUE;

Bool_t plotNumberOfV0s                                     = kTRUE;
Bool_t plotNumberOfSurvivingV0s                            = kTRUE;

//  debug histograms
Bool_t plotV0MassDebugCut1                                 = kTRUE;
Bool_t plotV0MassDebugCut2                                 = kTRUE;
Bool_t plotV0MassDebugCut3                                 = kTRUE;
Bool_t plotV0MassDebugCut4                                 = kTRUE;
Bool_t plotV0MassDebugCut5                                 = kTRUE;
Bool_t plotV0MassDebugCut6                                 = kTRUE;
Bool_t plotV0MassDebugCut7                                 = kTRUE;
Bool_t plotV0MassDebugCut8                                 = kTRUE;

Bool_t plotPi0Spectra                                      = kTRUE;
Bool_t plotEtaSpectra                                      = kTRUE;


/** ----------------- end define which histograms to plot here -------------------------------*/



/** ----------- Define the binning for the different plot types here -------------------------*/
//R-plots
Int_t nXBinsR = 1000;
Double_t firstXBinR = 0.;
Double_t lastXBinR = 250.;

//ZR-plots
Int_t nXBinsZR = 2000;
Double_t firstXBinZR = -10.;
Double_t lastXBinZR = 10.;
Int_t nYBinsZR = 1000;
Double_t firstYBinZR = 0.;
Double_t lastYBinZR = 250.;

//XY-plots
Int_t nXBinsXY = 2000;
Double_t firstXBinXY = -250.;
Double_t lastXBinXY = 250.;
Int_t nYBinsXY = 2000;
Double_t firstYBinXY = -250.;
Double_t lastYBinXY = 250.;

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
//RESdPt
Int_t nXBinsResdPt=500;
Int_t firstXBinResdPt= 0;
Int_t lastXBinResdPt=5;
Int_t nYBinsResdPt=1000;
Int_t firstYBinResdPt= -5;
Int_t lastYBinResdPt=5;

//RESdR
Int_t nXBinsResdR=1000;
Int_t firstXBinResdR= 0;
Int_t lastXBinResdR=250;
Int_t nYBinsResdR=1000;
Int_t firstYBinResdR= -25;
Int_t lastYBinResdR=25;

//RESdZ
Int_t nXBinsResdZ=2000;
Int_t firstXBinResdZ= -20;
Int_t lastXBinResdZ=20;
Int_t nYBinsResdZ=1000;
Int_t firstYBinResdZ= -20;
Int_t lastYBinResdZ=20;

//RESdRdPt
Int_t nXBinsResdRdPt=1000;
Int_t firstXBinResdRdPt= -22;
Int_t lastXBinResdRdPt=22;
Int_t nYBinsResdRdPt=1000;
Int_t firstYBinResdRdPt= -5;
Int_t lastYBinResdRdPt=5;

//RESMCPt
Int_t nXBinsResPt=500;
Int_t firstXBinResPt= 0;
Int_t lastXBinResPt=5;

//RESMCR
Int_t nXBinsResR=1000;
Int_t firstXBinResR= 0;
Int_t lastXBinResR=250;

//RESMCZ
Int_t nXBinsResZ=1000;
Int_t firstXBinResZ= 0;
Int_t lastXBinResZ=250;

//GammaMass-plots
Int_t nXBinsGammaMass = 100;
Double_t firstXBinGammaMass = 0.;
Double_t lastXBinGammaMass = 1.;

//Pi0Mass-plots
Int_t nXBinsPi0Mass = 100;
Double_t firstXBinPi0Mass = 0.;
Double_t lastXBinPi0Mass = 1.;

//EtaMass-plots
Int_t nXBinsEtaMass = 100;
Double_t firstXBinEtaMass = 0.;
Double_t lastXBinEtaMass = 1.;

//GammaWidth-plots
Int_t nXBinsGammaWidth = 100;
Double_t firstXBinGammaWidth = 0.;
Double_t lastXBinGammaWidth = 1.;

//GammaChi2-plots
Int_t nXBinsGammaChi2 = 100;
Double_t firstXBinGammaChi2 = 0;
Double_t lastXBinGammaChi2 = 100.;

//GammaNDF-plots
Int_t nXBinsGammaNDF = 10;
Double_t firstXBinGammaNDF = 0.;
Double_t lastXBinGammaNDF = 10.;

//Spectra-plots
Int_t nXBinsSpectra = 100;
Double_t firstXBinSpectra = 0.;
Double_t lastXBinSpectra = 1.;
Int_t nYBinsSpectra = 500;
Double_t firstYBinSpectra = 0.;
Double_t lastYBinSpectra = 100.;

/** ---------- end Define the binning for the different plot types here ----------------------*/


/************************************************************************************************
 *                                                                                              *
 *                                                                                              *
 *                     EVERYTHING BELOW IS FOR DEVELOPERS ONLY                                  *
 *                                                                                              *
 *                                                                                              *
 ************************************************************************************************/
TString outputFileName = "histogramsGammaConversion";
TString outputFileAppendix = "";
TString dataList = "";
Bool_t writeNtuple = kFALSE;

Bool_t scanArguments(TString arguments){
  Bool_t iResult = kTRUE;

  TString allArgs=arguments;
  TString argument;
  int bMissingParam=0;

  TObjArray* pTokens=allArgs.Tokenize(" ");
  if (pTokens) {
    
    for(int i=0; i<pTokens->GetEntries() && iResult==kTRUE; i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();
      
      if(argument.IsNull()) continue;
      // -- deconvolute-time option
      if(argument.CompareTo("-data-list") == 0){
	if((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	dataList = ((TObjString*)pTokens->At(i))->GetString();
	if(dataList.IsNull()){
	  cout<<"-data-list is NULL"<<endl;
	  iResult=kFALSE;
	}
	else{
	  cout<<"Data list is set to: "<<dataList<<endl;
	}
      }
      else if(argument.CompareTo("-output-file-name") == 0){
	if((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	outputFileName = ((TObjString*)pTokens->At(i))->GetString();
	if(outputFileName.IsNull()){
	  cout<<"-output-file-name is NULL"<<endl;
	  iResult=kFALSE;
	}
	else{
	  cout<<"Setting output file name to: "<<outputFileName<<endl;
	}
      }
      else if (argument.CompareTo("-write-ntuple") == 0){
	cout<<"Writing ntuple to file."<<endl;
	writeNtuple = kTRUE;
      }
      else if(argument.CompareTo("-append-to-output-file") == 0){
	if((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	outputFileAppendix = "_"+((TObjString*)pTokens->At(i))->GetString();
	if(outputFileAppendix.IsNull()){
	  cout<<"-appending-to-output-file is NULL"<<endl;
	  iResult=kFALSE;
	}
	else{
	  cout<<"Appending to the output file: "<<outputFileAppendix<<endl;
	}
      }
    }

    delete pTokens;
  }
  if (bMissingParam) {
    cout<<"Missing parameter for argument "<< argument.Data()<<endl;
    iResult=kFALSE;
  }
  return iResult;
}

void ConfigGammaConversion(TString arguments){

  if(!scanArguments(arguments)){
    break;
  }

  if(numberOfFilesToAnalyze==0){
    ifstream dataInStream;
    dataInStream.open(dataList.Data());
    if ( !dataInStream ){
      cout<<"Data list file does not exist: "<<dataList.Data()<<endl;
      return 0;
    }
    string line;
    while ( !dataInStream.eof() )
      {
	getline(dataInStream, line);
	numberOfFilesToAnalyze++;
      }
  }
  cout<<"Number Of files to analyze: "<<numberOfFilesToAnalyze<<endl;

  build();//build (if necessary) and load the libraries needed

  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C"); // load the CreateChain macro

  AliLog::SetGlobalLogLevel(AliLog::kError);

  //-------------------------------- Creating the histograms -------------------------------
  AliGammaConversionHistograms * histograms = new AliGammaConversionHistograms();

  if(plotMCEPR == kTRUE){ histograms->AddHistogram("MC_EP_R","Radius of gamma conversion points",nXBinsR, firstXBinR, lastXBinR,"counts","cm");}
  if(plotMCEPZR == kTRUE){ histograms->AddHistogram("MC_EP_ZR","Radius of gamma conversion points vs Z",nXBinsZR, firstXBinZR, lastXBinZR, nYBinsZR, firstYBinZR, lastYBinZR, "cm", "cm");}
  if(plotMCEPXY == kTRUE){ histograms->AddHistogram("MC_EP_XY","Gamma XY converison point.",nXBinsXY, firstXBinXY, lastXBinXY, nYBinsXY, firstYBinXY, lastYBinXY, "cm", "cm");}
  if(plotMCEPOpeningAngle == kTRUE){ histograms->AddHistogram("MC_EP_OpeningAngle","Opening angle of e+e- pairs from gamma conversion",nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "counts", "cm");}

  if(plotMCEEnergy == kTRUE){ histograms->AddHistogram("MC_E_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotMCEPt == kTRUE){ histograms->AddHistogram("MC_E_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotMCEEta == kTRUE){ histograms->AddHistogram("MC_E_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMCEPhi == kTRUE){ histograms->AddHistogram("MC_E_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}

  if(plotMCPEnergy == kTRUE){ histograms->AddHistogram("MC_P_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotMCPPt == kTRUE){ histograms->AddHistogram("MC_P_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotMCPEta == kTRUE){ histograms->AddHistogram("MC_P_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMCPPhi == kTRUE){ histograms->AddHistogram("MC_P_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}

  if(plotMCGammaEnergy == kTRUE){ histograms->AddHistogram("MC_Gamma_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotMCGammaPt == kTRUE){ histograms->AddHistogram("MC_Gamma_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotMCGammaEta == kTRUE){ histograms->AddHistogram("MC_Gamma_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMCGammaPhi == kTRUE){ histograms->AddHistogram("MC_Gamma_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}

  if(plotMCDirectGammaEnergy == kTRUE){ histograms->AddHistogram("MC_DirectGamma_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotMCDirectGammaPt == kTRUE){ histograms->AddHistogram("MC_DirectGamma_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotMCDirectGammaEta == kTRUE){ histograms->AddHistogram("MC_DirectGamma_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMCDirectGammaPhi == kTRUE){ histograms->AddHistogram("MC_DirectGamma_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}

  if(plotMCMatchGammaEta == kTRUE){ histograms->AddHistogram("MC_Match_Gamma_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMCMatchGammaPhi == kTRUE){ histograms->AddHistogram("MC_Match_Gamma_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotMCMatchGammaPt == kTRUE){ histograms->AddHistogram("MC_Match_Gamma_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotMCMatchGammaEnergy == kTRUE){ histograms->AddHistogram("MC_Match_Gamma_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotMCMatchGammaMass == kTRUE){ histograms->AddHistogram("MC_Match_Gamma_Mass" ,"" , nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass, "", "");}
  if(plotMCMatchGammaOpeningAngle == kTRUE){ histograms->AddHistogram("MC_Match_Gamma_OpeningAngle" ,"" , nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}
  if(plotMCMatchGammaR == kTRUE){ histograms->AddHistogram("MC_Match_Gamma_R" ,"" , nXBinsR, firstXBinR, lastXBinR, "", "");}
  if(plotMCMatchGammaZR == kTRUE){ histograms->AddHistogram("MC_Match_Gamma_ZR" ,"" , nXBinsZR, firstXBinZR, lastXBinZR, nYBinsZR, firstYBinZR, lastYBinZR, "", "");}
  if(plotMCMatchGammaXY == kTRUE){ histograms->AddHistogram("MC_Match_Gamma_XY" ,"" , nXBinsXY, firstXBinXY, lastXBinXY, nYBinsXY, firstYBinXY, lastYBinXY, "", "");}

  if(plotMCPi0Eta == kTRUE){ histograms->AddHistogram("MC_Pi0_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMCPi0Phi == kTRUE){ histograms->AddHistogram("MC_Pi0_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotMCPi0Pt == kTRUE){ histograms->AddHistogram("MC_Pi0_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotMCPi0Energy == kTRUE){ histograms->AddHistogram("MC_Pi0_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotMCPi0Mass == kTRUE){ histograms->AddHistogram("MC_Pi0_Mass" ,"" , nXBinsPi0Mass, firstXBinPi0Mass, lastXBinPi0Mass, "", "");}
  if(plotMCPi0OpeningAngle == kTRUE){ histograms->AddHistogram("MC_Pi0_GammaDaughter_OpeningAngle" ,"" , nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}
  if(plotMCPi0R == kTRUE){ histograms->AddHistogram("MC_Pi0_R" ,"" , nXBinsR, firstXBinR, lastXBinR, "", "");}
  if(plotMCPi0ZR == kTRUE){ histograms->AddHistogram("MC_Pi0_ZR" ,"" , nXBinsZR, firstXBinZR, lastXBinZR, nYBinsZR, firstYBinZR, lastYBinZR, "", "");}
  if(plotMCPi0XY == kTRUE){ histograms->AddHistogram("MC_Pi0_XY" ,"" , nXBinsXY, firstXBinXY, lastXBinXY, nYBinsXY, firstYBinXY, lastYBinXY, "", "");}

  if(plotMCPi0Eta == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMCPi0Phi == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotMCPi0Pt == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotMCPi0Energy == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotMCPi0Mass == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_Mass" ,"" , nXBinsPi0Mass, firstXBinPi0Mass, lastXBinPi0Mass, "", "");}
  if(plotMCPi0OpeningAngle == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_GammaDaughter_OpeningAngle" ,"" , nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}
  if(plotMCPi0R == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_R" ,"" , nXBinsR, firstXBinR, lastXBinR, "", "");}
  if(plotMCPi0ZR == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_ZR" ,"" , nXBinsZR, firstXBinZR, lastXBinZR, nYBinsZR, firstYBinZR, lastYBinZR, "", "");}
  if(plotMCPi0XY == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_XY" ,"" , nXBinsXY, firstXBinXY, lastXBinXY, nYBinsXY, firstYBinXY, lastYBinXY, "", "");}

  if(plotMCEtaEta == kTRUE){ histograms->AddHistogram("MC_Eta_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMCEtaPhi == kTRUE){ histograms->AddHistogram("MC_Eta_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotMCEtaPt == kTRUE){ histograms->AddHistogram("MC_Eta_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotMCEtaEnergy == kTRUE){ histograms->AddHistogram("MC_Eta_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotMCEtaMass == kTRUE){ histograms->AddHistogram("MC_Eta_Mass" ,"" , nXBinsEtaMass, firstXBinEtaMass, lastXBinEtaMass, "", "");}
  if(plotMCEtaOpeningAngleGamma == kTRUE){ histograms->AddHistogram("MC_Eta_GammaDaughter_OpeningAngle" ,"" , nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}
  if(plotMCEtaR == kTRUE){ histograms->AddHistogram("MC_Eta_R" ,"" , nXBinsR, firstXBinR, lastXBinR, "", "");}
  if(plotMCEtaZR == kTRUE){ histograms->AddHistogram("MC_Eta_ZR" ,"" , nXBinsZR, firstXBinZR, lastXBinZR, nYBinsZR, firstYBinZR, lastYBinZR, "", "");}
  if(plotMCEtaXY == kTRUE){ histograms->AddHistogram("MC_Eta_XY" ,"" , nXBinsXY, firstXBinXY, lastXBinXY, nYBinsXY, firstYBinXY, lastYBinXY, "", "");}
    
  // Histograms from esd tracks
  if(plotESDEPR == kTRUE){ histograms->AddHistogram("ESD_EP_R" ,"" , nXBinsR, firstXBinR, lastXBinR, "", "");}
  if(plotESDEPZR == kTRUE){ histograms->AddHistogram("ESD_EP_ZR" ,"" , nXBinsZR, firstXBinZR, lastXBinZR, nYBinsZR, firstYBinZR, lastYBinZR, "", "");}
  if(plotESDEPXY == kTRUE){ histograms->AddHistogram("ESD_EP_XY" ,"" , nXBinsXY, firstXBinXY, lastXBinXY, nYBinsXY, firstYBinXY, lastYBinXY, "", "");}
  if(plotESDEPOpeningAngle == kTRUE){ histograms->AddHistogram("ESD_EP_OpeningAngle" ,"" , nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}

  if(plotESDEEnergy == kTRUE){ histograms->AddHistogram("ESD_E_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotESDEPt == kTRUE){ histograms->AddHistogram("ESD_E_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotESDEEta == kTRUE){ histograms->AddHistogram("ESD_E_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotESDEPhi == kTRUE){ histograms->AddHistogram("ESD_E_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}

  if(plotESDPEnergy == kTRUE){ histograms->AddHistogram("ESD_P_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotESDPPt == kTRUE){ histograms->AddHistogram("ESD_P_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotESDPEta == kTRUE){ histograms->AddHistogram("ESD_P_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotESDPPhi == kTRUE){ histograms->AddHistogram("ESD_P_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}

  if(plotESDGammaEnergy == kTRUE){ histograms->AddHistogram("ESD_Gamma_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotESDGammaPt == kTRUE){ histograms->AddHistogram("ESD_Gamma_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotESDGammaEta == kTRUE){ histograms->AddHistogram("ESD_Gamma_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotESDGammaPhi == kTRUE){ histograms->AddHistogram("ESD_Gamma_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}

  if(plotESDMatchGammaOpeningAngle == kTRUE){ histograms->AddHistogram("ESD_Match_Gamma_OpeningAngle" ,"" , nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}
  if(plotESDMatchGammaEnergy == kTRUE){ histograms->AddHistogram("ESD_Match_Gamma_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotESDMatchGammaPt == kTRUE){ histograms->AddHistogram("ESD_Match_Gamma_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotESDMatchGammaEta == kTRUE){ histograms->AddHistogram("ESD_Match_Gamma_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotESDMatchGammaPhi == kTRUE){ histograms->AddHistogram("ESD_Match_Gamma_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotESDMatchGammaMass == kTRUE){ histograms->AddHistogram("ESD_Match_Gamma_Mass" ,"" , nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass, "", "");}
  if(plotESDMatchGammaWidth == kTRUE){ histograms->AddHistogram("ESD_Match_Gamma_Width" ,"" , nXBinsGammaWidth, firstXBinGammaWidth, lastXBinGammaWidth, "", "");}
  if(plotESDMatchGammaChi2 == kTRUE){ histograms->AddHistogram("ESD_Match_Gamma_Chi2" ,"" , nXBinsGammaChi2, firstXBinGammaChi2, lastXBinGammaChi2, "", "");}
  if(plotESDMatchGammaNDF == kTRUE){ histograms->AddHistogram("ESD_Match_Gamma_NDF" ,"" , nXBinsGammaNDF, firstXBinGammaNDF, lastXBinGammaNDF, "", "");}
  if(plotESDMatchGammaR == kTRUE){ histograms->AddHistogram("ESD_Match_Gamma_R" ,"" , nXBinsR, firstXBinR, lastXBinR, "", "");}
  if(plotESDMatchGammaZR == kTRUE){ histograms->AddHistogram("ESD_Match_Gamma_ZR" ,"" , nXBinsZR, firstXBinZR, lastXBinZR, nYBinsZR, firstYBinZR, lastYBinZR, "", "");}
  if(plotESDMatchGammaXY == kTRUE){ histograms->AddHistogram("ESD_Match_Gamma_XY" ,"" , nXBinsXY, firstXBinXY, lastXBinXY, nYBinsXY, firstYBinXY, lastYBinXY, "", "");}

  if(plotESDTwoGammaCombinationOpeningAngleGamma == kTRUE){ histograms->AddHistogram("ESD_TwoGammaCombination_GammaDaughter_OpeningAngle" ,"" , nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}
  if(plotESDTwoGammaCombinationEnergy == kTRUE){ histograms->AddHistogram("ESD_TwoGammaCombination_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotESDTwoGammaCombinationPt == kTRUE){ histograms->AddHistogram("ESD_TwoGammaCombination_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotESDTwoGammaCombinationEta == kTRUE){ histograms->AddHistogram("ESD_TwoGammaCombination_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotESDTwoGammaCombinationPhi == kTRUE){ histograms->AddHistogram("ESD_TwoGammaCombination_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotESDTwoGammaCombinationMass == kTRUE){ histograms->AddHistogram("ESD_TwoGammaCombination_Mass" ,"" , nXBinsPi0Mass, firstXBinPi0Mass, lastXBinPi0Mass, "", "");}
  if(plotESDTwoGammaCombinationR == kTRUE){ histograms->AddHistogram("ESD_TwoGammaCombination_R" ,"" , nXBinsR, firstXBinR, lastXBinR, "", "");}
  if(plotESDTwoGammaCombinationZR == kTRUE){ histograms->AddHistogram("ESD_TwoGammaCombination_ZR" ,"" , nXBinsZR, firstXBinZR, lastXBinZR, nYBinsZR, firstYBinZR, lastYBinZR, "", "");}
  if(plotESDTwoGammaCombinationXY == kTRUE){ histograms->AddHistogram("ESD_TwoGammaCombination_XY" ,"" , nXBinsXY, firstXBinXY, lastXBinXY, nYBinsXY, firstYBinXY, lastYBinXY, "", "");}

  if(plotESDBackgroundOpeningAngleGamma == kTRUE){ histograms->AddHistogram("ESD_Background_GammaDaughter_OpeningAngle" ,"" , nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}
  if(plotESDBackgroundEnergy == kTRUE){ histograms->AddHistogram("ESD_Background_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotESDBackgroundPt == kTRUE){ histograms->AddHistogram("ESD_Background_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotESDBackgroundEta == kTRUE){ histograms->AddHistogram("ESD_Background_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotESDBackgroundPhi == kTRUE){ histograms->AddHistogram("ESD_Background_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotESDBackgroundMass == kTRUE){ histograms->AddHistogram("ESD_Background_Mass" ,"" , nXBinsEtaMass, firstXBinEtaMass, lastXBinEtaMass, "", "");}
  if(plotESDBackgroundR == kTRUE){ histograms->AddHistogram("ESD_Background_R" ,"" , nXBinsR, firstXBinR, lastXBinR, "", "");}
  if(plotESDBackgroundZR == kTRUE){ histograms->AddHistogram("ESD_Background_ZR" ,"" , nXBinsZR, firstXBinZR, lastXBinZR, nYBinsZR, firstYBinZR, lastYBinZR, "", "");}
  if(plotESDBackgroundXY == kTRUE){ histograms->AddHistogram("ESD_Background_XY" ,"" , nXBinsXY, firstXBinXY, lastXBinXY, nYBinsXY, firstYBinXY, lastYBinXY, "", "");}

  if(plotMapping == kTRUE){
    histograms->InitializeMappingValues(nPhiIndex,nRIndex,nXBinsMapping,minRadius,maxRadius,nYBinsMapping,minPhi,maxPhi);
    histograms->AddMappingHistograms(nPhiIndex,nRIndex,nXBinsMapping,minRadius,maxRadius,nYBinsMapping,minPhi,maxPhi);
  }

  if(plotResolutiondPt == kTRUE){histograms->AddHistogram("Resolution_dPt" ,"" , nXBinsResdPt, firstXBinResdPt, lastXBinResdPt, nYBinsResdPt, firstYBinResdPt, lastYBinResdPt, "", "");}
  if(plotResolutiondR == kTRUE){histograms->AddHistogram("Resolution_dR" ,"" , nXBinsResdR, firstXBinResdR, lastXBinResdR, nYBinsResdR, firstYBinResdR, lastYBinResdR, "", "");}
  if(plotResolutiondZ == kTRUE){histograms->AddHistogram("Resolution_dZ" ,"" , nXBinsResdZ, firstXBinResdZ, lastXBinResdZ, nYBinsResdZ, firstYBinResdZ, lastYBinResdZ, "", "");}
  
  if(plotResolutiondRdPt == kTRUE){histograms->AddHistogram("Resolution_dR_dPt" ,"" , nXBinsResdRdPt, firstXBinResdRdPt, lastXBinResdRdPt, nYBinsResdRdPt, firstYBinResdRdPt, lastYBinResdRdPt, "", "");}
  
  if(plotResolutionMCPt == kTRUE){histograms->AddHistogram("Resolution_MC_Pt" ,"" , nXBinsResPt, firstXBinResPt, lastXBinResPt,"","");}
  if(plotResolutionMCR == kTRUE){histograms->AddHistogram("Resolution_MC_R" ,"" , nXBinsResR, firstXBinResR, lastXBinResR,"","");}
  if(plotResolutionMCZ == kTRUE){histograms->AddHistogram("Resolution_MC_Z" ,"" , nXBinsResZ, firstXBinResZ, lastXBinResZ,"","");}
  
  if(plotResolutionESDPt == kTRUE){histograms->AddHistogram("Resolution_ESD_Pt" ,"" , nXBinsResPt, firstXBinResPt, lastXBinResPt,"","");}
  if(plotResolutionESDR == kTRUE){histograms->AddHistogram("Resolution_ESD_R" ,"" , nXBinsResR, firstXBinResR, lastXBinResR,"","");}
  if(plotResolutionESDZ == kTRUE){histograms->AddHistogram("Resolution_ESD_Z" ,"" , nXBinsResZ, firstXBinResZ, lastXBinResZ,"","");}
  
  if(plotNumberOfV0s == kTRUE){histograms->AddHistogram("NumberOfV0s","Number of v0s",100, 0, 100,"","");}
  if(plotNumberOfSurvivingV0s == kTRUE){histograms->AddHistogram("NumberOfSurvivingV0s","Number of surviving v0s",100, 0, 100,"","");}

  //  debug histograms
  if(plotV0MassDebugCut1 == kTRUE){histograms->AddHistogram("V0MassDebugCut1" ,"debug1" , nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass,"","");}
  if(plotV0MassDebugCut2 == kTRUE){histograms->AddHistogram("V0MassDebugCut2" ,"debug2" , nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass,"","");}
  if(plotV0MassDebugCut3 == kTRUE){histograms->AddHistogram("V0MassDebugCut3" ,"debug3" , nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass,"","");}
  if(plotV0MassDebugCut4 == kTRUE){histograms->AddHistogram("V0MassDebugCut4" ,"debug4" , nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass,"","");}
  if(plotV0MassDebugCut5 == kTRUE){histograms->AddHistogram("V0MassDebugCut5" ,"debug5" , nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass,"","");}
  if(plotV0MassDebugCut6 == kTRUE){histograms->AddHistogram("V0MassDebugCut6" ,"debug6" , nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass,"","");}
  if(plotV0MassDebugCut7 == kTRUE){histograms->AddHistogram("V0MassDebugCut7" ,"debug7" , nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass,"","");}
  if(plotV0MassDebugCut8 == kTRUE){histograms->AddHistogram("V0MassDebugCut8" ,"debug8" , nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass,"","");}


  if(plotPi0Spectra == kTRUE){histograms->AddHistogram("InvMass_vs_Pt_Spectra" ,"Invariant Mass vs Pt" , nXBinsSpectra, firstXBinSpectra, lastXBinSpectra,nYBinsSpectra, firstYBinSpectra, lastYBinSpectra,"InvMass [GeV]","Pt [GeV]");}

  if(plotPi0Spectra == kTRUE && calculateBackground == kTRUE){histograms->AddHistogram("Background_InvMass_vs_Pt_Spectra" ,"Background Invariant Mass vs Pt" , nXBinsSpectra, firstXBinSpectra, lastXBinSpectra,nYBinsSpectra, firstYBinSpectra, lastYBinSpectra,"InvMass [GeV]","Pt [GeV]");}

  

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
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();

  // Common Output Tree in common âdefaultâ output file
  AliAnalysisDataContainer *coutput1 = mgr->GetCommonOutputContainer();

  // Private output objects
  outputFileName.ReplaceAll(".root","");
  outputFileAppendix..ReplaceAll(".root","");
  TString fileOut = outputFileName + outputFileAppendix + ".root";

  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer("histogramsAliGammaConversion", TList::Class(),AliAnalysisManager::kOutputContainer, fileOut);



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
  v0Reader->SetChi2CutConversion(chi2CutConversion);
  v0Reader->SetChi2CutMeson(chi2CutMeson);
  v0Reader->SetPIDProbability(probElectron);
  v0Reader->SetXVertexCut(xVertexCut);
  v0Reader->SetYVertexCut(yVertexCut);
  v0Reader->SetZVertexCut(zVertexCut);
  v0Reader->SetSigmaMass(sigmaCutGammaMass);
  v0Reader->SetUseImprovedVertex(useImprovedVertex);
  v0Reader->SetDoMCTruth(doMCTruth);

  // Create the GammaConversionTask
  AliAnalysisTaskGammaConversion *gammaconversion = new AliAnalysisTaskGammaConversion("GammaConversionTask");
  gammaconversion->SetDebugLevel(10);

  gammaconversion->SetWriteNtuple(writeNtuple);

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

  if(dataList.IsNull()){
    cout<<"Data list is not set, aborting."<<endl;
    return;
  }

  TChain* chain= CreateESDChain(dataList,numberOfFilesToAnalyze);
  
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

  ////
  //Setting up ESD.par//
  ////
  cout<<"compiling ESD"<<endl;
  setupPar("ESD");
  gSystem->Load("libVMC.so");
  gSystem->Load("libESD.so");

  ////
  ////
  //Setting up STEERBase.par//
  ////
  cout<<"compiling STEERBase"<<endl;
  setupPar("STEERBase");
  gSystem->Load("libSTEERBase.so");

  ////
  //Setting up AOD.par//
  ////
  cout<<"compiling AOD"<<endl;
  setupPar("AOD");
  gSystem->Load("libAOD.so");
                                                                
  ////
  //Setting up ANALYSIS.par//
  ////
  cout<<"compiling ANALYSIS"<<endl;
  setupPar("ANALYSIS");
  gSystem->Load("libANALYSIS.so");

  ////
  //Setting up ANALYSISalice.par//
  ////
  cout<<"compiling ANALUSISalice"<<endl;
  setupPar("ANALYSISalice");
  gSystem->Load("libANALYSISalice.so");
                                                                                                                                  
  ////
  //Setting up PWG4Gamma.par//
  ////
  //  cout<<"compiling GammaConv"<<endl;

  if(usePWG4PartCorr == kTRUE){
    cout<<"Using PWG4PartCorr library"<<endl;
    setupPar("PWG4PartCorr");
    gSystem->Load("libPWG4PartCorr.so");  
  }
  else{
    setupPar("PWG4GammaConv");
    gSystem->Load("libPWG4GammaConv.so");
  }
  //if head:: use PWG4PartCorr
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
