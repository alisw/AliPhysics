/** VERSION NUMBER 0 */
/** new Version Kenneth */

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

Double_t LineCutZRSlope = 0.662487;
Double_t LineCutZValue = 7.;

Double_t maxRCut   = 160.;
Double_t etaCut    = 1.2;
Double_t ptCut     = 0.02;
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
Double_t probElectron = 0.002;

Double_t minOpeningAngleGhostCut = 0.01;

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
Bool_t plotMCConversionR             = kTRUE;
Bool_t plotMCConversionZR            = kTRUE;
Bool_t plotMCConversionXY            = kTRUE;
Bool_t plotMCConversionOpeningAngle  = kTRUE;

Bool_t plotMCEEnergy  = kTRUE;
Bool_t plotMCEPt      = kTRUE;
Bool_t plotMCEEta     = kTRUE;
Bool_t plotMCEPhi     = kTRUE;

Bool_t plotMCPEnergy  = kTRUE;
Bool_t plotMCPPt      = kTRUE;
Bool_t plotMCPEta     = kTRUE;
Bool_t plotMCPPhi     = kTRUE;

Bool_t plotMCallGammaEnergy = kTRUE;
Bool_t plotMCallGammaPt     = kTRUE;
Bool_t plotMCallGammaEta    = kTRUE;
Bool_t plotMCallGammaPhi    = kTRUE;
Bool_t plotMCallGammaRapid  = kTRUE;


Bool_t plotMCConvGammaEnergy  = kTRUE;
Bool_t plotMCConvGammaPt      = kTRUE;
Bool_t plotMCConvGammaEta     = kTRUE;
Bool_t plotMCConvGammaPhi     = kTRUE;
Bool_t plotMCConvGammaRapid   = kTRUE;
Bool_t plotMCConvGammaPtvsEta = kTRUE;

Bool_t plotMCallDirectGammaEnergy  = kTRUE;
Bool_t plotMCallDirectGammaPt      = kTRUE;
Bool_t plotMCallDirectGammaEta     = kTRUE;
Bool_t plotMCallDirectGammaPhi     = kTRUE;
Bool_t plotMCallDirectGammaRapid   = kTRUE;


Bool_t plotMCConvDirectGammaEnergy  = kTRUE;
Bool_t plotMCConvDirectGammaPt      = kTRUE;
Bool_t plotMCConvDirectGammaEta     = kTRUE;
Bool_t plotMCConvDirectGammaPhi     = kTRUE;
Bool_t plotMCConvDirectGammaRapid   = kTRUE;

Bool_t plotMCMotherEta					= kTRUE;
Bool_t plotMCMotherRapid                                = kTRUE;
Bool_t plotMCMotherPhi					= kTRUE;
Bool_t plotMCMotherPt					= kTRUE;
Bool_t plotMCMotherEnergy				= kTRUE;
Bool_t plotMCMotherMass					= kTRUE;
Bool_t plotMCMotherOpeningAngle				= kTRUE;
Bool_t plotMCMotherR					= kTRUE;
Bool_t plotMCMotherZR					= kTRUE;
Bool_t plotMCMotherXY	       				= kTRUE;
Bool_t plotMCMotherPtvsEtaWithinAcceptance              = kTRUE;
Bool_t plotMCMotherPtvsRapidWithinAcceptance            = kTRUE;
Bool_t plotMCMotherPtvsEtaConvGammaWithinAcceptance     = kTRUE;
Bool_t plotMCMotherPtvsRapidConvGammaWithinAcceptance   = kTRUE;
Bool_t plotMCMotherSpectra				= kTRUE;

Bool_t plotMCPi0Eta				   = kTRUE;
Bool_t plotMCPi0Rapid                              = kTRUE;
Bool_t plotMCPi0Phi                                = kTRUE;
Bool_t plotMCPi0Pt                                 = kTRUE;
Bool_t plotMCPi0Energy                             = kTRUE;
Bool_t plotMCPi0Mass                               = kTRUE;
Bool_t plotMCPi0OpeningAngle                       = kTRUE;
Bool_t plotMCPi0R                                  = kTRUE;
Bool_t plotMCPi0ZR                                 = kTRUE;
Bool_t plotMCPi0XY                                 = kTRUE;
Bool_t plotMCPi0PtvsEtaWithinAcceptance            = kTRUE;
Bool_t plotMCPi0PtvsRapidWithinAcceptance          = kTRUE;
Bool_t plotMCPi0PtvsEtaConvGammaWithinAcceptance   = kTRUE;
Bool_t plotMCPi0PtvsRapidConvGammaWithinAcceptance = kTRUE;


Bool_t plotMCPi0SecondaryEta                                = kTRUE;
Bool_t plotMCPi0SecondaryRapid                              = kTRUE;
Bool_t plotMCPi0SecondaryPhi                                = kTRUE;
Bool_t plotMCPi0SecondaryPt                                 = kTRUE;
Bool_t plotMCPi0SecondaryEnergy                             = kTRUE;
Bool_t plotMCPi0SecondaryMass                               = kTRUE;
Bool_t plotMCPi0SecondaryOpeningAngle                       = kTRUE;
Bool_t plotMCPi0SecondaryR                                  = kTRUE;
Bool_t plotMCPi0SecondaryZR                                 = kTRUE;
Bool_t plotMCPi0SecondaryXY                                 = kTRUE;
Bool_t plotMCPi0SecondaryPtvsEtaWithinAcceptance            = kTRUE;
Bool_t plotMCPi0SecondaryPtvsRapidWithinAcceptance          = kTRUE;
Bool_t plotMCPi0SecondaryPtvsEtaConvGammaWithinAcceptance   = kTRUE;
Bool_t plotMCPi0SecondaryPtvsRapidConvGammaWithinAcceptance = kTRUE;


Bool_t plotMCEtaEta                                = kTRUE;
Bool_t plotMCEtaRapid                              = kTRUE;
Bool_t plotMCEtaPhi                                = kTRUE;
Bool_t plotMCEtaPt                                 = kTRUE;
Bool_t plotMCEtaEnergy                             = kTRUE;
Bool_t plotMCEtaMass                               = kTRUE;
Bool_t plotMCEtaOpeningAngleGamma                  = kTRUE;
Bool_t plotMCEtaR                                  = kTRUE;
Bool_t plotMCEtaZR                                 = kTRUE;
Bool_t plotMCEtaXY                                 = kTRUE;
Bool_t plotMCEtaPtvsEtaWithinAcceptance		   = kTRUE;
Bool_t plotMCEtaPtvsRapidWithinAcceptance	   = kTRUE;
Bool_t plotMCEtaPtvsEtaConvGammaWithinAcceptance   = kTRUE;
Bool_t plotMCEtaPtvsRapidConvGammaWithinAcceptance = kTRUE;


// Histograms from esd tracks
Bool_t plotESDConversionR            = kTRUE;
Bool_t plotESDConversionZR           = kTRUE;
Bool_t plotESDConversionXY           = kTRUE;
Bool_t plotESDConversionOpeningAngle = kTRUE;

Bool_t plotESDEEnergy = kTRUE;
Bool_t plotESDEPt     = kTRUE;
Bool_t plotESDEEta    = kTRUE;
Bool_t plotESDEPhi    = kTRUE;

Bool_t plotESDPEnergy = kTRUE;
Bool_t plotESDPPt     = kTRUE;
Bool_t plotESDPEta    = kTRUE;
Bool_t plotESDPPhi    = kTRUE;

Bool_t plotESDConvGammaEnergy = kTRUE;
Bool_t plotESDConvGammaPt     = kTRUE;
Bool_t plotESDConvGammaEta    = kTRUE;
Bool_t plotESDConvGammaPhi    = kTRUE;
Bool_t plotESDConvGammaMass   = kTRUE;
Bool_t plotESDConvGammaWidth  = kTRUE;
Bool_t plotESDConvGammaChi2   = kTRUE;
Bool_t plotESDConvGammaNDF    = kTRUE;
Bool_t plotESDConvGammaRapid  = kTRUE;
Bool_t plotESDConvGammaPtvsEta = kTRUE;

Bool_t plotESDTrueConvGammaEnergy         = kTRUE;
Bool_t plotESDTrueConvGammaPt             = kTRUE;
Bool_t plotESDTrueConvGammaEta            = kTRUE;
Bool_t plotESDTrueConvGammaPhi            = kTRUE;
Bool_t plotESDTrueConvGammaMass           = kTRUE;
Bool_t plotESDTrueConvGammaWidth          = kTRUE;
Bool_t plotESDTrueConvGammaChi2           = kTRUE;
Bool_t plotESDTrueConvGammaNDF            = kTRUE;
Bool_t plotESDTrueConvGammaRapid          = kTRUE;
Bool_t plotESDTrueConvGammaPtvsEta        = kTRUE;
Bool_t plotESDTrueConversionR             = kTRUE;
Bool_t plotESDTrueConversionZR            = kTRUE;
Bool_t plotESDTrueConversionXY            = kTRUE;
Bool_t plotESDTrueConversionOpeningAngle  = kTRUE;

Bool_t plotESDNoCutConvGammaEnergy         = kTRUE;
Bool_t plotESDNoCutConvGammaPt             = kTRUE;
Bool_t plotESDNoCutConvGammaEta            = kTRUE;
Bool_t plotESDNoCutConvGammaPhi            = kTRUE;
Bool_t plotESDNoCutConvGammaMass           = kTRUE;
Bool_t plotESDNoCutConvGammaWidth          = kTRUE;
Bool_t plotESDNoCutConvGammaChi2           = kTRUE;
Bool_t plotESDNoCutConvGammaNDF            = kTRUE;
Bool_t plotESDNoCutConvGammaRapid          = kTRUE;
Bool_t plotESDNoCutConvGammaPtvsEta        = kTRUE;
Bool_t plotESDNoCutConversionR             = kTRUE;
Bool_t plotESDNoCutConversionZR            = kTRUE;
Bool_t plotESDNoCutConversionXY            = kTRUE;
Bool_t plotESDNoCutConversionOpeningAngle  = kTRUE;

Bool_t plotESDMotherOpeningAngleGamma = kTRUE;
Bool_t plotESDMotherEnergy            = kTRUE;
Bool_t plotESDMotherPt                = kTRUE;
Bool_t plotESDMotherEta               = kTRUE;
Bool_t plotESDMotherPhi               = kTRUE;
Bool_t plotESDMotherMass              = kTRUE;
Bool_t plotESDMotherR                 = kTRUE;
Bool_t plotESDMotherZR                = kTRUE;
Bool_t plotESDMotherXY                = kTRUE;
Bool_t plotESDMotherRapid             = kTRUE;


Bool_t plotESDBackgroundOpeningAngleGamma = kTRUE;
Bool_t plotESDBackgroundEnergy            = kTRUE;
Bool_t plotESDBackgroundPt                = kTRUE;
Bool_t plotESDBackgroundEta               = kTRUE;
Bool_t plotESDBackgroundPhi               = kTRUE;
Bool_t plotESDBackgroundMass              = kTRUE;
Bool_t plotESDBackgroundR                 = kTRUE;
Bool_t plotESDBackgroundZR                = kTRUE;
Bool_t plotESDBackgroundXY                = kTRUE;
Bool_t plotESDBackgroundRapid             = kTRUE;



Bool_t plotMapping = kFALSE;       

Bool_t plotResolutiondPt = kTRUE;
Bool_t plotResolutiondR  = kTRUE;
Bool_t plotResolutiondZ  = kTRUE;

Bool_t plotResolutiondRdPt = kTRUE;

Bool_t plotResolutionMCPt = kTRUE;
Bool_t plotResolutionMCR  = kTRUE;
Bool_t plotResolutionMCZ  = kTRUE;

Bool_t plotResolutionESDPt = kTRUE;
Bool_t plotResolutionESDR  = kTRUE;
Bool_t plotResolutionESDZ  = kTRUE;

Bool_t plotESDNumberOfV0s          = kTRUE;
Bool_t plotESDNumberOfSurvivingV0s = kTRUE;

//  debug histograms
Bool_t plotESDCutGetOnFly      = kTRUE;
Bool_t plotESDCutNContributors = kTRUE;
Bool_t plotESDCutLikeSign      = kTRUE;
Bool_t plotESDCutRefit         = kTRUE;
Bool_t plotESDCutKink          = kTRUE;
Bool_t plotESDCutPIDProb       = kTRUE;
Bool_t plotESDCutR             = kTRUE;
Bool_t plotESDCutLine          = kTRUE;
Bool_t plotESDCutNDF           = kTRUE;
Bool_t plotESDCutChi2          = kTRUE;
Bool_t plotESDCutEta           = kTRUE;
Bool_t plotESDCutPt            = kTRUE;
Bool_t plotESDTrueConvGammaTrackLength =kTRUE;
Bool_t plotESDTrueConvGammaTrackLengthVSInvMass =kTRUE;


Bool_t plotPi0Spectra = kTRUE;
Bool_t plotEtaSpectra = kTRUE;


/** ----------------- end define which histograms to plot here -------------------------------*/



/** ----------- Define the binning for the different plot types here -------------------------*/
//R-plots
Int_t nXBinsR = 500;
Double_t firstXBinR = 0.;
Double_t lastXBinR = 250.;

//ZR-plots
Int_t nXBinsZR = 1200;
Double_t firstXBinZR = -300.;
Double_t lastXBinZR = 300.;
Int_t nYBinsZR = 500;
Double_t firstYBinZR = 0.;
Double_t lastYBinZR = 250.;

//XY-plots
Int_t nXBinsXY = 1000;
Double_t firstXBinXY = -250.;
Double_t lastXBinXY = 250.;
Int_t nYBinsXY = 1000;
Double_t firstYBinXY = -250.;
Double_t lastYBinXY = 250.;

//OpenAngle-plots
Int_t nXBinsOpeningAngle = 400;
Double_t firstXBinOpeningAngle = 0.;
Double_t lastXBinOpeningAngle = TMath::Pi();

//Energy-plots
Int_t nXBinsEnergy = 200;
Double_t firstXBinEnergy = 0.;
Double_t lastXBinEnergy = 50.;

//Pt-plots
Int_t nXBinsPt = 200;
Double_t firstXBinPt = 0.;
Double_t lastXBinPt = 50.;

//Eta-plots
Int_t nXBinsEta = 40;
Double_t firstXBinEta = -2.;
Double_t lastXBinEta = 2.;

//Rapidity
Int_t nXBinsRapid = 200;
Double_t firstXBinRapid = -10.;
Double_t lastXBinRapid = 10.;

//Phi-plots
Int_t nXBinsPhi = 72;
Double_t firstXBinPhi = -TMath::Pi();
Double_t lastXBinPhi = TMath::Pi();

//Mapping-plots
Int_t nXBinsMapping = 400;
Double_t firstXBinMapping = -100.;
Double_t lastXBinMapping = 100.;
Int_t nYBinsMapping = 40;
Double_t firstYBinMapping = -2;
Double_t lastYBinMapping = 2;

//ResolutionPlots
//RESdPt
Int_t nXBinsResdPt=500;
Int_t firstXBinResdPt= 0;
Int_t lastXBinResdPt=5;
Int_t nYBinsResdPt=1000;
Int_t firstYBinResdPt= -5;
Int_t lastYBinResdPt=5;

//RESdR
Int_t nXBinsResdR=500;
Int_t firstXBinResdR= 0;
Int_t lastXBinResdR=250;
Int_t nYBinsResdR=100;
Int_t firstYBinResdR= -25;
Int_t lastYBinResdR=25;

//RESdZ
Int_t nXBinsResdZ=80;
Int_t firstXBinResdZ= -20;
Int_t lastXBinResdZ=20;
Int_t nYBinsResdZ=80;
Int_t firstYBinResdZ= -20;
Int_t lastYBinResdZ=20;

//RESdRdPt
Int_t nXBinsResdRdPt=440;
Int_t firstXBinResdRdPt= -22;
Int_t lastXBinResdRdPt=22;
Int_t nYBinsResdRdPt=100;
Int_t firstYBinResdRdPt= -5;
Int_t lastYBinResdRdPt=5;

//RESMCPt
Int_t nXBinsResPt=100;
Int_t firstXBinResPt= 0;
Int_t lastXBinResPt=5;

//RESMCR
Int_t nXBinsResR=500;
Int_t firstXBinResR= 0;
Int_t lastXBinResR=250;

//RESMCZ
Int_t nXBinsResZ=500;
Int_t firstXBinResZ= 0;
Int_t lastXBinResZ=250;

//GammaMass-plots
Int_t nXBinsGammaMass = 4000;
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
Int_t nXBinsSpectra = 500;
Double_t firstXBinSpectra = 0.;
Double_t lastXBinSpectra = 1.;
Int_t nYBinsSpectra = 100;
Double_t firstYBinSpectra = 0.;
Double_t lastYBinSpectra = 50.;

//track length plots
Int_t nXBinsTrackLength = 1000;
Double_t firstXBinTrackLength = 0;
Double_t lastXBinTrackLength = 500;

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
	
  //  cout<<"All arguments: "<<arguments<<endl;
	
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
	if(line.compare("") != 0){//checks if there is an empty line in the data list
	  numberOfFilesToAnalyze++;
	}
      }
  }
  cout<<"Number Of files to analyze: "<<numberOfFilesToAnalyze<<endl;
	
  build();//build (if necessary) and load the libraries needed
	
  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C"); // load the CreateChain macro
	
  AliLog::SetGlobalLogLevel(AliLog::kError);
	
  //-------------------------------- Creating the histograms -------------------------------
  AliGammaConversionHistograms * histograms = new AliGammaConversionHistograms();
	
  if(plotMCConversionR == kTRUE){ histograms->AddHistogram("MC_Conversion_R","Radius of gamma conversion points",nXBinsR, firstXBinR, lastXBinR,"counts","cm");}
  if(plotMCConversionZR == kTRUE){ histograms->AddHistogram("MC_Conversion_ZR","Radius of gamma conversion points vs Z",nXBinsZR, firstXBinZR, lastXBinZR, nYBinsZR, firstYBinZR, lastYBinZR, "cm", "cm");}
  if(plotMCConversionXY == kTRUE){ histograms->AddHistogram("MC_Conversion_XY","Gamma XY converison point.",nXBinsXY, firstXBinXY, lastXBinXY, nYBinsXY, firstYBinXY, lastYBinXY, "cm", "cm");}
  if(plotMCConversionOpeningAngle == kTRUE){ histograms->AddHistogram("MC_Conversion_OpeningAngle","Opening angle of e+e- pairs from gamma conversion",nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "counts", "cm");}
	
  if(plotMCEEnergy == kTRUE){ histograms->AddHistogram("MC_E_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotMCEPt == kTRUE){ histograms->AddHistogram("MC_E_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotMCEEta == kTRUE){ histograms->AddHistogram("MC_E_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMCEPhi == kTRUE){ histograms->AddHistogram("MC_E_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
	
  if(plotMCPEnergy == kTRUE){ histograms->AddHistogram("MC_P_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotMCPPt == kTRUE){ histograms->AddHistogram("MC_P_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotMCPEta == kTRUE){ histograms->AddHistogram("MC_P_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMCPPhi == kTRUE){ histograms->AddHistogram("MC_P_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
	
  if(plotMCallGammaEnergy == kTRUE){ histograms->AddHistogram("MC_allGamma_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotMCallGammaPt == kTRUE){ histograms->AddHistogram("MC_allGamma_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotMCallGammaEta == kTRUE){ histograms->AddHistogram("MC_allGamma_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMCallGammaPhi == kTRUE){ histograms->AddHistogram("MC_allGamma_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotMCallGammaRapid == kTRUE){ histograms->AddHistogram("MC_allGamma_Rapid" ,"" , nXBinsRapid, firstXBinRapid, lastXBinRapid, "", "");}
	
  if(plotMCConvGammaEnergy == kTRUE){ histograms->AddHistogram("MC_ConvGamma_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotMCConvGammaPt == kTRUE){ histograms->AddHistogram("MC_ConvGamma_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotMCConvGammaEta == kTRUE){ histograms->AddHistogram("MC_ConvGamma_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMCConvGammaPhi == kTRUE){ histograms->AddHistogram("MC_ConvGamma_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotMCConvGammaRapid == kTRUE){ histograms->AddHistogram("MC_ConvGamma_Rapid" ,"" , nXBinsRapid, firstXBinRapid, lastXBinRapid, "", "");}
  if(plotMCConvGammaPtvsEta == kTRUE){ histograms->AddHistogram("MC_ConvGamma_Pt_Eta","", nXBinsPt, firstXBinPt, lastXBinPt,nXBinsEta, firstXBinEta, lastXBinEta,"","");}
	
  if(plotMCallDirectGammaEnergy == kTRUE){ histograms->AddHistogram("MC_allDirectGamma_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotMCallDirectGammaPt == kTRUE){ histograms->AddHistogram("MC_allDirectGamma_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotMCallDirectGammaEta == kTRUE){ histograms->AddHistogram("MC_allDirectGamma_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMCallDirectGammaPhi == kTRUE){ histograms->AddHistogram("MC_allDirectGamma_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotMCallDirectGammaRapid == kTRUE){ histograms->AddHistogram("MC_allDirectGamma_Rapid" ,"" , nXBinsRapid, firstXBinRapid, lastXBinRapid, "", "");}
	
  if(plotMCConvDirectGammaEnergy == kTRUE){ histograms->AddHistogram("MC_ConvDirectGamma_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotMCConvDirectGammaPt == kTRUE){ histograms->AddHistogram("MC_ConvDirectGamma_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotMCConvDirectGammaEta == kTRUE){ histograms->AddHistogram("MC_ConvDirectGamma_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMCConvDirectGammaPhi == kTRUE){ histograms->AddHistogram("MC_ConvDirectGamma_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotMCConvDirectGammaRapid == kTRUE){ histograms->AddHistogram("MC_ConvDirectGamma_Rapid" ,"" , nXBinsRapid, firstXBinRapid, lastXBinRapid, "", "");}
	
  if(plotMCMotherEta == kTRUE){ histograms->AddHistogram("MC_Mother_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMCMotherPhi == kTRUE){ histograms->AddHistogram("MC_Mother_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotMCMotherRapid == kTRUE){ histograms->AddHistogram("MC_Mother_Rapid" ,"" , nXBinsRapid, firstXBinRapid, lastXBinRapid, "", "");}
  if(plotMCMotherPt == kTRUE){ histograms->AddHistogram("MC_Mother_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotMCMotherEnergy == kTRUE){ histograms->AddHistogram("MC_Mother_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotMCMotherMass == kTRUE){ histograms->AddHistogram("MC_Mother_Mass" ,"" , nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass, "", "");}
  if(plotMCMotherOpeningAngle == kTRUE){ histograms->AddHistogram("MC_Mother_GammaDaughter_OpeningAngle" ,"" , nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}
  if(plotMCMotherR == kTRUE){ histograms->AddHistogram("MC_Mother_R" ,"" , nXBinsR, firstXBinR, lastXBinR, "", "");}
  if(plotMCMotherZR == kTRUE){ histograms->AddHistogram("MC_Mother_ZR" ,"" , nXBinsZR, firstXBinZR, lastXBinZR, nYBinsZR, firstYBinZR, lastYBinZR, "", "");}
  if(plotMCMotherXY == kTRUE){ histograms->AddHistogram("MC_Mother_XY" ,"" , nXBinsXY, firstXBinXY, lastXBinXY, nYBinsXY, firstYBinXY, lastYBinXY, "", "");}
  if(plotMCMotherPtvsEtaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Mother_Pt_Eta_withinAcceptance" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMCMotherPtvsRapidWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Mother_Pt_Rapid_withinAcceptance" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, nXBinsRapid, firstXBinRapid, lastXBinRapid, "", "");}
  if(plotMCMotherPtvsEtaConvGammaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Mother_Pt_Eta_ConvGamma_withinAcceptance" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMCMotherPtvsRapidConvGammaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Mother_Pt_Rapid_ConvGamma_withinAcceptance" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, nXBinsRapid, firstXBinRapid, lastXBinRapid, "", "");}

  if(plotMCMotherSpectra == kTRUE){ 
    histograms->AddHistogram("MC_Mother_InvMass_vs_Pt" ,"" ,nXBinsSpectra, firstXBinSpectra, lastXBinSpectra, nYBinsSpectra, firstYBinSpectra, lastYBinSpectra, "", "");
    histograms->AddHistogram("MC_Mother_InvMass_vs_Pt_withinAcceptance" ,"" ,nXBinsSpectra, firstXBinSpectra, lastXBinSpectra, nYBinsSpectra, firstYBinSpectra, lastYBinSpectra, "", "");
    histograms->AddHistogram("MC_Mother_InvMass_vs_Pt_ConvGamma_withinAcceptance" ,"" ,nXBinsSpectra, firstXBinSpectra, lastXBinSpectra, nYBinsSpectra, firstYBinSpectra, lastYBinSpectra, "", "");
  }
	
	
  if(plotMCPi0Eta == kTRUE){ histograms->AddHistogram("MC_Pi0_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}	
  if(plotMCPi0Rapid == kTRUE){ histograms->AddHistogram("MC_Pi0_Rapid" ,"" , nXBinsRapid, firstXBinRapid, lastXBinRapid, "", "");}	
  if(plotMCPi0Phi == kTRUE){ histograms->AddHistogram("MC_Pi0_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotMCPi0Pt == kTRUE){ histograms->AddHistogram("MC_Pi0_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotMCPi0Energy == kTRUE){ histograms->AddHistogram("MC_Pi0_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotMCPi0Mass == kTRUE){ histograms->AddHistogram("MC_Pi0_Mass" ,"" , nXBinsPi0Mass, firstXBinPi0Mass, lastXBinPi0Mass, "", "");}
  if(plotMCPi0OpeningAngle == kTRUE){ histograms->AddHistogram("MC_Pi0_GammaDaughter_OpeningAngle" ,"" , nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}
  if(plotMCPi0R == kTRUE){ histograms->AddHistogram("MC_Pi0_R" ,"" , nXBinsR, firstXBinR, lastXBinR, "", "");}
  if(plotMCPi0ZR == kTRUE){ histograms->AddHistogram("MC_Pi0_ZR" ,"" , nXBinsZR, firstXBinZR, lastXBinZR, nYBinsZR, firstYBinZR, lastYBinZR, "", "");}
  if(plotMCPi0XY == kTRUE){ histograms->AddHistogram("MC_Pi0_XY" ,"" , nXBinsXY, firstXBinXY, lastXBinXY, nYBinsXY, firstYBinXY, lastYBinXY, "", "");}
  if(plotMCPi0PtvsEtaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Pi0_Pt_Eta_withinAcceptance" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMCPi0PtvsRapidWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Pi0_Pt_Rapid_withinAcceptance" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, nXBinsRapid, firstXBinRapid, lastXBinRapid, "", "");}
  if(plotMCPi0PtvsEtaConvGammaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Pi0_Pt_Eta_ConvGamma_withinAcceptance" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMCPi0PtvsRapidConvGammaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Pi0_Pt_Rapid_ConvGamma_withinAcceptance" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, nXBinsRapid, firstXBinRapid, lastXBinRapid, "", "");}
	
	
  if(plotMCPi0SecondaryEta == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMCPi0SecondaryRapid == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_Rapid" ,"" , nXBinsRapid, firstXBinRapid, lastXBinRapid, "", "");}
  if(plotMCPi0SecondaryPhi == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotMCPi0SecondaryPt == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotMCPi0SecondaryEnergy == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotMCPi0SecondaryMass == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_Mass" ,"" , nXBinsPi0Mass, firstXBinPi0Mass, lastXBinPi0Mass, "", "");}
  if(plotMCPi0SecondaryOpeningAngle == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_GammaDaughter_OpeningAngle" ,"" , nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}
  if(plotMCPi0SecondaryR == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_R" ,"" , nXBinsR, firstXBinR, lastXBinR, "", "");}
  if(plotMCPi0SecondaryZR == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_ZR" ,"" , nXBinsZR, firstXBinZR, lastXBinZR, nYBinsZR, firstYBinZR, lastYBinZR, "", "");}
  if(plotMCPi0SecondaryXY == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_XY" ,"" , nXBinsXY, firstXBinXY, lastXBinXY, nYBinsXY, firstYBinXY, lastYBinXY, "", "");}
  if(plotMCPi0SecondaryPtvsEtaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_Pt_Eta_withinAcceptance" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMCPi0SecondaryPtvsRapidWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_Pt_Rapid_withinAcceptance" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, nXBinsRapid, firstXBinRapid, lastXBinRapid, "", "");}
  if(plotMCPi0SecondaryPtvsEtaConvGammaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_Pt_Eta_ConvGamma_withinAcceptance" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMCPi0SecondaryPtvsRapidConvGammaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Pi0_Secondaries_Pt_Rapid_ConvGamma_withinAcceptance" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, nXBinsRapid, firstXBinRapid, lastXBinRapid, "", "");}
	
	
	
  if(plotMCEtaEta == kTRUE){ histograms->AddHistogram("MC_Eta_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMCEtaRapid == kTRUE){ histograms->AddHistogram("MC_Eta_Rapid" ,"" , nXBinsRapid, firstXBinRapid, lastXBinRapid, "", "");}
  if(plotMCEtaPhi == kTRUE){ histograms->AddHistogram("MC_Eta_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotMCEtaPt == kTRUE){ histograms->AddHistogram("MC_Eta_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotMCEtaEnergy == kTRUE){ histograms->AddHistogram("MC_Eta_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotMCEtaMass == kTRUE){ histograms->AddHistogram("MC_Eta_Mass" ,"" , nXBinsEtaMass, firstXBinEtaMass, lastXBinEtaMass, "", "");}
  if(plotMCEtaOpeningAngleGamma == kTRUE){ histograms->AddHistogram("MC_Eta_GammaDaughter_OpeningAngle" ,"" , nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}
  if(plotMCEtaR == kTRUE){ histograms->AddHistogram("MC_Eta_R" ,"" , nXBinsR, firstXBinR, lastXBinR, "", "");}
  if(plotMCEtaZR == kTRUE){ histograms->AddHistogram("MC_Eta_ZR" ,"" , nXBinsZR, firstXBinZR, lastXBinZR, nYBinsZR, firstYBinZR, lastYBinZR, "", "");}
  if(plotMCEtaXY == kTRUE){ histograms->AddHistogram("MC_Eta_XY" ,"" , nXBinsXY, firstXBinXY, lastXBinXY, nYBinsXY, firstYBinXY, lastYBinXY, "", "");}
  if(plotMCEtaPtvsEtaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Eta_Pt_Eta_withinAcceptance" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMCEtaPtvsRapidWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Eta_Pt_Rapid_withinAcceptance" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, nXBinsRapid, firstXBinRapid, lastXBinRapid, "", "");}
  if(plotMCEtaPtvsEtaConvGammaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Eta_Pt_Eta_ConvGamma_withinAcceptance" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotMCEtaPtvsRapidConvGammaWithinAcceptance == kTRUE){ histograms->AddHistogram("MC_Eta_Pt_Rapid_ConvGamma_withinAcceptance" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, nXBinsRapid, firstXBinRapid, lastXBinRapid, "", "");}
	
	
  // Histograms from esd tracks	
  if(plotESDEEnergy == kTRUE){ histograms->AddHistogram("ESD_E_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotESDEPt == kTRUE){ histograms->AddHistogram("ESD_E_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotESDEEta == kTRUE){ histograms->AddHistogram("ESD_E_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotESDEPhi == kTRUE){ histograms->AddHistogram("ESD_E_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
	
  if(plotESDPEnergy == kTRUE){ histograms->AddHistogram("ESD_P_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotESDPPt == kTRUE){ histograms->AddHistogram("ESD_P_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotESDPEta == kTRUE){ histograms->AddHistogram("ESD_P_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotESDPPhi == kTRUE){ histograms->AddHistogram("ESD_P_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
	
  if(plotESDConvGammaEnergy == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotESDConvGammaPt == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotESDConvGammaEta == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotESDConvGammaPhi == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotESDConvGammaMass == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_Mass" ,"" ,  nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass, "", "");}
  if(plotESDConvGammaWidth == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_Width" ,"" , nXBinsGammaWidth, firstXBinGammaWidth, lastXBinGammaWidth, "", "");}
  if(plotESDConvGammaChi2 == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_Chi2" ,"" , nXBinsGammaChi2, firstXBinGammaChi2, lastXBinGammaChi2, "", "");}
  if(plotESDConvGammaNDF == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_NDF" ,"" , nXBinsGammaNDF, firstXBinGammaNDF, lastXBinGammaNDF, "", "");}
  if(plotESDConvGammaRapid == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_Rapid" ,"" , nXBinsRapid, firstXBinRapid, lastXBinRapid, "", "");}
  if(plotESDConvGammaPtvsEta == kTRUE){ histograms->AddHistogram("ESD_ConvGamma_Pt_Eta","", nXBinsPt, firstXBinPt, lastXBinPt,nXBinsEta, firstXBinEta, lastXBinEta,"","" );}

  if(plotESDConversionR == kTRUE){ histograms->AddHistogram("ESD_Conversion_R" ,"" , nXBinsR, firstXBinR, lastXBinR, "", "");}
  if(plotESDConversionZR == kTRUE){ histograms->AddHistogram("ESD_Conversion_ZR" ,"" , nXBinsZR, firstXBinZR, lastXBinZR, nYBinsZR, firstYBinZR, lastYBinZR, "", "");}
  if(plotESDConversionXY == kTRUE){ histograms->AddHistogram("ESD_Conversion_XY" ,"" , nXBinsXY, firstXBinXY, lastXBinXY, nYBinsXY, firstYBinXY, lastYBinXY, "", "");}
  if(plotESDConversionOpeningAngle == kTRUE){ histograms->AddHistogram("ESD_Conversion_OpeningAngle" ,"" , nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}


  if(plotESDTrueConvGammaEnergy == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotESDTrueConvGammaPt == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotESDTrueConvGammaEta == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotESDTrueConvGammaPhi == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotESDTrueConvGammaMass == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_Mass" ,"" ,  nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass, "", "");}
  if(plotESDTrueConvGammaWidth == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_Width" ,"" , nXBinsGammaWidth, firstXBinGammaWidth, lastXBinGammaWidth, "", "");}
  if(plotESDTrueConvGammaChi2 == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_Chi2" ,"" , nXBinsGammaChi2, firstXBinGammaChi2, lastXBinGammaChi2, "", "");}
  if(plotESDTrueConvGammaNDF == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_NDF" ,"" , nXBinsGammaNDF, firstXBinGammaNDF, lastXBinGammaNDF, "", "");}
  if(plotESDTrueConvGammaRapid == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_Rapid" ,"" , nXBinsRapid, firstXBinRapid, lastXBinRapid, "", "");}
  if(plotESDTrueConvGammaPtvsEta == kTRUE){ histograms->AddHistogram("ESD_TrueConvGamma_Pt_Eta" ,"" , nXBinsPt, firstXBinPt, lastXBinPt,nXBinsEta, firstXBinEta, lastXBinEta, "", "");}

  if(plotESDTrueConversionR == kTRUE){ histograms->AddHistogram("ESD_TrueConversion_R" ,"" , nXBinsR, firstXBinR, lastXBinR, "", "");}
  if(plotESDTrueConversionZR == kTRUE){ histograms->AddHistogram("ESD_TrueConversion_ZR" ,"" , nXBinsZR, firstXBinZR, lastXBinZR, nYBinsZR, firstYBinZR, lastYBinZR, "", "");}
  if(plotESDTrueConversionXY == kTRUE){ histograms->AddHistogram("ESD_TrueConversion_XY" ,"" , nXBinsXY, firstXBinXY, lastXBinXY, nYBinsXY, firstYBinXY, lastYBinXY, "", "");}
  if(plotESDTrueConversionOpeningAngle == kTRUE){ histograms->AddHistogram("ESD_TrueConversion_OpeningAngle" ,"" , nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}



  if(plotESDNoCutConvGammaEnergy == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotESDNoCutConvGammaPt == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotESDNoCutConvGammaEta == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotESDNoCutConvGammaPhi == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotESDNoCutConvGammaMass == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_Mass" ,"" ,  nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass, "", "");}
  if(plotESDNoCutConvGammaWidth == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_Width" ,"" , nXBinsGammaWidth, firstXBinGammaWidth, lastXBinGammaWidth, "", "");}
  if(plotESDNoCutConvGammaChi2 == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_Chi2" ,"" , nXBinsGammaChi2, firstXBinGammaChi2, lastXBinGammaChi2, "", "");}
  if(plotESDNoCutConvGammaNDF == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_NDF" ,"" , nXBinsGammaNDF, firstXBinGammaNDF, lastXBinGammaNDF, "", "");}
  if(plotESDNoCutConvGammaRapid == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_Rapid" ,"" , nXBinsRapid, firstXBinRapid, lastXBinRapid, "", "");}
  if(plotESDNoCutConvGammaPtvsEta == kTRUE){ histograms->AddHistogram("ESD_NoCutConvGamma_Pt_Eta" ,"" , nXBinsPt, firstXBinPt, lastXBinPt,nXBinsEta, firstXBinEta, lastXBinEta, "", "");}

  if(plotESDNoCutConversionR == kTRUE){ histograms->AddHistogram("ESD_NoCutConversion_R" ,"" , nXBinsR, firstXBinR, lastXBinR, "", "");}
  if(plotESDNoCutConversionZR == kTRUE){ histograms->AddHistogram("ESD_NoCutConversion_ZR" ,"" , nXBinsZR, firstXBinZR, lastXBinZR, nYBinsZR, firstYBinZR, lastYBinZR, "", "");}
  if(plotESDNoCutConversionXY == kTRUE){ histograms->AddHistogram("ESD_NoCutConversion_XY" ,"" , nXBinsXY, firstXBinXY, lastXBinXY, nYBinsXY, firstYBinXY, lastYBinXY, "", "");}
  if(plotESDNoCutConversionOpeningAngle == kTRUE){ histograms->AddHistogram("ESD_NoCutConversion_OpeningAngle" ,"" , nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}




  if(plotESDMotherOpeningAngleGamma == kTRUE){ histograms->AddHistogram("ESD_Mother_GammaDaughter_OpeningAngle" ,"" , nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}
  if(plotESDMotherEnergy == kTRUE){ histograms->AddHistogram("ESD_Mother_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotESDMotherPt == kTRUE){ histograms->AddHistogram("ESD_Mother_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotESDMotherEta == kTRUE){ histograms->AddHistogram("ESD_Mother_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotESDMotherPhi == kTRUE){ histograms->AddHistogram("ESD_Mother_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotESDMotherMass == kTRUE){ histograms->AddHistogram("ESD_Mother_Mass" ,"" , nXBinsPi0Mass, firstXBinPi0Mass, lastXBinPi0Mass, "", "");}
  if(plotESDMotherR == kTRUE){ histograms->AddHistogram("ESD_Mother_R" ,"" , nXBinsR, firstXBinR, lastXBinR, "", "");}
  if(plotESDMotherZR == kTRUE){ histograms->AddHistogram("ESD_Mother_ZR" ,"" , nXBinsZR, firstXBinZR, lastXBinZR, nYBinsZR, firstYBinZR, lastYBinZR, "", "");}
  if(plotESDMotherXY == kTRUE){ histograms->AddHistogram("ESD_Mother_XY" ,"" , nXBinsXY, firstXBinXY, lastXBinXY, nYBinsXY, firstYBinXY, lastYBinXY, "", "");}
  if(plotESDMotherRapid == kTRUE){ histograms->AddHistogram("ESD_Mother_Rapid" ,"" , nXBinsRapid, firstXBinRapid, lastXBinRapid, "", "");}

	
  if(plotESDBackgroundOpeningAngleGamma == kTRUE){ histograms->AddHistogram("ESD_Background_GammaDaughter_OpeningAngle" ,"" , nXBinsOpeningAngle, firstXBinOpeningAngle, lastXBinOpeningAngle, "", "");}
  if(plotESDBackgroundEnergy == kTRUE){ histograms->AddHistogram("ESD_Background_Energy" ,"" , nXBinsEnergy, firstXBinEnergy, lastXBinEnergy, "", "");}
  if(plotESDBackgroundPt == kTRUE){ histograms->AddHistogram("ESD_Background_Pt" ,"" , nXBinsPt, firstXBinPt, lastXBinPt, "", "");}
  if(plotESDBackgroundEta == kTRUE){ histograms->AddHistogram("ESD_Background_Eta" ,"" , nXBinsEta, firstXBinEta, lastXBinEta, "", "");}
  if(plotESDBackgroundPhi == kTRUE){ histograms->AddHistogram("ESD_Background_Phi" ,"" , nXBinsPhi, firstXBinPhi, lastXBinPhi, "", "");}
  if(plotESDBackgroundMass == kTRUE){ histograms->AddHistogram("ESD_Background_Mass" ,"" , nXBinsEtaMass, firstXBinEtaMass, lastXBinEtaMass, "", "");}
  if(plotESDBackgroundR == kTRUE){ histograms->AddHistogram("ESD_Background_R" ,"" , nXBinsR, firstXBinR, lastXBinR, "", "");}
  if(plotESDBackgroundZR == kTRUE){ histograms->AddHistogram("ESD_Background_ZR" ,"" , nXBinsZR, firstXBinZR, lastXBinZR, nYBinsZR, firstYBinZR, lastYBinZR, "", "");}
  if(plotESDBackgroundXY == kTRUE){ histograms->AddHistogram("ESD_Background_XY" ,"" , nXBinsXY, firstXBinXY, lastXBinXY, nYBinsXY, firstYBinXY, lastYBinXY, "", "");}
  if(plotESDBackgroundRapid == kTRUE){ histograms->AddHistogram("ESD_Background_Rapid" ,"" , nXBinsRapid, firstXBinRapid, lastXBinRapid, "", "");}

	
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
	
  if(plotESDNumberOfV0s == kTRUE){histograms->AddHistogram("ESD_NumberOfV0s","Number of v0s",100, 0, 100,"","");}
  if(plotESDNumberOfSurvivingV0s == kTRUE){histograms->AddHistogram("ESD_NumberOfSurvivingV0s","Number of surviving v0s",100, 0, 100,"","");}
	
  //  debug histograms
  if(plotESDCutGetOnFly == kTRUE){histograms->AddHistogram("ESD_CutGetOnFly_InvMass" ,"Not GetOnFly" , nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass,"","");}
  if(plotESDCutNContributors == kTRUE){histograms->AddHistogram("ESD_CutNContributors_InvMass" ,"NContributors <= 0" , nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass,"","");}
  if(plotESDCutLikeSign == kTRUE){histograms->AddHistogram("ESD_CutLikeSign_InvMass" ,"LikeSign" , nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass,"","");}
  if(plotESDCutRefit == kTRUE){histograms->AddHistogram("ESD_CutRefit_InvMass" ,"No TPC refit" , nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass,"","");}
  if(plotESDCutKink == kTRUE){histograms->AddHistogram("ESD_CutKink_InvMass" ,"Kinks" , nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass,"","");}
  if(plotESDCutPIDProb == kTRUE){histograms->AddHistogram("ESD_CutPIDProb_InvMass" ,"wrong TPC PID" , nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass,"","");}
  if(plotESDCutR == kTRUE){histograms->AddHistogram("ESD_CutR_InvMass" ,"Above RMax" , nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass,"","");}
  if(plotESDCutNDF == kTRUE){histograms->AddHistogram("ESD_CutNDF_InvMass" ,"NDF <= 0" , nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass,"","");}
  if(plotESDCutChi2 == kTRUE){histograms->AddHistogram("ESD_CutChi2_InvMass" ,"#chi^{2} > Max" , nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass,"","");}
  if(plotESDCutEta == kTRUE){histograms->AddHistogram("ESD_CutEta_InvMass" ,"Above #eta max" , nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass,"","");}
  if(plotESDCutPt == kTRUE){histograms->AddHistogram("ESD_CutPt_InvMass" ,"Below p_{t} min" , nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass,"","");}
  if(plotESDCutLine == kTRUE){histograms->AddHistogram("ESD_CutLine_InvMass" ,"Out of reconstruction area" , nXBinsGammaMass, firstXBinGammaMass, lastXBinGammaMass,"","");}
  if(plotESDTrueConvGammaTrackLength == kTRUE){histograms->AddHistogram("ESD_TrueConvGamma_TrackLength","Track length of TrueConvGamma",nXBinsTrackLength,firstXBinTrackLength,lastXBinTrackLength,"","");}
  if(plotESDTrueConvGammaTrackLengthVSInvMass == kTRUE){histograms->AddHistogram("ESD_TrueConvGamma_TrackLengthVSInvMass","Track length of TrueConvGamma vs Inv mass",nXBinsTrackLength,firstXBinTrackLength,lastXBinTrackLength,nXBinsPt, firstXBinPt, lastXBinPt,"","");}


  if(plotPi0Spectra == kTRUE){
histograms->AddHistogram("ESD_Mother_InvMass_vs_Pt" ,"Invariant Mass vs Pt" , nXBinsSpectra, firstXBinSpectra, lastXBinSpectra,nYBinsSpectra, firstYBinSpectra, lastYBinSpectra,"InvMass [GeV]","Pt [GeV]");
}
  if(plotPi0Spectra == kTRUE && calculateBackground == kTRUE){
histograms->AddHistogram("ESD_Background_InvMass_vs_Pt" ,"Background Invariant Mass vs Pt" , nXBinsSpectra, firstXBinSpectra, lastXBinSpectra,nYBinsSpectra, firstYBinSpectra, lastYBinSpectra,"InvMass [GeV]","Pt [GeV]");
}
	
	
	
  //------------------------------ end Creating the histograms -----------------------------
	
  // Create the Analysis manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager", "My Analysis");
	
  // Define Input Event Handler 
  AliESDInputHandler* inpHandler = new AliESDInputHandler();
	
  // Define Output Event Handler
  AliAODHandler* aodHandler = new AliAODHandler();
  TString fileOutAOD = "AOD_"+ outputFileName + outputFileAppendix + ".root";
  aodHandler->SetOutputFileName(fileOutAOD);
  //  aodHandler->SetOutputFileName("aodAliGammaConversion.root");
	
  // Define MC Truth Event Handler
  AliMCEventHandler* mcHandler = new AliMCEventHandler();
	
  // Add Handlers to the Task Manager
  mgr->SetInputEventHandler  (inpHandler);
  mgr->SetOutputEventHandler (aodHandler);
  mgr->SetMCtruthEventHandler(mcHandler);
	
  // Be sure you are told what you are doing
  //  mgr->SetDebugLevel(10);
	
  // Declare Common Input Tchain
  AliAnalysisDataContainer *cinput1 = NULL;
  if(usePWG4PartCorr){
    cinput1 = mgr->CreateContainer("Chain",TChain::Class(),AliAnalysisManager::kInputContainer);
  }
  else{
    cinput1 = mgr->GetCommonInputContainer();
  }
	
  // Common Output Tree in common âdefaultâ output file
  AliAnalysisDataContainer *coutput1 = NULL;
  if(usePWG4PartCorr){
    coutput1 = mgr->CreateContainer("tree",TTree::Class(),AliAnalysisManager::kOutputContainer, "default");
  }
  else{
    coutput1 = mgr->GetCommonOutputContainer();
  }
	
  // Private output objects
  if(outputFileName.Contains(".root")){
    outputFileName.ReplaceAll(".root","");
  }
  if(outputFileAppendix.Contains(".root")){
    outputFileAppendix.ReplaceAll(".root","");
  }
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
  v0Reader->SetLineCutZRSlope(LineCutZRSlope);
  v0Reader->SetLineCutZValue(LineCutZValue);	
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

  gammaconversion->SetMinOpeningAngleGhostCut(minOpeningAngleGhostCut);
	
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
