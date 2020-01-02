#include <iostream>
#include <fstream>

runCascadeAnalysis(TString       lCascType = "XiMinus",
                   TString lWhichEstimator = "V0M",
                   Double_t   lMultBoundLo = 0.00,
                   Double_t   lMultBoundHi = 0.01,
                   TString   lWhichSystVar = "",
                   Bool_t       lSweepFlag = kFALSE){

  TString lWhichMult = Form("%05.0lfto%05.0lf", 100.*lMultBoundLo, 100.*lMultBoundHi); 
  TString lCascTypeBoth = lCascType; lCascTypeBoth.ReplaceAll("Minus", "").ReplaceAll("Plus", "");

  // Load common libraries
  gSystem->Load("libCore.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libOADB");
  gSystem->Load("libPWGLFSTRANGENESS");
  // Use AliRoot includes to compile our task
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");

  //Argument: System analysed - one of "XiMinus", "XiPlus", "OmegaMinus"
  cout<<"Macro to test Cascade analysis module"<<endl;
  cout<<"----------------------------------------------------"<<endl;
  cout<<"               Cascade Analysis Macro "<<endl;
  cout<<"----------------------------------------------------"<<endl;
  cout<<endl;
  cout<<"----------------------------------------------------"<<endl;
  cout<<" ---> Compiling needed class, please wait... "<<endl;
  //Compile Macro
  Int_t workedornot = gSystem->CompileMacro("AliCascadeModule.cxx","-kfo");
  cout<<"----------------------------------------------------"<<endl;
  cout<<endl;
  if( workedornot == 0 ){ 
      cout<<"*************************************"<<endl;
      cout<<" AliCascadeModule.cxx compilation failed! "<<endl;
      cout<<"*************************************"<<endl;
      return;
  }

  //Load Class
  gSystem->Load("AliCascadeModule_cxx.so");

  //Initialize Analysis Object
  AliCascadeModule *casc = new AliCascadeModule(lCascType);
  
  if(lCascType.Contains("Xi") ){
    /// Binning from Fiorella's analysis (2016-11-30) ///////////////////////////
    //
    // // default
    //Double_t ptbinlimits[] = { 0.6, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 2.9, 3.4, 4.0, 5.0, 6.5, 9.0 };
    //
    // // wide
    // Double_t ptbinlimits[] = { 0.6, 2.0, 4.0, 6.5 };
    //
    /////////////////////////////////////////////////////////////////////////////
    //
    // my binning
    //Double_t ptbinlimits[] = { 0.60, 0.80, 1.00, 1.20, 1.40, 1.60, 2.00, 2.60, 3.40, 4.5, 6.0, 8.0, 10.0 }; 
    Double_t ptbinlimits[] = { 0.60, 0.80, 1.00, 1.20, 1.40, 1.60, 1.80, 2.00, 2.20, 2.50, 2.90, 3.40, 4.00, 5.00, 6.50, 8.50, 12.00 };
  }
  if(lCascType.Contains("Omega") ){
    /// Binning from Fiorella's analysis (2016-11-30) ///////////////////////////
    //
    // // default
    //Double_t ptbinlimits[] = { 0.9, 1.6, 2.2, 2.6, 3.0, 3.8, 5.5, 8.0 };
    //
    // // wide
    // Double_t ptbinlimits[] = { 0.9, 3.0, 5.5 };
    //
    /////////////////////////////////////////////////////////////////////////////
    //
    // my binning
    //Double_t ptbinlimits[] = { 0.80, 1.20, 1.60, 1.90, 2.20, 2.50, 3.00, 3.80, 5.00, 8.00 }; 
    Double_t ptbinlimits[]   = { 0.90, 1.20, 1.60, 1.90, 2.20, 2.60, 3.00, 3.80, 5.50, 8.00, 12.0 };
  }

  Long_t ptbinnumb = sizeof(ptbinlimits)/sizeof(Double_t) - 1;

  //Set Pt Bins Used in Analysis
  casc->SetPtBinLimits( ptbinnumb, ptbinlimits );

  //Set Normalization strategy to analysis level
  casc->SetNormalizationStrategy("FinalAnalysis");
  
  //Set Rapidity Window
  //casc->SetRapidityType("CMS");
  //Shift: has to be consistent with the value used in
  //       grid extraction task!
  //casc->SetRapidityShift(0.465);
  casc->SetRapidityWindow(-0.5, 0.5);

  // Configuration to compare with Fiorella's results ---------------------------------------------
  casc->SetDefaultCuts();
  //
  // topological variables
  casc->SetCutCascRadius                                   ( (lCascType.Contains("Xi")) ? 0.6 : 0.5 );
  casc->SetCutV0Radius                                     ( (lCascType.Contains("Xi")) ? 1.2 : 1.1 );
  casc->SetCutDCABachToPV                                  ( 0.04 );
  casc->SetCutDCAV0ToPV                                    ( 0.06 );
  casc->SetCutDCANegToPV                                   ( (lCascType.Contains("Minus")) ? 0.04 : 0.03 );
  casc->SetCutDCAPosToPV                                   ( (lCascType.Contains("Minus")) ? 0.03 : 0.04 );
  casc->SetCutDCAV0Daughters                               ( 1.5 );
  casc->SetCutDCACascDaughters                             ( 1.3 );
  casc->SetCutCascCosPA                                    ( 0.97 );
  casc->SetCutV0CosPA                                      ( 0.97 );
  casc->SetCutV0Mass                                       ( 0.008 );
  //
  // selections
  // Vertex Type.........................................: Offline Vertexer (default)
  casc->SetRapidityWindow                                  ( -0.5, 0.5 ); // MC value used for MC 
  //casc->SetRapidityType                                    ("CMS");
  //Shift: has to be consistent with the value used in
  //       grid extraction task!
  //casc->SetRapidityShift                                   (0.465);
  casc->SetCutTPCPIDNSigmas                                ( 5.0 );
  casc->SetCutProperLifetime                               ( (lCascType.Contains("Xi")) ? (3.*4.91) : (3.*2.461) );
  casc->SetCutCompetingSpecies                             ( (lCascType.Contains("Xi")) ? -1.0 : 0.008);
  casc->SetCutDaughterEta                                  ( 0.8 );
  // Primary Selection (MC only).........................: AliStack::IsPhysicalPrimary()
  // MC Association (MC only)............................: PDG code association for V0s and daughter tracks
  // Tracking flags for daughters........................: kTPCrefit
  casc->SetCutLeastNumberOfClusters                        ( 70 );
  //
  // signal extraction
  casc->SetCutSigmaForSignalExtraction                     ( 5.0 ); // 4.0 );
  //-----------------------------------------------------------------------------------------------
  // Additional cuts for testing the new vertexer
  //casc->SetCutMinTrackLength                               ( 80.0 );
  //-----------------------------------------------------------------------------------------------
  // TOF match
  //casc->SetUseTOFMatchBach( kTRUE );
  //casc->SetUseTOFMatchNeg( kTRUE );
  //casc->SetUseTOFMatchPos( kTRUE );
  //casc->SetUseTOFMatchOne( kTRUE );
  //-----------------------------------------------------------------------------------------------
  // ITS refit
  //casc->SetUseITSrefitBach( kTRUE );
  //casc->SetUseITSrefitNeg( kTRUE );
  //casc->SetUseITSrefitPos( kTRUE );
  //casc->SetUseITSrefitOne( kTRUE );
  //-----------------------------------------------------------------------------------------------
  // TOF match or ITS refit
  //casc->SetUseITSTOFBach( kTRUE );
  //casc->SetUseITSTOFNeg( kTRUE );
  //casc->SetUseITSTOFPos( kTRUE );
  casc->SetUseITSTOFOne( kTRUE );
  //-----------------------------------------------------------------------------------------------
  // Pt-dependent cuts for high-multiplicity analysis
  TFile* fPtDepCuts = new TFile(Form("../compute_ptdep_cuts/PtDepCuts_%s_%s_00000to00005_SL1.root", lCascType.Data(), lWhichEstimator.Data()), "READ");
  //
  if(lCascType.Contains("Xi")) {
      TH1F* hCutBBCosPA = (TH1F*)fPtDepCuts->Get(Form("hCutBBCosPA_SL1_%s", lCascType.Data()));
      casc->SetCutBBCosPA( hCutBBCosPA );
  }
  //
  //TH1F* hCutV0CosPA = (TH1F*)fPtDepCuts->Get(Form("hCutV0CosPA_SL5_%s", lCascType.Data()));
  //casc->SetCutV0CosPA( hCutV0CosPA );
  ////
  //if(lCascType.Contains("Minus")) {
  //    TH1F* hCutDCAzPosToPV = (TH1F*)fPtDepCuts->Get(Form("hCutDCAzPosToPV_SL5_%s", lCascType.Data()));
  //    casc->SetCutDCAzPosToPV( hCutDCAzPosToPV );
  //}
  //if(lCascType.Contains("Plus")) {
  //    TH1F* hCutDCAzNegToPV = (TH1F*)fPtDepCuts->Get(Form("hCutDCAzNegToPV_SL5_%s", lCascType.Data()));
  //    casc->SetCutDCAzNegToPV( hCutDCAzNegToPV );
  //}
  ////
  //TH1F* hCutDCAzBachToPV = (TH1F*)fPtDepCuts->Get(Form("hCutDCAzBachToPV_SL5_%s", lCascType.Data()));
  //casc->SetCutDCAzBachToPV( hCutDCAzBachToPV );
  ////
  ////TH1F* hCutCascCosPA = (TH1F*)fPtDepCuts->Get(Form("hCutCascCosPA_SL5_%s", lCascType.Data()));
  ////casc->SetCutCascCosPA( hCutCascCosPA );
  //-----------------------------------------------------------------------------------------------
  // Additional cuts used for 2016 data only (not applied for 2015 data) --------------------------
  casc->SetMVPileupRejection( kTRUE );
  //casc->SetMinDistToClosestNonEmptyBC( 2 );
  //-----------------------------------------------------------------------------------------------

  //Set CINT1B/INEL Ratio for normalization
  casc->SetCINT1BoverINEL                ( 1.0 );

  //Set further configurations
  casc->SetPerformMultiplicityStudy( kTRUE );
  casc->SetCentralityEstimator( lWhichEstimator );
  casc->SetLowMultValue ( lMultBoundLo );
  casc->SetHighMultValue( lMultBoundHi );
  //
  casc->SetUseIntegratedEfficiency( kTRUE );

  // Geant-Fluka correction for anti-protons ------------------------------------------------------
  if( lCascType.Contains("Plus") ) {
      TFile* fGeantFluka = new TFile("../../G3G4Correction/antiPrCorrFunc.root", "READ");
      TF1* funcGeantFlukaCorr = fGeantFluka->Get("funcCorrAntiProtonGEANT4");
      //TF1* funcGeantFlukaCorr = fGeantFluka->Get("funcCorrAntiProtonFLUKA"); // used for systematics
      casc->SetGeantFlukaCorrection ( funcGeantFlukaCorr );
  }
  //-----------------------------------------------------------------------------------------------

//  //Provide peak position and width from fit ------------------------------------------------------
//  TString lDataSetLow = lDataSet; lDataSetLow.ToLower();
//  TString lSigExtParamsPath     = Form(".", lDataSetLow.Data());
//  TString lSigExtParamsFilename = Form("%s/SigExtParams-%s-13TeV-%s-00000to10000.root", lSigExtParamsPath.Data(), lCascTypeBoth.Data(), lWhichEstimator.Data());
//  cout << "[runCascadeAnalysis] Using user provided peak position and width for signal extraction:" << endl;
//  cout << "                     " << lSigExtParamsFilename << endl;
//  TFile *fSigExtParams = TFile::Open(lSigExtParamsFilename.Data(), "READ");
//  TF1* fPeakPositionFit = (TF1*)fSigExtParams->Get("fMeanFit");
//  TF1* fPeakWidthFit    = (TF1*)fSigExtParams->Get("fSigmaFit");
//  //
//  casc->SetPeakPositionAndWidthFromFit( fPeakPositionFit, fPeakWidthFit );
//  //-----------------------------------------------------------------------------------------------

  //Set data files
  TString lDataFilename = Form("../AnalysisResultsData_V0M_%s.root", lWhichMult.Data());//0x0;
  TString lMCFilename   = "../AnalysisResultsMC.root";//0x0;
  // 
  casc->SetRealDataFile( lDataFilename.Data() );
  casc->SetMCDataFile  ( lMCFilename.Data()   );

  //Set output file name
  cout << "[runCascadeAnalysis] Performing multiplicity selection: [" << lMultBoundLo << "," << lMultBoundHi << "[ (" << lWhichEstimator << ")" << endl;
  if( lWhichEstimator.Contains("V0") ) casc->SetOutputFile( Form("Results-%s-13TeV-%s-%.5ito%.5i.root", lCascType.Data(), lWhichEstimator.Data(), (Int_t)(100.*lMultBoundLo), (Int_t)(100.*lMultBoundHi)) );
  else                                 casc->SetOutputFile( Form("Results-%s-13TeV-%s-%.3ito%.3i.root", lCascType.Data(), lWhichEstimator.Data(), (Int_t)(lMultBoundLo),      (Int_t)(lMultBoundHi))      );

    //Run the analysis 
    if( lWhichSystVar == "" )  casc->DoAnalysis();
    else {

        TString lSystFile = "Results-Systematics";
        if( lWhichEstimator.Contains("V0") ) lSystFile.Append( Form("-%s-13TeV-%s-%.5ito%.5i-",lCascType.Data(),lWhichEstimator.Data(),(Int_t)(100.*lMultBoundLo),(Int_t)(100.*lMultBoundHi)) );
        else                                 lSystFile.Append( Form("-%s-13TeV-%s-%.3ito%.3i-",lCascType.Data(),lWhichEstimator.Data(),(Int_t)(lMultBoundLo),(Int_t)(lMultBoundHi)) );

        //---------------------------------------------------------------------------------------------
        // Topological Selection Variables Systematics

        ////////////////////////////////////////////////////////////////////////////
        // V0RADIUS ////////////////////////////////////////////////////////////////
        if( lWhichSystVar == "V0Radius" ) {
            cout<<endl<<"---> Performing Systematics Studies: V0 Radius"<<endl;
            if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
                casc->SetCutV0Radius(1.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutV0Radius(1.10); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutV0Radius(2.50); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                casc->SetCutV0Radius(5.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
            }
            if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
                casc->SetCutV0Radius(1.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutV0Radius(1.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutV0Radius(2.50); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                casc->SetCutV0Radius(6.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
            }
        }
        ////////////////////////////////////////////////////////////////////////////
        // CASCRADIUS //////////////////////////////////////////////////////////////
        if( lWhichSystVar == "CascRadius" ) {
            cout<<endl<<"---> Performing Systematics Studies: Cascade Radius"<<endl;
            if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
                casc->SetCutCascRadius(0.40); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutCascRadius(0.50); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutCascRadius(0.80); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                casc->SetCutCascRadius(1.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
            }
            if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
                casc->SetCutCascRadius(0.30); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutCascRadius(0.40); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutCascRadius(0.60); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                casc->SetCutCascRadius(0.70); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
            }
        }
        ////////////////////////////////////////////////////////////////////////////
        // DCABACHTOPV /////////////////////////////////////////////////////////////
        if( lWhichSystVar == "DCABachToPV" ) {
            cout<<endl<<"---> Performing Systematics Studies: DCA Bachelor to Primary Vertex"<<endl;
            if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
                casc->SetCutDCABachToPV(0.03); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutDCABachToPV(0.03); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutDCABachToPV(0.10); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                casc->SetCutDCABachToPV(0.17); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
            }
            if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
                casc->SetCutDCABachToPV(0.03); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutDCABachToPV(0.03); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutDCABachToPV(0.07); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                casc->SetCutDCABachToPV(0.10); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
            }
        }
        ////////////////////////////////////////////////////////////////////////////
        // DCAV0TOPV ///////////////////////////////////////////////////////////////
        if( lWhichSystVar == "DCAV0ToPV" ) {
            cout<<endl<<"---> Performing Systematics Studies: DCA V0 to Primary Vertex"<<endl;
            if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
                casc->SetCutDCAV0ToPV(0.05); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutDCAV0ToPV(0.05); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutDCAV0ToPV(0.10); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                casc->SetCutDCAV0ToPV(0.15); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
            }
            if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
                casc->SetCutDCAV0ToPV(0.05); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutDCAV0ToPV(0.05); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutDCAV0ToPV(0.08); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                casc->SetCutDCAV0ToPV(0.10); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
            }
        }
        ////////////////////////////////////////////////////////////////////////////
        // DCANEGTOPV //////////////////////////////////////////////////////////////
        if( lWhichSystVar == "DCANegToPV" ) {
            cout<<endl<<"---> Performing Systematics Studies: DCA Negative to Primary Vertex"<<endl;
            if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
                casc->SetCutDCANegToPV( (lCascType.Contains("Minus")) ? 0.02 : 0.02 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutDCANegToPV( (lCascType.Contains("Minus")) ? 0.03 : 0.02 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutDCANegToPV( (lCascType.Contains("Minus")) ? 0.15 : 0.09 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                casc->SetCutDCANegToPV( (lCascType.Contains("Minus")) ? 0.30 : 0.11 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
            }
            if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
                casc->SetCutDCANegToPV( (lCascType.Contains("Minus")) ? 0.02 : 0.02 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutDCANegToPV( (lCascType.Contains("Minus")) ? 0.03 : 0.02 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutDCANegToPV( (lCascType.Contains("Minus")) ? 0.10 : 0.05 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                casc->SetCutDCANegToPV( (lCascType.Contains("Minus")) ? 0.30 : 0.10 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
            }
        }
        ////////////////////////////////////////////////////////////////////////////
        // DCAPOSTOPV //////////////////////////////////////////////////////////////
        if( lWhichSystVar == "DCAPosToPV" ) {
            cout<<endl<<"---> Performing Systematics Studies: DCA Positive to Primary Vertex"<<endl;
            if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
                casc->SetCutDCAPosToPV( (lCascType.Contains("Minus")) ? 0.02 : 0.02 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutDCAPosToPV( (lCascType.Contains("Minus")) ? 0.02 : 0.03 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutDCAPosToPV( (lCascType.Contains("Minus")) ? 0.09 : 0.15 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                casc->SetCutDCAPosToPV( (lCascType.Contains("Minus")) ? 0.11 : 0.30 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
            }
            if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
                casc->SetCutDCAPosToPV( (lCascType.Contains("Minus")) ? 0.02 : 0.02 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutDCAPosToPV( (lCascType.Contains("Minus")) ? 0.02 : 0.03 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutDCAPosToPV( (lCascType.Contains("Minus")) ? 0.05 : 0.10 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                casc->SetCutDCAPosToPV( (lCascType.Contains("Minus")) ? 0.10 : 0.30 ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
            }
        }
        ////////////////////////////////////////////////////////////////////////////
        // DCAV0DAUGHTERS //////////////////////////////////////////////////////////
        if( lWhichSystVar == "DCAV0Daughters" ) {
            cout<<endl<<"---> Performing Systematics Studies: DCA V0 Daughters"<<endl;
            if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
                casc->SetCutDCAV0Daughters(2.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutDCAV0Daughters(1.80); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutDCAV0Daughters(1.20); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                casc->SetCutDCAV0Daughters(1.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
            }
            if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
                casc->SetCutDCAV0Daughters(2.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutDCAV0Daughters(1.80); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutDCAV0Daughters(1.30); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                casc->SetCutDCAV0Daughters(1.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
            }
        }
        ////////////////////////////////////////////////////////////////////////////
        // DCACASCDAUGHTERS ////////////////////////////////////////////////////////
        if( lWhichSystVar == "DCACascDaughters" ) {
            cout<<endl<<"---> Performing Systematics Studies: DCA Cascade Daughters"<<endl;
            if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
                casc->SetCutDCACascDaughters(2.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutDCACascDaughters(1.80); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutDCACascDaughters(1.20); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                casc->SetCutDCACascDaughters(1.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
            }
            if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
                casc->SetCutDCACascDaughters(2.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutDCACascDaughters(1.80); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutDCACascDaughters(1.00); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                casc->SetCutDCACascDaughters(0.60); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
            }
        }
        ////////////////////////////////////////////////////////////////////////////
        // V0COSPA /////////////////////////////////////////////////////////////////
        if( lWhichSystVar == "V0CosPA" ) {
            cout<<endl<<"---> Performing Systematics Studies: Cosine of Pointing Angle of the V0"<<endl;
            if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
                casc->SetCutV0CosPA(0.95); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutV0CosPA(0.96); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutV0CosPA(0.98); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                casc->SetCutV0CosPA(0.99); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
            }
            if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
                casc->SetCutV0CosPA(0.95); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutV0CosPA(0.96); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutV0CosPA(0.98); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                casc->SetCutV0CosPA(0.99); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
            }
        }
        ////////////////////////////////////////////////////////////////////////////
        // CASCCOSPA //////////////////////////////////////////////////////////////
        if( lWhichSystVar == "CascCosPA" ) {
            cout<<endl<<"---> Performing Systematics Studies: Cosine of Pointing Angle of the Cascade"<<endl;
            if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
                casc->SetCutCascCosPA(0.95); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutCascCosPA(0.96); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutCascCosPA(0.98); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                casc->SetCutCascCosPA(0.99); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
            }
            if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
                casc->SetCutCascCosPA(0.95); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutCascCosPA(0.96); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutCascCosPA(0.98); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                casc->SetCutCascCosPA(0.99); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
            }
        }
        ////////////////////////////////////////////////////////////////////////////
        // V0MASS //////////////////////////////////////////////////////////////////
        if( lWhichSystVar == "V0Mass" ) {
            cout<<endl<<"---> Performing Systematics Studies: V0 Mass Window"<<endl;
            if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
                casc->SetCutV0Mass(0.010); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutV0Mass(0.009); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutV0Mass(0.007); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                casc->SetCutV0Mass(1.006); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
            }
            if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
                casc->SetCutV0Mass(0.010); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutV0Mass(0.009); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutV0Mass(0.007); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                casc->SetCutV0Mass(0.006); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
            }
        }
        ////////////////////////////////////////////////////////////////////////////
//        // DCAZBACHTOPV //////////////////////////////////////////////////////////////
//        if( lWhichSystVar == "DCAzBachToPV" ) {
//            cout<<endl<<"---> Performing Systematics Studies: DCAz Bachelor to Primary Vertex"<<endl;
//            if( lCascType == "XiMinus" || lCascType == "XiPlus" || lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
//                //
//                TFile* fPtDepCuts_loose = new TFile(Form("PtDepCuts_%s_%s_00000to00005_SL2.root", lCascType.Data(), lWhichEstimator.Data()), "READ");
//                TH1F* hCutDCAzBachToPV_loose = (TH1F*)fPtDepCuts_loose->Get(Form("hCutDCAzBachToPV_SL2_%s", lCascType.Data()));
//                casc->SetCutDCAzBachToPV( hCutDCAzBachToPV_loose );
//                casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); 
//                casc->DoAnalysis();
//                casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); 
//                casc->DoAnalysis();
//                //
//                TFile* fPtDepCuts_tight = new TFile(Form("PtDepCuts_%s_%s_00000to00005_SL10.root", lCascType.Data(), lWhichEstimator.Data()), "READ");
//                TH1F* hCutDCAzBachToPV_tight = (TH1F*)fPtDepCuts_tight->Get(Form("hCutDCAzBachToPV_SL10_%s", lCascType.Data()));
//                casc->SetCutDCAzBachToPV( hCutDCAzBachToPV_tight );
//                casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); 
//                casc->DoAnalysis();
//                casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); 
//                casc->DoAnalysis();
//            }
//        }
//        ////////////////////////////////////////////////////////////////////////////
        // BBCOSPA /////////////////////////////////////////////////////////////////
        if( lWhichSystVar == "BBCosPA" ) {
            cout<<endl<<"---> Performing Systematics Studies: Bachelor-Baryon Cosine of Pointing Angle"<<endl;
            if( lCascType == "XiMinus" || lCascType == "XiPlus") {
                //
                TFile* fPtDepCuts_systematics = new TFile(Form("../PtDepCuts_%s_%s_00000to00005_SL1.root", lCascType.Data(), lWhichEstimator.Data()), "READ");
                TH1F* hCutBBCosPA_loose = (TH1F*)((TH1F*)fPtDepCuts_systematics->Get(Form("hCutBBCosPA_SL1_%s", lCascType.Data())))->Clone("hCutBBCosPA_loose");
                TH1F* hCutBBCosPA_tight = (TH1F*)((TH1F*)fPtDepCuts_systematics->Get(Form("hCutBBCosPA_SL1_%s", lCascType.Data())))->Clone("hCutBBCosPA_tight");
                //
                TH1F* hUnity = (TH1F*)hCutBBCosPA_loose->Clone("hUnity");
                hUnity->Reset(); for(Int_t i=0; i<hUnity->GetNbinsX(); ++i) hUnity->SetBinContent(i+1, 1.);
                hCutBBCosPA_loose->Add(hUnity);
                hCutBBCosPA_loose->Scale(0.5);
                casc->SetCutBBCosPA( hCutBBCosPA_loose );
                casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); 
                casc->DoAnalysis();
                //casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); 
                //casc->DoAnalysis();
                TString systname = lSystFile + lWhichSystVar;
                gSystem->Exec(Form("cp %s-1.root %s-2.root", systname.Data(), systname.Data()));
                //
                hCutBBCosPA_tight->Scale(1.0); // the default is already the tightest
                casc->SetCutBBCosPA( hCutBBCosPA_tight );
                casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); 
                casc->DoAnalysis();
                //casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); 
                //casc->DoAnalysis();
                gSystem->Exec(Form("cp %s-3.root %s-4.root", systname.Data(), systname.Data()));
            }
        }
        ////////////////////////////////////////////////////////////////////////////

        //---------------------------------------------------------------------------------------------
        // Other Selection Variables Systematics

        ////////////////////////////////////////////////////////////////////////////
        // TPCPIDNSIGMAS ///////////////////////////////////////////////////////////
        if( lWhichSystVar == "TPCPIDNSigmas" ) {
            cout<<endl<<"---> Performing Systematics Studies: TPC PID Nsigmas"<<endl;
            if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
                casc->SetCutTPCPIDNSigmas(7.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutTPCPIDNSigmas(6.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutTPCPIDNSigmas(4.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                //casc->SetCutTPCPIDNSigmas(4.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
                //
                TString systname = lSystFile + lWhichSystVar;
                gSystem->Exec(Form("cp %s-3.root %s-4.root", systname.Data(), systname.Data()));
            }
            if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
                casc->SetCutTPCPIDNSigmas(7.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutTPCPIDNSigmas(6.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutTPCPIDNSigmas(4.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                //casc->SetCutTPCPIDNSigmas(4.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
                //
                TString systname = lSystFile + lWhichSystVar;
                gSystem->Exec(Form("cp %s-3.root %s-4.root", systname.Data(), systname.Data()));
            }
        }
        ////////////////////////////////////////////////////////////////////////////
        // TPCNCLUSTERS ///////////////////////////////////////////////////////////
        if( lWhichSystVar == "TPCNClusters" ) {
            cout<<endl<<"---> Performing Systematics Studies: TPC Least Number of Clusters"<<endl;
            if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
                casc->SetCutLeastNumberOfClusters(70); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                //casc->SetCutLeastNumberOfClusters(70); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutLeastNumberOfClusters(75); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                casc->SetCutLeastNumberOfClusters(80); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
                //
                TString systname = lSystFile + lWhichSystVar;
                gSystem->Exec(Form("cp %s-1.root %s-2.root", systname.Data(), systname.Data()));
            }
            if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
                casc->SetCutLeastNumberOfClusters(70); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                //casc->SetCutLeastNumberOfClusters(70); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutLeastNumberOfClusters(73); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                casc->SetCutLeastNumberOfClusters(76); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
                //
                TString systname = lSystFile + lWhichSystVar;
                gSystem->Exec(Form("cp %s-1.root %s-2.root", systname.Data(), systname.Data()));
            }
        }
        ////////////////////////////////////////////////////////////////////////////
        // PROPERLIFETIME //////////////////////////////////////////////////////////
        if( lWhichSystVar == "ProperLifetime" ) {
            cout<<endl<<"---> Performing Systematics Studies: Proper Lifetime"<<endl;
            if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
                casc->SetCutProperLifetime(5.0*4.91); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutProperLifetime(4.0*4.91); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutProperLifetime(2.5*4.91); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                //casc->SetCutProperLifetime(2.5*4.91); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
                //
                TString systname = lSystFile + lWhichSystVar;
                gSystem->Exec(Form("cp %s-3.root %s-4.root", systname.Data(), systname.Data()));
            }
            if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
                casc->SetCutProperLifetime(5.0*2.461); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                casc->SetCutProperLifetime(4.0*2.461); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutProperLifetime(3.0*2.461); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                //casc->SetCutProperLifetime(3.0*2.461); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
                //
                TString systname = lSystFile + lWhichSystVar;
                gSystem->Exec(Form("cp %s-3.root %s-4.root", systname.Data(), systname.Data()));
            }
        }
        ////////////////////////////////////////////////////////////////////////////
        // COMPETINGSPECIES ////////////////////////////////////////////////////////
        if( lWhichSystVar == "CompetingSpecies" ) {
            cout<<endl<<"---> Performing Systematics Studies: Competing Species"<<endl;
            if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
                casc->SetCutCompetingSpecies(-1.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                //casc->SetCutCompetingSpecies(-1.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                //casc->SetCutCompetingSpecies(-1.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                //casc->SetCutCompetingSpecies(-1.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
                //
                TString systname = lSystFile + lWhichSystVar;
                gSystem->Exec(Form("cp %s-1.root %s-2.root", systname.Data(), systname.Data()));
                gSystem->Exec(Form("cp %s-1.root %s-3.root", systname.Data(), systname.Data()));
                gSystem->Exec(Form("cp %s-1.root %s-4.root", systname.Data(), systname.Data()));
            }
            if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
                casc->SetCutCompetingSpecies( -1.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                //casc->SetCutCompetingSpecies( -1.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutCompetingSpecies(0.008); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                //casc->SetCutCompetingSpecies(0.008); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
                //
                TString systname = lSystFile + lWhichSystVar;
                gSystem->Exec(Form("cp %s-1.root %s-2.root", systname.Data(), systname.Data()));
                gSystem->Exec(Form("cp %s-3.root %s-4.root", systname.Data(), systname.Data()));
            }
        }
        ////////////////////////////////////////////////////////////////////////////
        // SIGEXTBINCOUNT //////////////////////////////////////////////////////////
        if( lWhichSystVar == "SigExtBinCount" ) {
            cout<<endl<<"---> Performing Systematics Studies: Sigma Cut for Signal Extraction (bin count)"<<endl;
            if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
                casc->SetCutSigmaForSignalExtraction(5.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                //casc->SetCutSigmaForSignalExtraction(5.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutSigmaForSignalExtraction(3.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                //casc->SetCutSigmaForSignalExtraction(3.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
                //
                TString systname = lSystFile + lWhichSystVar;
                gSystem->Exec(Form("cp %s-1.root %s-2.root", systname.Data(), systname.Data()));
                gSystem->Exec(Form("cp %s-3.root %s-4.root", systname.Data(), systname.Data()));
            }
            if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
                casc->SetCutSigmaForSignalExtraction(4.5); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                //casc->SetCutSigmaForSignalExtraction(4.5); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutSigmaForSignalExtraction(3.5); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                //casc->SetCutSigmaForSignalExtraction(3.5); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
                //
                TString systname = lSystFile + lWhichSystVar;
                gSystem->Exec(Form("cp %s-1.root %s-2.root", systname.Data(), systname.Data()));
                gSystem->Exec(Form("cp %s-3.root %s-4.root", systname.Data(), systname.Data()));
            }
        }
        ////////////////////////////////////////////////////////////////////////////
        // SIGEXTFITTING ///////////////////////////////////////////////////////////
        if( lWhichSystVar == "SigExtFitting" ) {
            cout<<endl<<"---> Performing Systematics Studies: Sigma Cut for Signal Extraction (fitting)"<<endl;
            casc->SetFitBackground( kTRUE );
            if( lCascType == "XiMinus" || lCascType == "XiPlus" ) {
                casc->SetCutSigmaForSignalExtraction(5.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                //casc->SetCutSigmaForSignalExtraction(5.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutSigmaForSignalExtraction(3.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                //casc->SetCutSigmaForSignalExtraction(3.0); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
                //
                TString systname = lSystFile + lWhichSystVar;
                gSystem->Exec(Form("cp %s-1.root %s-2.root", systname.Data(), systname.Data()));
                gSystem->Exec(Form("cp %s-3.root %s-4.root", systname.Data(), systname.Data()));
            }
            if( lCascType == "OmegaMinus" || lCascType == "OmegaPlus" ) {
                casc->SetCutSigmaForSignalExtraction(4.5); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                //casc->SetCutSigmaForSignalExtraction(4.5); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                casc->SetCutSigmaForSignalExtraction(3.5); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                //casc->SetCutSigmaForSignalExtraction(3.5); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
                //
                TString systname = lSystFile + lWhichSystVar;
                gSystem->Exec(Form("cp %s-1.root %s-2.root", systname.Data(), systname.Data()));
                gSystem->Exec(Form("cp %s-3.root %s-4.root", systname.Data(), systname.Data()));
            }
        }
        ////////////////////////////////////////////////////////////////////////////
        // GEANTFLUKACORRECTION ////////////////////////////////////////////////////
        if( lWhichSystVar == "GeantFlukaCorrection" ) {
            cout<<endl<<"---> Performing Systematics Studies: Geant-Fluka correction"<<endl;
            if( lCascType == "XiPlus" || lCascType == "OmegaPlus" ) {
                TFile* fGeantFluka = new TFile("antiPrCorrFunc.root", "READ");
                TF1* funcGeantFlukaCorr = (TF1*)fGeantFluka->Get("funcCorrAntiProtonFLUKA"); // used for systematics
                casc->SetGeantFlukaCorrection( funcGeantFlukaCorr ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-1.root" ); casc->DoAnalysis();
                //casc->SetGeantFlukaCorrection( funcGeantFlukaCorr ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-2.root" ); casc->DoAnalysis();
                //casc->SetGeantFlukaCorrection( funcGeantFlukaCorr ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-3.root" ); casc->DoAnalysis();
                //casc->SetGeantFlukaCorrection( funcGeantFlukaCorr ); casc->SetOutputFile( lSystFile + lWhichSystVar + "-4.root" ); casc->DoAnalysis();
                //
                TString systname = lSystFile + lWhichSystVar;
                gSystem->Exec(Form("cp %s-1.root %s-2.root", systname.Data(), systname.Data()));
                gSystem->Exec(Form("cp %s-1.root %s-3.root", systname.Data(), systname.Data()));
                gSystem->Exec(Form("cp %s-1.root %s-4.root", systname.Data(), systname.Data()));
            }
        }
        ////////////////////////////////////////////////////////////////////////////
    }
    //---------------------------------------------------------------------------------------------





}

