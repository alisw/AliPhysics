runV0Analysis(TString lV0Type = "Lambda", Bool_t doSystematics = kFALSE){
  //Argument: System analysed - one of "Lambda", "AntiLambda", "K0Short"
  cout<<"Macro to test V0 analysis module"<<endl;

  cout<<"----------------------------------------------------"<<endl;
  cout<<"               V0 Analysis Macro "<<endl;
  cout<<"----------------------------------------------------"<<endl;
  cout<<endl;
  cout<<"----------------------------------------------------"<<endl;
  cout<<" ---> Compiling needed class, please wait... "<<endl;
  //Compile Macro
  Int_t workedornot = gSystem->CompileMacro("AliV0Module.cxx","-kfo");
  cout<<"----------------------------------------------------"<<endl;
  cout<<endl;
   if( workedornot == 0 ){ 
      cout<<"*************************************"<<endl;
      cout<<" AliV0Module.cxx compilation failed! "<<endl;
      cout<<"*************************************"<<endl;
      return;
   }

  //Load Class
  gSystem->Load("AliV0Module_cxx.so");

  //Initialize Analysis Object
  AliV0Module *v0 = new AliV0Module(lV0Type);

  if ( lV0Type == "Lambda" || lV0Type == "AntiLambda" ){
    //Low-Stat binning
    Double_t ptbinlimits[] = {0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2.0,2.3,2.6,2.9,3.4,4.0,5.0};
    //High-Stat binning
    //Double_t ptbinlimits[] = {0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4,4.4,4.9,5.5,6.0,7,8,10};
  }
  if ( lV0Type == "K0Short"){
    //Low-Stat binning
    //Double_t ptbinlimits[] = {
    //  0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4,
    //  1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 5.0 };
    //High-Stat binning
    Double_t ptbinlimits[] ={ 0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8, 0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.,2.2,2.4,2.6,2.8,3.0,
    3.3,3.6,3.9,4.2,4.6,5,5.4,5.9, 6.5,7,7.5,8,8.5,9.2,10,11,12,13.5,15};
  }

  Long_t ptbinnumb = sizeof(ptbinlimits)/sizeof(Double_t) - 1;

  //Set Pt Bins Used in Analysis
  v0->SetPtBinLimits( ptbinnumb, ptbinlimits );

  //Set Default Cuts - note: Particle dependent
  v0->SetDefaultCuts(); 

  //Feeddown treatment switch (applies only to Lambda, AntiLambda) 
  v0->SetFeeddownTreatment             ("UseMCRatio");

  //Set CINT1B/INEL Ratio for normalization
  v0->SetCINT1BoverINEL                ( 0.851 );

  //Set Filenames 
  v0->SetRealDataFile    ( "MergedAnalysisResults-125085-LHC10d.root" );
  v0->SetMCDataFile      ( "MergedAnalysisResults-125085-LHC10f6a.root"   );
  v0->SetFeedDownDataFile( "MergedAnalysisResults-125085-LHC10f6a.root"   );

  if ( lV0Type == "Lambda"     ) v0->SetOutputFile( "Results-Lambda-7TeV.root"     );
  if ( lV0Type == "AntiLambda" ) v0->SetOutputFile( "Results-AntiLambda-7TeV.root" );
  if ( lV0Type == "K0Short"    ) v0->SetOutputFile( "Results-K0Short-7TeV.root"    );

  //Run the analysis 
  v0->DoAnalysis();

  //Attempt DoubleChargedXi feeddown and see what it looks like
  if(lV0Type == "Lambda" || lV0Type == "AntiLambda" ){   
    v0->SetFeeddownTreatment("DoubleChargedXi");
    if ( lV0Type == "Lambda"     ) v0->SetOutputFile( "Results-Lambda-7TeV-DoubleChargedXiFD.root"     );
    if ( lV0Type == "AntiLambda" ) v0->SetOutputFile( "Results-AntiLambda-7TeV-DoubleChargedXiFD.root" );
    v0->DoAnalysis();
  }

  //Perform Systematics if demanded
  if ( doSystematics ){ 
    //------------------------------------------------------------
    // Topological Selection Variables Systematics
    TString lSystFile = "Results-Systematics"; 
    if(lV0Type == "Lambda"    ) lSystFile.Append("-Lambda-");
    if(lV0Type == "AntiLambda") lSystFile.Append("-AntiLambda-");
    if(lV0Type == "K0Short"   ) lSystFile.Append("-K0Short-");

    cout<<endl<<"---> Performing Systematics Studies: V0 Decay Radius"<<endl;
    v0->SetDefaultCuts(); 
    v0->SetCutV0Radius(0.300); v0->SetOutputFile( lSystFile + "V0Radius-1.root" ); v0->DoAnalysis(); 
    v0->SetCutV0Radius(0.400); v0->SetOutputFile( lSystFile + "V0Radius-2.root" ); v0->DoAnalysis(); 
    v0->SetCutV0Radius(0.600); v0->SetOutputFile( lSystFile + "V0Radius-3.root" ); v0->DoAnalysis(); 
    v0->SetCutV0Radius(0.700); v0->SetOutputFile( lSystFile + "V0Radius-4.root" ); v0->DoAnalysis(); 

    cout<<endl<<"---> Performing Systematics Studies: DCA Neg track to PV"<<endl;
    v0->SetDefaultCuts(); 
    v0->SetCutDCANegToPV(0.050); v0->SetOutputFile( lSystFile + "DCANegToPV-1.root" ); v0->DoAnalysis(); 
    v0->SetCutDCANegToPV(0.055); v0->SetOutputFile( lSystFile + "DCANegToPV-2.root" ); v0->DoAnalysis(); 
    v0->SetCutDCANegToPV(0.070); v0->SetOutputFile( lSystFile + "DCANegToPV-3.root" ); v0->DoAnalysis(); 
    v0->SetCutDCANegToPV(0.080); v0->SetOutputFile( lSystFile + "DCANegToPV-4.root" ); v0->DoAnalysis(); 

    cout<<endl<<"---> Performing Systematics Studies: DCA Pos track to PV"<<endl;
    v0->SetDefaultCuts(); 
    v0->SetCutDCAPosToPV(0.050); v0->SetOutputFile( lSystFile + "DCAPosToPV-1.root" ); v0->DoAnalysis(); 
    v0->SetCutDCAPosToPV(0.055); v0->SetOutputFile( lSystFile + "DCAPosToPV-2.root" ); v0->DoAnalysis(); 
    v0->SetCutDCAPosToPV(0.070); v0->SetOutputFile( lSystFile + "DCAPosToPV-3.root" ); v0->DoAnalysis(); 
    v0->SetCutDCAPosToPV(0.080); v0->SetOutputFile( lSystFile + "DCAPosToPV-4.root" ); v0->DoAnalysis(); 

    cout<<endl<<"---> Performing Systematics Studies: DCA Pos track to PV"<<endl;
    v0->SetDefaultCuts(); 
    if( lV0Type == "Lambda" || lV0Type == "AntiLambda" ){ 
      v0->SetCutV0CosPA(0.993); v0->SetOutputFile( lSystFile + "V0CosPA-1.root" ); v0->DoAnalysis(); 
      v0->SetCutV0CosPA(0.994); v0->SetOutputFile( lSystFile + "V0CosPA-2.root" ); v0->DoAnalysis(); 
      v0->SetCutV0CosPA(0.996); v0->SetOutputFile( lSystFile + "V0CosPA-3.root" ); v0->DoAnalysis(); 
      v0->SetCutV0CosPA(0.997); v0->SetOutputFile( lSystFile + "V0CosPA-4.root" ); v0->DoAnalysis(); 
    } else {
      v0->SetCutV0CosPA(0.95); v0->SetOutputFile( lSystFile + "V0CosPA-1.root" ); v0->DoAnalysis(); 
      v0->SetCutV0CosPA(0.96); v0->SetOutputFile( lSystFile + "V0CosPA-2.root" ); v0->DoAnalysis(); 
      v0->SetCutV0CosPA(0.98); v0->SetOutputFile( lSystFile + "V0CosPA-3.root" ); v0->DoAnalysis(); 
      v0->SetCutV0CosPA(0.99); v0->SetOutputFile( lSystFile + "V0CosPA-4.root" ); v0->DoAnalysis(); 
    }
    cout<<endl<<"---> Performing Systematics Studies: DCA V0 Daughters"<<endl;
    v0->SetDefaultCuts(); 
    v0->SetCutDCAPosToPV(0.050); v0->SetOutputFile( lSystFile + "DCAV0Daughters-1.root" ); v0->DoAnalysis(); 
    v0->SetCutDCAPosToPV(0.055); v0->SetOutputFile( lSystFile + "DCAV0Daughters-2.root" ); v0->DoAnalysis(); 
    v0->SetCutDCAPosToPV(0.070); v0->SetOutputFile( lSystFile + "DCAV0Daughters-3.root" ); v0->DoAnalysis(); 
    v0->SetCutDCAPosToPV(0.080); v0->SetOutputFile( lSystFile + "DCAV0Daughters-4.root" ); v0->DoAnalysis(); 
    //End topological selection variables systematics
    //------------------------------------------------------------

    //other systematics: could be placed here in the future...
  }

}

