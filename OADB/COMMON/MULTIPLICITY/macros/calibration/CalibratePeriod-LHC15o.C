

CalibratePeriod_LHC15o(TString lPeriodName = "LHC15o"){

  //Load ALICE stuff
  TString gLibs[] =   {"STEER",
      "ANALYSIS", "ANALYSISalice", "ANALYSIScalib"};
  TString thislib = "lib";
  for(Int_t ilib = 0; ilib<4; ilib++){
    thislib="lib";
    thislib.Append(gLibs[ilib].Data());
    cout<<"Will load "<<thislib.Data()<<endl;
    gSystem->Load(thislib.Data());
  }
  gSystem->SetIncludePath("-I$ROOTSYS/include  -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

  cout<<"Alive! "<<endl;

  //All fine, let's try the calibrator
  AliMultSelectionCalibrator *lCalib = new AliMultSelectionCalibrator("lCalib");
  AliMultSelection *lMultSelDefault = new AliMultSelection();

  //============================================================
  // --- Definition of Boundaries ---
  //============================================================

  //Set Adaptive Percentile Boundaries, adjust if finer selection desired
  Double_t lDesiredBoundaries[1000];
  Long_t   lNDesiredBoundaries=0;
  lDesiredBoundaries[0] = 100;
  /*
  //From Low To High Multiplicity
  for( Int_t ib = 1; ib < 91; ib++) {
      lNDesiredBoundaries++;
      lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] - 1.0;
  }
  for( Int_t ib = 1; ib < 91; ib++) {
      lNDesiredBoundaries++;
      lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] - 0.1;
  }
  for( Int_t ib = 1; ib < 91; ib++) {
      lNDesiredBoundaries++;
      lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] - 0.01;
  }
  for( Int_t ib = 1; ib < 101; ib++) {
      lNDesiredBoundaries++;
      lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] - 0.001;
  }
  lNDesiredBoundaries++;
    */

  //Very simple 1%-wide bins all the way
  for( Int_t ib = 1; ib < 101; ib++) {
      lNDesiredBoundaries++;
      lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] - 1.0;
  }
  lNDesiredBoundaries++;
  lDesiredBoundaries[lNDesiredBoundaries] = 0;

  //cout<< "Dump for debug: "<<endl;
  //for( Int_t ib=0;ib<101; ib++)
  //    cout<<"Boundary #"<<ib<<" at "<<lDesiredBoundaries[ib]<<endl;

  lCalib->SetBoundaries( lNDesiredBoundaries, lDesiredBoundaries );
  cout<<"Boundaries set. Will attempt calibration now... "<<endl;

  //// ANCHOR POINT /////
  Double_t lDefaultV0MAnchor     = 0;
  Double_t lDefaultV0MPercentile = 0;

  Double_t lDefaultCL0Anchor     = 0;
  Double_t lDefaultCL0Percentile = 0;

  Double_t lDefaultCL1Anchor     = 0;
  Double_t lDefaultCL1Percentile = 0;

  if ( lPeriodName.Contains("LHC15o") ){
    lDefaultV0MAnchor     = 133.5;
    lDefaultV0MPercentile = 90.007;
    lDefaultCL0Anchor     = 35.5;
    lDefaultCL0Percentile = 90.04;
    lDefaultCL1Anchor     = 30.5;
    lDefaultCL1Percentile = 90.09;
  }
  if ( lPeriodName.Contains("LHC15m") ){
    lDefaultV0MAnchor     = 115.0;
    lDefaultV0MPercentile = 87.5;
    lDefaultCL0Anchor     = 39.5;
    lDefaultCL0Percentile = 88.9;
    lDefaultCL1Anchor     = 40.5;
    lDefaultCL1Percentile = 88.1;
  }
  if ( lPeriodName.Contains("LHC10h") ){
    lDefaultV0MAnchor     = 81.0;
    lDefaultV0MPercentile = 90.00;
    lDefaultCL0Anchor     = 23.5;
    lDefaultCL0Percentile = 90.20;
    lDefaultCL1Anchor     = 21.5;
    lDefaultCL1Percentile = 90.19;
  }

  ///// EVENT SELECTION /////
  if ( lPeriodName.Contains("LHC10h") ){
    cout<<"Setting event selection criteria for Pb-Pb..."<<endl;
    lCalib->GetEventCuts()->SetVzCut(10.0);
    lCalib->GetEventCuts()->SetTriggerCut                (kTRUE );
    lCalib->GetEventCuts()->SetINELgtZEROCut             (kFALSE);
    lCalib->GetEventCuts()->SetTrackletsVsClustersCut    (kFALSE );
    lCalib->GetEventCuts()->SetRejectPileupInMultBinsCut (kFALSE);
    lCalib->GetEventCuts()->SetVertexConsistencyCut      (kFALSE);
    lCalib->GetEventCuts()->SetNonZeroNContribs          (kTRUE);
  }

  if ( lPeriodName.Contains("LHC15m") ||lPeriodName.Contains("LHC15o") ){
    cout<<"Setting event selection criteria for Pb-Pb..."<<endl;
    lCalib->GetEventCuts()->SetVzCut(10.0);
    lCalib->GetEventCuts()->SetTriggerCut                (kTRUE );
    lCalib->GetEventCuts()->SetINELgtZEROCut             (kFALSE);
    lCalib->GetEventCuts()->SetTrackletsVsClustersCut    (kFALSE );
    lCalib->GetEventCuts()->SetRejectPileupInMultBinsCut (kFALSE);
    lCalib->GetEventCuts()->SetVertexConsistencyCut      (kFALSE);
    lCalib->GetEventCuts()->SetNonZeroNContribs          (kTRUE );
  }

  //============================================================
  // --- Definition of Input Variables ---
  //============================================================

  lCalib->SetupStandardInput();

  //============================================================
  // --- Definition of Estimators ---
  //============================================================
  //AliMultEstimator *fEstV0A = new AliMultEstimator("V0A", "", "(fAmplitude_V0A)");
  //AliMultEstimator *fEstV0C = new AliMultEstimator("V0C", "", "(fAmplitude_V0C)");

  // AliMultEstimator *fEstV0M = new AliMultEstimator("V0M", "", "(((fAmplitude_V0A)+(fAmplitude_V0C)) / (1 + ((fEvSel_VtxZ)-1.0)*(0.003 + ((fEvSel_VtxZ)-1.0)*(-0.0004)))");
  // AliMultEstimator *fEstV0Mplus05 = new AliMultEstimator("V0Mplus05", "", "(((fAmplitude_V0A)+(fAmplitude_V0C)) / (1 + ((fEvSel_VtxZ)-1.0)*(0.003 + ((fEvSel_VtxZ)-1.0)*(-0.0004)))");
  // AliMultEstimator *fEstV0Mplus10 = new AliMultEstimator("V0Mplus10", "", "(((fAmplitude_V0A)+(fAmplitude_V0C)) / (1 + ((fEvSel_VtxZ)-1.0)*(0.003 + ((fEvSel_VtxZ)-1.0)*(-0.0004)))");
  // AliMultEstimator *fEstV0Mminus05 = new AliMultEstimator("V0Mminus05", "", "(((fAmplitude_V0A)+(fAmplitude_V0C)) / (1 + ((fEvSel_VtxZ)-1.0)*(0.003 + ((fEvSel_VtxZ)-1.0)*(-0.0004)))");
  // AliMultEstimator *fEstV0Mminus10 = new AliMultEstimator("V0Mminus10", "", "(((fAmplitude_V0A)+(fAmplitude_V0C)) / (1 + ((fEvSel_VtxZ)-1.0)*(0.003 + ((fEvSel_VtxZ)-1.0)*(-0.0004)))");

  //WARNING: DO NOT FORGET PARENTHESIS!
  //TString lV0MDefinition = "((fAmplitude_V0A)+(fAmplitude_V0C)) / (1 + ((fEvSel_VtxZ)-1.0)*(0.001 + ((fEvSel_VtxZ)-1.0)*(0.0001) + ((fEvSel_VtxZ)-1.0)*((fEvSel_VtxZ)-1.0)*0.000008) + ((fEvSel_VtxZ)-1.0)*((fEvSel_VtxZ)-1.0)*((fEvSel_VtxZ)-1.0)*(-0.000001))";

  //New string, from the 29th January 2016
  //TString lV0MDefinition = "((fAmplitude_V0A)+(fAmplitude_V0C)) / (1 + ((fEvSel_VtxZ)-1.0)*(0.001 + ((fEvSel_VtxZ)-1.0)*(-1.0e-04) + TMath::Power(((fEvSel_VtxZ)-1.0),2)*(4.0e-06) + TMath::Power(((fEvSel_VtxZ)-1.0),3)*(7.0e-07) + TMath::Power(((fEvSel_VtxZ)-1.0),4)*(8.0e-09) + TMath::Power(((fEvSel_VtxZ)-1.0),5)*(-4.0e-09)))";


  TString formulaA ="[1]+[2]*(x-[0])+[3]*TMath::Power(x-[0],2)+[4]*TMath::Power(x-[0],3)+[5]*TMath::Power(x-[0],4)";

  TString formulaB = "[1]+[2]*(x-[0])+[3]*TMath::Power(x-[0],2)";
  TString formulaC = "[1]+[2]*(x-[0])+[3]*TMath::Power(x-[0],2)";

  Double_t lParametersA[6] ={3.50488e-04,9.99677e-01,1.47706e-03,-4.44979e-06,5.85466e-06,-6.43362e-07};
  Double_t lParametersB[4]= {1.68730e+00,-1.56165e+00,-3.05809e-01,-9.40820e-03 };
  Double_t lParametersC[4]= {4.08987e+00,4.46269e-01,1.05592e-01,-4.94802e-03};

  for( Int_t i=0; i<6; i++) formulaA.ReplaceAll(Form("[%i]",i),Form("(%.10f)",lParametersA[i]));
  for( Int_t i=0; i<4; i++) formulaB.ReplaceAll(Form("[%i]",i),Form("(%.10f)",lParametersB[i]));
  for( Int_t i=0; i<4; i++) formulaC.ReplaceAll(Form("[%i]",i),Form("(%.10f)",lParametersC[i]));


  TString formulaTot;
  formulaTot = Form("(x<-14.5)*(%s) + (TMath::Abs(x)<14.5)*(%s)+(x>14.5)*(%s)", formulaB.Data(), formulaA.Data(), formulaC.Data());
  formulaTot.ReplaceAll("x","(fEvSel_VtxZ)");

  TString lV0MDefinition = "";
  lV0MDefinition.Append(Form("((fAmplitude_V0A)+(fAmplitude_V0C))/(%s)",formulaTot.Data()));

  cout<<"Constructed V0M Definition from scratch: "<<endl;
  cout<<lV0MDefinition.Data()<<endl;

  AliMultEstimator *fEstV0M = new AliMultEstimator("V0M", "", lV0MDefinition.Data());
  AliMultEstimator *fEstV0Mplus05 = new AliMultEstimator("V0Mplus05", "", lV0MDefinition.Data());
  AliMultEstimator *fEstV0Mplus10 = new AliMultEstimator("V0Mplus10", "", lV0MDefinition.Data());
  AliMultEstimator *fEstV0Mminus05 = new AliMultEstimator("V0Mminus05", "", lV0MDefinition.Data());
  AliMultEstimator *fEstV0Mminus10 = new AliMultEstimator("V0Mminus10", "", lV0MDefinition.Data());

  if( lPeriodName.Contains("LHC10h") ) {
    fEstV0M -> SetUseAnchor        ( kTRUE   ) ;
    fEstV0M -> SetAnchorPoint      ( lDefaultV0MAnchor    ) ;
    fEstV0M -> SetAnchorPercentile ( lDefaultV0MPercentile     ) ;

    fEstV0Mplus05 -> SetUseAnchor        ( kTRUE   ) ;
    fEstV0Mplus05 -> SetAnchorPoint      ( lDefaultV0MAnchor   ) ;
    fEstV0Mplus05 -> SetAnchorPercentile ( lDefaultV0MPercentile +0.5  ) ;

    fEstV0Mplus10 -> SetUseAnchor        ( kTRUE   ) ;
    fEstV0Mplus10 -> SetAnchorPoint      ( lDefaultV0MAnchor   ) ;
    fEstV0Mplus10 -> SetAnchorPercentile ( lDefaultV0MPercentile +1.  ) ;

    fEstV0Mminus05 -> SetUseAnchor        ( kTRUE   ) ;
    fEstV0Mminus05 -> SetAnchorPoint      ( lDefaultV0MAnchor   ) ;
    fEstV0Mminus05 -> SetAnchorPercentile ( lDefaultV0MPercentile -0.5  ) ;

    fEstV0Mminus10 -> SetUseAnchor        ( kTRUE   ) ;
    fEstV0Mminus10 -> SetAnchorPoint      ( lDefaultV0MAnchor   ) ;
    fEstV0Mminus10 -> SetAnchorPercentile ( lDefaultV0MPercentile -1.  ) ;

  }
  if( lPeriodName.Contains("LHC15m") ) {
    fEstV0M -> SetUseAnchor        ( kTRUE   ) ;
    fEstV0M -> SetAnchorPoint      ( lDefaultV0MAnchor     ) ;
    fEstV0M -> SetAnchorPercentile ( lDefaultV0MPercentile ) ;
  }
  if( lPeriodName.Contains("LHC15o") ) {
    AliMultEstimator *fEstV0A = new AliMultEstimator("V0A", "", "((fAmplitude_V0A)/(1 + ((fEvSel_VtxZ)-1.11974)*(-0.00529266 - ((fEvSel_VtxZ)-1.11974)*0.000153883)))");
    fEstV0A -> SetUseAnchor        ( kTRUE   ) ;
    fEstV0A -> SetAnchorPoint      ( 49.0854 ) ;
    fEstV0A -> SetAnchorPercentile ( 90.007  ) ;

    AliMultEstimator *fEstV0C = new AliMultEstimator("V0C", "", "((fAmplitude_V0C)/(1 + ((fEvSel_VtxZ)-1.12997)*(4.04273e-07 - ((fEvSel_VtxZ)-1.12997)* 2.94567e-05 )))");
    fEstV0C -> SetUseAnchor        ( kTRUE   ) ;
    fEstV0C -> SetAnchorPoint      ( 83.7402 ) ;
    fEstV0C -> SetAnchorPercentile ( 90.007  ) ;

    fEstV0M -> SetUseAnchor        ( kTRUE   ) ;
    fEstV0M -> SetAnchorPoint      ( lDefaultV0MAnchor   ) ;
    fEstV0M -> SetAnchorPercentile ( lDefaultV0MPercentile   ) ;

    fEstV0Mplus05 -> SetUseAnchor        ( kTRUE   ) ;
    fEstV0Mplus05 -> SetAnchorPoint      ( lDefaultV0MAnchor   ) ;
    fEstV0Mplus05 -> SetAnchorPercentile ( lDefaultV0MPercentile +0.5  ) ;

    fEstV0Mplus10 -> SetUseAnchor        ( kTRUE   ) ;
    fEstV0Mplus10 -> SetAnchorPoint      ( lDefaultV0MAnchor   ) ;
    fEstV0Mplus10 -> SetAnchorPercentile ( lDefaultV0MPercentile +1.  ) ;

    fEstV0Mminus05 -> SetUseAnchor        ( kTRUE   ) ;
    fEstV0Mminus05 -> SetAnchorPoint      ( lDefaultV0MAnchor   ) ;
    fEstV0Mminus05 -> SetAnchorPercentile ( lDefaultV0MPercentile -0.5  ) ;

    fEstV0Mminus10 -> SetUseAnchor        ( kTRUE   ) ;
    fEstV0Mminus10 -> SetAnchorPoint      ( lDefaultV0MAnchor   ) ;
    fEstV0Mminus10 -> SetAnchorPercentile ( lDefaultV0MPercentile -1.  ) ;
  }

  AliMultEstimator *fEstOnlineV0M = new AliMultEstimator("OnlineV0M", "", "(fAmplitude_OnlineV0A)+(fAmplitude_OnlineV0C)");
  AliMultEstimator *fEstOnlineV0A = new AliMultEstimator("OnlineV0A", "", "(fAmplitude_OnlineV0A)");
  AliMultEstimator *fEstOnlineV0C = new AliMultEstimator("OnlineV0C", "", "(fAmplitude_OnlineV0C)");

  AliMultEstimator *fEstADM = new AliMultEstimator("ADM", "", "(fMultiplicity_ADA)+(fMultiplicity_ADC)");
  AliMultEstimator *fEstADA = new AliMultEstimator("ADA", "", "(fMultiplicity_ADA)");
  AliMultEstimator *fEstADC = new AliMultEstimator("ADC", "", "(fMultiplicity_ADC)");

  //Integer estimators
  AliMultEstimator *fEstnSPDClusters = new AliMultEstimator("SPDClusters", "", "(fnSPDClusters)");
  fEstnSPDClusters->SetIsInteger(kTRUE);
  AliMultEstimator *fEstnSPDTracklets = new AliMultEstimator("SPDTracklets", "", "(fnTracklets)");
  fEstnSPDTracklets->SetIsInteger(kTRUE);
  AliMultEstimator *fEstRefMultEta5 = new AliMultEstimator("RefMult05", "", "(fRefMultEta5)");
  fEstRefMultEta5->SetIsInteger(kTRUE);
  AliMultEstimator *fEstRefMultEta8 = new AliMultEstimator("RefMult08", "", "(fRefMultEta8)");
  fEstRefMultEta8->SetIsInteger(kTRUE);




  //Only do this in run 2, AD didn't exist in Run 1
  //Will also save space in the OADB for old datasets!
  //if( lPeriodName.Contains("LHC15") ){
    //lCalib->GetMultSelection() -> AddEstimator( fEstOnlineV0M );
    //lCalib->GetMultSelection() -> AddEstimator( fEstOnlineV0A );
    //lCalib->GetMultSelection() -> AddEstimator( fEstOnlineV0C );
    //lCalib->GetMultSelection() -> AddEstimator( fEstADM );
    //lCalib->GetMultSelection() -> AddEstimator( fEstADA );
    //lCalib->GetMultSelection() -> AddEstimator( fEstADC );
  //}

  //Universal: Tracking, etc
  //lCalib->GetMultSelection() -> AddEstimator( fEstnSPDClusters  );
  //lCalib->GetMultSelection() -> AddEstimator( fEstnSPDTracklets );
  //lCalib->GetMultSelection() -> AddEstimator( fEstRefMultEta5 );
  //lCalib->GetMultSelection() -> AddEstimator( fEstRefMultEta8 );

  if( lPeriodName.Contains("LHC10h") ) {
    AliMultEstimator *fEstnSPDClustersCorr = new AliMultEstimator("SPDClustersCorr", "", "(fnSPDClusters)/(1 + ((fEvSel_VtxZ)-1.08384)*(0.00131615 - ((fEvSel_VtxZ)-1.08384)*(-0.000659206)))");
    fEstnSPDClustersCorr -> SetUseAnchor        ( kTRUE    ) ;
    fEstnSPDClustersCorr -> SetAnchorPercentile ( 90.0     ) ;
    fEstnSPDClustersCorr -> SetAnchorPoint      ( 48.22    ) ;

    AliMultEstimator *fEstCL0 = new AliMultEstimator("CL0", "", "(fnSPDClusters0)/(1+((fEvSel_VtxZ)-1.08424)*(0.000662816-((fEvSel_VtxZ)-1.08424)*(-0.00154626)))");
    fEstCL0 -> SetUseAnchor        ( kTRUE  ) ;
    fEstCL0 -> SetAnchorPoint      ( lDefaultCL0Anchor    ) ;
    fEstCL0 -> SetAnchorPercentile ( lDefaultCL0Percentile  );

    AliMultEstimator *fEstCL1 = new AliMultEstimator("CL1", "", "(fnSPDClusters1)/(1+((fEvSel_VtxZ)-1.07316)*(-0.000665863-((fEvSel_VtxZ)-1.07316)*(-0.0011758)))");
    fEstCL1 -> SetUseAnchor        ( kTRUE  ) ;
    fEstCL1 -> SetAnchorPoint      ( lDefaultCL1Anchor     ) ;
    fEstCL1 -> SetAnchorPercentile ( lDefaultCL1Percentile ) ;
  }

  if( lPeriodName.Contains("LHC15o") ) {

    AliMultEstimator *fEstnSPDClustersCorr = new AliMultEstimator("SPDClustersCorr", "", "(fnSPDClusters)/(1 + ((fEvSel_VtxZ)-1.14951)*(0.000249854 - ((fEvSel_VtxZ)-1.14951)*(-0.000176228)))");
    fEstnSPDClustersCorr -> SetUseAnchor        ( kTRUE    ) ;
    fEstnSPDClustersCorr -> SetAnchorPercentile ( 90.0    ) ;
    fEstnSPDClustersCorr -> SetAnchorPoint      ( 71.5     ) ;


    //Config central CL0
    AliMultEstimator *fEstCL0 = new AliMultEstimator("CL0", "", "(fnSPDClusters0)/(1+((fEvSel_VtxZ)-1.0)*(0.004+((fEvSel_VtxZ)-1.0)*(-0.001)))");
    fEstCL0 -> SetUseAnchor        ( kTRUE  ) ;
    fEstCL0 -> SetAnchorPoint      ( lDefaultCL0Anchor   ) ;
    fEstCL0 -> SetAnchorPercentile ( lDefaultCL0Percentile  ) ;




    //Config central CL1
    AliMultEstimator *fEstCL1 = new AliMultEstimator("CL1", "", "(fnSPDClusters1)/(1+((fEvSel_VtxZ)-1.0)*(0.0008+((fEvSel_VtxZ)-1.0)*(-0.0009)))");
    fEstCL1 -> SetUseAnchor        ( kTRUE   ) ;
    fEstCL1 -> SetAnchorPoint      ( lDefaultCL1Anchor    ) ;
    fEstCL1 -> SetAnchorPercentile ( lDefaultCL1Percentile  ) ;

    //Plus and Minus
    AliMultEstimator *fEstCL1plus10 = new AliMultEstimator("CL1plus10", "", "(fnSPDClusters1)/(1+((fEvSel_VtxZ)-1.0)*(0.0008+((fEvSel_VtxZ)-1.0)*(-0.0009)))");
    fEstCL1plus10 -> SetUseAnchor        ( kTRUE   ) ;
    fEstCL1plus10 -> SetAnchorPoint      ( lDefaultCL1Anchor    ) ;
    fEstCL1plus10 -> SetAnchorPercentile ( lDefaultCL1Percentile + 1.0 ) ;
    AliMultEstimator *fEstCL1plus05 = new AliMultEstimator("CL1plus05", "", "(fnSPDClusters1)/(1+((fEvSel_VtxZ)-1.0)*(0.0008+((fEvSel_VtxZ)-1.0)*(-0.0009)))");
    fEstCL1plus05 -> SetUseAnchor        ( kTRUE   ) ;
    fEstCL1plus05 -> SetAnchorPoint      ( lDefaultCL1Anchor    ) ;
    fEstCL1plus05 -> SetAnchorPercentile ( lDefaultCL1Percentile + 0.5 ) ;
    AliMultEstimator *fEstCL1minus05 = new AliMultEstimator("CL1minus05", "", "(fnSPDClusters1)/(1+((fEvSel_VtxZ)-1.0)*(0.0008+((fEvSel_VtxZ)-1.0)*(-0.0009)))");
    fEstCL1minus05 -> SetUseAnchor        ( kTRUE   ) ;
    fEstCL1minus05 -> SetAnchorPoint      ( lDefaultCL1Anchor ) ;
    fEstCL1minus05 -> SetAnchorPercentile ( lDefaultCL1Percentile - 0.5 ) ;
    AliMultEstimator *fEstCL1minus10 = new AliMultEstimator("CL1minus10", "", "(fnSPDClusters1)/(1+((fEvSel_VtxZ)-1.0)*(0.0008+((fEvSel_VtxZ)-1.0)*(-0.0009)))");
    fEstCL1minus10 -> SetUseAnchor        ( kTRUE   ) ;
    fEstCL1minus10 -> SetAnchorPoint      ( lDefaultCL1Anchor ) ;
    fEstCL1minus10 -> SetAnchorPercentile ( lDefaultCL1Percentile - 1.0 ) ;


    //////
    AliMultEstimator *fEstV0I = new AliMultEstimator("V0I", "", "((fAmplitude_V0A1)+(fAmplitude_V0A2)+(fAmplitude_V0C1)+(fAmplitude_V0C2))");
    fEstV0I -> SetUseAnchor        ( kTRUE   ) ;
    fEstV0I -> SetAnchorPoint      ( 56.9147 ) ;
    fEstV0I -> SetAnchorPercentile ( 90.0  ) ;


    AliMultEstimator *fEstV0O = new AliMultEstimator("V0O", "", "((fAmplitude_V0A3)+(fAmplitude_V0A4)+(fAmplitude_V0C3)+(fAmplitude_V0C4))");
    fEstV0O -> SetUseAnchor        ( kTRUE   ) ;
    fEstV0O -> SetAnchorPoint      ( 75.8783 ) ;
    fEstV0O -> SetAnchorPercentile ( 90.0  ) ;

  }


  //Univeral: V0
  lMultSelDefault -> AddEstimator( fEstV0M );
  lMultSelDefault -> AddEstimator( fEstCL0 );
  lMultSelDefault -> AddEstimator( fEstCL1 );

  lMultSelDefault -> AddEstimator( fEstV0Mplus05  );
  lMultSelDefault -> AddEstimator( fEstV0Mplus10  );
  lMultSelDefault -> AddEstimator( fEstV0Mminus05 );
  lMultSelDefault -> AddEstimator( fEstV0Mminus10 );

  //lCalib->GetMultSelection() -> AddEstimator( fEstCL0plus05 );
  //lCalib->GetMultSelection() -> AddEstimator( fEstCL0plus10 );
  //lCalib->GetMultSelection() -> AddEstimator( fEstCL0minus05 );
  //lCalib->GetMultSelection() -> AddEstimator( fEstCL0minus10 );

  //lCalib->GetMultSelection() -> AddEstimator( fEstCL1plus05 );
  //lCalib->GetMultSelection() -> AddEstimator( fEstCL1plus10 );
  //lCalib->GetMultSelection() -> AddEstimator( fEstCL1minus05 );
  //lCalib->GetMultSelection() -> AddEstimator( fEstCL1minus10 );

  //lCalib->GetMultSelection() -> AddEstimator( fEstV0A );
  //lCalib->GetMultSelection() -> AddEstimator( fEstV0C );
  //lCalib->GetMultSelection() -> AddEstimator( fEstV0I );
  //lCalib->GetMultSelection() -> AddEstimator( fEstV0O );

  lMultSelDefault -> AddEstimator( fEstnSPDClustersCorr );

  //Needed as basis of MC calibrator -- keep this here, please!
  AliMultEstimator *fEstnSPDTracklets = new AliMultEstimator("SPDTracklets", "", "(fnTracklets)");
  fEstnSPDTracklets->SetIsInteger(kTRUE);
  lMultSelDefault -> AddEstimator( fEstnSPDTracklets );

  cout<<"Inspect default object:"<<endl;
  lMultSelDefault->PrintInfo();

  //Run-wise changes

  //Run List from Processing - 23rd December 2015
  Long_t lRuns[] = {
    244824 , 244827 , 244889 , 244911 , 244917 , 244918 , 244947 , 244972 , 244975 , 244980 , 244982 , 244983 , 245061 , 245064 , 245066 , 245068 , 245096 , 245099 , 245100 , 245101 , 245103 , 245145 , 245146 , 245148 , 245151 , 245152 , 245230 , 245231 , 245232 , 245233 , 245253 , 245256 , 245259 , 245341 , 245343 , 245345 , 245346 , 245347 , 245349 , 245353 , 245396 , 245397 , 245401 , 245407 , 245409 , 245410 , 245411 , 245439 , 245441 , 245446 , 245450 , 245452 , 245453 , 245454 , 245496 , 245497 , 245501 , 245504 , 245505 , 245507 , 245535 , 245540 , 245542 , 245543 , 245544 , 245545 , 245554 , 245683 , 245692 , 245700 , 245702 , 245705 , 245729 , 245731 , 245738 , 245752 , 245759 , 245766 , 245775 , 245785 , 245793 , 245829 , 245831 , 245833 , 245923 , 245949 , 245952 , 245954 , 245963 , 245996 , 246001 , 246003 , 246012 , 246036 , 246037 , 246042 , 246048 , 246049 , 246052 , 246053 , 246087 , 246089 , 246113 , 246115 , 246148 , 246151 , 246152 , 246153 , 246178 , 246180 , 246181 , 246182 , 246185 , 246187 , 246217 , 246220 , 246222 , 246225 , 246271 , 246272 , 246275 , 246276 , 246390 , 246391 , 246392 , 246424 , 246428 , 246431 , 246433 , 246434 , 246469 , 246487 , 246488 , 246493 , 246495 , 246540 , 246543 , 246553 , 246567 , 246568 , 246575 , 246583 , 246639 , 246648 , 246671 , 246675 , 246676 , 246750 , 246751 , 246755 , 246757 , 246758 , 246759 , 246760 , 246763 , 246765 , 246766 , 246804 , 246805 , 246806 , 246807 , 246808 , 246809 , 246810 , 246844 , 246845 , 246846 , 246847 , 246851 , 246855 , 246858 , 246859 , 246864 , 246865 , 246867 , 246870 , 246871 , 246928 , 246930 , 246937 , 246942 , 246945 , 246948 , 246949 , 246980 , 246982 , 246984 , 246989 , 246991 , 246994
  };

  Float_t lV0MAnchors[] = {
    132.5,  132.5,  132.5,  132.5,  132.5,  132.5,  132.5,  132.5, 132.5,  132.5,  132.5,  132.5,  132.5,  132.5,  132.5,  132.5, 132.5,  132.5,  132.5,  132.5,  132.5,  132.5, 132.414 , 132.378 , 132.342 , 132.306 , 132.269 , 132.233 , 132.197 , 132.161 , 132.125 , 132.088 , 132.052 , 132.016 , 131.98 , 131.944 , 131.907 , 131.871 , 131.835 , 131.799 , 131.763 , 131.726 , 131.69 , 131.654 , 131.618 , 131.582 , 131.545 , 131.509 , 131.473 , 131.437 , 131.401 , 131.364 , 131.328 , 131.292 , 131.256 , 131.22 , 131.183 , 131.147 , 131.111 , 131.075 , 131.038 , 131.002 , 130.966 , 130.93 , 130.894 , 130.857 , 130.821 , 130.785 , 130.749 , 130.713 , 130.676 , 130.64 , 130.604 , 130.568 , 130.532 , 130.495 , 130.459 , 130.423 , 130.387 , 130.351 , 130.314 , 130.278 , 130.242 , 130.206 , 130.17 , 130.133 , 130.097 , 130.061 , 130.025 , 129.989 , 129.952 , 129.916 , 129.88 , 129.844 , 129.807 , 129.771 , 129.735 , 129.699 , 129.663 , 129.626 , 129.59 , 129.554 , 129.518 , 129.482 , 129.445 , 129.409 , 129.373 , 129.337 , 129.301 , 129.264 , 129.228 , 129.192 , 129.156 , 129.12 , 129.083 , 129.047 , 129.011 , 128.975 , 128.939 , 128.902 , 128.866 , 128.83 , 128.794 , 128.758 , 128.721 , 128.685 , 128.649 , 128.613 , 128.576 , 128.54 , 128.504 , 128.468 , 128.432 , 128.395 , 128.359 , 128.323 , 128.287 , 128.251 , 128.214 , 128.178 , 128.142 , 128.106 , 128.07 , 128.033 , 127.997 , 127.961 , 127.925 , 127.889 , 127.852 , 127.816 , 127.78 , 127.744 , 127.708 , 127.671 , 127.635 , 127.599 , 127.563 , 127.527 , 127.49 , 127.454 , 127.418 , 127.382 , 127.345 , 127.309 , 127.273 , 127.237 , 127.201 , 127.164 , 127.128 , 127.092 , 127.056 , 127.02 , 126.983 , 126.947 , 126.911 , 126.875 , 126.839 , 126.802 , 126.766 , 126.73 , 126.694 , 126.658 , 126.621 , 126.585 , 126.549 , 126.513 , 126.477 , 126.44 , 126.404 , 126.368
  };


  const Long_t lNRuns = sizeof(lRuns) / sizeof(Long_t);
  AliMultSelection *lMultSelArray[500]; //We need one AliMultSelection for each run above!

  //Sweep array
  for(Long_t iRun = 0; iRun<lNRuns; iRun++) {
    //cout<<"Instantiating AliMultSelection for run "<<lRuns[iRun]<<"... "<<endl;
    lMultSelArray[iRun] = new AliMultSelection(lMultSelDefault);

    //Adjust V0M Anchors, please
    lMultSelArray[iRun]->GetEstimator("V0M")->SetAnchorPoint( lV0MAnchors[iRun] );
    lMultSelArray[iRun]->GetEstimator("V0Mplus05")->SetAnchorPoint( lV0MAnchors[iRun] );
    lMultSelArray[iRun]->GetEstimator("V0Mplus10")->SetAnchorPoint( lV0MAnchors[iRun] );
    lMultSelArray[iRun]->GetEstimator("V0Mminus05")->SetAnchorPoint( lV0MAnchors[iRun] );
    lMultSelArray[iRun]->GetEstimator("V0Mminus10")->SetAnchorPoint( lV0MAnchors[iRun] );

    //Debug: check if changed!
    //lMultSelArray[iRun] -> PrintInfo();

    //cout<<"Adding run range..."<<endl;
    //Add run to run ranges (still RbyR for now)
    lCalib->AddRunRange( lRuns[iRun], lRuns[iRun], lMultSelArray[iRun] );
  }
  lCalib->SetRunToUseAsDefault(246434);

  //============================================================
  // --- Definition of Input/Output ---
  //============================================================
  //lCalib -> SetInputFile  ( "~/Dropbox/MultSelCalib/LHC15o/MergedLHC15o.root");
  //lCalib -> SetInputFile  ( "/hera/alice/alberica/centrality/Trees/Singles/LHC15o/muon_calo_pass1_1301/files/AnalysisResults_244918.root");
  lCalib -> SetInputFile  ( "~/Dropbox/MultSelCalib/LHC15o/MergedLHC15o.root");
  //lCalib -> SetInputFile  ( "~/Dropbox/MultSelCalib/LHC15o/files/AnalysisResults_245064.root");

  lCalib -> SetBufferFile ( "buffer.root" );
  lCalib -> SetOutputFile ( "OADB-LHC15o-SuperCorrection.root" );
  lCalib -> Calibrate     ();

}
