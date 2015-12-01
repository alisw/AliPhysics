

CalibratePeriodPbPb(TString lPeriodName = "LHC10h"){
    
    //Load ALICE stuff
    TString gLibs[] =   {"STEER",
        "ANALYSIS", "ANALYSISalice", "ANALYSIScalib","OADB"};
    TString thislib = "lib";
    for(Int_t ilib = 0; ilib<5; ilib++){
        thislib="lib";
        thislib.Append(gLibs[ilib].Data());
        cout<<"Will load "<<thislib.Data()<<endl;
        gSystem->Load(thislib.Data());
    }
    gSystem->SetIncludePath("-I$ROOTSYS/include  -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
    
    cout<<"Alive! "<<endl;
    
    //All fine, let's try the calibrator
    AliMultSelectionCalibrator *lCalib = new AliMultSelectionCalibrator("lCalib");
    
    lCalib->SetRunToUseAsDefault( 244918 );
    
    //============================================================
    // --- Definition of Boundaries ---
    //============================================================
    
    //Set Adaptive Percentile Boundaries, adjust if finer selection desired
    Double_t lDesiredBoundaries[1000];
    Long_t   lNDesiredBoundaries=0;
    lDesiredBoundaries[0] = 100;
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
    lDesiredBoundaries[lNDesiredBoundaries] = 0;
    
    lCalib->SetBoundaries( lNDesiredBoundaries, lDesiredBoundaries );
    cout<<"Boundaries set. Will attempt calibration now... "<<endl;
    
    if ( lPeriodName.Contains("LHC10h") ){
        cout<<"Setting event selection criteria for Pb-Pb..."<<endl;
        lCalib->GetEventCuts()->SetVzCut(10.0);
        lCalib->GetEventCuts()->SetTriggerCut                (kTRUE );
        lCalib->GetEventCuts()->SetINELgtZEROCut             (kFALSE);
        lCalib->GetEventCuts()->SetTrackletsVsClustersCut    (kFALSE);
        lCalib->GetEventCuts()->SetRejectPileupInMultBinsCut (kFALSE);
        lCalib->GetEventCuts()->SetVertexConsistencyCut      (kFALSE);
        lCalib->GetEventCuts()->SetNonZeroNContribs          (kTRUE);
    }
    
    if ( lPeriodName.Contains("LHC15m") || lPeriodName.Contains("LHC15o") ){
        cout<<"Setting event selection criteria for Pb-Pb..."<<endl;
        lCalib->GetEventCuts()->SetVzCut(10.0);
        lCalib->GetEventCuts()->SetTriggerCut                (kTRUE );
        lCalib->GetEventCuts()->SetINELgtZEROCut             (kFALSE);
        lCalib->GetEventCuts()->SetTrackletsVsClustersCut    (kFALSE);
        lCalib->GetEventCuts()->SetRejectPileupInMultBinsCut (kFALSE);
        lCalib->GetEventCuts()->SetVertexConsistencyCut      (kFALSE);
        lCalib->GetEventCuts()->SetNonZeroNContribs          (kTRUE);
    }
    
    //============================================================
    // --- Definition of Input Variables ---
    //============================================================
    
    lCalib->SetupStandardInput();
    
    //Changes in new version: create AliMultSelection here
    AliMultSelection *lMultSel = new AliMultSelection();
    lCalib->SetMultSelection(lMultSel);
    
    //============================================================
    // --- Definition of Estimators ---
    //============================================================
    
    
    Double_t lDefaultV0MAnchor     = 0;
    Double_t lDefaultV0MPercentile = 0;
    
    Double_t lDefaultCL0Anchor     = 0;
    Double_t lDefaultCL0Percentile = 0;
    
    Double_t lDefaultCL1Anchor     = 0;
    Double_t lDefaultCL1Percentile = 0;
    
    if ( lPeriodName.Contains("LHC15o") ){
        lDefaultV0MAnchor     = 133.5;
        lDefaultV0MPercentile = 90.007;
        lDefaultCL0Anchor     = 33.5;
        lDefaultCL0Percentile = 90.64;
        lDefaultCL1Anchor     = 30.5;
        lDefaultCL1Percentile = 90.485;
    }
    if ( lPeriodName.Contains("LHC15m") ){
        lDefaultV0MAnchor     = 115.0;
        lDefaultV0MPercentile = 87.5;
        lDefaultCL0Anchor     = 39.5;
        lDefaultCL0Percentile = 88.9;
        lDefaultCL1Anchor     = 40.5;
        lDefaultCL1Percentile = 88.1;
    }
    
    AliMultEstimator *fEstV0M = new AliMultEstimator("V0M", "", "(fAmplitude_V0A)+(fAmplitude_V0C)");
    //Config central V0M
    fEstV0M -> SetUseAnchor        ( kTRUE   ) ;
    fEstV0M -> SetAnchorPoint      ( lDefaultV0MAnchor   ) ;
    fEstV0M -> SetAnchorPercentile ( lDefaultV0MPercentile  ) ;
    
    //Plus and Minus
    AliMultEstimator *fEstV0Mplus10 = new AliMultEstimator("V0Mplus10", "", "(fAmplitude_V0A)+(fAmplitude_V0C)");
    fEstV0Mplus10 -> SetUseAnchor        ( kTRUE   ) ;
    fEstV0Mplus10 -> SetAnchorPoint      ( lDefaultV0MAnchor   ) ;
    fEstV0Mplus10 -> SetAnchorPercentile ( lDefaultV0MPercentile + 1.0    ) ;
    AliMultEstimator *fEstV0Mplus05 = new AliMultEstimator("V0Mplus05", "", "(fAmplitude_V0A)+(fAmplitude_V0C)");
    fEstV0Mplus05 -> SetUseAnchor        ( kTRUE   ) ;
    fEstV0Mplus05 -> SetAnchorPoint      ( lDefaultV0MAnchor   ) ;
    fEstV0Mplus05 -> SetAnchorPercentile ( lDefaultV0MPercentile + 0.5    ) ;
    AliMultEstimator *fEstV0Mminus05 = new AliMultEstimator("V0Mminus05", "", "(fAmplitude_V0A)+(fAmplitude_V0C)");
    fEstV0Mminus05 -> SetUseAnchor        ( kTRUE   ) ;
    fEstV0Mminus05 -> SetAnchorPoint      ( lDefaultV0MAnchor   ) ;
    fEstV0Mminus05 -> SetAnchorPercentile ( lDefaultV0MPercentile - 0.5    ) ;
    AliMultEstimator *fEstV0Mminus10 = new AliMultEstimator("V0Mminus10", "", "(fAmplitude_V0A)+(fAmplitude_V0C)");
    fEstV0Mminus10 -> SetUseAnchor        ( kTRUE   ) ;
    fEstV0Mminus10 -> SetAnchorPoint      ( lDefaultV0MAnchor   ) ;
    fEstV0Mminus10 -> SetAnchorPercentile ( lDefaultV0MPercentile - 1.0    ) ;
    
    //Config central CL0
    AliMultEstimator *fEstCL0 = new AliMultEstimator("CL0", "", "(fnSPDClusters0)/(1+((fEvSel_VtxZ)-1.13374)*(0.00258567-((fEvSel_VtxZ)-1.13374)*(-0.00131249)))");
    fEstCL0 -> SetUseAnchor        ( kTRUE  ) ;
    fEstCL0 -> SetAnchorPoint      ( lDefaultCL0Anchor   ) ;
    fEstCL0 -> SetAnchorPercentile ( lDefaultCL0Percentile  ) ;
    
    //Plus and Minus
    AliMultEstimator *fEstCL0plus10 = new AliMultEstimator("CL0plus10", "", "(fnSPDClusters0)/(1+((fEvSel_VtxZ)-1.13374)*(0.00258567-((fEvSel_VtxZ)-1.13374)*(-0.00131249)))");
    fEstCL0plus10 -> SetUseAnchor        ( kTRUE  ) ;
    fEstCL0plus10 -> SetAnchorPoint      ( lDefaultCL0Anchor   ) ;
    fEstCL0plus10 -> SetAnchorPercentile ( lDefaultCL0Percentile + 1.0 ) ;
    AliMultEstimator *fEstCL0plus05 = new AliMultEstimator("CL0plus05", "", "(fnSPDClusters0)/(1+((fEvSel_VtxZ)-1.13374)*(0.00258567-((fEvSel_VtxZ)-1.13374)*(-0.00131249)))");
    fEstCL0plus05 -> SetUseAnchor        ( kTRUE  ) ;
    fEstCL0plus05 -> SetAnchorPoint      ( lDefaultCL0Anchor   ) ;
    fEstCL0plus05 -> SetAnchorPercentile ( lDefaultCL0Percentile + 0.5 ) ;
    AliMultEstimator *fEstCL0minus05 = new AliMultEstimator("CL0minus05", "", "(fnSPDClusters0)/(1+((fEvSel_VtxZ)-1.13374)*(0.00258567-((fEvSel_VtxZ)-1.13374)*(-0.00131249)))");
    fEstCL0minus05 -> SetUseAnchor        ( kTRUE  ) ;
    fEstCL0minus05 -> SetAnchorPoint      ( lDefaultCL0Anchor   ) ;
    fEstCL0minus05 -> SetAnchorPercentile ( lDefaultCL0Percentile - 0.5 ) ;
    AliMultEstimator *fEstCL0minus10 = new AliMultEstimator("CL0minus10", "", "(fnSPDClusters0)/(1+((fEvSel_VtxZ)-1.13374)*(0.00258567-((fEvSel_VtxZ)-1.13374)*(-0.00131249)))");
    fEstCL0minus10 -> SetUseAnchor        ( kTRUE  ) ;
    fEstCL0minus10 -> SetAnchorPoint      ( lDefaultCL0Anchor   ) ;
    fEstCL0minus10 -> SetAnchorPercentile ( lDefaultCL0Percentile - 1.0 ) ;
    
    //Config central CL1
    AliMultEstimator *fEstCL1 = new AliMultEstimator("CL1", "", "(fnSPDClusters1)/(1+((fEvSel_VtxZ)-1.12)*(-0.000346139-((fEvSel_VtxZ)-1.12)*(-0.000849613)))");
    fEstCL1 -> SetUseAnchor        ( kTRUE   ) ;
    fEstCL1 -> SetAnchorPoint      ( lDefaultCL1Anchor    ) ;
    fEstCL1 -> SetAnchorPercentile ( lDefaultCL1Percentile  ) ;
    
    //Plus and Minus
    AliMultEstimator *fEstCL1plus10 = new AliMultEstimator("CL1plus10", "", "(fnSPDClusters1)/(1+((fEvSel_VtxZ)-1.12)*(-0.000346139-((fEvSel_VtxZ)-1.12)*(-0.000849613)))");
    fEstCL1plus10 -> SetUseAnchor        ( kTRUE   ) ;
    fEstCL1plus10 -> SetAnchorPoint      ( lDefaultCL1Anchor    ) ;
    fEstCL1plus10 -> SetAnchorPercentile ( lDefaultCL1Percentile + 1.0 ) ;
    AliMultEstimator *fEstCL1plus05 = new AliMultEstimator("CL1plus05", "", "(fnSPDClusters1)/(1+((fEvSel_VtxZ)-1.12)*(-0.000346139-((fEvSel_VtxZ)-1.12)*(-0.000849613)))");
    fEstCL1plus05 -> SetUseAnchor        ( kTRUE   ) ;
    fEstCL1plus05 -> SetAnchorPoint      ( lDefaultCL1Anchor    ) ;
    fEstCL1plus05 -> SetAnchorPercentile ( lDefaultCL1Percentile + 0.5 ) ;
    AliMultEstimator *fEstCL1minus05 = new AliMultEstimator("CL1minus05", "", "(fnSPDClusters1)/(1+((fEvSel_VtxZ)-1.12)*(-0.000346139-((fEvSel_VtxZ)-1.12)*(-0.000849613)))");
    fEstCL1minus05 -> SetUseAnchor        ( kTRUE   ) ;
    fEstCL1minus05 -> SetAnchorPoint      ( lDefaultCL1Anchor    ) ;
    fEstCL1minus05 -> SetAnchorPercentile ( lDefaultCL1Percentile - 0.5 ) ;
    AliMultEstimator *fEstCL1minus10 = new AliMultEstimator("CL1minus10", "", "(fnSPDClusters1)/(1+((fEvSel_VtxZ)-1.12)*(-0.000346139-((fEvSel_VtxZ)-1.12)*(-0.000849613)))");
    fEstCL1minus10 -> SetUseAnchor        ( kTRUE   ) ;
    fEstCL1minus10 -> SetAnchorPoint      ( lDefaultCL1Anchor    ) ;
    fEstCL1minus10 -> SetAnchorPercentile ( lDefaultCL1Percentile - 1.0 ) ;
    
    //Needed as basis of MC calibrator
    AliMultEstimator *fEstnSPDTracklets = new AliMultEstimator("SPDTracklets", "", "(fnTracklets)");
    fEstnSPDTracklets->SetIsInteger(kTRUE);
    
    //Univeral: V0
    lCalib->GetMultSelection() -> AddEstimator( fEstV0M );
    lCalib->GetMultSelection() -> AddEstimator( fEstCL0 );
    lCalib->GetMultSelection() -> AddEstimator( fEstCL1 );
    
    //Add Plus and Minus estimators
    lCalib->GetMultSelection() -> AddEstimator( fEstV0Mplus10 );
    lCalib->GetMultSelection() -> AddEstimator( fEstV0Mplus05 );
    lCalib->GetMultSelection() -> AddEstimator( fEstV0Mminus05 );
    lCalib->GetMultSelection() -> AddEstimator( fEstV0Mminus10 );
    lCalib->GetMultSelection() -> AddEstimator( fEstCL0plus10 );
    lCalib->GetMultSelection() -> AddEstimator( fEstCL0plus05 );
    lCalib->GetMultSelection() -> AddEstimator( fEstCL0minus05 );
    lCalib->GetMultSelection() -> AddEstimator( fEstCL0minus10 );
    lCalib->GetMultSelection() -> AddEstimator( fEstCL1plus10 );
    lCalib->GetMultSelection() -> AddEstimator( fEstCL1plus05 );
    lCalib->GetMultSelection() -> AddEstimator( fEstCL1minus05 );
    lCalib->GetMultSelection() -> AddEstimator( fEstCL1minus10 );
    
    lCalib->GetMultSelection() -> AddEstimator( fEstnSPDTracklets );

    
    //============================================================
    // --- Definition of Input/Output ---
    //============================================================
    
    if( !lPeriodName.Contains("test") ){
        //Per Period calibration: standard locations...
        lCalib -> SetInputFile  ( Form("~/work/calibs/Merged%s.root",lPeriodName.Data() ) );
        lCalib -> SetBufferFile ( Form("~/work/fast/buffer-%s.root", lPeriodName.Data() ) );
        
        //Local running please
        lCalib -> SetInputFile  ( Form("~/Dropbox/MultSelCalib/%s/Merged%s.root",lPeriodName.Data(), lPeriodName.Data() ) );
        lCalib -> SetBufferFile ( "buffer.root" );
        
        lCalib -> SetOutputFile ( Form("OADB-%s.root", lPeriodName.Data() ) );
        lCalib -> Calibrate     ();
    }else{
        lCalib -> SetInputFile  ( "../MultSelCalib/LHC10h/files/AnalysisResults_137161.root");
        lCalib -> SetBufferFile ( "buffer-test.root" );
        lCalib -> SetOutputFile ( "OADB-testing.root" );
        lCalib -> Calibrate     ();
    }
}
