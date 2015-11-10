

CalibratePeriod(TString lPeriodName = "LHC10h"){
    
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
        lCalib->GetEventCuts()->SetVzCut(1e+6);
        lCalib->GetEventCuts()->SetTriggerCut                (kTRUE );
        lCalib->GetEventCuts()->SetINELgtZEROCut             (kFALSE);
        lCalib->GetEventCuts()->SetTrackletsVsClustersCut    (kTRUE );
        lCalib->GetEventCuts()->SetRejectPileupInMultBinsCut (kFALSE);
        lCalib->GetEventCuts()->SetVertexConsistencyCut      (kFALSE);
    }

    //============================================================
    // --- Definition of Input Variables ---
    //============================================================
    
    lCalib->SetupStandardInput();
    
    //============================================================
    // --- Definition of Estimators ---
    //============================================================
    
    AliMultEstimator *fEstV0M = new AliMultEstimator("V0M", "", "(fAmplitude_V0A)+(fAmplitude_V0C)");
    AliMultEstimator *fEstV0A = new AliMultEstimator("V0A", "", "(fAmplitude_V0A)");
    AliMultEstimator *fEstV0C = new AliMultEstimator("V0C", "", "(fAmplitude_V0C)");
    
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
    
    //Univeral: V0
    lCalib->GetMultSelection() -> AddEstimator( fEstV0M );
    lCalib->GetMultSelection() -> AddEstimator( fEstV0A );
    lCalib->GetMultSelection() -> AddEstimator( fEstV0C );
    
    //Testing!
    AliMultEstimator *fEstV0MAnchored = new AliMultEstimator("V0MAnchored", "", "(fAmplitude_V0A)+(fAmplitude_V0C)");
    fEstV0MAnchored -> SetUseAnchor        ( kTRUE   ) ;
    fEstV0MAnchored -> SetAnchorPoint      ( 81.0    ) ;
    fEstV0MAnchored -> SetAnchorPercentile ( 90.0    ) ;
    lCalib->GetMultSelection() -> AddEstimator( fEstV0MAnchored );
    
    //Anchoring test
    
    //Only do this in run 2, AD didn't exist in Run 1
    //Will also save space in the OADB for old datasets!
    if( lPeriodName.Contains("LHC15") ){
        lCalib->GetMultSelection() -> AddEstimator( fEstOnlineV0M );
        lCalib->GetMultSelection() -> AddEstimator( fEstOnlineV0A );
        lCalib->GetMultSelection() -> AddEstimator( fEstOnlineV0C );
        lCalib->GetMultSelection() -> AddEstimator( fEstADM );
        lCalib->GetMultSelection() -> AddEstimator( fEstADA );
        lCalib->GetMultSelection() -> AddEstimator( fEstADC );
    }
    
    //Universal: Tracking, etc
    lCalib->GetMultSelection() -> AddEstimator( fEstnSPDClusters  );
    lCalib->GetMultSelection() -> AddEstimator( fEstnSPDTracklets );
    lCalib->GetMultSelection() -> AddEstimator( fEstRefMultEta5 );
    lCalib->GetMultSelection() -> AddEstimator( fEstRefMultEta8 );
    
    //============================================================
    // --- Definition of Input/Output ---
    //============================================================
    
    if( !lPeriodName.Contains("test") ){
        //Per Period calibration: standard locations...
        lCalib -> SetInputFile  ( Form("/home/daviddc/Dropbox/MultSelCalib/%s/Merged%s.root",lPeriodName.Data(), lPeriodName.Data() ) );
        lCalib -> SetBufferFile ( Form("/home/daviddc/work/fast/buffer-%s.root", lPeriodName.Data() ) );
        lCalib -> SetOutputFile ( Form("OADB-%s.root", lPeriodName.Data() ) );
        lCalib -> Calibrate     ();
    }else{
        lCalib -> SetInputFile  ( "../MultSelCalib/LHC10h/files/AnalysisResults_137161.root");
        lCalib -> SetBufferFile ( "buffer-test.root" );
        lCalib -> SetOutputFile ( "OADB-testing.root" );
        lCalib -> Calibrate     ();
    }
}
