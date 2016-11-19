

CalibratePeriodPP(TString lPeriodName = "LHC15f") {

    //Load ALICE stuff
    TString gLibs[] =    {"STEER",
                            "ANALYSIS", "ANALYSISalice", "ANALYSIScalib","OADB"
                           };
    TString thislib = "lib";
    for(Int_t ilib = 0; ilib<5; ilib++) {
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

    if ( lPeriodName.Contains("LHC10h") ) {
        cout<<"Setting event selection criteria for Pb-Pb..."<<endl;
        lCalib->GetEventCuts()->SetVzCut(10.0);
        lCalib->GetEventCuts()->SetTriggerCut                (kTRUE );
        lCalib->GetEventCuts()->SetINELgtZEROCut             (kFALSE);
        lCalib->GetEventCuts()->SetTrackletsVsClustersCut    (kFALSE);
        lCalib->GetEventCuts()->SetRejectPileupInMultBinsCut (kFALSE);
        lCalib->GetEventCuts()->SetVertexConsistencyCut      (kFALSE);
        lCalib->GetEventCuts()->SetNonZeroNContribs          (kTRUE);
    }

    if ( lPeriodName.Contains("LHC15f") ) {
        cout<<"Extra cleanup in LHC15f"<<endl;
        lCalib->GetEventCuts()->SetIsNotAsymmetricInVZERO (kTRUE );
        lCalib->GetEventCuts()->SetIsNotIncompleteDAQ     (kTRUE );
    }

    
    if ( lPeriodName.Contains("LHC15m") || lPeriodName.Contains("LHC15o") ) {
        cout<<"Setting event selection criteria for Pb-Pb..."<<endl;
        lCalib->GetEventCuts()->SetVzCut(10.0);
        lCalib->GetEventCuts()->SetTriggerCut                (kTRUE );
        lCalib->GetEventCuts()->SetINELgtZEROCut             (kFALSE);
        lCalib->GetEventCuts()->SetTrackletsVsClustersCut    (kFALSE);
        lCalib->GetEventCuts()->SetRejectPileupInMultBinsCut (kFALSE);
        lCalib->GetEventCuts()->SetVertexConsistencyCut      (kFALSE);
        lCalib->GetEventCuts()->SetNonZeroNContribs          (kTRUE);
    }

    //Additional selections for pp: incompleteDAQ and asymmetric vzero 



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

    if ( lPeriodName.Contains("LHC15o") ) {
        lDefaultV0MAnchor     = 133.5;
        lDefaultV0MPercentile = 90.007;
        lDefaultCL0Anchor     = 33.5;
        lDefaultCL0Percentile = 90.64;
        lDefaultCL1Anchor     = 30.5;
        lDefaultCL1Percentile = 90.485;
    }
    if ( lPeriodName.Contains("LHC15m") ) {
        lDefaultV0MAnchor     = 115.0;
        lDefaultV0MPercentile = 87.5;
        lDefaultCL0Anchor     = 39.5;
        lDefaultCL0Percentile = 88.9;
        lDefaultCL1Anchor     = 40.5;
        lDefaultCL1Percentile = 88.1;
    }

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
    AliMultEstimator *fEstCL0 = new AliMultEstimator("CL0", "", "(fnSPDClusters0)");
    fEstCL0->SetIsInteger(kTRUE);
    AliMultEstimator *fEstCL1 = new AliMultEstimator("CL1", "", "(fnSPDClusters1)");
    fEstCL1->SetIsInteger(kTRUE);
    
    //Tracklet-based
    AliMultEstimator *fEstnSPDTracklets = new AliMultEstimator("SPDTracklets", "", "(fnTracklets)");
    fEstnSPDTracklets->SetIsInteger(kTRUE);
    AliMultEstimator *fEstnSPDTracklets08 = new AliMultEstimator("SPDTracklets08", "", "(fnTracklets08)");
    fEstnSPDTracklets08->SetIsInteger(kTRUE);
    AliMultEstimator *fEstnSPDTracklets08to15 = new AliMultEstimator("SPD08to15", "", "(fnTracklets08)<-0.5?(fnTracklets08):(fnTracklets15)-(fnTracklets08)");
    fEstnSPDTracklets08to15->SetIsInteger(kTRUE);
    
    //Ref-mult-based
    AliMultEstimator *fEstRefMultEta5 = new AliMultEstimator("RefMult05", "", "(fRefMultEta5)");
    fEstRefMultEta5->SetIsInteger(kTRUE);
    AliMultEstimator *fEstRefMultEta8 = new AliMultEstimator("RefMult08", "", "(fRefMultEta8)");
    fEstRefMultEta8->SetIsInteger(kTRUE);

    //ZDC-based estimators
    /* From mail exchange with Cvetan
     -> ZNApp
    "-fZnaFired * fZnaTower + !fZnaFired * 1e6"

    -> ZNCpp
    "-fZncFired * fZncTower + !fZncFired * 1e6"

    -> ZNACpp
    "-0.89 * fZnaFired * fZnaTower - fZncFired * fZncTower + !fZnaFired * !fZncFired * 1e6"
    */

    AliMultEstimator *fEstZNApp = new AliMultEstimator("ZNApp","", "-(fZnaFired) * (fZnaTower) + !(fZnaFired) * 1e6");
    AliMultEstimator *fEstZNCpp = new AliMultEstimator("ZNCpp","", "-(fZncFired) * (fZncTower) + !(fZncFired) * 1e6");
    AliMultEstimator *fEstZNACpp = new AliMultEstimator("ZNACpp","", "-0.89 * (fZnaFired) * (fZnaTower) - (fZncFired) * (fZncTower) + !(fZnaFired) * !(fZncFired) * 1e6");

    //Universal: V0
    lCalib->GetMultSelection() -> AddEstimator( fEstV0M );
    lCalib->GetMultSelection() -> AddEstimator( fEstV0A );
    lCalib->GetMultSelection() -> AddEstimator( fEstV0C );
    //Only do this in run 2, AD didn't exist in Run 1
    //Will also save space in the OADB for old datasets!
    if( lPeriodName.Contains("LHC15") ) {
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
    lCalib->GetMultSelection() -> AddEstimator( fEstnSPDTracklets08 );
    lCalib->GetMultSelection() -> AddEstimator( fEstnSPDTracklets08to15 );
    lCalib->GetMultSelection() -> AddEstimator( fEstCL0 );
    lCalib->GetMultSelection() -> AddEstimator( fEstCL1 );
    lCalib->GetMultSelection() -> AddEstimator( fEstRefMultEta5 );
    lCalib->GetMultSelection() -> AddEstimator( fEstRefMultEta8 );

    lCalib->GetMultSelection() -> AddEstimator( fEstZNApp );
    lCalib->GetMultSelection() -> AddEstimator( fEstZNCpp );
    lCalib->GetMultSelection() -> AddEstimator( fEstZNACpp );


    //============================================================
    // --- Definition of Input/Output ---
    //============================================================

    if( !lPeriodName.Contains("test") ) {
        //Per Period calibration: standard locations...
        lCalib -> SetInputFile  ( Form("~/work/calibs/Merged%s.root",lPeriodName.Data() ) );
        lCalib -> SetBufferFile ( Form("~/work/fast/buffer-%s.root", lPeriodName.Data() ) );

        //Local running please
        lCalib -> SetInputFile  ( Form("~/Dropbox/MultSelCalib/%s/Merged%s.root",lPeriodName.Data(), lPeriodName.Data() ) );
        lCalib -> SetBufferFile ( "buffer.root" );

        lCalib -> SetOutputFile ( Form("OADB-%s.root", lPeriodName.Data() ) );
        lCalib -> Calibrate     ();
    } else {
        lCalib -> SetInputFile  ( "../MultSelCalib/LHC10h/files/AnalysisResults_137161.root");
        lCalib -> SetBufferFile ( "buffer-test.root" );
        lCalib -> SetOutputFile ( "OADB-testing.root" );
        lCalib -> Calibrate     ();
    }
}
