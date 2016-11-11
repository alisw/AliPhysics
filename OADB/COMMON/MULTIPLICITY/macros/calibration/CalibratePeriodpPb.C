

CalibratePeriodpPb(TString lPeriodName = "LHC16q") {

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
    for( Int_t ib = 1; ib < 91; ib++) {
        lNDesiredBoundaries++;
        lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] - 1.0;
    }
        for( Int_t ib = 1; ib < 101; ib++) {
        lNDesiredBoundaries++;
        lDesiredBoundaries[lNDesiredBoundaries] = lDesiredBoundaries[lNDesiredBoundaries-1] - 0.1;
    }

    lNDesiredBoundaries++;
    lDesiredBoundaries[lNDesiredBoundaries] = 0;
    
    //cout<< "Dump for debug: "<<endl;
    //for( Int_t ib=0;ib<101; ib++)
    //    cout<<"Boundary #"<<ib<<" at "<<lDesiredBoundaries[ib]<<endl;
    
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

    if ( lPeriodName.Contains("LHC15m") || lPeriodName.Contains("LHC15o")  ) {
        cout<<"Setting event selection criteria for Pb-Pb..."<<endl;
        lCalib->GetEventCuts()->SetVzCut(10.0);
        lCalib->GetEventCuts()->SetTriggerCut                (kTRUE );
        lCalib->GetEventCuts()->SetINELgtZEROCut             (kFALSE);
        lCalib->GetEventCuts()->SetTrackletsVsClustersCut    (kFALSE);
        lCalib->GetEventCuts()->SetRejectPileupInMultBinsCut (kFALSE);
        lCalib->GetEventCuts()->SetVertexConsistencyCut      (kFALSE);
        lCalib->GetEventCuts()->SetNonZeroNContribs          (kTRUE);
    }

        if ( lPeriodName.Contains("LHC16q")   ) {
        cout<<"Setting event selection criteria for p-Pb..."<<endl;
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
    lCalib->SetRunToUseAsDefault( 265309 );
    //============================================================
    // --- Definition of Estimators ---
    //============================================================


    Double_t lDefaultV0MAnchor     = 0;
    Double_t lDefaultV0MPercentile = 0;

    Double_t lDefaultCL0Anchor     = 0;
    Double_t lDefaultCL0Percentile = 0;

    Double_t lDefaultCL1Anchor     = 0;
    Double_t lDefaultCL1Percentile = 0;

    Double_t lDefaultZNAAnchor     = 0;
    Double_t lDefaultZNAPercentile = 0;

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
    if ( lPeriodName.Contains("LHC16q") ) {
        lDefaultZNAAnchor     = 0;
        lDefaultZNAPercentile = 96.0;
    }

    
    AliMultEstimator *fEstV0M = new AliMultEstimator("V0M", "", "(fAmplitude_V0A)+(fAmplitude_V0C)");
    AliMultEstimator *fEstV0A = new AliMultEstimator("V0A", "", "((fAmplitude_V0A)/(1+((fEvSel_VtxZ)-1.04)*(-0.004)))");
    AliMultEstimator *fEstV0C = new AliMultEstimator("V0C", "", "(fAmplitude_V0C)");

    AliMultEstimator *fEstOnlineV0M = new AliMultEstimator("OnlineV0M", "", "(fAmplitude_OnlineV0A)+(fAmplitude_OnlineV0C)");
    AliMultEstimator *fEstOnlineV0A = new AliMultEstimator("OnlineV0A", "", "(fAmplitude_OnlineV0A)");
    AliMultEstimator *fEstOnlineV0C = new AliMultEstimator("OnlineV0C", "", "(fAmplitude_OnlineV0C)");

    AliMultEstimator *fEstADM = new AliMultEstimator("ADM", "", "(fMultiplicity_ADA)+(fMultiplicity_ADC)");
    AliMultEstimator *fEstADA = new AliMultEstimator("ADA", "", "(fMultiplicity_ADA)");
    AliMultEstimator *fEstADC = new AliMultEstimator("ADC", "", "(fMultiplicity_ADC)");
    

    AliMultEstimator *fEstCL0 = new AliMultEstimator("CL0", "", "(fnSPDClusters0)/(1+((fEvSel_VtxZ)-1.18)*((-0.0005)+((fEvSel_VtxZ)-1.18)*(-0.001)))");

    AliMultEstimator *fEstCL1 = new AliMultEstimator("CL1", "", "(fnSPDClusters1)/(1+((fEvSel_VtxZ)-1.16)*((-0.0026)+((fEvSel_VtxZ)-1.16)*(-0.0009)))");

    //Integer estimators
    AliMultEstimator *fEstnSPDClusters = new AliMultEstimator("SPDClusters", "", "(fnSPDClusters)");
    fEstnSPDClusters->SetIsInteger(kTRUE);
    AliMultEstimator *fEstnSPDTracklets = new AliMultEstimator("SPDTracklets", "", "(fnTracklets)");
    fEstnSPDTracklets->SetIsInteger(kTRUE);
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

    AliMultEstimator *fEstZNA = new AliMultEstimator("ZNA", "", "(fZnaTower)*((fZnaTower)>0)");
    fEstZNA -> SetUseAnchor        ( kTRUE   ) ;
    fEstZNA -> SetAnchorPoint      ( lDefaultZNAAnchor    ) ;
    fEstZNA -> SetAnchorPercentile ( lDefaultZNAPercentile  ) ;


    
    //Universal: V0
    lCalib->GetMultSelection() -> AddEstimator( fEstV0M );
    lCalib->GetMultSelection() -> AddEstimator( fEstV0A );
    //lCalib->GetMultSelection() -> AddEstimator( fEstV0C );

    
    lCalib ->GetMultSelection() -> AddEstimator( fEstCL0 );
    lCalib ->GetMultSelection() -> AddEstimator( fEstCL1 );        

    lCalib->GetMultSelection() -> AddEstimator( fEstZNA );
    lCalib->GetMultSelection() -> AddEstimator( fEstnSPDTracklets );

    //============================================================
    // --- Definition of Input/Output ---
    //============================================================

    //lCalib -> SetInputFile  ( "/lustre/nyx/alice/users/alberica/centrality/Trees/Singles/LHC16q/muon_calo_pass1_1111/files/AnalysisResults_265309.root");
    lCalib -> SetInputFile  ( "/lustre/nyx/alice/users/alberica/centrality/Trees/Merged/pPb/Merged_LHC16q.root");
    lCalib -> SetBufferFile ( "/lustre/nyx/alice/users/alberica/centrality/OADB/LHC16q/buffer.root" );
    lCalib -> SetOutputFile ( "/lustre/nyx/alice/users/alberica/centrality/OADB/LHC16q/OADB-LHC16q.root" );
    lCalib -> Calibrate     ();

    
}
