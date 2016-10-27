////////////////////////////////////////////////////////////
//
// Macro for calibrating high-multiplicity-triggered pp
// sample. Uses the same functionality that is used for
// Pb-Pb anchoring.
//
// --- FOR TESTING PURPOSES
//
////////////////////////////////////////////////////////////

CalibratePeriodPP-LHC16l-HM(TString lPeriodName = "LHC16l") {

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

    //============================================================
    // --- Definition of Input Variables ---
    //============================================================

    lCalib->SetupStandardInput();

    //Changes in new version: create AliMultSelection here
    AliMultSelection *lMultSel = new AliMultSelection();
    
    //Commented: we'll work on this later!
    //lCalib->SetMultSelection(lMultSel);

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

    //Universal: V0
    lMultSel -> AddEstimator( fEstV0M );
    lMultSel -> AddEstimator( fEstV0A );
    lMultSel -> AddEstimator( fEstV0C );
    //Only do this in run 2, AD didn't exist in Run 1
    //Will also save space in the OADB for old datasets!
    if( lPeriodName.Contains("LHC15") ) {
        lMultSel -> AddEstimator( fEstOnlineV0M );
        lMultSel -> AddEstimator( fEstOnlineV0A );
        lMultSel -> AddEstimator( fEstOnlineV0C );
        lMultSel -> AddEstimator( fEstADM );
        lMultSel -> AddEstimator( fEstADA );
        lMultSel -> AddEstimator( fEstADC );
    }

    //Universal: Tracking, etc
    lMultSel -> AddEstimator( fEstnSPDClusters  );
    lMultSel -> AddEstimator( fEstnSPDTracklets );
    lMultSel -> AddEstimator( fEstRefMultEta5 );
    lMultSel -> AddEstimator( fEstRefMultEta8 );

    lMultSel -> AddEstimator( fEstZNApp );
    lMultSel -> AddEstimator( fEstZNCpp );
    lMultSel -> AddEstimator( fEstZNACpp );

    
    //============================================================
    // --- Definition of Anchor points ---
    //============================================================
    
    Long_t lRuns[] = {258883, 258884, 258885, 258886, 258889, 258890, 258919, 258920, 258921, 258923, 258926, 258931, 258962, 258964, 259086, 259088, 259090, 259091, 259096, 259099, 259117, 259118, 259162, 259164, 259204, 259257, 259261, 259263, 259264, 259269, 259270, 259271, 259272, 259273, 259274, 259302, 259303, 259305, 259307, 259334, 259336, 259339, 259340, 259341, 259342, 259378, 259381, 259382, 259388, 259389, 259394, 259395, 259396, 259469, 259473, 259477, 259649, 259650, 259668, 259697, 259700, 259703, 259704, 259705, 259711, 259713, 259747, 259748, 259750, 259751, 259752, 259756, 259781, 259788, 259789, 259792, 259822, 259841, 259842, 259860, 259866, 259867, 259868, 259888, 259961, 259979, 260010, 260011, 260014};
    
    Float_t lAPs[] = {493.879, 504.509, 500.107, 495.953, 498.561, 448.549, 444.875, 445.283, 449.538, 447.706, 448.11, 446.215, 444.947, 441.085, 443.288, 440.476, 439.291, 438.114, 444.897, 441.659, 444.007, 448.376, 441.823, 440.908, 440.401, 441.003, 439.717, 438.799, 438.091, 436.965, 437.072, 436.919, 439.53, 440.268, 435.717, 446.6, 439.121, 435.649, 430.636, 437.265, 438.004, 437.388, 432.947, 435.139, 431.794, 440.464, 436.485, 440.164, 424.503, 436.452, 436.673, 435.215, 438.357, 442.597, 438.327, 437.191, 437.562, 436.487, 438.799, 435.908, 435.639, 435.279, 436.202, 434.805, 437.008, 432.496, 437.726, 444.541, 437.591, 438.887, 435.292, 434.204, 437.705, 437.221, 435.314, 429.356, 434.347, 434.321, 439.071, 437.444, 434.344, 428.718, 434.926, 441.239, 430.507, 432.405, 441.327, 442.942, 439.949};
    
    const Long_t lNRuns = sizeof(lRuns) / sizeof(Long_t);
    AliMultSelection *lMultSelArray[500]; //We need one AliMultSelection for each run above!
    
    //Set Anchor point (common)
    lMultSel->GetEstimator("V0M")->SetUseAnchor( kTRUE ) ;
    lMultSel->GetEstimator("V0M")->SetAnchorPercentile( 0.1 );
    
    //Sweep array
    for(Long_t iRun = 0; iRun<lNRuns; iRun++) {
        //cout<<"Instantiating AliMultSelection for run "<<lRuns[iRun]<<"... "<<endl;
        lMultSelArray[iRun] = new AliMultSelection(lMultSel);
        
        //Adjust V0M Anchors, please
        lMultSelArray[iRun]->GetEstimator("V0M")->SetAnchorPoint( lAPs[iRun] );
        
        //Debug: check if changed!
        //lMultSelArray[iRun] -> PrintInfo();
        
        //cout<<"Adding run range..."<<endl;
        //Add run to run ranges (still RbyR for now)
        lCalib->AddRunRange( lRuns[iRun], lRuns[iRun], lMultSelArray[iRun] );
    }
    lCalib->SetRunToUseAsDefault(259668);
    
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
