////////////////////////////////////////////////////////////
//
// Macro for calibrating high-multiplicity-triggered pp
// sample. Uses the same functionality that is used for
// Pb-Pb anchoring.
//
////////////////////////////////////////////////////////////

ReCalibratePeriodPP(const Char_t* inputDir, Int_t runNo, 
                    TString lPeriodName = "LHC16k",
                    TString lWhichData = "VHM",
                    Long_t lRunToUseAsDefault = 257630 ) {

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
    lCalib->SetSelectedTriggerClass(AliVEvent::kHighMultV0);
    
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
    lMultSel->AddEstimator( fEstV0M );
    lMultSel->AddEstimator( fEstV0A );
    lMultSel->AddEstimator( fEstV0C );
    //Only do this in run 2, AD didn't exist in Run 1
    //Will also save space in the OADB for old datasets!
    if( lPeriodName.Contains("LHC15") ) {
        lMultSel->AddEstimator( fEstOnlineV0M );
        lMultSel->AddEstimator( fEstOnlineV0A );
        lMultSel->AddEstimator( fEstOnlineV0C );
        lMultSel->AddEstimator( fEstADM );
        lMultSel->AddEstimator( fEstADA );
        lMultSel->AddEstimator( fEstADC );
    }

    //Universal: Tracking, etc
    lMultSel->AddEstimator( fEstnSPDClusters  );
    lMultSel->AddEstimator( fEstnSPDTracklets );
    lMultSel->AddEstimator( fEstRefMultEta5 );
    lMultSel->AddEstimator( fEstRefMultEta8 );

    lMultSel->AddEstimator( fEstZNApp );
    lMultSel->AddEstimator( fEstZNCpp );
    lMultSel->AddEstimator( fEstZNACpp );

    
    //============================================================
    // --- Definition of Anchor points ---
    //============================================================
    vector<int> lRuns;
    vector<float> lAnchorPercentile;
    vector<float> lAnchorPoint;
    FILE *fAP = fopen(Form("Anchors_%s_%s.txt", lPeriodName.Data(), lWhichData.Data()), "r");
    Char_t buffer[500];
    while(fgets(buffer, 500, fAP)) {
      Int_t run_a;
      Int_t run_b;
      Float_t ancPer;
      Float_t ancPoint;
      sscanf(buffer, "%d %d %f %f", &run_a, &run_b, &ancPer, &ancPoint);
      if(run_a!=run_b) { cout << "Error: runs do not match" << endl; return -1; }
      if(run_a!=runNo) continue;
      //Debug; print info retrieved from file
      printf("Run: %d %d \t %.2f%%   %f\n", run_a, run_b, ancPer, ancPoint);
      if(ancPer<0.0 || ancPoint<0.0) return;
      lRuns.push_back(run_a);
      lAnchorPercentile.push_back(ancPer);
      lAnchorPoint.push_back(ancPoint);
    }
    
    const Long_t lNRuns = lRuns.size();
    AliMultSelection *lMultSelArray[500]; //We need one AliMultSelection for each run above!
    
    //Set Anchor point (common)
    lMultSel->GetEstimator("V0M")->SetUseAnchor( kTRUE ) ;
    //lMultSel->GetEstimator("V0M")->SetAnchorPercentile( 0.1 ); // will be done inside the loop
    
    //Sweep array
    for(Long_t iRun = 0; iRun<lNRuns; iRun++) {
        cout<<"Instantiating AliMultSelection for run "<<lRuns[iRun]<<"... "<<endl;
        lMultSelArray[iRun] = new AliMultSelection(lMultSel);
        
        //Adjust V0M Anchors, please
        cout << "Setting Anchors........: " << lAnchorPercentile[iRun] << "% --> " << lAnchorPoint[iRun] << endl;
        lMultSelArray[iRun]->GetEstimator("V0M")->SetAnchorPercentile( lAnchorPercentile[iRun] );
        lMultSelArray[iRun]->GetEstimator("V0M")->SetAnchorPoint     ( lAnchorPoint[iRun]      );
        
        //Debug: check if changed!
        lMultSelArray[iRun] -> PrintInfo();
        
        //cout<<"Adding run range..."<<endl;
        //Add run to run ranges (still RbyR for now)
        lCalib->AddRunRange( lRuns[iRun], lRuns[iRun], lMultSelArray[iRun] );
    }
    lCalib->SetRunToUseAsDefault( lRunToUseAsDefault );
    
    //============================================================
    // --- Definition of Input/Output ---
    //============================================================
    
    lCalib -> SetInputFile  ( Form("%s/AnalysisResults_%d.root", inputDir, runNo) );
    lCalib -> SetBufferFile ( Form("temp/buffers/buffer-%s-%d-%s.root", lPeriodName.Data(), runNo, lWhichData.Data()) );
    lCalib -> SetOutputFile ( Form("temp/partialOADBs/OADB-%s-%d-%s.root", lPeriodName.Data(), runNo, lWhichData.Data()) );
    lCalib -> Calibrate     ();
}
