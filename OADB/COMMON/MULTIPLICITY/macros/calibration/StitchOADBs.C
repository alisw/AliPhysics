////////////////////////////////////////////////////////////
//
// Macro for stitching calibration histograms obtained with
// minimum bias and high-multiplicity triggered data to 
// provide better V0M percentile estimation for pp high-
// multiplicity analyses 
//
////////////////////////////////////////////////////////////

StitchOADBs(TString       lPeriodName = "LHC16k", 
            Long_t lRunToUseAsDefault = 257630, 
            Bool_t          lDrawFlag = kTRUE ) {

    //Read run ranges and anchor percentiles/points from file
    vector<Int_t> lRuns;
    vector<float> lAnchorPercentile;
    vector<float> lAnchorPoint;
    FILE *fAP = fopen(Form("Anchors_%s_VHM.txt", lPeriodName.Data()), "r");
    Char_t buffer[500];
    while(fgets(buffer, 500, fAP)) {
        Int_t run_a;
        Int_t run_b;
        Float_t ancPer;
        Float_t ancPoint;
        sscanf(buffer, "%d %d %f %f", &run_a, &run_b, &ancPer, &ancPoint);
        //Debug; print info retrieved from file
        printf("Run: %d %d \t %.2f%%   %f\n", run_a, run_b, ancPer, ancPoint);
        if(run_a!=run_b) {
            cout << "Error: runs don't match!" << endl;
            return;
        }
        lRuns.push_back(run_a);
        lAnchorPercentile.push_back(ancPer);
        lAnchorPoint.push_back(ancPoint);
    }
    
    Long_t lNRuns = lRuns.size();
    cout<<"Registered "<<lNRuns<<" runs. "<<endl;
    
    //Stitch information into these arrays
    Float_t lBins [1000];
    Float_t lPerc [1000];
    Long_t lNBins = 0;
    TString lTag  [1000]; // just for cross-checking
    
    for(Long_t ibin=0; ibin<1000; ibin++){
        lBins[ibin] = 0;
    }
    
    //ORIGINAL
    TString lOADBfile1 = Form("OADB-%s-MB.root", lPeriodName.Data());
    
    //EXTRA HIGH-MULTIPLICITY
    TString lOADBfile2 = Form("OADB-%s-VHM.root", lPeriodName.Data());
    
    //STITCHED
    TString lOADBfileout = Form("OADB-%s.root", lPeriodName.Data());
    
    cout<<"Opening minimum bias info ... "<<endl;
    TFile *file1 = new TFile( lOADBfile1.Data(), "READ" );
    AliOADBContainer *c1 = file1->Get("MultSel");
    c1->SetName("MultSel1") ;

    cout<<"Opening high-multiplicity info..."<<endl;
    TFile *file2 = new TFile( lOADBfile2.Data(), "READ" );
    AliOADBContainer *c2 = file2->Get("MultSel");
    c2->SetName("MultSel2") ;
    
    //Open output OADB file, generate everything within loop
    TFile * fileout = new TFile (lOADBfileout.Data(), "recreate");
    AliOADBContainer * oadbContMS = new AliOADBContainer("MultSel");
    
    AliOADBMultSelection * oadbMultSelection = 0x0;
    AliMultSelectionCuts * cuts = 0x0;
    AliMultSelection     * fsels = 0x0;

    AliOADBMultSelection* oadbMultSelectionDefault = 0x0; 
    AliMultSelectionCuts * cutsDefault = 0x0;
    AliMultSelection     * fselsDefault = 0x0;
    
    TH1F *lThisCalibHisto1 = 0x0;
    TH1F *lThisCalibHisto2 = 0x0;

    for( Long_t iRun = 0; iRun < lNRuns; iRun++){

        cout << endl;
        cout << "Run " << iRun << ": " << lRuns[iRun] << endl; 
        cout << " - Anchor Percentile: " << lAnchorPercentile[iRun] << endl;
        cout << " - Anchor Point.....: " << lAnchorPoint[iRun] << endl;

        lNBins = 0;
        
        //______________________________________________________________________________
        //Step 1: Process low-multiplicity boundaries all the way up to the anchor point
        lThisCalibHisto1 = 0x0;
        //get calib histo corresponding to the first run in the range (should be the same for the others)
        AliOADBMultSelection *lOADB1 = (AliOADBMultSelection*)c1->GetObject( lRuns[iRun], "Default");
        if( (Int_t)c1->GetIndexForRun( lRuns[iRun] )<0 ) {
            cout << "Warning: no MB calibration histo found for run " << lRuns[iRun] << " - skipping..." << endl;
            continue;
        }
        lThisCalibHisto1 = lOADB1->GetCalibHisto( "hCalib_V0M" );
        lThisCalibHisto1->SetName( "hCalib_V0M_1" );
        lBins[0] = 0.; //starting at zero.
        for(Long_t ibin=1; ibin<lThisCalibHisto1->GetNbinsX()+1; ibin++) {
            //Check if still above anchor point: if yes, add; if no, break
            if( lThisCalibHisto1->GetBinContent(ibin) > lAnchorPercentile[iRun] ) {
                //Beware of bin number offsets!
                lTag [lNBins     ] = "MB";
                lPerc[lNBins     ] = lThisCalibHisto1->GetBinContent(ibin  );
                lBins[lNBins + 1 ] = lThisCalibHisto1->GetBinLowEdge(ibin+1);
                lNBins ++ ;
            } 
            else break;
        }
        cout << " - Retrieved calibration histogram from MB..." << endl;
        
        //______________________________________________________________________________
        //Step 2: Process high-multiplicity boundaries starting from the anchor point
        lThisCalibHisto2 = 0x0;
        AliOADBMultSelection *lOADB2 = (AliOADBMultSelection*)c2->GetObject( lRuns[iRun], "Default");
        if( (Int_t)c2->GetIndexForRun( lRuns[iRun] )<0 ) {
            cout << "Warning: no VHM calibration histo found for range " << lRuns[iRun] << " - using calibration from MB..." << endl;
            for(Long_t ibin=1; ibin<lThisCalibHisto1->GetNbinsX()+1; ibin++) {
                //Simply continue adding for percentiles below anchor point
                if( lThisCalibHisto1->GetBinContent(ibin) <= lAnchorPercentile[iRun] ) {
                    //Beware of bin number offsets!
                    lTag [lNBins     ] = "MB";
                    lPerc[lNBins     ] = lThisCalibHisto1->GetBinContent(ibin  );
                    lBins[lNBins + 1 ] = lThisCalibHisto1->GetBinLowEdge(ibin+1);
                    lNBins ++ ;
                }
            }
        }
        else {
            lThisCalibHisto2 = (TH1F*)lOADB2->GetCalibHisto( "hCalib_V0M" );
            lThisCalibHisto2->SetName("hCalib_V0M_2");
            for(Long_t ibin=1; ibin<lThisCalibHisto2->GetNbinsX()+1; ibin++) {
                //Check if below anchor point: if yes, add; if no, do nothing
                if( lThisCalibHisto2->GetBinContent(ibin) <= lAnchorPercentile[iRun] || lThisCalibHisto2->GetBinLowEdge(ibin) >= lAnchorPoint[iRun] ) {
                //if( lThisCalibHisto2->GetBinContent(ibin) <= lAnchorPercentile[iRun] ) {
                    //Beware of bin number offsets!
                    lTag [lNBins     ] = "VHM";
                    lPerc[lNBins     ] = lThisCalibHisto2->GetBinContent(ibin  );
                    lBins[lNBins + 1 ] = lThisCalibHisto2->GetBinLowEdge(ibin+1);
                    lNBins ++ ;
                }
            }
        }
        cout << " - Retrieved calibration histogram from VHM..." << endl;
        AliMultEstimator * multEstimator = (AliMultEstimator*)((AliMultSelection*)lOADB2->GetMultSelection())->GetEstimator("V0M");
        cout << " - Anchor Percentile: " << (Double_t)multEstimator->GetAnchorPercentile() << endl;
        cout << " - Anchor Point.....: " << (Double_t)multEstimator->GetAnchorPoint() << endl;

        //______________________________________________________________________________
        //Step 3: Dump boundaries for cross-checking, if requested
        //for(Long_t ibin=0; ibin<lNBins; ibin++) {
        //    cout<<"ibin = "<<ibin<<" low: "<<lBins[ibin]<<" up: "<<lBins[ibin+1]<<", perc = "<<lPerc[ibin]<< " (" << lTag[ibin] << ")" <<endl;
        //}


        //______________________________________________________________________________
        //Step 4: Create and store output into output AliOADBContainer
        
        //Basic properties
        oadbMultSelection = new AliOADBMultSelection();
        cuts              = new AliMultSelectionCuts();
        cuts = lOADB1->GetEventCuts();
        fsels             = new AliMultSelection( lOADB1->GetMultSelection() );
        
        oadbMultSelection->SetEventCuts        ( cuts  );
        oadbMultSelection->SetMultSelection    ( fsels );

        // default run
        if( lRuns[iRun] == lRunToUseAsDefault ) {
            cout << " - Default Run Calibration" << endl;
            oadbMultSelectionDefault = new AliOADBMultSelection("Default"); 
            oadbMultSelectionDefault->SetEventCuts( cuts );
            oadbMultSelectionDefault->SetMultSelection( fsels );
        }

        //Determine Quantiles from calibration histogram
        TH1F *lThisCalibHisto = 0x0;
        TString lThisCalibHistoName;
        TString lV0M = "V0M";
        
        for(Long_t iEst=0; iEst<fsels->GetNEstimators(); iEst++) {
            //is this V0M? Otherwise, I don't care
            //cout<<"Processing estimator named: "<<fsels->GetEstimator(iEst)->GetName()<<endl;
            if( fsels->GetEstimator(iEst)->GetName() != lV0M ){
                lThisCalibHistoName = Form("hCalib_%s",fsels->GetEstimator(iEst)->GetName());
                lThisCalibHisto = 0x0;
                lThisCalibHisto = lOADB1->GetCalibHisto( lThisCalibHistoName );
                oadbMultSelection->AddCalibHisto( lThisCalibHisto );
                if( lRuns[iRun] == lRunToUseAsDefault ) oadbMultSelectionDefault->AddCalibHisto( lThisCalibHisto );
            }
            else {
                //Create special calibration histo for V0M to supercede the one in the OADB
                TH1F *hCalib_V0M = new TH1F("hCalib_V0M","",lNBins,lBins);
                for(Long_t ibin=1; ibin<hCalib_V0M->GetNbinsX()+1; ibin++){
                    hCalib_V0M -> SetBinContent(ibin, lPerc[ibin-1]);
                }
                //Debug: draw calib histos for checking
                if(lDrawFlag && !gROOT->IsBatch()) {
                    TCanvas *cc = new TCanvas("cc", "", 10, 10, 1000, 800);
                    TH1D* h1 = (TH1D*)lThisCalibHisto1->Clone(Form("h1_%d", iRun));
                    h1->GetXaxis()->SetRangeUser(lAnchorPoint[iRun]-20., lAnchorPoint[iRun]+20.);
                    h1->SetLineWidth(2);
                    h1->Draw();
                    if(lThisCalibHisto2) {
                        TH1D* h2 = (TH1D*)lThisCalibHisto2->Clone(Form("h2_%d", iRun));
                        h2->SetLineColor(2);
                        h2->SetLineWidth(2);
                        h2->Draw("same");
                    }
                    TH1D* h3 = (TH1D*)hCalib_V0M->Clone(Form("h3_%d", iRun));
                    h3->SetLineColor(3);
                    h3->Draw("same");
                    cc->Update();
                    cc->WaitPrimitive();
                    //cin.get(); // wait
                }

                oadbMultSelection->AddCalibHisto( hCalib_V0M );
                if( lRuns[iRun] == lRunToUseAsDefault ) oadbMultSelectionDefault->AddCalibHisto( hCalib_V0M );
                //
                // store anchor percentile for high multiplicity triggers
                fsels->GetEstimator(iEst)->SetAnchorPercentile( lAnchorPercentile[iRun] );
            }
        }
        oadbContMS->AppendObject(oadbMultSelection, lRuns[iRun], lRuns[iRun] );
    }

    // add default calibration
    if( oadbMultSelectionDefault ) oadbContMS->AddDefaultObject( oadbMultSelectionDefault );

    oadbContMS->Write();
}

