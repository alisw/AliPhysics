////////////////////////////////////////////////////////////
//
// Macro for determining which raw V0M value corresponds
// to a certain high-multiplicity-percentile of a min-bias
// data sample for anchoring high-multiplicity-triggered
// data
//
// --- FOR TESTING PURPOSES
//
////////////////////////////////////////////////////////////

StitchOADBs( ){
    //Macro to determine anchor points
    
    //Long_t lRuns[] = {258883, 258884, 258885, 258886, 258889, 258890, 258919, 258920, 258921, 258923, 258926, 258931, 258962, 258964, 259086, 259088, 259090, 259091, 259096, 259099, 259117, 259118, 259162, 259164, 259204, 259257, 259261, 259263, 259264, 259269, 259270, 259271, 259272, 259273, 259274, 259302, 259303, 259305, 259307, 259334, 259336, 259339, 259340, 259341, 259342, 259378, 259381, 259382, 259388, 259389, 259394, 259395, 259396, 259469, 259473, 259477, 259649, 259650, 259668, 259697, 259700, 259703, 259704, 259705, 259711, 259713, 259747, 259748, 259750, 259751, 259752, 259756, 259781, 259788, 259789, 259792, 259822, 259841, 259842, 259860, 259866, 259867, 259868, 259888, 259961, 259979, 260010, 260011, 260014};
    Long_t lRuns[] = { 259668 };
    Long_t lNRuns = sizeof( lRuns ) / sizeof( Long_t );
    cout<<"Registered "<<lNRuns<<" runs. "<<endl;
    
    //Stitch information into these arrays
    Float_t lBins [1000];
    Float_t lPerc [1000];
    Long_t lNBins = 0;
    
    for(Long_t ibin=0; ibin<1000; ibin++){
        lBins[ibin] = 0;
    }
    
    //ORIGINAL
    TString lOADBfile1 = "OADB-LHC16l.root";
    
    //EXTRA HIGH-MULTIPLICITY
    TString lOADBfile2 = "OADB-testing.root";
    
    //STITCHED
    TString lOADBfileout = "OADB-out.root";
    
    cout<<"Opening low-multiplicity info ... "<<endl;
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
    
    TH1F *lThisCalibHisto1 = 0x0;
    TH1F *lThisCalibHisto2 = 0x0;
    for( Long_t iRun = 0; iRun < lNRuns; iRun++){
        AliOADBMultSelection *lOADB1 = c1->GetObject(lRuns[iRun], "Default");
        lThisCalibHisto1 = lOADB1->GetCalibHisto( "hCalib_V0M" );
        lThisCalibHisto1->SetName("hCalib_V0M_1");
        
        AliOADBMultSelection *lOADB2 = c2->GetObject(lRuns[iRun], "Default");
        lThisCalibHisto2 = lOADB2->GetCalibHisto( "hCalib_V0M" );
        lThisCalibHisto2->SetName("hCalib_V0M_2");
        
        //______________________________________________________________________________
        //Step 1: Process low-multiplicity boundaries all the way up to the anchor point
        lBins[0] = 0.; //starting at zero.
        lNBins = 0;
        for(Long_t ibin=1; ibin<lThisCalibHisto1->GetNbinsX()+1; ibin++) {
            //Check if still above anchor point: if yes, add; if no, break
            if( lThisCalibHisto1->GetBinContent(ibin) > 0.1 ){
                //Beware of bin number offsets!
                lPerc[lNBins     ] = lThisCalibHisto1->GetBinContent(ibin  );
                lBins[lNBins + 1 ] = lThisCalibHisto1->GetBinLowEdge(ibin+1);
                lNBins ++ ;
            }else{
                break;
            }
        }
        //______________________________________________________________________________
        //Step 2: Process high-multiplicity boundaries starting from the anchor point
        for(Long_t ibin=1; ibin<lThisCalibHisto2->GetNbinsX()+1; ibin++) {
            //Check if still above anchor point: if yes, add; if no, do nothing
            if( lThisCalibHisto1->GetBinContent(ibin) < 0.1 ){
                //Beware of bin number offsets!
                lPerc[lNBins     ] = lThisCalibHisto1->GetBinContent(ibin  );
                lBins[lNBins + 1 ] = lThisCalibHisto2->GetBinLowEdge(ibin+1);
                lNBins ++ ;
            }
        }
        //______________________________________________________________________________
        //Step 3: Dump boundaries for cross-checking, if requested
        for(Long_t ibin=0; ibin<lNBins; ibin++) {
        //    cout<<"ibin = "<<ibin<<" low: "<<lBins[ibin]<<" up: "<<lBins[ibin+1]<<", perc = "<<lPerc[ibin]<<endl;
        }

        //______________________________________________________________________________
        //Step 4: Create and store output into output AliOADBContainer
        
        //Basic properties
        oadbMultSelection = new AliOADBMultSelection();
        cuts              = new AliMultSelectionCuts();
        cuts = lOADB1->GetEventCuts();
        fsels             = new AliMultSelection( lOADB1->GetMultSelection() );
        
        oadbMultSelection->SetEventCuts        ( cuts  );
        oadbMultSelection->SetMultSelection    ( fsels );
        
        //Determine Quantiles from calibration histogram
        TH1F *lThisCalibHisto = 0x0;
        TString lThisCalibHistoName;
        TString lV0M = "V0M";
        
        for(Long_t iEst=0; iEst<fsels->GetNEstimators(); iEst++) {
            //is this V0M? Otherwise, I don't care
            cout<<"Processing estimator named: "<<fsels->GetEstimator(iEst)->GetName()<<endl;
            if( fsels->GetEstimator(iEst)->GetName() != lV0M ){
                lThisCalibHistoName = Form("hCalib_%s",fsels->GetEstimator(iEst)->GetName());
                lThisCalibHisto = 0x0;
                lThisCalibHisto = lOADB1->GetCalibHisto( lThisCalibHistoName );
                oadbMultSelection->AddCalibHisto( lThisCalibHisto );
            }else{
                //Create special calibration histo for V0M to supercede the one in the OADB
                TH1F *hCalib_V0M = new TH1F("hCalib_V0M","",lNBins,lBins);
                for(Long_t ibin=1; ibin<hCalib_V0M->GetNbinsX()+1; ibin++){
                    hCalib_V0M -> SetBinContent(ibin, lPerc[ibin-1]);
                }
                oadbMultSelection->AddCalibHisto( hCalib_V0M );
            }
        }
        oadbContMS->AppendObject(oadbMultSelection, lRuns[iRun], lRuns[iRun] );
    }
    
    oadbContMS->Write();
}
