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

DetermineAnchorsPP( TString lData = "LHC16l" ){
    //Macro to determine anchor points
    
    Long_t lRuns[] = {258883, 258884, 258885, 258886, 258889, 258890, 258919, 258920, 258921, 258923, 258926, 258931, 258962, 258964, 259086, 259088, 259090, 259091, 259096, 259099, 259117, 259118, 259162, 259164, 259204, 259257, 259261, 259263, 259264, 259269, 259270, 259271, 259272, 259273, 259274, 259302, 259303, 259305, 259307, 259334, 259336, 259339, 259340, 259341, 259342, 259378, 259381, 259382, 259388, 259389, 259394, 259395, 259396, 259469, 259473, 259477, 259649, 259650, 259668, 259697, 259700, 259703, 259704, 259705, 259711, 259713, 259747, 259748, 259750, 259751, 259752, 259756, 259781, 259788, 259789, 259792, 259822, 259841, 259842, 259860, 259866, 259867, 259868, 259888, 259961, 259979, 260010, 260011, 260014};
    Long_t lNRuns = sizeof( lRuns ) / sizeof( Long_t );
    cout<<"Registered "<<lNRuns<<" runs. "<<endl;
    
    Float_t lAPs[200];
    
    cout<<"Opening OADB ... "<<endl;
    TFile *file = new TFile(Form("../../data/OADB-%s.root",lData.Data() ) );
    AliOADBContainer *c = file->Get("MultSel");
    
    TH1F *lThisCalibHisto = 0x0;
    for( Long_t iRun = 0; iRun < lNRuns; iRun++){
        AliOADBMultSelection *lOADB = c->GetObject(lRuns[iRun], "Default");
        lThisCalibHisto = lOADB->GetCalibHisto( "hCalib_V0M" );
        for (Long_t ibin=1; ibin<lThisCalibHisto->GetNbinsX()+1; ibin++){
            if( lThisCalibHisto->GetBinContent(ibin)<0.1){
                lAPs[iRun] = lThisCalibHisto->GetBinLowEdge(ibin);
                break;
            }
        }
    }
    
    cout<<"Printout: "<<endl;
    for (Long_t iRun = 0; iRun<lNRuns; iRun++) {
        cout<<"Run #"<<iRun<<": "<<lRuns[iRun]<<" AP = "<<lAPs[iRun]<<endl;
    }
    
    cout<<"For copy-pasting:"<<endl;
    cout<<" Float_t lAPs[] = {"<<flush;
    for (Long_t iRun = 0; iRun<lNRuns; iRun++){
        cout<<lAPs[iRun]<<flush;
        if (iRun!=lNRuns-1) cout<<", "<<flush;
    }
    cout<<"};"<<endl;
    
    
    
} 
