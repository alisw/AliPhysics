

CalibratePeriodMC(TString lPeriodName = "LHC15k1"){
    cout<<"Run!"<<endl;
    
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

    //All fine, let's try the new MC calibrator
    AliMultSelectionCalibratorMC *lCalib = new AliMultSelectionCalibratorMC("lCalib");
    
    lCalib->SetupStandardInput();
    
    //Actual Input files 
    //lCalib -> SetInputFileData  ( "~/Dropbox/MultSelCalib/LHC10h/files/AnalysisResults_137161.root") ;
    lCalib -> SetInputFileData  ( "~/Dropbox/MultSelCalib/LHC15o/MergedLHC15o.root" ) ;
    //lCalib -> SetInputFileData  ( "~/Dropbox/MultSelCalib/LHC15o/files/AnalysisResults_244918.root") ;
    lCalib -> SetInputFileOADB  ( "../OADB-LHC15o-Complete.root") ;
    //lCalib -> SetInputFileMC    ( "~/Dropbox/MultSelCalib/LHC11a10a_bis/files/AnalysisResults_137161.root") ;
    lCalib -> SetInputFileMC    ( Form("~/Dropbox/MultSelCalib/%s/Merged%s.root",lPeriodName.Data(),lPeriodName.Data()) ) ;
    
    lCalib->SetRunToUseAsDefault(245064);
    
    //Buffer files 
    //lCalib -> SetBufferFileData ( "/home/daviddc/work/fast/buffer-data.root" );
    //lCalib -> SetBufferFileMC   ( "/home/daviddc/work/fast/buffer-mc.root" );
    //lCalib -> SetBufferFileData ( "buffer-data.root" );
    //lCalib -> SetBufferFileMC   ( "buffer-mc.root" );

    
    //Output OADB 
    lCalib -> SetOutputFile     ( Form("OADB-%s.root",lPeriodName.Data()) );
    
    lCalib -> Calibrate     ();
    
}
