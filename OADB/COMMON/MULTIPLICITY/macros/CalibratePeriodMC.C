CalibratePeriodMC(){
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
    
    //Setup the usual standard inputs for the calibration
    lCalib->SetupStandardInput();
    
    //Actual Input files
    lCalib -> SetInputFileData  ( "~/Dropbox/MultSelCalib/LHC15mo/MergedLHC15mo.root") ;
    lCalib -> SetInputFileOADB  ( "$ALICE_PHYSICS/../src/OADB/COMMON/MULTIPLICITY/data/OADB-LHC15m.root") ;
    lCalib -> SetInputFileMC    ( "~/Dropbox/MultSelCalib/LHC15k1_plus/MergedLHC15k1_plus.root") ;
    
    //Default run: this is the run whose scaling factor will be saved as default
    //OADB object. Note that this is only a good guess and using this default will print
    //out a warning message
    lCalib->SetRunToUseAsDefault(244918);
    
    //Buffer files (these may become large, position there whereever convenient)
    lCalib -> SetBufferFileData ( "buffer-data.root" );
    lCalib -> SetBufferFileMC   ( "buffer-mc.root" );

    //Output micro-OADB
    lCalib -> SetOutputFile     ( "OADB-LHC15k1_plus.root" );
    lCalib -> Calibrate     ();
    
}
