void LoadLibraries() {
    //
    // load libraries needed for Pass0
    //
    gSystem->Load("libSTAT");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libANALYSIScalib");
    //
    // detector libraries
    //    
    gSystem->Load("libTPCcalib");
    gSystem->Load("libTRDcalib");
    gSystem->Load("libT0calib");
    gSystem->Load("libTOFcalib");
    //
    // PWGPP libraries
    //    
    gSystem->Load("libANALYSISalice.so");
    gSystem->Load("libANALYSIScalib.so");
    gSystem->Load("libTENDER.so");
    gSystem->Load("libPWGPP.so");
    //gSystem->Load("libPWG4PartCorrBase.so");
    //gSystem->Load("libHMPIDbase.so");
}
