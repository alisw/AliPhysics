void LoadLibraries() {
    //
    // load libraries needed for CPass1
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
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libANALYSIScalib");
    gSystem->Load("libTender");
    gSystem->Load("libPWGPP");
    //gSystem->Load("libPWG4PartCorrBase");
    //gSystem->Load("libHMPIDbase");
}
