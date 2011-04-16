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
}
