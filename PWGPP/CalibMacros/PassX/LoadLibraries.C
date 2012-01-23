void LoadLibraries(){
    gSystem->Load("libTree");
    gSystem->Load("libGeom");
    gSystem->Load("libPhysics");
    gSystem->Load("libVMC");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");

    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libANALYSIScalib");
    gSystem->Load("libCORRFW");
    gSystem->Load("libPWGmuon");
    //
    // detector libraries
    //    
    gSystem->Load("libTPCcalib");
}
