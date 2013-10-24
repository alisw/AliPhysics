void aod(){

    gSystem->Load("libTree");
    gSystem->Load("libGeom");
    gSystem->Load("libPhysics");
    gSystem->Load("libVMC");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libCORRFW");
    gSystem->Load("libPWGHFbase");
    gSystem->Load("libPWGmuon");

    gROOT->Macro("${ALICE_ROOT}/STEER/CreateAODfromESD.C(\"AliESDs.root\",\"AliAODs.root\",kFALSE)");
}
