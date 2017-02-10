void aod(){

    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");

    if( gSystem->Load("libCORRFW") <0) {
      cerr << "Error: AliPhysics is not installed, no AOD test is possible!" << endl;
      return;
    }
    gSystem->Load("libPWGHFbase");
    gSystem->Load("libPWGmuon");
    gSystem->Load("libESDfilter");
    gSystem->Load("libTender");
    gSystem->Load("libPWGPP");

    AliEMCALGeometry::GetInstance("EMCAL_COMPLETEV1");

    gROOT->Macro("${ALICE_ROOT}/STEER/macros/CreateAODfromESD.C(\"AliESDs.root\",\"AliAOD.root\",\"local://$ALICE_ROOT/OCDB\",\"local://.\")");
}
