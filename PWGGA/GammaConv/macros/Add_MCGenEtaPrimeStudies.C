AliGenerator* CreatePythia8Gen( Float_t e_cms,
                                Int_t tune,
                                TString specialSetting,
                                Double_t etaPrimeSup,
                                Int_t kProcess,
                                Double_t ptHardMin,
                                Double_t ptHardMax
                            );

AliGenerator* Add_MCGenEtaPrimeStudies(     Float_t e_cms           = 13000., 
                                            Int_t tune              = 5, 
                                            TString specialSetting  = "EtaPrimeBiased",
                                            Double_t etaPrimeSup    = 0.12,
                                            Int_t kProcess          = 0, 
                                            Double_t ptHardMin      = 0, 
                                            Double_t ptHardMax      = 1.
                                            ) {
    // Add Pythia 8 generator: 
    //    -kProcess=0  MB generation
    //    -kProcess=1  Jet production, pthard generation
    //
    // Special setting:
    //      EtaPrimeBiased: require at least one eta' in each event, all eta' decays into pi+ pi- eta, all eta into gamma gamma
    //      EtaPrimeNotSuppressedAndFixedDecays: etaPrimeSup factor in Pythia set to given value (default 0.12); all eta' decay into pi+ pi- eta, all eta into gamma gamma
    //      EtaPrimeNotSuppressed: etaPrimeSup factor in Pythia set to given value (default 0.12)

    gSystem->Load("liblhapdf");
    
    AliGenerator *genP  = NULL;
    genP                = CreatePythia8Gen(e_cms, tune, specialSetting.Data(), etaPrimeSup, kProcess, ptHardMin, ptHardMax);
    
    return genP;
}

AliGenerator* CreatePythia8Gen( Float_t e_cms, 
                                Int_t tune, 
                                TString specialSetting,
                                Double_t etaPrimeSup,
                                Int_t kProcess, 
                                Double_t ptHardMin, 
                                Double_t ptHardMax
                            ) {
    
    gSystem->Load("libpythia6");
    gSystem->Load("libEGPythia6");
    gSystem->Load("libAliPythia6");
    gSystem->Load("libpythia8");
    gSystem->Load("libAliPythia8");
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
    gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));


    AliPythia8 *pythia = AliPythia8::Instance();            // For ROOT6 needs to be created before AliGenPythiaPlus object, otherwise ending in "illegal instruction"
    AliGenPythiaPlus* gener = new AliGenPythiaPlus(pythia);

    std::cout << "*****************************************************************" << std::endl;
    std::cout << "Process: "<< kProcess << "\t Tune: " << tune << "\t Special eta' setting: "<< specialSetting.Data() << "\t Eta suppression factor: "<< etaPrimeSup <<  std::endl;
    if (kProcess == 1){
        std::cout << "pTHardMin: "<< ptHardMin << "\t ptHardMax: "<< ptHardMax <<  std::endl;
    }	
    std::cout << "*****************************************************************" << std::endl;
    // set process (MB)
    if(kProcess==0){
        std::cout << "running in min bias mode"<< std::endl;
        gener->SetProcess(kPyMbDefault);
    }
    
    if(kProcess==1) {
        gener->SetProcess(kPyJets);
        std::cout << "went into jet loop" << std::endl;
        if(ptHardMin > 0.){ 
            std::cout << "Setting pTHardMin: "<< ptHardMin << "\t ptHardMax: "<< ptHardMax <<  std::endl;
            gener->SetPtHard(ptHardMin,ptHardMax);
        }
    } 

    //Centre of mass energy 
    gener->SetEnergyCMS(e_cms); // in GeV

    // Event list
    gener->SetEventListRange(-1, -1);

    //random seed based on time
    (AliPythia8::Instance())->ReadString("Random:setSeed = on");
    (AliPythia8::Instance())->ReadString("Random:seed = 0");

    if (specialSetting.EqualTo("EtaPrimeBiased")) {
        std::cout << "special setting " << specialSetting.Data() << std::endl;
        gener->SetTriggerParticle(331, 1.2);
        (AliPythia8::Instance())->ReadString("331:onMode = off");
        (AliPythia8::Instance())->ReadString("331:onIfMatch = 211 -211 221");
        (AliPythia8::Instance())->ReadString("221:onMode = off");
        (AliPythia8::Instance())->ReadString("221:onIfMatch = 22 22");
    } else if (specialSetting.EqualTo("EtaPrimeNotSuppressedAndFixedDecays")) {
        std::cout << "special setting " << specialSetting.Data() << std::endl;
        (AliPythia8::Instance())->ReadString( Form("StringFlav:etaPrimeSup = %f",etaPrimeSup) );
        (AliPythia8::Instance())->ReadString("331:onMode = off");
        (AliPythia8::Instance())->ReadString("331:onIfMatch = 211 -211 221");
        (AliPythia8::Instance())->ReadString("221:onMode = off");
        (AliPythia8::Instance())->ReadString("221:onIfMatch = 22 22");
    } else if (specialSetting.EqualTo("EtaPrimeNotSuppressed")) {
        std::cout << "special setting " << specialSetting.Data() << std::endl;
        (AliPythia8::Instance())->ReadString( Form("StringFlav:etaPrimeSup = %f",etaPrimeSup) );
    } else {
        std::cout << Form("specialSetting %s not found... returning\n", specialSetting.Data());
        return nullptr;
    }
    

    return gener;
}
