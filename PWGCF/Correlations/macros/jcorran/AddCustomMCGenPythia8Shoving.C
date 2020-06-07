AliGenerator* CreatePythia8Gen( TString lTune          = "pp",
                               Float_t e_cms       = 13000.,
                               Float_t fgAmplitude = 4.0
                               );

AliGenerator* AddCustomMCGenPythia8Shoving(   TString lTune          = "pp-shoving",
                                       Float_t e_cms       = 13000.,
                                       Float_t fgAmplitude = 4.0
                                       ) {
    
    
    gSystem->Load("liblhapdf");
    
    AliGenerator *genP  = NULL;
    genP                = CreatePythia8Gen(lTune.Data(), e_cms, fgAmplitude);
    
    return genP;
}

AliGenerator* CreatePythia8Gen( TString lTune,
                               Float_t e_cms,
                               Float_t fgAmplitude = 4.0
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
    std::cout << " Desired PYTHIA configuration: "<< lTune.Data()<< std::endl;
    std::cout << "*****************************************************************" << std::endl;
    // set process (MB)
    gener->SetProcess(kPyMbDefault);
    
    //Centre of mass energy
    gener->SetEnergyCMS(e_cms); // in GeV
    
    //random seed based on time
    (AliPythia8::Instance())->ReadString("Beams:idA = 2212");
    (AliPythia8::Instance())->ReadString("Beams:idB = 2212");
    (AliPythia8::Instance())->ReadString("PhaseSpace:pTHatMax = -1.0"); //this should be fixed in the constructor
    (AliPythia8::Instance())->ReadString("Main:timesAllowErrors = 50000");
    
    if ( lTune.EqualTo("pp") ){
        // Specific settings go here
        // default: do nothing, Monash 2013 will do its thing
        (AliPythia8::Instance())->ReadString(Form("Tune:pp = %d",14));
    }
    if ( lTune.EqualTo("pp-shoving") ){
        // This is a shoving setting acquired from Bierlich, Christian 1901.07447 and 1710.09725
        //===========================================================================
        (AliPythia8::Instance())->ReadString("SoftQCD:nonDiffractive = on");
        (AliPythia8::Instance())->ReadString("SoftQCD:singleDiffractive = on");
        (AliPythia8::Instance())->ReadString("SoftQCD:doubleDiffractive = on");
        (AliPythia8::Instance())->ReadString("Ropewalk:RopeHadronization = on"); //
        (AliPythia8::Instance())->ReadString("Ropewalk:doShoving = on"); //
        (AliPythia8::Instance())->ReadString("Ropewalk:doFlavour = off"); //
        (AliPythia8::Instance())->ReadString("Ropewalk:rCutOff = 10.0"); //
        (AliPythia8::Instance())->ReadString("Ropewalk:limitMom = on"); //
        (AliPythia8::Instance())->ReadString(Form("Ropewalk:pTcut = %lf",2.0)); //
        (AliPythia8::Instance())->ReadString("Ropewalk:r0 = 0.41"); //
        (AliPythia8::Instance())->ReadString("Ropewalk:m0 = 0.2"); //
        (AliPythia8::Instance())->ReadString(Form("Ropewalk:gAmplitude = %lf",fgAmplitude)); //
        (AliPythia8::Instance())->ReadString("Ropewalk:gExponent = 1.0"); //
        (AliPythia8::Instance())->ReadString("Ropewalk:deltat = 0.1"); //
        (AliPythia8::Instance())->ReadString("Ropewalk:tShove = 1.");  //
        (AliPythia8::Instance())->ReadString("Ropewalk:deltay = 0.1"); //
        (AliPythia8::Instance())->ReadString("Ropewalk:tInit = 1.5");  //
        // Enabling setting of vertex information.
        (AliPythia8::Instance())->ReadString("PartonVertex:setVertex = on"); //
        (AliPythia8::Instance())->ReadString("PartonVertex:protonRadius = 0.7"); //
        (AliPythia8::Instance())->ReadString("PartonVertex:emissionWidth = 0.1"); //
        //===========================================================================
    }
    if ( lTune.EqualTo("pp-default") ){
        (AliPythia8::Instance())->ReadString("SoftQCD:nonDiffractive = on");
        (AliPythia8::Instance())->ReadString("SoftQCD:singleDiffractive = on");
        (AliPythia8::Instance())->ReadString("SoftQCD:doubleDiffractive = on");
        (AliPythia8::Instance())->ReadString("PhaseSpace:pTHatMax = 0");
        (AliPythia8::Instance())->ReadString("PhaseSpace:pTHatMax = 13000");
        (AliPythia8::Instance())->ReadString("ParticleDecays:limitTau0 = On");
        (AliPythia8::Instance())->ReadString("ParticleDecays:tau0Max = 10.0");
    }
    
    (AliPythia8::Instance())->SetDecayLonglived();
    return gener;
}
