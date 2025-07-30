AliGenerator* CreatePythia8Gen( Float_t e_cms,
                                Int_t tune,
                                Bool_t kCR,
                                Int_t kF,
                                Int_t kProcess,
                                Double_t ptHardMin,
                                Double_t ptHardMax,
                                Bool_t longlived,
                                TString specialTune
                            );

AliGenerator* Add_MCGenPythia8_TuneX(   Float_t e_cms       = 2760., 
                                        Int_t tune          = 5, 
                                        Bool_t kCR          = kTRUE, 
                                        Int_t kF            = 1, 
                                        Int_t kProcess      = 0, 
                                        Double_t ptHardMin  = 0, 
                                        Double_t ptHardMax  = 1.,
                                        Bool_t longlived    = kTRUE,
                                        TString specialTune = ""

                                    ) {
    // Add Pythia 8 generator: 
    //    -kProcess=0  MB generation
    //    -kProcess=1  Jet production, pthard generation
    //    - Color reconnection = ON/OFF
    //    - Set k factor, default = 1; range of possible values in xmldoc/CouplingsAndScales.xml

    gSystem->Load("liblhapdf");
    
    AliGenerator *genP  = NULL;
    genP                = CreatePythia8Gen(e_cms, tune, kCR, kF, kProcess, ptHardMin, ptHardMax,longlived, specialTune);
    
    return genP;
}

AliGenerator* CreatePythia8Gen( Float_t e_cms, 
                                Int_t tune, 
                                Bool_t kCR, 
                                Int_t kF, 
                                Int_t kProcess, 
                                Double_t ptHardMin, 
                                Double_t ptHardMax,
                                Bool_t longlived,
                                TString specialTune
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
    std::cout << "Process: "<< kProcess << "\t Color reconnection: "<< kCR << "\t Tune: " << tune << "\t kFactor: "<< kF <<  std::endl;
    if(!(specialTune.EqualTo("") || specialTune.EqualTo("Monash"))){
        std::cout << ">>> Running with special settings (not standard Monash tune): " << specialTune << std::endl;
    }
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

    if(kProcess==2){
        std::cout << "running direct photon mode"<< std::endl;
        gener->SetProcess(kPyDirectGamma);
        if(ptHardMin > 0.){ 
            std::cout << "Setting pTHardMin: "<< ptHardMin << "\t ptHardMax: "<< ptHardMax <<  std::endl;
            gener->SetPtHard(ptHardMin,ptHardMax);
        }
    }

    if(kProcess==3){
        std::cout << "running direct photon mode + min bias"<< std::endl;
        gener->SetProcess(kPyMbWithDirectPhoton);
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

    if(specialTune.EqualTo("Junctions_Mode2")){
        //Paper reference: https://arxiv.org/pdf/1505.01681.pdf ("mode 0")
        //===========================================================================
        cout << "running with junctions\n";
        (AliPythia8::Instance())->ReadString("StringPT:sigma = 0.335");
        (AliPythia8::Instance())->ReadString("StringZ:aLund = 0.36");
        (AliPythia8::Instance())->ReadString("StringZ:bLund = 0.56");
        (AliPythia8::Instance())->ReadString("StringFlav:probQQtoQ = 0.078");
        (AliPythia8::Instance())->ReadString("StringFlav:ProbStoUD = 0.2");
        (AliPythia8::Instance())->ReadString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        //===========================================================================
        (AliPythia8::Instance())->ReadString("MultiPartonInteractions:pT0Ref = 2.15");
        //===========================================================================
        (AliPythia8::Instance())->ReadString("BeamRemnants:remnantMode = 1");
        (AliPythia8::Instance())->ReadString("BeamRemnants:saturation = 5");
        //===========================================================================
        (AliPythia8::Instance())->ReadString("ColourReconnection:mode = 1");
        (AliPythia8::Instance())->ReadString("ColourReconnection:allowDoubleJunRem = off");
        (AliPythia8::Instance())->ReadString("ColourReconnection:m0 = 0.3");
        (AliPythia8::Instance())->ReadString("ColourReconnection:allowJunctions = on");
        (AliPythia8::Instance())->ReadString("ColourReconnection:junctionCorrection = 1.20");
        (AliPythia8::Instance())->ReadString("ColourReconnection:timeDilationMode = 2");
        (AliPythia8::Instance())->ReadString("ColourReconnection:timeDilationPar = 0.18");

    //==============================================================================
    // Standard Monash tune
    } else if (specialTune.EqualTo("Ropes")) {
        //===========================================================================
        (AliPythia8::Instance())->ReadString("MultiPartonInteractions:pT0Ref = 2.15");
        //===========================================================================
        (AliPythia8::Instance())->ReadString("BeamRemnants:remnantMode = 1");
        (AliPythia8::Instance())->ReadString("BeamRemnants:saturation = 5");
        //===========================================================================
        (AliPythia8::Instance())->ReadString("ColourReconnection:mode = 1");
        (AliPythia8::Instance())->ReadString("ColourReconnection:allowDoubleJunRem = off");
        (AliPythia8::Instance())->ReadString("ColourReconnection:m0 = 0.3");
        (AliPythia8::Instance())->ReadString("ColourReconnection:allowJunctions = on");
        (AliPythia8::Instance())->ReadString("ColourReconnection:junctionCorrection = 1.2");
        (AliPythia8::Instance())->ReadString("ColourReconnection:timeDilationMode = 2");
        (AliPythia8::Instance())->ReadString("ColourReconnection:timeDilationPar = 0.18");
        //===========================================================================
        (AliPythia8::Instance())->ReadString("Ropewalk:RopeHadronization = on");
        (AliPythia8::Instance())->ReadString("Ropewalk:doShoving = on");
        (AliPythia8::Instance())->ReadString("Ropewalk:tInit = 1.5"); // Propagation time
        (AliPythia8::Instance())->ReadString("Ropewalk:deltat = 0.05");
        (AliPythia8::Instance())->ReadString("Ropewalk:tShove 0.1");
        (AliPythia8::Instance())->ReadString("Ropewalk:gAmplitude = 0."); // Set shoving strength to 0 explicitly
        (AliPythia8::Instance())->ReadString("Ropewalk:doFlavour = on");
        (AliPythia8::Instance())->ReadString("Ropewalk:r0 = 0.5");
        (AliPythia8::Instance())->ReadString("Ropewalk:m0 = 0.2");
        (AliPythia8::Instance())->ReadString("Ropewalk:beta = 0.1");
        //===========================================================================
        // Enabling setting of vertex information.
        (AliPythia8::Instance())->ReadString("PartonVertex:setVertex = on");
        (AliPythia8::Instance())->ReadString("PartonVertex:protonRadius = 0.7");
        (AliPythia8::Instance())->ReadString("PartonVertex:emissionWidth = 0.1");
    } else if (specialTune.EqualTo("Shoving")){
        // This is a shoving setting acquired from Bierlich, Christian 1901.07447 and 1710.09725
        //===========================================================================
        (AliPythia8::Instance())->ReadString("SoftQCD:nonDiffractive = on");
        (AliPythia8::Instance())->ReadString("SoftQCD:singleDiffractive = on");
        (AliPythia8::Instance())->ReadString("SoftQCD:doubleDiffractive = on");
        (AliPythia8::Instance())->ReadString("Ropewalk:RopeHadronization = on");
        (AliPythia8::Instance())->ReadString("Ropewalk:doShoving = on");
        (AliPythia8::Instance())->ReadString("Ropewalk:doFlavour = off");
        (AliPythia8::Instance())->ReadString("Ropewalk:rCutOff = 10.0");
        (AliPythia8::Instance())->ReadString("Ropewalk:limitMom = on");
        (AliPythia8::Instance())->ReadString("Ropewalk:pTcut = 2");
        (AliPythia8::Instance())->ReadString("Ropewalk:r0 = 0.41");
        (AliPythia8::Instance())->ReadString("Ropewalk:m0 = 0.2");
        (AliPythia8::Instance())->ReadString("Ropewalk:gAmplitude = 10");
        (AliPythia8::Instance())->ReadString("Ropewalk:gExponent = 1.0");
        (AliPythia8::Instance())->ReadString("Ropewalk:deltat = 0.1");
        (AliPythia8::Instance())->ReadString("Ropewalk:tShove = 1.");
        (AliPythia8::Instance())->ReadString("Ropewalk:deltay = 0.1");
        (AliPythia8::Instance())->ReadString("Ropewalk:tInit = 1.5");
        // Enabling setting of vertex information.
        (AliPythia8::Instance())->ReadString("PartonVertex:setVertex = on");
        (AliPythia8::Instance())->ReadString("PartonVertex:protonRadius = 0.7");
        (AliPythia8::Instance())->ReadString("PartonVertex:emissionWidth = 0.1");
    } else if (specialTune.EqualTo("JetSetting")){
        // For comparison to Pythia+POWHEG from Markus Fasel and Hadi Hassan
        std::cout << "Running with specialTune JetSetting" << std::endl;
        (AliPythia8::Instance())->ReadString("Next:numberShowInfo = 1");
        (AliPythia8::Instance())->ReadString("Next:numberShowProcess = 1");
        (AliPythia8::Instance())->ReadString("Next:numberShowEvent = 1");
        (AliPythia8::Instance())->ReadString("Main:timesAllowErrors = 10");

        (AliPythia8::Instance())->ReadString("Init:showChangedSettings = on");
        (AliPythia8::Instance())->ReadString("Init:showChangedParticleData = off");

        (AliPythia8::Instance())->ReadString("PartonLevel:MPI = on");
        // Switch On Pi0 Decay
        (AliPythia8::Instance())->ReadString("111:mayDecay  = on");
        (AliPythia8::Instance())->ReadString("310:mayDecay  = off");
        (AliPythia8::Instance())->ReadString("3122:mayDecay = off");
        (AliPythia8::Instance())->ReadString("3112:mayDecay = off");
        (AliPythia8::Instance())->ReadString("3212:mayDecay = off");
        (AliPythia8::Instance())->ReadString("3222:mayDecay = off");
        (AliPythia8::Instance())->ReadString("3312:mayDecay = off");
        (AliPythia8::Instance())->ReadString("3322:mayDecay = off");
        (AliPythia8::Instance())->ReadString("3334:mayDecay = off");

        // Tune Parameters
        (AliPythia8::Instance())->ReadString("Tune:preferLHAPDF = 2");
        (AliPythia8::Instance())->ReadString("Tune:pp = 21"); // Tune A14 NNPDF3.2

        // PDF Selection
        (AliPythia8::Instance())->ReadString("PDF:pSet = 14");

    } else if(specialTune.EqualTo("MonashWithFSR") || specialTune.EqualTo("MonashWithoutFSR")){
        // setting tune
        std::cout << "setting tune: " << tune << std::endl;
        gener->SetTune(tune);
        if(kCR){            
            (AliPythia8::Instance())->ReadString("ColourReconnection:reconnect = on");
        } else {
            (AliPythia8::Instance())->ReadString("ColourReconnection:reconnect = off");
        }
        (AliPythia8::Instance())->ReadString(Form("MultipartonInteractions:kFactor = %i", kF));

        // Setting FSR on
        if(specialTune.EqualTo("MonashWithFSR")) {
            (AliPythia8::Instance())->ReadString("TimeShower:QEDshowerByQ = on");
        } else {
            (AliPythia8::Instance())->ReadString("TimeShower:QEDshowerByQ = off");
        }

    } else if(specialTune.EqualTo("") || specialTune.EqualTo("Monash")){
        // setting tune
        std::cout << "setting tune: " << tune << std::endl;
        gener->SetTune(tune);
        if(kCR){            
            (AliPythia8::Instance())->ReadString("ColourReconnection:reconnect = on");
        } else {
            (AliPythia8::Instance())->ReadString("ColourReconnection:reconnect = off");
        }
            
        (AliPythia8::Instance())->ReadString(Form("MultipartonInteractions:kFactor = %i", kF));
    } else {
        std::cout << Form("specialTune %s not found... returning\n", specialTune.Data());
        return nullptr;
    }
    
    (AliPythia8::Instance())->SetDecayLonglived(longlived);

    return gener;
}
