
AliGenerator* AddMCGenPythia8ShovingAndColorRopes( Float_t e_cms       = 13000.,
                                                   Float_t fgAmplitude = 5.0,
                                                   Float_t fpTcut      = 2.0,
                                                   TString flimitMom   = "on"
                                                 ) {
    
    gSystem->Load("liblhapdf");

    gSystem->Load("libpythia6");
    gSystem->Load("libEGPythia6");
    gSystem->Load("libAliPythia6");
    gSystem->Load("libpythia8");
    gSystem->Load("libAliPythia8");
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8243/xmldoc"));
    gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));
    
    
    AliPythia8 *pythia = AliPythia8::Instance();            // For ROOT6 needs to be created before AliGenPythiaPlus object, otherwise ending in "illegal instruction"
    AliGenPythiaPlus* gener = new AliGenPythiaPlus(pythia);

    std::cout << "*****************************************************************" << std::endl;
    std::cout << " Desired PYTHIA configuration: pp-shoving"<< std::endl;
    std::cout << "*****************************************************************" << std::endl;
    // set process (MB)
    gener->SetProcess(kPyMbDefault);
    
    //Centre of mass energy
    gener->SetEnergyCMS(e_cms); // in GeV
    
    //random seed based on time
    pythia->ReadString("Beams:idA = 2212");
    pythia->ReadString("Beams:idB = 2212");
    pythia->ReadString("PhaseSpace:pTHatMax = -1.0"); //this should be fixed in the constructor
    pythia->ReadString("Main:timesAllowErrors = 50000");

    // This is a shoving setting acquired from Valentina Zacolo and Christian Bierlich
    //===========================================================================
    pythia->ReadString("SoftQCD:nonDiffractive = on");
    pythia->ReadString("SoftQCD:singleDiffractive = on");
    pythia->ReadString("SoftQCD:doubleDiffractive = on");

    // ALICE primary particles
    pythia->ReadString("ParticleDecays:limitTau0 = on");
    pythia->ReadString("ParticleDecays:tau0Max = 10.");

    //QCD based 
    pythia->ReadString("MultiPartonInteractions:pT0Ref = 2.15");
    pythia->ReadString("BeamRemnants:remnantMode = 1"); //
    pythia->ReadString("BeamRemnants:saturation = 5"); //
    pythia->ReadString("ColourReconnection:mode = 1"); //
    pythia->ReadString("ColourReconnection:allowDoubleJunRem = off"); //
    pythia->ReadString("ColourReconnection:m0 = 0.3"); //
    pythia->ReadString("ColourReconnection:allowJunctions = on"); //
    pythia->ReadString("ColourReconnection:junctionCorrection = 1.2"); //
    pythia->ReadString("ColourReconnection:timeDilationMode = 2"); //
    pythia->ReadString("ColourReconnection:timeDilationPar = 0.18"); //
    pythia->ReadString("Ropewalk:RopeHadronization = on"); //

    //Shoving
    pythia->ReadString("Ropewalk:doShoving = on"); //
    pythia->ReadString("Ropewalk:tInit = 1.5");  //   
    pythia->ReadString("Ropewalk:deltat = 0.1"); //
    pythia->ReadString("Ropewalk:tShove = 1.");  //
    pythia->ReadString(Form("Ropewalk:gAmplitude = %lf",fgAmplitude)); //

    pythia->ReadString(Form("Ropewalk:limitMom = %s",flimitMom.Data())); //
    pythia->ReadString(Form("Ropewalk:pTcut = %lf",fpTcut)); //

    //Ropes
    pythia->ReadString("Ropewalk:doFlavour = on"); //
    pythia->ReadString("Ropewalk:r0 = 0.41"); //
    pythia->ReadString("Ropewalk:m0 = 0.2"); //
    pythia->ReadString("Ropewalk:beta = 0.1"); //
    pythia->ReadString("Ropewalk:rCutOff = 10.0"); //
    pythia->ReadString("Ropewalk:gExponent = 1.0"); //
    pythia->ReadString("Ropewalk:deltay = 0.1"); //
        
    // Enabling setting of vertex information.
    pythia->ReadString("PartonVertex:setVertex = on"); //
    pythia->ReadString("PartonVertex:protonRadius = 0.7"); //
    pythia->ReadString("PartonVertex:emissionWidth = 0.1"); //*/
    //===========================================================================
    
    return gener;
}
