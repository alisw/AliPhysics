AliGenerator* CreatePythia8Gen( TString lTune          = "gluon",
                                Float_t e_cms       = 13000.,
                                Float_t mHatMin     = 5.
                                );

AliGenerator* AddCustomMCGenPythia8QuarkGluonJet(   TString lTune          = "gluon",
                                       Float_t e_cms       = 13000.,
                                       Float_t mHatMin     = 5.
                                       ) {
    
    
    gSystem->Load("liblhapdf");
    
    AliGenerator *genP  = NULL;
    genP                = CreatePythia8Gen(lTune.Data(), e_cms,mHatMin);
    
    return genP;
}

AliGenerator* CreatePythia8Gen( TString lTune,
                               Float_t e_cms, 
                               Float_t mHatMin) {
    
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
    pythia->ReadString("Beams:idA = 2212");
    pythia->ReadString("Beams:idB = 2212");
    pythia->ReadString("PhaseSpace:pTHatMax = -1.0"); //this should be fixed in the constructor
    pythia->ReadString(Form("PhaseSpace:mHatMin=%f ",mHatMin));
    pythia->ReadString("Main:timesAllowErrors = 50000");
    if( lTune.EqualTo("gluon") ) pythia->ReadString("HardQCD:gg2gg = on");
    if( lTune.EqualTo("quark") ) pythia->ReadString("HardQCD:gg2qqbar = on");
    
    pythia->SetDecayLonglived();
    return gener;
}
