AliGenerator* CreatePythia8Gen( Float_t e_cms,
                                Int_t tune,
                                Bool_t kCR,
                                Int_t kF,
                                Int_t kProcess,
                                Double_t ptHardMin,
                                Double_t ptHardMax,
                                Bool_t longlived
                            );

AliGenerator* Add_MCGenPythia8_TuneX(   Float_t e_cms       = 2760., 
                                        Int_t tune          = 5, 
                                        Bool_t kCR          = kTRUE, 
                                        Int_t kF            = 1, 
                                        Int_t kProcess      = 0, 
                                        Double_t ptHardMin  = 0, 
                                        Double_t ptHardMax  = 1.,
                                        Bool_t longlived = kFALSE

                                    ) {
    // Add Pythia 8 generator: 
    //    -kProcess=0  MB generation
    //    -kProcess=1  Jet production, pthard generation
    //    - Color reconnection = ON/OFF
    //    - Set k factor, default = 1; range of possible values in xmldoc/CouplingsAndScales.xml

    gSystem->Load("liblhapdf");
    
    AliGenerator *genP  = NULL;
    genP                = CreatePythia8Gen(e_cms, tune, kCR, kF, kProcess, ptHardMin, ptHardMax,longlived);
    
    return genP;
}

AliGenerator* CreatePythia8Gen( Float_t e_cms, 
                                Int_t tune, 
                                Bool_t kCR, 
                                Int_t kF, 
                                Int_t kProcess, 
                                Double_t ptHardMin, 
                                Double_t ptHardMax,
                                Bool_t longlived
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

    // setting tune
    std::cout << "setting tune: " << tune << std::endl;
    gener->SetTune(tune);
// 	(AliPythia8::Instance())->ReadString(Form("Tune:pp = %i", tune));//CR

    //random seed based on time
    (AliPythia8::Instance())->ReadString("Random:setSeed = on");
    (AliPythia8::Instance())->ReadString("Random:seed = 0");

    if(kCR)             
        (AliPythia8::Instance())->ReadString("ColourReconnection:reconnect = on");
    else
        (AliPythia8::Instance())->ReadString("ColourReconnection:reconnect = off");
        
    (AliPythia8::Instance())->ReadString(Form("MultipartonInteractions:kFactor = %i", kF));
    
    (AliPythia8::Instance())->SetDecayLonglived();

    return gener;
}
