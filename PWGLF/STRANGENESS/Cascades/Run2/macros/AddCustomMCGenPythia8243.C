#include "AliPythia8.h"
#include "AliGenPythiaPlus.h"

AliGenerator* CreatePythia8Gen( TString lTune          = "pp",
                               Float_t e_cms       = 13000.
                               );

AliGenerator* AddCustomMCGenPythia8243(   TString lTune          = "pp",
                                       Float_t e_cms       = 13000.
                                       
                                       ) {
    // Add Pythia 8 generator:
    //    -kProcess=0  MB generation
    //    -kProcess=1  Jet production, pthard generation
    //    - Color reconnection = ON/OFF
    //    - Set k factor, default = 1; range of possible values in xmldoc/CouplingsAndScales.xml
    
    gSystem->Load("liblhapdf");
    
    AliGenerator *genP  = NULL;
    genP                = CreatePythia8Gen(lTune.Data(), e_cms);
    
    return genP;
}

AliGenerator* CreatePythia8Gen( TString lTune,
                               Float_t e_cms
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
    if ( lTune.EqualTo("pp-experimental") ){
        // This is me testing a few things
        TRandom3 gRand3;
        gRand3.SetSeed(0);
        Float_t lRandpT0Ref = 2.0+(((Float_t)((Int_t)(10.*gRand3.Uniform())))/10.);
        (AliPythia8::Instance())->ReadString(Form("MultipartonInteractions:pT0Ref = %.2f",lRandpT0Ref));
        std::cout << " Random pT0Ref: "<< lRandpT0Ref << std::endl;
    }
    if ( lTune.EqualTo("pp-moreqcd") ){
        std::cout << " Setting pp-moreqcd parameters..." << std::endl;
        //Paper reference: https://arxiv.org/pdf/1505.01681.pdf ("mode 0")
        //===========================================================================
        (AliPythia8::Instance())->ReadString("StringPT:sigma = 0.335");
        (AliPythia8::Instance())->ReadString("StringZ:aLund = 0.36");
        (AliPythia8::Instance())->ReadString("StringZ:bLund = 0.56");
        (AliPythia8::Instance())->ReadString("StringFlav:probQQtoQ = 0.078");
        (AliPythia8::Instance())->ReadString("StringFlav:ProbStoUD = 0.2");
        (AliPythia8::Instance())->ReadString("StringFlav:probQQ1toQQ0join = 0.0275,0.0275,0.0275,0.0275");
        //===========================================================================
        (AliPythia8::Instance())->ReadString("MultiPartonInteractions:pT0Ref = 2.12");
        //===========================================================================
        (AliPythia8::Instance())->ReadString("BeamRemnants:remnantMode = 1");
        (AliPythia8::Instance())->ReadString("BeamRemnants:saturation = 5");
        //===========================================================================
        (AliPythia8::Instance())->ReadString("ColourReconnection:mode = 1");
        (AliPythia8::Instance())->ReadString("ColourReconnection:allowDoubleJunRem = off");
        (AliPythia8::Instance())->ReadString("ColourReconnection:m0 = 2.9");
        (AliPythia8::Instance())->ReadString("ColourReconnection:allowJunctions = on");
        (AliPythia8::Instance())->ReadString("ColourReconnection:junctionCorrection = 1.43");
        (AliPythia8::Instance())->ReadString("ColourReconnection:timeDilationMode = 0");
        //===========================================================================
    }
    if ( lTune.EqualTo("pp-nocr") ){
        (AliPythia8::Instance())->ReadString("ColourReconnection:reconnect = off");
        (AliPythia8::Instance())->ReadString("PartonLevel:earlyResDec = off");
        (AliPythia8::Instance())->ReadString("MultipartonInteractions:pT0Ref = 2.30");
    }
    if ( lTune.EqualTo("pp-ropes") ){
        // This is a ropes setting acquired from Peter Christiansen on the Lund thing
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
        //===========================================================================
    }
    if ( lTune.EqualTo("pp-shoving") ){
        // Enabling flavour ropes, setting model parameters.
        // The model is still untuned. These parameter values
        // are choosen for illustrative purposes.
        // This is a shoving setting acquired from Peter Christiansen on the Lund thing
        //===========================================================================
        (AliPythia8::Instance())->ReadString("Ropewalk:RopeHadronization = on");
        (AliPythia8::Instance())->ReadString("Ropewalk:doShoving = on");
        (AliPythia8::Instance())->ReadString("Ropewalk:doFlavour = off");
        (AliPythia8::Instance())->ReadString("Ropewalk:rCutOff = 10.0");
        (AliPythia8::Instance())->ReadString("Ropewalk:limitMom = on");
        (AliPythia8::Instance())->ReadString("Ropewalk:pTcut = 2.0");
        (AliPythia8::Instance())->ReadString("Ropewalk:r0 = 0.41");
        (AliPythia8::Instance())->ReadString("Ropewalk:m0 = 0.2");
        (AliPythia8::Instance())->ReadString("Ropewalk:gAmplitude = 10.0");
        (AliPythia8::Instance())->ReadString("Ropewalk:gExponent = 1.0");
        (AliPythia8::Instance())->ReadString("Ropewalk:deltat = 0.1");
        (AliPythia8::Instance())->ReadString("Ropewalk:tShove = 1.");
        (AliPythia8::Instance())->ReadString("Ropewalk:deltay = 0.1");
        (AliPythia8::Instance())->ReadString("Ropewalk:tInit = 1.5");
        //===========================================================================
        // Enabling setting of vertex information.
        (AliPythia8::Instance())->ReadString("PartonVertex:setVertex = on");
        (AliPythia8::Instance())->ReadString("PartonVertex:protonRadius = 0.7");
        (AliPythia8::Instance())->ReadString("PartonVertex:emissionWidth = 0.1");
        //===========================================================================
    }
    
    (AliPythia8::Instance())->SetDecayLonglived();
    return gener;
}
