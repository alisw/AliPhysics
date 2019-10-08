R__LOAD_LIBRARY(liblhapdf)
R__LOAD_LIBRARY(libpythia6)

#include "AliGenerator.h"
#include "AliGenPythia.h"

AliGenerator* CreatePythiaMonash(Float_t e_cms);
AliGenerator* CreatePythiaMonashExperimental(Float_t e_cms);
AliGenerator* CreatePythiaMonashMoreQCD(Float_t e_cms);
AliGenerator* CreatePythiaMonashRopes(Float_t e_cms);
AliGenerator* CreatePythiaMonashShoving(Float_t e_cms);

AliGenerator* AddCustomMCGenPythia8(TString lSystem = "pp", TString lConfig = "", Float_t e_cms = 13000.)
{
    gSystem->Load("libpythia6");
    gSystem->Load("libEGPythia6");
    gSystem->Load("libAliPythia6");
    gSystem->Load("libpythia8");
    gSystem->Load("libAliPythia8");
    gSystem->Load("liblhapdf");
    AliGenerator *genP = NULL;
    if( lSystem.EqualTo("pp") )
        genP = CreatePythiaMonash(e_cms);
    if( lSystem.EqualTo("pp-experimental") )
        genP = CreatePythiaMonashExperimental(e_cms);
    if( lSystem.EqualTo("pp-nocr") )
        genP = CreatePythiaMonashNoCR(e_cms);
    if( lSystem.EqualTo("pp-moreqcd") )
        genP = CreatePythiaMonashMoreQCD(e_cms);
    if( lSystem.EqualTo("pp-ropes") )
        genP = CreatePythiaMonashRopes(e_cms);
    if( lSystem.EqualTo("pp-shoving") )
        genP = CreatePythiaMonashShoving(e_cms);
    
    return genP;
}


AliGenerator* CreatePythiaMonash(Float_t e_cms)
{
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
    gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));
    
    AliPythia8 *pythia = AliPythia8::Instance();            // For ROOT6 needs to be created before AliGenPythiaPlus object, otherwise ending in "illegal instruction"
    AliGenPythiaPlus* gener = new AliGenPythiaPlus(AliPythia8::Instance());
    
    //Standard setting setup
    //Set process (min-bias)
    gener->SetProcess(kPyMbDefault);
    //Centre of mass energy
    gener->SetEnergyCMS(e_cms); // in GeV
    //random seed based on time
    AliPythia8::Instance()->ReadString("Random:setSeed = on");
    AliPythia8::Instance()->ReadString("Random:seed = 0");
    
    //============================================================
    // Specific settings go here
    (AliPythia8::Instance())->ReadString("Beams:idA = 2212");
    (AliPythia8::Instance())->ReadString("Beams:idB = 2212");
    (AliPythia8::Instance())->ReadString(Form("Tune:pp = %d",14));
    //============================================================
    (AliPythia8::Instance())->SetDecayLonglived();

    return gener;
}

AliGenerator* CreatePythiaMonashExperimental(Float_t e_cms)
{
    gSystem->Load("libpythia6");
    gSystem->Load("libEGPythia6");
    gSystem->Load("libAliPythia6");
    gSystem->Load("libpythia8");
    gSystem->Load("libAliPythia8");
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
    gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));
    
    AliPythia8 *pythia = AliPythia8::Instance();            // For ROOT6 needs to be created before AliGenPythiaPlus object, otherwise ending in "illegal instruction"
    AliGenPythiaPlus* gener = new AliGenPythiaPlus(AliPythia8::Instance());
    
    //Standard setting setup
    //Set process (min-bias)
    gener->SetProcess(kPyMbDefault);
    //Centre of mass energy
    gener->SetEnergyCMS(e_cms); // in GeV
    //random seed based on time
    AliPythia8::Instance()->ReadString("Random:setSeed = on");
    AliPythia8::Instance()->ReadString("Random:seed = 0");
    
    //============================================================
    // Specific settings go here
    (AliPythia8::Instance())->ReadString("Beams:idA = 2212");
    (AliPythia8::Instance())->ReadString("Beams:idB = 2212");
    (AliPythia8::Instance())->ReadString(Form("Tune:pp = %d",14));
    (AliPythia8::Instance())->ReadString(Form("MultipartonInteractions:pT0Ref = %.2f",2.0+gRandom->Uniform()));
    //============================================================
    (AliPythia8::Instance())->SetDecayLonglived();

    return gener;
}

AliGenerator* CreatePythiaMonashMoreQCD(Float_t e_cms)
{
    gSystem->Load("libpythia6");
    gSystem->Load("libEGPythia6");
    gSystem->Load("libAliPythia6");
    gSystem->Load("libpythia8");
    gSystem->Load("libAliPythia8");
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
    gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));
    
    AliPythia8 *pythia = AliPythia8::Instance();            // For ROOT6 needs to be created before AliGenPythiaPlus object, otherwise ending in "illegal instruction"
    AliGenPythiaPlus* gener = new AliGenPythiaPlus(AliPythia8::Instance());
    
    //Standard setting setup
    //Set process (min-bias)
    gener->SetProcess(kPyMbDefault);
    //Centre of mass energy
    gener->SetEnergyCMS(e_cms); // in GeV
    //random seed based on time
    AliPythia8::Instance()->ReadString("Random:setSeed = on");
    AliPythia8::Instance()->ReadString("Random:seed = 0");
    
    //============================================================
    // Specific settings go here: More-QCD scheme / extra junctions
    (AliPythia8::Instance())->ReadString("Beams:idA = 2212");
    (AliPythia8::Instance())->ReadString("Beams:idB = 2212");
    (AliPythia8::Instance())->ReadString(Form("Tune:pp = %d",14));
    (AliPythia8::Instance())->ReadString("ColourReconnection:reconnect = on");
    (AliPythia8::Instance())->ReadString("ColourReconnection:mode = 1");
    (AliPythia8::Instance())->ReadString("ColourReconnection:timeDilationPar = 0.18");
    (AliPythia8::Instance())->ReadString("ColourReconnection:junctionCorrection = 1.2");
    (AliPythia8::Instance())->ReadString("MultiPartonInteractions:pT0Ref = 2.15");
    //============================================================
    (AliPythia8::Instance())->SetDecayLonglived();

    return gener;
}

AliGenerator* CreatePythiaMonashNoCR(Float_t e_cms)
{
    gSystem->Load("libpythia6");
    gSystem->Load("libEGPythia6");
    gSystem->Load("libAliPythia6");
    gSystem->Load("libpythia8");
    gSystem->Load("libAliPythia8");
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
    gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));
    
    AliPythia8 *pythia = AliPythia8::Instance();            // For ROOT6 needs to be created before AliGenPythiaPlus object, otherwise ending in "illegal instruction"
    AliGenPythiaPlus* gener = new AliGenPythiaPlus(AliPythia8::Instance());
    
    //Standard setting setup
    //Set process (min-bias)
    gener->SetProcess(kPyMbDefault);
    //Centre of mass energy
    gener->SetEnergyCMS(e_cms); // in GeV
    //random seed based on time
    AliPythia8::Instance()->ReadString("Random:setSeed = on");
    AliPythia8::Instance()->ReadString("Random:seed = 0");
    
    //============================================================
    // Specific settings go here: no color reconnection
    // (warning: NOT RETUNED)
    (AliPythia8::Instance())->ReadString("Beams:idA = 2212");
    (AliPythia8::Instance())->ReadString("Beams:idB = 2212");
    (AliPythia8::Instance())->ReadString(Form("Tune:pp = %d",14));
    (AliPythia8::Instance())->ReadString("ColourReconnection:reconnect = off");
    (AliPythia8::Instance())->ReadString("PartonLevel:earlyResDec = off");
    (AliPythia8::Instance())->ReadString("MultipartonInteractions:pT0Ref = 2.30");
    //============================================================
    (AliPythia8::Instance())->SetDecayLonglived();

    return gener;
}

AliGenerator* CreatePythiaMonashRopes(Float_t e_cms)
{
    gSystem->Load("libpythia6");
    gSystem->Load("libEGPythia6");
    gSystem->Load("libAliPythia6");
    gSystem->Load("libpythia8");
    gSystem->Load("libAliPythia8");
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
    gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));

    AliPythia8 *pythia = AliPythia8::Instance();            // For ROOT6 needs to be created before AliGenPythiaPlus object, otherwise ending in "illegal instruction"
    AliGenPythiaPlus* gener = new AliGenPythiaPlus(AliPythia8::Instance());
    
    //Standard setting setup
    //Set process (min-bias)
    gener->SetProcess(kPyMbDefault);
    //Centre of mass energy
    gener->SetEnergyCMS(e_cms); // in GeV
    //random seed based on time
    AliPythia8::Instance()->ReadString("Random:setSeed = on");
    AliPythia8::Instance()->ReadString("Random:seed = 0");
    
    //============================================================
    // Specific settings go here: color ropes
    // updated 26th August 2019
    // brand new configuration from Christian  Bierlich
    //============================================================
    (AliPythia8::Instance())->ReadString("Beams:idA = 2212");
    (AliPythia8::Instance())->ReadString("Beams:idB = 2212");
    (AliPythia8::Instance())->ReadString(Form("Tune:pp = %d",14)); //should be default, but ok

    // QCD based CR
    (AliPythia8::Instance())->ReadString("MultiPartonInteractions:pT0Ref = 2.15");
    (AliPythia8::Instance())->ReadString("BeamRemnants:remnantMode = 1");
    (AliPythia8::Instance())->ReadString("BeamRemnants:saturation = 5");
    (AliPythia8::Instance())->ReadString("ColourReconnection:mode = 1");
    (AliPythia8::Instance())->ReadString("ColourReconnection:allowDoubleJunRem = off");
    (AliPythia8::Instance())->ReadString("ColourReconnection:m0 = 0.3");
    (AliPythia8::Instance())->ReadString("ColourReconnection:allowJunctions = on");
    (AliPythia8::Instance())->ReadString("ColourReconnection:junctionCorrection = 1.2");
    (AliPythia8::Instance())->ReadString("ColourReconnection:timeDilationMode = 2");
    (AliPythia8::Instance())->ReadString("ColourReconnection:timeDilationPar = 0.18");
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
    
    // Enabling setting of vertex information.
    (AliPythia8::Instance())->ReadString("PartonVertex:setVertex = on");
    (AliPythia8::Instance())->ReadString("PartonVertex:protonRadius = 0.7");
    (AliPythia8::Instance())->ReadString("PartonVertex:emissionWidth = 0.1");
    //============================================================
    (AliPythia8::Instance())->SetDecayLonglived();

    return gener;
}

AliGenerator* CreatePythiaMonashShoving(Float_t e_cms)
{
    gSystem->Load("libpythia6");
    gSystem->Load("libEGPythia6");
    gSystem->Load("libAliPythia6");
    gSystem->Load("libpythia8");
    gSystem->Load("libAliPythia8");
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
    gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));
    
    AliPythia8 *pythia = AliPythia8::Instance();            // For ROOT6 needs to be created before AliGenPythiaPlus object, otherwise ending in "illegal instruction"
    AliGenPythiaPlus* gener = new AliGenPythiaPlus(AliPythia8::Instance());
    
    //Standard setting setup
    //Set process (min-bias)
    gener->SetProcess(kPyMbDefault);
    //Centre of mass energy
    gener->SetEnergyCMS(e_cms); // in GeV
    //random seed based on time
    AliPythia8::Instance()->ReadString("Random:setSeed = on");
    AliPythia8::Instance()->ReadString("Random:seed = 0");
    
    //============================================================
    // Specific settings go here: shoving
    // updated 26th August 2019
    // brand new configuration from Christian  Bierlich
    //============================================================
    (AliPythia8::Instance())->ReadString("Beams:idA = 2212");
    (AliPythia8::Instance())->ReadString("Beams:idB = 2212");
    (AliPythia8::Instance())->ReadString(Form("Tune:pp = %d",14)); //should be default, but ok
    
    // Enabling flavour ropes, setting model parameters.
    // The model is still untuned. These parameter values
    // are choosen for illustrative purposes.
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
    
    // Enabling setting of vertex information.
    (AliPythia8::Instance())->ReadString("PartonVertex:setVertex = on");
    (AliPythia8::Instance())->ReadString("PartonVertex:protonRadius = 0.7");
    (AliPythia8::Instance())->ReadString("PartonVertex:emissionWidth = 0.1");
    //============================================================
    (AliPythia8::Instance())->SetDecayLonglived();

    return gener;
}

