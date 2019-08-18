R__LOAD_LIBRARY(liblhapdf)
R__LOAD_LIBRARY(libpythia6)

#include "AliGenerator.h"
#include "AliGenPythia.h"

AliGenerator* AddMCGenPythia8(TString lSystem = "pp", TString lConfig = "", Float_t e_cms = 13000.)
{
    AliGenerator *genP = NULL;
    if( lSystem.EqualTo("pp") )
        genP = CreatePythiaMonash(e_cms);
    if( lSystem.EqualTo("pp-nocr") )
        genP = CreatePythiaMonashNoCR(e_cms);
    if( lSystem.EqualTo("pp-moreqcd") )
        genP = CreatePythiaMonashMoreQCD(e_cms);
    if( lSystem.EqualTo("pp-ropes") )
        genP = CreatePythiaMonashRopes(e_cms);
    
    return gener;
}

AliGenerator* CreatePythiaMonash(Float_t e_cms)
{
    
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
    gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));
    
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
    
    return gener;
}

AliGenerator* CreatePythiaMonashMoreQCD(Float_t e_cms)
{
    
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
    gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));
    
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
    
    return gener;
}

AliGenerator* CreatePythiaMonashNoCR(Float_t e_cms)
{
    
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
    gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));
    
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
    
    return gener;
}

AliGenerator* CreatePythiaMonashRopes(Float_t e_cms)
{
    
    gSystem->Setenv("PYTHIA8DATA", gSystem->ExpandPathName("$ALICE_ROOT/PYTHIA8/pythia8/xmldoc"));
    gSystem->Setenv("LHAPDF",      gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF"));
    gSystem->Setenv("LHAPATH",     gSystem->ExpandPathName("$ALICE_ROOT/LHAPDF/PDFsets"));
    
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
    // Specific settings go here: color ropes - CHECK ME !
    (AliPythia8::Instance())->ReadString("Beams:idA = 2212");
    (AliPythia8::Instance())->ReadString("Beams:idB = 2212");
    (AliPythia8::Instance())->ReadString(Form("Tune:pp = %d",14));
    (AliPythia8::Instance())->ReadString("Ropewalk:RopeHadronization = on"); //! Rope Hadronization framework
    (AliPythia8::Instance())->ReadString("Ropewalk:doShoving = off"); //! Enable the string shoving mechanism
    (AliPythia8::Instance())->ReadString("Ropewalk:doFlavour = on"); //! Enable the flavour ropes mechanism
    (AliPythia8::Instance())->ReadString("Ropewalk:r0 = 0.5"); //! The transverse radius of a string, in units of fm
    (AliPythia8::Instance())->ReadString("Ropewalk:m0 = 0.2"); //! Imposed lower mass cutoff
    (AliPythia8::Instance())->ReadString("Ropewalk:beta = 0.1"); //! This parameter controls how large a fraction of the parameter will scale with string tension
    (AliPythia8::Instance())->ReadString("PartonVertex:setVertex = on"); //! Enabling setting of vertex information.
    //============================================================
    
    return gener;
}
