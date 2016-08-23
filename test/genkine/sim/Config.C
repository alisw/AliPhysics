// Configuration of simulation

void Config()
{
  // The common part contains the description of geometry,
  // magnetic field and generator. The initialization of
  // the random generator seed is also there

  gROOT->Macro("../commonConfig.C");


  //=======================================================================
  // Create transport
  new     TGeant3TGeo("C++ Interface to Geant3");

  //=======================================================================
  // Set External decayer
  AliDecayer *decayer = new AliDecayerPythia();
  decayer->SetForceDecay(kAll);
  decayer->Init();

  //forbid some decays
  AliPythia * py= AliPythia::Instance();
  py->SetMDME(737,1,0); //forbid D*+->D+ + pi0
  py->SetMDME(738,1,0); //forbid D*+->D+ + gamma
  
  for(Int_t d=747; d<=762; d++){ 
    py->SetMDME(d,1,0);
  }
  
  for(Int_t d=764; d<=807; d++){ 
    py->SetMDME(d,1,0);
  }

  TVirtualMC * vmc = TVirtualMC::GetMC();
  vmc->SetExternalDecayer(decayer);
  
  //=======================================================================
  // ************* STEERING parameters FOR ALICE SIMULATION **************
  // --- Specify event type to be tracked through the ALICE setup
  // --- All positions are in cm, angles in degrees, and P and E in GeV

  vmc->SetProcess("DCAY",1);
  vmc->SetProcess("PAIR",1);
  vmc->SetProcess("COMP",1);
  vmc->SetProcess("PHOT",1);
  vmc->SetProcess("PFIS",0);
  vmc->SetProcess("DRAY",0);
  vmc->SetProcess("ANNI",1);
  vmc->SetProcess("BREM",1);
  vmc->SetProcess("MUNU",1);
  vmc->SetProcess("CKOV",1);
  vmc->SetProcess("HADR",1);
  vmc->SetProcess("LOSS",2);
  vmc->SetProcess("MULS",1);
  vmc->SetProcess("RAYL",1);
  
  Float_t cut = 1.e-3;        // 1MeV cut by default
  Float_t tofmax = 1.e10;
  
  vmc->SetCut("CUTGAM", cut);
  vmc->SetCut("CUTELE", cut);
  vmc->SetCut("CUTNEU", cut);
  vmc->SetCut("CUTHAD", cut);
  vmc->SetCut("CUTMUO", cut);
  vmc->SetCut("BCUTE",  cut); 
  vmc->SetCut("BCUTM",  cut); 
  vmc->SetCut("DCUTE",  cut); 
  vmc->SetCut("DCUTM",  cut); 
  vmc->SetCut("PPCUTM", cut);
  vmc->SetCut("TOFMAX", tofmax); 
  
}
