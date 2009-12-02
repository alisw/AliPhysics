// $Id$
//
// AliRoot Configuration for running aliroot with Monte Carlo.
// ConfigCommon2() includes the common setting for all MCs
// which has to be called after MC is instantiated.
// Called from MC specific configs (g3Config.C, g4Config.C).
//
// Extracted from G3 specific Config.C 
// by I. Hrivnacova, IPN Orsay

// Functions
Float_t EtaToTheta(Float_t arg);
AliGenerator* GeneratorFactory();

void genExtFileConfig()
{
  cout << "Running genExtFileConfig.C ... " << endl;

  //=======================================================================
  // Steering parameters for ALICE simulation
  //=======================================================================

  gMC->SetProcess("DCAY",1);
  gMC->SetProcess("PAIR",1);
  gMC->SetProcess("COMP",1);
  gMC->SetProcess("PHOT",1);
  gMC->SetProcess("PFIS",0);
  gMC->SetProcess("DRAY",0);
  gMC->SetProcess("ANNI",1);
  gMC->SetProcess("BREM",1);
  gMC->SetProcess("MUNU",1);
  gMC->SetProcess("CKOV",1);
  gMC->SetProcess("HADR",1);
  gMC->SetProcess("LOSS",2);
  gMC->SetProcess("MULS",1);
  //gMC->SetProcess("RAYL",1);

  Float_t cut = 1.e-3;        // 1MeV cut by default
  Float_t tofmax = 1.e10;

  gMC->SetCut("CUTGAM", cut);
  gMC->SetCut("CUTELE", cut);
  gMC->SetCut("CUTNEU", cut);
  gMC->SetCut("CUTHAD", cut);
  gMC->SetCut("CUTMUO", cut);
  gMC->SetCut("BCUTE",  cut); 
  gMC->SetCut("BCUTM",  cut); 
  gMC->SetCut("DCUTE",  cut); 
  gMC->SetCut("DCUTM",  cut); 
  gMC->SetCut("PPCUTM", cut);
  gMC->SetCut("TOFMAX", tofmax); 

  //=======================================================================
  // External decayer
  //=======================================================================

  TVirtualMCDecayer *decayer = new AliDecayerPythia();
  decayer->SetForceDecay(kAll);
  decayer->Init();
/*
  //forbid some decays
  AliPythia * py= AliPythia::Instance();
  py->SetMDME(737,1,0); //forbid D*+->D+ + pi0
  py->SetMDME(738,1,0);//forbid D*+->D+ + gamma

  for(Int_t d=747; d<=762; d++){ 
    py->SetMDME(d,1,0);
  }

  for(Int_t d=764; d<=807; d++){ 
    py->SetMDME(d,1,0);
  }
*/
  gMC->SetExternalDecayer(decayer);
  
  //=======================================================================
  // Event generator
  //=======================================================================

  // External generator configuration
  AliGenerator* gener = GeneratorFactory();
  gener->SetOrigin(0, 0, 0);    // vertex position
  //gener->SetSigma(0, 0, 5.3);   // Sigma in (X,Y,Z) (cm) on IP position
  //gener->SetCutVertexZ(1.);     // Truncate at 1 sigma
  //gener->SetVertexSmear(kPerEvent); 
  gener->SetTrackingFlag(1);
  gener->Init();
    
  cout << "Running genExtFileConfig.C finished ... " << endl;
}

Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}

AliGenerator* GeneratorFactory() {

  AliGenExtFile *gener = new AliGenExtFile(-1);
  AliGenReaderTreeK * reader = new AliGenReaderTreeK();

  reader->SetFileName("galice.root");
  reader->AddDir("gen");
  gener->SetReader(reader);
     
  return gener; 
}

