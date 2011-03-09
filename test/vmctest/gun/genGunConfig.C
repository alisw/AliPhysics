/// $Id$
//
// AliRoot Configuration for gun vmctest
//
// Extracted from G3 specific Config.C in test/gun.
// by I. Hrivnacova, IPN Orsay

Float_t EtaToTheta(Float_t arg);

AliGenerator* genGunConfig()
{
  cout << "Running genGunConfig.C ... " << endl;

  //=======================================================================
  // Event generator
  //=======================================================================

  // The cocktail itself

  AliGenCocktail *gener = new AliGenCocktail();
  gener->SetPhiRange(0, 360);
  // Set pseudorapidity range from -8 to 8.
  Float_t thmin = EtaToTheta(8);   // theta min. <---> eta max
  Float_t thmax = EtaToTheta(-8);  // theta max. <---> eta min 
  gener->SetThetaRange(thmin,thmax);
  gener->SetOrigin(0, 0, 0);  //vertex position
  gener->SetSigma(0, 0, 0);   //Sigma in (X,Y,Z) (cm) on IP position


  // Particle guns for the barrel part (taken from RichConfig)

  AliGenFixed *pG1=new AliGenFixed(1);
  pG1->SetPart(kProton);
  pG1->SetMomentum(2.5);
  pG1->SetTheta(109.5-3);
  pG1->SetPhi(10);
  gener->AddGenerator(pG1,"g1",1);
  
  AliGenFixed *pG2=new AliGenFixed(1);
  pG2->SetPart(kPiPlus);
  pG2->SetMomentum(1.0);
  pG2->SetTheta( 90.0-3);
  pG2->SetPhi(10);
  gener->AddGenerator(pG2,"g2",1);

  AliGenFixed *pG3=new AliGenFixed(1);
  pG3->SetPart(kPiMinus);
  pG3->SetMomentum(1.5);
  pG3->SetTheta(109.5-3);
  pG3->SetPhi(30);
  gener->AddGenerator(pG3,"g3",1);
  
  AliGenFixed *pG4=new AliGenFixed(1);
  pG4->SetPart(kKPlus);
  pG4->SetMomentum(0.7);
  pG4->SetTheta( 90.0-3);
  pG4->SetPhi(30);
  gener->AddGenerator(pG4,"g4",1);
  
  AliGenFixed *pG5=new AliGenFixed(1);
  pG5->SetPart(kKMinus);
  pG5->SetMomentum(1.0);
  pG5->SetTheta( 70.0-3);
  pG5->SetPhi(30);
  gener->AddGenerator(pG5,"g5",1);
  
  AliGenFixed *pG6=new AliGenFixed(1);
  pG6->SetPart(kProtonBar);
  pG6->SetMomentum(2.5);
  pG6->SetTheta( 90.0-3);
  pG6->SetPhi(50);
  gener->AddGenerator(pG6,"g6",1);
  
  AliGenFixed *pG7=new AliGenFixed(1);
  pG7->SetPart(kPiMinus);
  pG7->SetMomentum(0.7);
  pG7->SetTheta( 70.0-3);
  pG7->SetPhi(50);
  gener->AddGenerator(pG7,"g7",1);

  // Electrons for TRD

  AliGenFixed *pG8=new AliGenFixed(1);
  pG8->SetPart(kElectron);
  pG8->SetMomentum(1.2);
  pG8->SetTheta( 95.0);
  pG8->SetPhi(190);
  gener->AddGenerator(pG8,"g8",1);

  AliGenFixed *pG9=new AliGenFixed(1);
  pG9->SetPart(kPositron);
  pG9->SetMomentum(1.2);
  pG9->SetTheta( 85.0);
  pG9->SetPhi(190);
  gener->AddGenerator(pG9,"g9",1);

  // PHOS

  AliGenBox *gphos = new AliGenBox(1);
  gphos->SetMomentumRange(10,11.);
  gphos->SetPhiRange(270.5,270.7);
  gphos->SetThetaRange(90.5,90.7);
  gphos->SetPart(kGamma);
  gener->AddGenerator(gphos,"GENBOX GAMMA for PHOS",1);

  // EMCAL

  AliGenBox *gemcal = new AliGenBox(1);
  gemcal->SetMomentumRange(10,11.);
  gemcal->SetPhiRange(90.5,199.5);
  gemcal->SetThetaRange(90.5,90.7);
  gemcal->SetPart(kGamma);
  gener->AddGenerator(gemcal,"GENBOX GAMMA for EMCAL",1);

  // MUON
  AliGenBox * gmuon1 = new AliGenBox(1);
  gmuon1->SetMomentumRange(20.,20.1);
  gmuon1->SetPhiRange(0., 360.);         
  gmuon1->SetThetaRange(171.000,178.001);
  gmuon1->SetPart(kMuonMinus);           // Muons
  gener->AddGenerator(gmuon1,"GENBOX MUON1",1);

  AliGenBox * gmuon2 = new AliGenBox(1);
  gmuon2->SetMomentumRange(20.,20.1);
  gmuon2->SetPhiRange(0., 360.);         
  gmuon2->SetThetaRange(171.000,178.001);
  gmuon2->SetPart(kMuonPlus);           // Muons
  gener->AddGenerator(gmuon2,"GENBOX MUON1",1);

  //TOF
  AliGenFixed *gtof=new AliGenFixed(1);
  gtof->SetPart(kProton);
  gtof->SetMomentum(2.5);
  gtof->SetTheta(95);
  gtof->SetPhi(340);
  gener->AddGenerator(gtof,"Proton for TOF",1);

  //FMD1
  AliGenFixed *gfmd1=new AliGenFixed(1);
  gfmd1->SetPart(kGamma);
  gfmd1->SetMomentum(25);
  gfmd1->SetTheta(1.8);
  gfmd1->SetPhi(10);
  gener->AddGenerator(gfmd1,"Gamma for FMD1",1);
  
  //FMD2i
  AliGenFixed *gfmd2i=new AliGenFixed(1);
  gfmd2i->SetPart(kPiPlus);
  gfmd2i->SetMomentum(1.5);
  gfmd2i->SetTheta(7.3);
  gfmd2i->SetPhi(20);
  gener->AddGenerator(gfmd2i,"Pi+ for FMD2i",1);
  
  //FMD2o
  AliGenFixed *gfmd2o=new AliGenFixed(1);
  gfmd2o->SetPart(kPiMinus);
  gfmd2o->SetMomentum(1.5);
  gfmd2o->SetTheta(16.1);
  gfmd2o->SetPhi(30);
  gener->AddGenerator(gfmd2o,"Pi- for FMD2o",1);
  
  //FMD3o
  AliGenFixed *gfmd3o=new AliGenFixed(1);
  gfmd3o->SetPart(kPiPlus);
  gfmd3o->SetMomentum(1.5);
  gfmd3o->SetTheta(163.9);
  gfmd3o->SetPhi(40);
  gener->AddGenerator(gfmd3o,"Pi+ for FMD3o",1);
  
  //FMD3i
  AliGenFixed *gfmd3i=new AliGenFixed(1);
  gfmd3i->SetPart(kPiMinus);
  gfmd3i->SetMomentum(1.5);
  gfmd3i->SetTheta(170.5);
  gfmd3i->SetPhi(50);
  gener->AddGenerator(gfmd3i,"Pi- for FMD3i",1);
  
  //VZERO C
  AliGenFixed *gv0c=new AliGenFixed(1);
  gv0c->SetPart(kPiPlus);
  gv0c->SetMomentum(1.5);
  gv0c->SetTheta(170);
  gv0c->SetPhi(50);
  gener->AddGenerator(gv0c,"Pi+ for V0C",1);
  
  //VZERO A
  AliGenFixed *gv0a=new AliGenFixed(1);
  gv0a->SetPart(kPiMinus);
  gv0a->SetMomentum(1.5);
  gv0a->SetTheta(1.5);
  gv0a->SetPhi(70);
  gener->AddGenerator(gv0a,"Pi- for V0A",1);


  //PMD
  AliGenFixed *gpmd=new AliGenFixed(1);
  gpmd->SetPart(kGamma);
  gpmd->SetMomentum(2);
  gpmd->SetTheta(12.6);
  gpmd->SetPhi(60);
  gener->AddGenerator(gpmd,"Gamma for PMD",1);

  //ZDC
  AliGenFixed *gzdc1=new AliGenFixed(1);
  gzdc1->SetPart(kProton);
  gzdc1->SetMomentum(700);
  gzdc1->SetTheta(0.6);
  gzdc1->SetPhi(60);
  gener->AddGenerator(gzdc1,"Proton for ZDC",1);

  AliGenFixed *gzdc2=new AliGenFixed(1);
  gzdc2->SetPart(kNeutron);
  gzdc2->SetMomentum(500);
  gzdc2->SetTheta(0.6);
  gzdc2->SetPhi(60);
  gener->AddGenerator(gzdc2,"Neutron for ZDC",1);

  //T0
  AliGenFixed *gt0=new AliGenFixed(1);
  gt0->SetPart(kPiPlus);
  gt0->SetMomentum(2);
  gt0->SetTheta(5.1);
  gt0->SetPhi(60);
  gener->AddGenerator(gt0,"Pi+ for T0",1);

  AliGenFixed *gt01=new AliGenFixed(1);
  gt01->SetPart(kPiMinus);
  gt01->SetMomentum(2);
  gt01->SetTheta(5.1);
  gt01->SetPhi(60);
  gener->AddGenerator(gt01,"Pi- for T0",1);


  //ACORDE
  AliGenFixed *gacorde=new AliGenFixed(1);
  gacorde->SetPart(kMuonPlus);
  gacorde->SetMomentum(20);
  gacorde->SetTheta(90.);
  gacorde->SetPhi(90);
  gener->AddGenerator(gacorde,"Muon+ for ACORDE",1);

  AliGenFixed *gacorde1=new AliGenFixed(1);
  gacorde1->SetPart(kMuonMinus);
  gacorde1->SetMomentum(20);
  gacorde1->SetTheta(90.);
  gacorde1->SetPhi(90);
  gener->AddGenerator(gacorde1,"Muon- for ACORDE",1);

  gener->Init();
  
  return gener;
  
  cout << "Running genGunConfig.C finished ... " << endl;
}

Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}
