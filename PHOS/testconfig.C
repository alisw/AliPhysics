void Config()
{
  // 7-DEC-2000 09:00
  // Switch on Transition Radiation simulation. 6/12/00 18:00
  // iZDC=1  7/12/00 09:00
  // ThetaRange is (0., 180.). It was (0.28,179.72) 7/12/00 09:00
  // Theta range given through pseudorapidity limits 22/6/2001
  
  // Set Random Number seed
  
  TDatime dat ; 
  gRandom->SetSeed(dat.GetTime());
  
  gSystem->Load("libgeant321.so") ; 
  new     TGeant3("C++ Interface to Geant3");
  
  TFile  *rootfile = new TFile("testPHOS.root", "recreate");
  rootfile->SetCompressionLevel(2);
  
  TGeant3 *geant3 = (TGeant3 *) gMC;
  
  //
  // Set External decayer
  AliDecayer *decayer = new AliDecayerPythia();
  
  decayer->SetForceDecay(kAll);
  decayer->Init();
  gMC->SetExternalDecayer(decayer);
  //
  //
  //=======================================================================
  // ******* GEANT STEERING parameters FOR ALICE SIMULATION *******
  geant3->SetTRIG(1);         //Number of events to be processed 
  geant3->SetSWIT(4, 10);
  geant3->SetDEBU(0, 0, 1);
  //geant3->SetSWIT(2,2);
  geant3->SetDCAY(1);
  geant3->SetPAIR(1);
  geant3->SetCOMP(1);
  geant3->SetPHOT(1);
  geant3->SetPFIS(0);
  geant3->SetDRAY(0);
  geant3->SetANNI(1);
  geant3->SetBREM(1);
  geant3->SetMUNU(1);
  geant3->SetCKOV(1);
  geant3->SetHADR(1);         //Select pure GEANH (HADR 1) or GEANH/NUCRIN (HADR 3)
  geant3->SetLOSS(2);
  geant3->SetMULS(1);
  geant3->SetRAYL(1);
  geant3->SetAUTO(1);         //Select automatic STMIN etc... calc. (AUTO 1) or manual (AUTO 0)
  geant3->SetABAN(0);         //Restore 3.16 behaviour for abandoned tracks
  geant3->SetOPTI(2);         //Select optimisation level for GEANT geometry searches (0,1,2)
  geant3->SetERAN(5.e-7);
  
  Float_t cut = 1.e-3;        // 1MeV cut by default
  Float_t tofmax = 1.e10;
  
  //             GAM ELEC NHAD CHAD MUON EBREM MUHAB EDEL MUDEL MUPA TOFMAX
  geant3->SetCUTS(cut, cut, cut, cut, cut, cut, cut, cut, cut, cut,
		  tofmax);
  //
  //=======================================================================
  // ************* STEERING parameters FOR ALICE SIMULATION **************
  // --- Specify event type to be tracked through the ALICE setup
  // --- All positions are in cm, angles in degrees, and P and E in GeV
  
  int     nParticles = 1;
  AliGenBox *gener = new AliGenBox(nParticles);
  
  gener->SetPart(22) ;
  gener->SetPtRange(9.99, 10.00);
  gener->SetPhiRange(220, 320);
  // Set pseudorapidity range from -8 to 8.
  Float_t thmin = EtaToTheta(0.12);   // 220 theta min. <---> eta max
  Float_t thmax = EtaToTheta(-0.12);  // 320 theta max. <---> eta min 
 
  gener->SetThetaRange(thmin, thmax);
  gener->SetOrigin(0, 0, 0);  //vertex position
  gener->SetSigma(0, 0, 0);   //Sigma in (X,Y,Z) (cm) on IP position
  gener->Init();
  // 
  // Activate this line if you want the vertex smearing to happen
  // track by track
  //
  //gener->SetVertexSmear(perTrack); 
  
  gAlice->SetField(0,2);  //Specify maximum magnetic field in Tesla (neg. ==> default field)
  
  //=================== Alice BODY parameters =============================
  AliBODY *BODY = new AliBODY("BODY", "Alice envelop");
  
  
  AliPHOS *PHOS = new AliPHOSv1("PHOS", "IHEP");
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance("aliroot") ;
  AliPHOSQAMeanChecker * hm  = static_cast<AliPHOSQAMeanChecker *>gime->QATasks("HitsMul");
  AliPHOSQAMeanChecker * te  = static_cast<AliPHOSQAMeanChecker *>gime->QATasks("TotEner");
  AliPHOSQAMeanChecker * hmB = static_cast<AliPHOSQAMeanChecker *>gime->QATasks("HitsMulB");
  AliPHOSQAMeanChecker * teB = static_cast<AliPHOSQAMeanChecker *>gime->QATasks("TotEnerB");
  hm->Set(62.18, 23.81) ;
  hm->Print() ; 
  te->Set(8.092, 3.06) ;
  hmB->Set(63.498, 24.348) ;
  teB->Set(8.363, 3.44) ;
  hm->Print("") ;      
}

Float_t EtaToTheta(Float_t arg){
  return (180./TMath::Pi())*2.*atan(exp(-arg));
}
