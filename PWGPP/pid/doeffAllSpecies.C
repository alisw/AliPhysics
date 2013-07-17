doeffAllSpecies(){
  performAllPi();
  performAllKa();
  performAllPr();
}

performAllPi(){
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libNetx.so");
  gSystem->Load("libPWGPPpid.so");

  gSystem->AddIncludePath("-I$ALICE_ROOT/PWGPP/pid");

  gROOT->LoadMacro("$ALICE_ROOT/PWGPP/pid/doeffPi.C++");

  // tune these parameters
  cmin = 1; // centrality min 1
  cmax = 10;// centrality max 10
  Float_t etamin = -0.8;
  Float_t etamax = -0.8;

  // reset all flags
  rebinsize = 1; // don't change this, not choice here!!!
  kGoodMatch=kFALSE;
  kSigma2vs3 = kFALSE;
  require5sigma = kFALSE;
  bayesVsigma = kFALSE;
  kTOFmatch = kFALSE;
  kOverAll = kFALSE;
  kOverAllTOFmatch = kFALSE;
  kOverAll2Sigma = kFALSE;
  isMC = kFALSE;
  selectTrue = kTRUE;
  keepTrue = kFALSE;

  if(! LoadLib()) return;

  // matching and PID eff.
  doeffPi(1,0.1,etamin,etamax);
  doeffPi(0,0.1,etamin,etamax);
  doeffPi(1,0.2,etamin,etamax);
  doeffPi(1,0.4,etamin,etamax);
  doeffPi(1,0.6,etamin,etamax);
  doeffPi(1,0.8,etamin,etamax);
  doeffPi(0,0.2,etamin,etamax);
  doeffPi(0,0.4,etamin,etamax);
  doeffPi(0,0.6,etamin,etamax);
  doeffPi(0,0.8,etamin,etamax);

  // Good matching eff (1 - TOF mism)
  kGoodMatch=kTRUE;
  doeffPi(1,0.1,etamin,etamax);
  doeffPi(0,0.1,etamin,etamax);
  kGoodMatch=kFALSE;
  
  // eff 2 sigma / eff 3 sigma
  kSigma2vs3=kTRUE;
  doeffPi(1,0.1,etamin,etamax);
  doeffPi(0,0.1,etamin,etamax);
  kSigma2vs3=kFALSE;

  kOverAll=kTRUE;
  // TPC|TOF overall eff
  doeffPi(1,0.2,etamin,etamax);
  doeffPi(1,0.4,etamin,etamax);
  doeffPi(1,0.6,etamin,etamax);
  doeffPi(1,0.8,etamin,etamax);
  doeffPi(0,0.2,etamin,etamax);
  doeffPi(0,0.4,etamin,etamax);
  doeffPi(0,0.6,etamin,etamax);
  doeffPi(0,0.8,etamin,etamax);

  kOverAllTOFmatch=kTRUE;
  // TPC&TOF overall eff
  doeffPi(1,0.2,etamin,etamax);
  doeffPi(1,0.4,etamin,etamax);
  doeffPi(1,0.6,etamin,etamax);
  doeffPi(1,0.8,etamin,etamax);
  doeffPi(0,0.2,etamin,etamax);
  doeffPi(0,0.4,etamin,etamax);
  doeffPi(0,0.6,etamin,etamax);
  doeffPi(0,0.8,etamin,etamax);
  kOverAllTOFmatch=kFALSE;

  kOverAll2Sigma=kTRUE;
  // TPC&TOF 2 TOF sigma cut
  doeffPi(1,0.1,etamin,etamax);
  doeffPi(0,0.1,etamin,etamax);
  kOverAll2Sigma=kFALSE;
  kOverAll=kFALSE;
  gSystem->Unload("$ALICE_ROOT/PWGPP/pid/doeffPi_C.so");
}

performAllPr(){
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libNetx.so");
  gSystem->Load("libPWGPPpid.so");

  gSystem->AddIncludePath("-I$ALICE_ROOT/PWGPP/pid");

  gROOT->LoadMacro("$ALICE_ROOT/PWGPP/pid/doeffPr.C++");

  // tune these parameters
  cmin = 1; // centrality min 1
  cmax = 10;// centrality max 10
  Float_t etamin = -0.8;
  Float_t etamax = -0.8;

  // reset all flags
  rebinsize = 1; // don't change this, not choice here!!!
  kGoodMatch=kFALSE;
  kSigma2vs3 = kFALSE;
  require5sigma = kFALSE;
  bayesVsigma = kFALSE;
  kTOFmatch = kFALSE;
  kOverAll = kFALSE;
  kOverAllTOFmatch = kFALSE;
  kOverAll2Sigma = kFALSE;
  isMC = kFALSE;
  selectTrue = kTRUE;
  keepTrue = kFALSE;

  if(! LoadLib()) return;

  // matching and PID eff.
  doeffPr(1,0.1,etamin,etamax);
  doeffPr(0,0.1,etamin,etamax);
  doeffPr(1,0.2,etamin,etamax);
  doeffPr(1,0.4,etamin,etamax);
  doeffPr(1,0.6,etamin,etamax);
  doeffPr(1,0.8,etamin,etamax);
  doeffPr(0,0.2,etamin,etamax);
  doeffPr(0,0.4,etamin,etamax);
  doeffPr(0,0.6,etamin,etamax);
  doeffPr(0,0.8,etamin,etamax);

  // Good matching eff (1 - TOF mism)
  kGoodMatch=kTRUE;
  doeffPr(1,0.1,etamin,etamax);
  doeffPr(0,0.1,etamin,etamax);
  kGoodMatch=kFALSE;
  
  // eff 2 sigma / eff 3 sigma
  kSigma2vs3=kTRUE;
  doeffPr(1,0.1,etamin,etamax);
  doeffPr(0,0.1,etamin,etamax);
  kSigma2vs3=kFALSE;

  kOverAll=kTRUE;
  // TPC|TOF overall eff
  doeffPr(1,0.2,etamin,etamax);
  doeffPr(1,0.4,etamin,etamax);
  doeffPr(1,0.6,etamin,etamax);
  doeffPr(1,0.8,etamin,etamax);
  doeffPr(0,0.2,etamin,etamax);
  doeffPr(0,0.4,etamin,etamax);
  doeffPr(0,0.6,etamin,etamax);
  doeffPr(0,0.8,etamin,etamax);

  kOverAllTOFmatch=kTRUE;
  // TPC&TOF overall eff
  doeffPr(1,0.2,etamin,etamax);
  doeffPr(1,0.4,etamin,etamax);
  doeffPr(1,0.6,etamin,etamax);
  doeffPr(1,0.8,etamin,etamax);
  doeffPr(0,0.2,etamin,etamax);
  doeffPr(0,0.4,etamin,etamax);
  doeffPr(0,0.6,etamin,etamax);
  doeffPr(0,0.8,etamin,etamax);
  kOverAllTOFmatch=kFALSE;

  kOverAll2Sigma=kTRUE;
  // TPC&TOF 2 TOF sigma cut
  doeffPr(1,0.1,etamin,etamax);
  doeffPr(0,0.1,etamin,etamax);
  kOverAll2Sigma=kFALSE;
  kOverAll=kFALSE;
  gSystem->Unload("$ALICE_ROOT/PWGPP/pid/doeffPr_C.so");
}

performAllKa(){
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libTree.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libCORRFW.so");
  gSystem->Load("libNetx.so");
  gSystem->Load("libPWGPPpid.so");

  gSystem->AddIncludePath("-I$ALICE_ROOT/PWGPP/pid");

  gROOT->LoadMacro("$ALICE_ROOT/PWGPP/pid/doeffKa.C++");

  // tune these parameters
  cmin = 1; // centrality min 1
  cmax = 10;// centrality max 10
  Float_t etamin = -0.8;
  Float_t etamax = -0.8;

  // reset all flags
  rebinsize = 1; // don't change this, not choice here!!!
  kGoodMatch=kFALSE;
  kSigma2vs3 = kFALSE;
  require5sigma = kFALSE;
  bayesVsigma = kFALSE;
  kTOFmatch = kFALSE;
  kOverAll = kFALSE;
  kOverAllTOFmatch = kFALSE;
  kOverAll2Sigma = kFALSE;
  isMC = kFALSE;
  selectTrue = kTRUE;
  keepTrue = kFALSE;

  if(! LoadLib()) return;

  // matching and PID eff.
  doeffKa(1,0.1,etamin,etamax);
  doeffKa(0,0.1,etamin,etamax);
  doeffKa(1,0.2,etamin,etamax);
  doeffKa(1,0.4,etamin,etamax);
  doeffKa(1,0.6,etamin,etamax);
  doeffKa(1,0.8,etamin,etamax);
  doeffKa(0,0.2,etamin,etamax);
  doeffKa(0,0.4,etamin,etamax);
  doeffKa(0,0.6,etamin,etamax);
  doeffKa(0,0.8,etamin,etamax);

  // Good matching eff (1 - TOF mism)
  kGoodMatch=kTRUE;
  doeffKa(1,0.1,etamin,etamax);
  doeffKa(0,0.1,etamin,etamax);
  kGoodMatch=kFALSE;
  
  // eff 2 sigma / eff 3 sigma
  kSigma2vs3=kTRUE;
  doeffKa(1,0.1,etamin,etamax);
  doeffKa(0,0.1,etamin,etamax);
  kSigma2vs3=kFALSE;

  kOverAll=kTRUE;
  // TPC|TOF overall eff
  doeffKa(1,0.2,etamin,etamax);
  doeffKa(1,0.4,etamin,etamax);
  doeffKa(1,0.6,etamin,etamax);
  doeffKa(1,0.8,etamin,etamax);
  doeffKa(0,0.2,etamin,etamax);
  doeffKa(0,0.4,etamin,etamax);
  doeffKa(0,0.6,etamin,etamax);
  doeffKa(0,0.8,etamin,etamax);

  kOverAllTOFmatch=kTRUE;
  // TPC&TOF overall eff
  doeffKa(1,0.2,etamin,etamax);
  doeffKa(1,0.4,etamin,etamax);
  doeffKa(1,0.6,etamin,etamax);
  doeffKa(1,0.8,etamin,etamax);
  doeffKa(0,0.2,etamin,etamax);
  doeffKa(0,0.4,etamin,etamax);
  doeffKa(0,0.6,etamin,etamax);
  doeffKa(0,0.8,etamin,etamax);
  kOverAllTOFmatch=kFALSE;

  kOverAll2Sigma=kTRUE;
  // TPC&TOF 2 TOF sigma cut
  doeffKa(1,0.1,etamin,etamax);
  doeffKa(0,0.1,etamin,etamax);
  kOverAll2Sigma=kFALSE;
  kOverAll=kFALSE;
  gSystem->Unload("$ALICE_ROOT/PWGPP/pid/doeffKa_C.so");
}
