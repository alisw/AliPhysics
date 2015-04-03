doeffAllSpecies(Int_t isp=1){
  if(isp==1)  performAllPi();
  else if(isp==2) performAllKa();
  else if(isp == 3) performAllPr();
}

performAllPi(){
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libTree");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libAOD");
  gSystem->Load("libESD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libNetx");
  gSystem->Load("libPWGPPpid");
  
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/PWGPP/pid");
  
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/pid/doeffPi.C++");
  
  // tune these parameters
  cmin = 1; // centrality min 1
  cmax = 10;// centrality max 10
  Float_t etamin = -0.8;
  Float_t etamax = 0.8;

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
  kPid2Sigma = kFALSE;
  kPid3Sigma = kFALSE;
 
  if(! LoadLib()) return;

  kPid2Sigma = kTRUE;
  doeffPi(1,0.,etamin,etamax);
  doeffPi(0,0.,etamin,etamax);
  kPid2Sigma = kFALSE;
  kPid3Sigma = kTRUE;
  doeffPi(1,0.,etamin,etamax);
  doeffPi(0,0.,etamin,etamax);
  kPid3Sigma = kFALSE;

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
  gSystem->Unload("$ALICE_PHYSICS/PWGPP/pid/doeffPi_C");
}

performAllPr(){
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libTree");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libAOD");
  gSystem->Load("libESD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libNetx");
  gSystem->Load("libPWGPPpid");
  
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/PWGPP/pid");
  
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/pid/doeffPr.C++");
  
  // tune these parameters
  cmin = 1; // centrality min 1
  cmax = 10;// centrality max 10
  Float_t etamin = -0.8;
  Float_t etamax = 0.8;
  
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
  kPid2Sigma = kFALSE;
  kPid3Sigma = kFALSE;
 
  if(! LoadLib()) return;

  kPid2Sigma = kTRUE;
  doeffPr(1,0.,etamin,etamax);
  doeffPr(0,0.,etamin,etamax);
  kPid2Sigma = kFALSE;
  kPid3Sigma = kTRUE;
  doeffPr(1,0.,etamin,etamax);
  doeffPr(0,0.,etamin,etamax);
  kPid3Sigma = kFALSE;
  
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
  gSystem->Unload("$ALICE_PHYSICS/PWGPP/pid/doeffPr_C");
}

performAllKa(){
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libTree");
  gSystem->Load("libMinuit");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libAOD");
  gSystem->Load("libESD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libNetx");
  gSystem->Load("libPWGPPpid");

  gSystem->AddIncludePath("-I$ALICE_PHYSICS/PWGPP/pid");

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/pid/doeffKa.C++");

  // tune these parameters
  cmin = 4; // centrality min 1
  cmax = 10;// centrality max 10
  Float_t etamin = -0.8;
  Float_t etamax = 0.8;

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
  kPid2Sigma = kFALSE;
  kPid3Sigma = kFALSE;
 
  if(! LoadLib()) return;

  kPid2Sigma = kTRUE;
  doeffKa(1,0.,etamin,etamax);
  doeffKa(0,0.,etamin,etamax);
  kPid2Sigma = kFALSE;
  kPid3Sigma = kTRUE;
  doeffKa(1,0.,etamin,etamax);
  doeffKa(0,0.,etamin,etamax);
  kPid3Sigma = kFALSE;

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
  gSystem->Unload("$ALICE_PHYSICS/PWGPP/pid/doeffKa_C");
}

TGraphErrors *MakeRatio(const char *nf1,const char *nf2,const char *nfo=""/*output file*/,const char *title=""/*title*/){
  TFile *f1 = new TFile(nf1);
  TFile *f2 = new TFile(nf2);
  
  TGraphErrors *g1 = (TGraphErrors *) f1->Get("Graph");
  TGraphErrors *g2 = (TGraphErrors *) f2->Get("Graph");
  if(!(g1 && g2)) return NULL;
  if(!(g1 && g2)) return NULL;
  
  Int_t n1= g1->GetN();
  Int_t n2= g2->GetN();
  
  if(n1 != n2) return NULL;
  
  if(n1 > 100) n1 = 100;
  
  Float_t x[100],y[100],ex[100],ey[100];
  
  for(Int_t i=0;i < n1;i++){
    x[i] = g1->GetX()[i];
    ex[i] = g1->GetEX()[i];
    if(g1->GetY()[i] > 0 && g2->GetY()[i] > 0){
      y[i] = g1->GetY()[i] / g2->GetY()[i];
      ey[i] = (g1->GetEY()[i] / g1->GetY()[i])**2 + (g2->GetEY()[i] / g2->GetY()[i])**2;
      ey[i] = sqrt(ey[i])*y[i];
    }
    else{
      y[i] = 0.5;
      ey[i] = 0.5;
    }
  }
  TGraphErrors *gr = new TGraphErrors(n1,x,y,ex,ey);
  gr->SetTitle(title);

  if(nfo[0] != '\0'){ // write output
    printf("written in %s\n",nfo);
    TFile *fo = new TFile(nfo,"RECREATE");
    gr->Write();
    fo->Close();
  }

  return gr;
}
