void makeCentralityFit(const char * run="167693",const char * system="V0M", int Rebin=100,int Nevt=1e5, int Ncell=0)
{
 //load libraries
  gSystem->SetBuildDir("/tmp/");
  gSystem->Load("libCore.so");  
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase.so");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->LoadMacro("AliCentralityGlauberFit.cxx+");

  const char *finnameGlau ="/home/atoia/GlauberNtuple/GlauberMC_PbPb_ntuple_sigma64_mind4_r662_a546.root";


  TString finname = Form("/home/atoia/analysis/data2011/multRef/EventStat_temp_%s.root",run);
  TString foutname = Form("/home/atoia/analysis/data2011/fit/%s_fitTEST_%s.root",system,run);
  TString foutnameGlau = Form("/home/atoia/analysis/data2011/fit/%s_ntupleTEST_%s.root",system,run);
  const char *histname=Form("fHOutMult%s",system);


  /*
  TString finname = Form("/home/atoia/analysis/EPVzero/VZEROEquaFactorStat.root");
  TString foutname = Form("/home/atoia/analysis/data2011/fit/%sCell%d_fit_%s.root",system,Ncell,run);
  TString foutnameGlau = Form("/home/atoia/analysis/data2011/fit/%sCell%d_ntuple_%s.root",system,Ncell,run);
  const char *histname=Form("fMultCell_%d",Ncell);
  */
  AliCentralityGlauberFit *mPM = new AliCentralityGlauberFit(finnameGlau);
  mPM->SetInputFile(finname);        
  mPM->SetInputNtuple(finnameGlau);     
  mPM->SetOutputFile(foutname);  
  mPM->SetOutputNtuple(foutnameGlau);
  mPM->AddHisto(histname);

  mPM->SetRebin(Rebin);
  mPM->SetNevents(Nevt);
  mPM->SetAncestorMode(2); // 1: Npart**alpha, 2: alpha*Npart + (1-alpha)*Ncoll
  mPM->SetFastFitMode(0);  // 1:NBD, 2:Gauss
  mPM->UseChi2(kTRUE);     // If TRUE minimize Chi2
  mPM->UseAverage(kFALSE); // If TRUE use Average
  mPM->SetNtrials(1);

  // ----------range to fit--------------
  if (strncmp (system,"V0M",1) == 0) {
    mPM->SetRangeToFit(100., 22000.);   // range to fit
    mPM->SetRangeToScale(100.); // range to scale
    mPM->SetGlauberParam(1,28.7,29., 1,1.601,1.5, 1,0.80,0.805); // fit parameters
  }
  else if (strncmp (system,"CL1",1) == 0) {
    mPM->SetRangeToFit(40., 5400.);   
    mPM->SetRangeToScale(40.); 
    mPM->SetGlauberParam(1,7.13,7.3, 1,1.217,1.9, 1,0.802,0.815); 
  }
  else if (strncmp (system,"TRK",1) == 0) {
    mPM->SetRangeToFit(10., 2600.);   
    mPM->SetRangeToScale(10.);  
    mPM->SetGlauberParam(1,3.9,4.2, 1,1.3,2.5, 1,0.801,0.81);
  }

  mPM->MakeFits();  

  // ----------for Minuit--------------
  //mPM->MakeFitsMinuitNBD(0.8,28.,1.29);          // initial parameters


  TFile * f = new TFile (foutname);
  TH1 * hd = (TH1*) gDirectory->Get(Form("fHOutMult%s",system));
  TH1 * hg = (TH1*) gDirectory->Get(Form("fHOutMult%s_GLAU",system));
  //TH1 * hd = (TH1*) gDirectory->Get(Form("fMultCell_%d",Ncell));
  //TH1 * hg = (TH1*) gDirectory->Get(Form("fMultCell_%d_GLAU",Ncell));
  hg->SetLineColor(kRed);
  hd->Draw("e");
  hg->Draw("same");



}

