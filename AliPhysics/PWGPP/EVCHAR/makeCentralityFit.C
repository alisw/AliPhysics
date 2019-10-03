void makeCentralityFit(const char * run="188362",const char * system="TKL", int Rebin=1,int Nevt=1e6)
{
 //load libraries
  gSystem->SetBuildDir("/tmp/");
  gSystem->Load("libCore");
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libSTEERBase");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->LoadMacro("AliCentralityGlauberFit.cxx+");

  const char *finnameGlau ="/home/atoia/GlauberNtuplePA/GlauberMC_pPb_ntuple_sigma70_mind4_r662_a546_Rpro4.root";
  TString finname = Form("/home/atoia/analysis/data2012/multRef/EventStat_temp_%s.root",run);
  TString foutname = Form("/home/atoia/analysis/data2012/fit/%s_fit_%s.root",system,run);
  TString foutnameGlau = Form("/home/atoia/analysis/data2012/fit/%s_ntuple_%s.root",system,run);
  const char *histname=Form("fHOutMult%s",system);

  AliCentralityGlauberFit *mPM = new AliCentralityGlauberFit(finnameGlau);
  mPM->SetInputFile(finname);        
  mPM->SetInputNtuple(finnameGlau);     
  mPM->SetOutputFile(foutname);  
  mPM->SetOutputNtuple(foutnameGlau);
  mPM->AddHisto(histname);

  mPM->SetRebin(Rebin);
  mPM->SetNevents(Nevt);
  mPM->SetAncestorMode(1); // 1: Npart**alpha, 2: alpha*Npart + (1-alpha)*Ncoll
  mPM->SetFastFitMode(0);  // 1:NBD, 2:Gauss
  mPM->UseChi2(kTRUE);     // If TRUE minimize Chi2
  mPM->UseAverage(kFALSE); // If TRUE use Average
  mPM->SetNtrials(1);

  // ----------range to fit--------------
  if (strncmp (system,"V0A",1) == 0) {
    mPM->SetRangeToFit(15., 600.);   // range to fit
    mPM->SetRangeToScale(15.); // range to scale
    mPM->SetGlauberParam(1,12.8,13, 20,0.5,2.5, 1,1,1); // fit parameters
  }
  else if (strncmp (system,"TKL",1) == 0) {
    mPM->SetRangeToFit(15., 200.);   // range to fit
    mPM->SetRangeToScale(15.); // range to scale
    mPM->SetGlauberParam(1,4.9,5.5, 1,0.61,0.62, 1,1,1); // fit parameters
  }
  else if (strncmp (system,"CL1",1) == 0) {
    mPM->SetRangeToFit(40., 400.);   
    mPM->SetRangeToScale(40.); 
    mPM->SetGlauberParam(1,7.9,8, 1,0.43,0.46, 1,1,1); 
  }
  else if (strncmp (system,"TRK",1) == 0) {
    mPM->SetRangeToFit(10., 2600.);   
    mPM->SetRangeToScale(10.);  
    mPM->SetGlauberParam(1,3.9,4.2, 1,1.3,2.5, 1,0.801,0.81);
  }

  mPM->MakeFits();  

  // ----------for Minuit--------------
  //mPM->MakeFitsMinuitNBD(1,8,0.5);          // initial parameters


  TFile * f = new TFile (foutname);
  TH1 * hd = (TH1*) gDirectory->Get(Form("fHOutMult%s",system));
  TH1 * hg = (TH1*) gDirectory->Get(Form("fHOutMult%s_GLAU",system));
  hg->SetLineColor(kRed);
  hd->Draw("e");
  hg->Draw("same");



}

