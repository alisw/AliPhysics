void makeCentralityFit(const char * run="188359",const char * system="ZNA", int Rebin=1,int Nevt=1e6)
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

  const char *finnameGlau ="GlauberMC_pPb_ntuple_sigma70_mind4_r662_a546_Rpro6.root";

  // data
  TString finname = Form("Histos%s.root",run);
  TString foutname = Form("%s_fit_%s.root",system,run);
  TString foutnameGlau = Form("%s_ntuple_%s.root",system,run);

  const char *histname=("hnSpecn");
  //const char *histname=("hnSpecp");
  

  AliCentralityGlauberFit *mPM = new AliCentralityGlauberFit(finnameGlau);
  mPM->SetInputFile(finname);        
  mPM->SetInputNtuple(finnameGlau);     
  mPM->SetOutputFile(foutname);  
  mPM->SetOutputNtuple(foutnameGlau);
  mPM->AddHisto(histname);

  mPM->SetRebin(Rebin);
  mPM->SetNevents(Nevt);
  mPM->UseChi2(kTRUE);     // If TRUE minimize Chi2

  if (strncmp (system,"ZNA",3) == 0) {
    mPM->SetRangeToFit(1., 70.);   
    mPM->SetRangeToScale(1.);  
    mPM->SetGlauberParam(1,0.00,1., 1,0.97,1., 1,0.95,0.3, 1,0.65,0.8, 1,0.585,0.6);
  }
  else if (strncmp (system,"ZPA",3) == 0) {
    mPM->SetRangeToFit(1., 25.);  2
    mPM->SetRangeToScale(1.);  
    mPM->SetGlauberParam(1,0.0,1., 1,0.6,1, 1,0.25,0.3, 1,0.65,0.8, 1,0.585,0.6);
  }
  mPM->MakeFits();  

  TFile * f = new TFile (foutname);
  TH1 * hd = (TH1*) gDirectory->Get(("hnSpecn"));
  TH1 * hg = (TH1*) gDirectory->Get(("hnSpecn_GLAU"));
  hg->SetLineColor(kRed);
  hd->Draw("hist");
  hg->Draw("histsame");



}

