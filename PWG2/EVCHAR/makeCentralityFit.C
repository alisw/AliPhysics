void makeCentralityFit()
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

  const char * run="167693";
  const char * system="TRK";
  TString finname = Form("/home/atoia/analysis/data2011/multRef/AnalysisResults_%s.root",run);
  TString foutname = Form("/home/atoia/analysis/data2011/fit/%s_fit_%s.root",system,run);
  //TString foutnameGlau = Form("/home/atoia/analysis/data2011/fit/%s_ntuple_%s.root",system,run);
  TString foutnameGlau = Form("test_%s.root",system,run);
  //const char *histname="fHOutMultV0M";
  //const char *histname="fHOutMultCL1";
  const char *histname="fHOutMultTRK";


  AliCentralityGlauberFit *mPM = new AliCentralityGlauberFit(finnameGlau);
  mPM->SetInputFile(finname);        
  mPM->SetInputNtuple(finnameGlau);     
  mPM->SetOutputFile(foutname);  
  mPM->SetOutputNtuple(foutnameGlau);
  mPM->AddHisto(histname);

  mPM->SetRebin(1);
  mPM->SetAncestorMode(2);   // If 1 use Npart**alpha, if 2 use alpha*Npart + (1-alpha)*Ncoll
  mPM->SetFastFitMode(0);    // If 1 or 2 use cruder approximation to compute curve faster 1:NBD, 2:Gauss
  mPM->UseChi2(kTRUE);       // If TRUE minimize Chi2
  mPM->UseAverage(kFALSE);    // If TRUE use Average
  mPM->SetNtrials(1);
  mPM->SetNevents(1e5);

  // ----------range to fit--------------
  //mPM->SetRangeToFit(100., 20000.);   // V0M
  //mPM->SetRangeToFit(40., 5200.);   // CL1
  mPM->SetRangeToFit(10., 2400.);   // TRK
  // ----------range to scale--------------
  //mPM->SetRangeToScale(100.); // V0M  
  //mPM->SetRangeToScale(40.); // CL1  
  mPM->SetRangeToScale(10.);  //TRK 
  // ----------initial parameters--------------
  //mPM->SetGlauberParam(10,29,31, 10,0.7,1.5, 10,0.85,0.87);
  //mPM->MakeFitsMinuitNBD(0.801,30.,1.2);          // initial parameters
  //mPM->SetGlauberParam(1,28,28.5, 1,1.291,1.9, 1,0.801,0.82); // V0M
  //mPM->SetGlauberParam(1,7.25,7.5, 1,1.291,1.9, 1,0.801,0.82); // Cl1
  mPM->SetGlauberParam(1,3.5,3.7, 1,1.291,1.5, 1,0.801,0.81); // TRK
  // ----------done--------------

  mPM->MakeFits();  
  
  TFile * f = new TFile (foutname);
  TH1 * hd = (TH1*) gDirectory->Get("fHOutMultTRK");
  TH1 * hg = (TH1*) gDirectory->Get("fHOutMultTRK_GLAU");
  hg->SetLineColor(kRed);
  hd->Draw("e");
  hg->Draw("same");



}

