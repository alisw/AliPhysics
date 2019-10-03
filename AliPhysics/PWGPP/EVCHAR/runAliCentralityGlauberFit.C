{
  //load libraries
  //gSystem->SetBuildDir("/tmp");
  gSystem->Load("libCore");
  gSystem->Load("libTree");
  gSystem->Load("libGeom");
  gSystem->Load("libVMC");
  gSystem->Load("libPhysics");
  gSystem->Load("libSTEERBase");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->LoadMacro("AliCentralityGlauberFit.cxx+");
  
  const char *finname ="../centrality/cvetan/VZero_366.root"; // name input file
  const char *foutname="test.root"; // name output file

  float percentXsec=100.0;
  
  AliCentralityGlauberFit *mPM = new AliCentralityGlauberFit();
  mPM.SetOutputFile(foutname);

  // latest V0 cvetan corrected 
  mPM.AddHisto("fCorrMult1d");
  //mPM.SetGlauberParam(1,27.550,28., 10,0.4,0.8, 1,1.116,0.825, 1,0.974,0.98); //minuit Npart**alpha
  mPM.SetGlauberParam(1,23.5,25., 1,0.3,1.5, 1,1.115,0.83, 1,0.99,0.99); 
  mPM.SetRebin(10);
  mPM.SetRangeToFit(100., 21000.);   
  mPM.MakeFitsMinuitNBD(finname);  
  //mPM->MakeFits(finname);  

  TFile * f = new TFile (foutname);
  gDirectory->Get("fCorrMult1d")->Draw();;
  TH1 * hg = (TH1*) gDirectory->Get("fCorrMult1d_GLAU");
  hg->SetLineColor(kRed);
  hg->Draw("same");
}

