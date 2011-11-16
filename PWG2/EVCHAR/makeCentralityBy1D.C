void makeCentralityBy1D 
//(int run =167693, const char *histname ="fHOutMultV0M", float percentXsec=90.05, float mult=77.4) //V0M
//(int run =167693, const char *histname ="fHOutMultCL1", float percentXsec=90.14, float mult=20) // CL1
(int run =167693, const char *histname ="fHOutMultTRK", float percentXsec=90.09, float mult=10) // TRK
{
 //load libraries
  // gSystem->SetBuildDir("/tmp/");
  // gSystem->Load("libCore.so");  
  // gSystem->Load("libTree.so");
  // gSystem->Load("libGeom.so");
  // gSystem->Load("libVMC.so");
  // gSystem->Load("libPhysics.so");
  // gSystem->Load("libSTEERBase.so");
  // gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->LoadMacro("AliCentralityBy1D.cxx+");
  AliCentralityBy1D *mPM = new AliCentralityBy1D();

  TString finname = Form("/home/atoia/analysis/data2011/multRef/AnalysisResults_%i.root",run);
  TString foutname = Form("/home/atoia/analysis/data2011/by1D/AliCentralityBy1D_%i_TRK.root",run);

  mPM->AddHisto(histname);  
  mPM->SetInputFile(finname);        
  mPM->SetOutputFile(foutname);  
  mPM->SetMultLowBound(mult);
  mPM->SetPercentileCrossSection(percentXsec);
  mPM->MakePercentiles();

}

