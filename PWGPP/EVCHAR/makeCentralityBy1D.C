void makeCentralityBy1D 
//(int run =167987, const char *system ="V0M", float percentXsec=90.06, float mult=77.4) //V0M
//(int run =170162, const char *system ="CL1", float percentXsec=90.04, float mult=20) // CL1
//(int run =167693, const char *system ="TRK", float percentXsec=90.20, float mult=11) // TRK

// for pA pilot run
// for MC use percentXsec=100 and mult=0
// for MC available also NPA
// for pA use V0M, V0A, V0C, CL1, TRK, CND, (tbd FMD: not available in 146079) 
(int run =195344, const char *system ="FMD", float percentXsec=100.0, float mult=0.0) 
{
 //load libraries
  // gSystem->SetBuildDir("/tmp/");
  // gSystem->Load("libCore");
  // gSystem->Load("libTree");
  // gSystem->Load("libGeom");
  // gSystem->Load("libVMC");
  // gSystem->Load("libPhysics");
  // gSystem->Load("libSTEERBase");
  // gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->LoadMacro("AliCentralityBy1D.cxx+");
  AliCentralityBy1D *mPM = new AliCentralityBy1D();

  TString finname = Form("/home/atoia/analysis/data2013/multRef/EventStat_temp_%i.root",run);
  TString foutname = Form("/home/atoia/analysis/data2013/by1D/AliCentralityBy1D_%i_%s.root",run,system);
  const char *histname=Form("fHOutMult%s",system);

  mPM->AddHisto(histname);  
  mPM->SetInputFile(finname);        
  mPM->SetOutputFile(foutname);  
  mPM->SetMultLowBound(mult);
  mPM->SetPercentileCrossSection(percentXsec);
  mPM->MakePercentiles();

}

