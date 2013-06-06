void loadlibs(const char *dir=".")
{
  //fix ld path for par files
  gSystem->Exec("export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH");
  
  //set inlcude paths
  gSystem->AddIncludePath("-I$ROOTSYS/include");
   Bool_t hasAR=!TString(gSystem->Getenv("ALICE_ROOT")).IsNull();
  // if (hasAR) gSystem->AddIncludePath("-I$ALICE_ROOT/ -I$ALICE_ROOT/include -I$ALICE_ROOT/STEER  -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/RAW");
  if (hasAR) gSystem->AddIncludePath("-I$ALICE_ROOT/ -I$ALICE_ROOT/include -I$ALICE_ROOT/STEER  -I$ALICE_ROOT/ANALYSIS -I$ALICE_ROOT/TPC -I$ALICE_ROOT/RAW -I$ALICE_ROOT/STEER/STEERBase/ -I$ALICE_ROOT/STEER/ESD/");

  gSystem->Load("libCore");
  gSystem->Load("libPhysics");
  gSystem->Load("libMinuit");
  gSystem->Load("libGui.so");

  gSystem->Load("libGeom");
  gSystem->Load("libVMC");

  gSystem->Load("libNet");
  gSystem->Load("libTree");
  gSystem->Load("libProof");

  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libCDB");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libRAWDatarec");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libSTEER");
  gSystem->Load("libSTAT");

  gSystem->Load("libTPCbase");
  gSystem->Load("libTPCrec");
  
  // gSystem->Load("libANALYSISalice");
 gSystem->Load("libTPCcalib");
 gSystem->Load("libThread");

  gSystem->AddIncludePath(Form("-I%s/",dir));                                                                                                                      
  //  gSystem->AddIncludePath(Form("-I%s/code/Event",dir));

  //gSystem->Exec(Form("cd %s/code; make",dir));

  // gSystem->Load(Form("%s/code/libGEMEvent.so",dir));
  //if (!gSystem->AccessPathName(Form("%s/code/libGEMtest.so",dir))){
  //   gSystem->Load(Form("%s/code/libGEMtest.so",dir));
  // }
  
//   gROOT->LoadMacro(Form("%s/code/RawReader/AliRawReaderGEM.cxx+g",dir));
//   gROOT->LoadMacro(Form("%s/code/AliTPCSimpleEventDisplay/AliTPCSimpleEventDisplay.cxx+g",dir));
  // if (hasAR) gROOT->LoadMacro(Form("%s/code/AliTPCSimpleEventDisplay/TestSimpleEvDisp.C+g",dir));
  // gROOT->LoadMacro(Form("%s/CRSIMSubTrack.cxx+g",dir));
  // gROOT->LoadMacro(Form("%s/CRSIMTrack.cxx+g",dir));
  // gROOT->LoadMacro(Form("%s/CRSIMEvent.cxx+g",dir));
  // gROOT->LoadMacro(Form("%s/CRSIMDrawer.cxx+g",dir));
  // gROOT->LoadMacro(Form("%s/CRSIMDisplay.cxx+g",dir));
  // gROOT->LoadMacro(Form("%s/CRSIMEventGenerator.cxx+g",dir));
  // gROOT->LoadMacro(Form("%s/CRSIMEventGeneratorSimple.cxx+g",dir));
  // gROOT->LoadMacro(Form("%s/CRSIMRunGenerator.cxx+g",dir));
  // CRSIMDrawer *drawer = new CRSIMDrawer("test.root");
  //CRSIMDisplay *disp = new CRSIMDisplay();
  gROOT->LoadMacro(Form("%s/ToyMCTrack.cxx+g",dir));
  gROOT->LoadMacro(Form("%s/ToyMCEvent.cxx+g",dir));
  gROOT->LoadMacro(Form("%s/ToyMCEventGenerator.cxx+g",dir));
  gROOT->LoadMacro(Form("%s/ToyMCEventGeneratorSimple.cxx+g",dir));
  //  gROOT->LoadMacro(Form("%s/ToyMCDrawer.cxx+g",dir));
  // disp->SetDrawer(drawer);
  
}
