void loadlibs () 
{
  if (gClassTable->GetID("AliRun") >= 0) return;
   cout<<"RICH private loadlibs.C ...";
  gSystem->Load("libminicern");
  if (gClassTable->GetID("TVector3") < 0)// additional check as libPhysics may be loaded in logon.C
      gSystem->Load("libPhysics");
  gSystem->Load("libEG");
  gSystem->Load("libSTEER");
  gSystem->Load("libTGeant3Dummy");
  gSystem->Load("libdummyhijing");
  gSystem->Load("libTHijing");
  gSystem->Load("libdummypythia6");
  gSystem->Load("libEGPythia6");
  gSystem->Load("libCONTAINERS");
  gSystem->Load("libdummyHBTP");
  gSystem->Load("libTHbtp");
  gSystem->Load("libdummymevsim");
  gSystem->Load("libTMevSim");
  gSystem->Load("libEVGEN");
  gSystem->Load("libHBTP"); 
  gSystem->Load("libRALICE");
  gSystem->Load("libFMD");
  gSystem->Load("libMUON");
  gSystem->Load("libPHOS");
  gSystem->Load("libPMD");
  gSystem->Load("libRICH");
  gSystem->Load("libSTRUCT");
  gSystem->Load("libTOF");
  gSystem->Load("libTPC");
  gSystem->Load("libTRD");
  gSystem->Load("libZDC");
  gSystem->Load("libITS");
  gSystem->Load("libCASTOR");
  gSystem->Load("libSTART");
  cout<<"done\n";
}
