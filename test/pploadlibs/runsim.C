void runsim(){
  gROOT->Macro("loadlibssim.C");

  gSystem->Load("libgeant321");

  gROOT->Macro("sim.C");
}
