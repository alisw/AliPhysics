void runsim(){
  gSystem->Load("liblhapdf");
  gSystem->Load("libpythia6");
  gSystem->Load("libgeant321");

  gROOT->Macro("sim.C");
}
