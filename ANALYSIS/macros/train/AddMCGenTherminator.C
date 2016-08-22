AliGenerator *AddMCGenTherminator()
{  
// User defined generator  
  gSystem->Load("liblhapdf");      // Parton density functions
  gSystem->Load("libEGPythia6");   // TGenerator interface
  gSystem->Load("libpythia6");     // Pythia
  gSystem->Load("libAliPythia6");  // ALICE specific implementations
  gSystem->Load("libhijing");       
  gSystem->Load("libTHijing");
  gSystem->Load("libTTherminator");


  AliGenTherminator* gener = new AliGenTherminator(-1);

  gener->SetModel("SingleFreezeOut"); // setting for Cracow model
  gener->SetEventNumberInFile(500);
  gener->SetTau(9.0);
  gener->SetRhoMax(11.4);
  gener->SetMiuB(0.0);
  gener->SetPtHardMin(2.3);
  
  return gener;
}
  
  
