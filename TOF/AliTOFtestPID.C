void AliTOFtestPID(TString fileNameMatch, TString fileNameCuts) 
{
  /////////////////////////////////////////////////////////////////////////
  //
  // Test macro for TOF PID 
  // Author: F. Pierella
  // Report problems to pierella@bo.infn.it
  // input filenames:
  // fileNameMatch -> file with the result of matching
  // fileNameCuts  -> file containing the graphical cuts
  //
  // Use case:
  // start root
  // // load the macro
  // root[0] .L AliTOFtestPID.C
  // root[1] AliTOFtestPID("match-6KevPYTHIA-0.2T.root","stdCutspp.root")
  /////////////////////////////////////////////////////////////////////////

  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
    cout << "Loaded shared libraries" << endl;
  }
  
  AliTOFPID* tofpid=new AliTOFPID(fileNameMatch.Data(),fileNameCuts.Data());
  // make a choice: uncomment one of these lines and try
  // tofpid->Exec("pp","visual","asC");
  // tofpid->Exec("pp","novisual","asC");
  // tofpid->Exec("Pb-Pb","visual","asC");
  // e.g. for p-p events
  tofpid->Exec("pp","visual","asC");
  //tofpid->Exec("pp","novisual","asEPS");
}
