void AliTOFtestProb(TString fileNameMatch) 
{
  /////////////////////////////////////////////////////////////////////////
  //
  // Test macro for TOF Prob
  // Author: F. Pierella
  // Report problems to pierella@bo.infn.it
  // input filenames:
  // fileNameMatch -> file with the result of matching
  //
  // Use case:
  // start root
  // // load the macro
  // root[0] .L AliTOFtestProb.C
  // root[1] AliTOFtestProb("match-6KevPYTHIA-0.2T.root")
  /////////////////////////////////////////////////////////////////////////

  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
    cout << "Loaded shared libraries" << endl;
  }
  
  AliTOFProb* tofprob=new AliTOFProb(fileNameMatch.Data());
  tofprob->Exec("");
}
