void AliTOFhits2sdigits(TString fileNameHits, Int_t firstEvent=0,Int_t nEvents=1) 
{

  /////////////////////////////////////////////////////////////////////////
  //
  // Creates TOF summable digits from the hit information. 
  //
  // Report problems to pierella@bo.infn.it
  //
  // Use case:
  // start root
  // // load the macro
  // root[0] .L AliTOFhits2sdigits.C
  // root[1] AliTOFhits2sdigits("galice.root",0,1)
  /////////////////////////////////////////////////////////////////////////

  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
    cout << "Loaded shared libraries" << endl;
  }


  // Create the TOF sdigitzer and sdigitize by default the first event
  // (in fact by default Int_t firstEvent=0,Int_t nEvents=1)
  AliTOFSDigitizer *sdigitizer = new AliTOFSDigitizer(fileNameHits.Data(),fileNameHits.Data(),firstEvent,nEvents); // it is the same nevents numbering
  // scheme used by STEER/AliHits2SDigits.C

  // Activate this line if you want to print the parameters
  // used in sdigitization
  // sdigitizer->PrintParameters();

  // performs sdigitization of the above events with "all" verbose option
  // "tim" option is also available for benchmarking only
  sdigitizer->Exec("all");  

  // N.B.: in order to maintain the functionality to sdigitize 
  // all events in current file add a second option
  // sdigitizer->Exec("all","all");
  // the second "all" option overrides the previous settings for 
  // lower and upper  bounds for event to sdigitize and allow
  // the sdigitization for ALL events in fileNameHits
}
