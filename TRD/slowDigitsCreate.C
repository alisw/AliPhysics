void slowDigitsCreate() {

/////////////////////////////////////////////////////////////////////////
//
// Creates the digits from the hit information. 
//
/////////////////////////////////////////////////////////////////////////

  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
    cout << "Loaded shared libraries" << endl;
  }

  // Input (and output) file name
  Char_t *alifile = "galice.root"; 

  // Create the TRD digitzer 
  AliTRDdigitizer *Digitizer = new AliTRDdigitizer("digitizer","Digitizer class");

  // Open the AliRoot file
  Digitizer->Open(alifile);

  // Set the parameter
  Digitizer->SetDiffusion(0);
  Digitizer->SetVerbose(1);
  //Digitizer->SetTimeResponse(0);
  //Digitizer->SetExB();
  //Digitizer->SetElAttach();
  //Digitizer->SetAttachProb();

  // Create the digits
  Digitizer->MakeDigits();

  // Write the digits into the input file
  Digitizer->WriteDigits();

  // Save the digitizer class in the AliROOT file
  Digitizer->Write();

}
