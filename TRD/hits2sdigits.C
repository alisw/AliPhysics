void hits2sdigits() 
{

  /////////////////////////////////////////////////////////////////////////
  //
  // Creates summable digits from the hit information. 
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
  AliTRDdigitizer *digitizer = new AliTRDdigitizer("TRDdigitizer","Digitizer class");

  // Set the parameter
  digitizer->SetVerbose(1);

  // For the summable digits
  digitizer->SetSDigits(kTRUE);

  // Open the AliRoot file
  digitizer->Open(alifile);

  // Create the digits
  digitizer->MakeDigits();

  // Write the digits into the input file
  digitizer->WriteDigits();

  // Save the digitizer class in the AliROOT file
  digitizer->Write();

}
