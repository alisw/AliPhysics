void hits2digits() 
{

/////////////////////////////////////////////////////////////////////////
//
// Creates digits from the hit information. 
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
  AliTRDdigitizer *digitizer = new AliTRDdigitizer("digitizer","Digitizer class");

  // Set the parameter
  digitizer->SetVerbose(1);

  // Open the AliRoot file
  digitizer->Open(alifile);

  // Create the digits
  digitizer->MakeDigits();

  // Write the digits into the input file
  digitizer->WriteDigits();

  // Save the digitizer class in the AliROOT file
  digitizer->Write();

}
