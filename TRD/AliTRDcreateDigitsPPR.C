void AliTRDcreateDigitsPPR() {

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

  // Set the parameter
  Digitizer->SetDiffusion();
  Digitizer->SetVerbose(1);
  //Digitizer->SetExB();

  // Take out ADC conversion for raw data digits
  Digitzer->SetADCinRange(1000000);
  Digitzer->SetADCoutRange(1000000);
  Digitzer->SetADCthreshold(0);

  // Open the AliRoot file
  Digitizer->Open(alifile);

  // Create the digits
  Digitizer->MakeDigits();

  // Write the digits into the input file
  Digitizer->WriteDigits();

  // Save the digitizer class in the AliROOT file
  Digitizer->Write();

}
