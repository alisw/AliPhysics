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
  Char_t *alifile = "galice_d_v1.root"; 

  // Create the TRD digitzer 
  AliTRDdigitizer *Digitizer = new AliTRDdigitizer("digitizer","Digitizer class");

  // Set the parameter
  Digitizer->SetDiffusion();
  //Digitizer->SetExB();
  //Digitizer->SetElAttach();
  //Digitizer->SetAttachProb();

  // Open the AliRoot file
  Digitizer->Open(alifile);

  // Create the digits
  Digitizer->MakeDigits();

  // Write the digits into the input file
  Digitizer->WriteDigits();

  // Save the digitizer class in the AliROOT file
  Digitizer->Write();

}
