void sdigits2digits()
{

  /////////////////////////////////////////////////////////////////////////
  //
  // Converts s-digits to normal digits
  //
  /////////////////////////////////////////////////////////////////////////

  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
    cout << "Loaded shared libraries" << endl;
  }

  Char_t *fileName = "galice.root";

  // Create the TRD digits merger
  AliTRDdigitizer *digitizer = new AliTRDdigitizer("digitizer","Digitizer class");  

  // Set the parameter
  digitizer->SetDebug(1);

  // Initialize the geometry 
  digitizer->Open(fileName);

  // Create the digits manager for the input s-digits
  AliTRDdigitsManager *sdigitsManager = new AliTRDdigitsManager();
  sdigitsManager->SetDebug(1);
  sdigitsManager->SetSDigits(kTRUE);
  sdigitsManager->ReadDigits();

  // Add the s-digits to the input list 
  digitizer->AddSDigitsManager(sdigitsManager);

  // Convert the s-digits to normal digits
  digitizer->SDigits2Digits();

  // Store the digits
  digitizer->WriteDigits();

}
