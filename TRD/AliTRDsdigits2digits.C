void AliTRDsdigits2digits()
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
  AliTRDdigitizer *digitizer = new AliTRDdigitizer("TRDdigitizer"
                                                  ,"TRD digitizer class");  

  // Set the parameter
  digitizer->SetDebug(1);

  // Initialize the geometry 
  digitizer->Open(fileName);

  AliRunLoader* rl = AliRunLoader::GetRunLoader(AliConfig::GetDefaultEventFolderName());
  AliLoader* loader = rl->GetLoader("TRDLoader");
  // Define the parameter object
  // If no external parameter object is defined, 
  // default parameter will be used
  AliTRDparameter *parameter = new AliTRDparameter("TRDparameter","TRD parameter class");
  digitizer->SetParameter(parameter);

  // Create the digits manager for the input s-digits
  AliTRDdigitsManager *sdigitsManager = new AliTRDdigitsManager();
  sdigitsManager->SetDebug(1);
  sdigitsManager->SetSDigits(kTRUE);
  if (loader->TreeS() == 0x0) loader->LoadSDigits();
  
  sdigitsManager->ReadDigits(loader->TreeS());
  // Add the s-digits to the input list 
  digitizer->AddSDigitsManager(sdigitsManager);

  // Convert the s-digits to normal digits
  digitizer->SDigits2Digits();

  // Store the digits
  digitizer->WriteDigits();

  // Save the parameter object in the AliROOT file

  rl->CdGAFile();
  parameter->Write();

}
