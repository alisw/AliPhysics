void AliTRDhits2sdigits() 
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
  AliTRDdigitizer *digitizer = new AliTRDdigitizer("TRDdigitizer","TRD digitizer class");

  // Set the parameter
  digitizer->SetDebug(1);

  // For the summable digits
  digitizer->SetSDigits(kTRUE);

  // Open the AliRoot file
  digitizer->Open(alifile);

  // Define the parameter object
  // If no external parameter object is defined, 
  // default parameter will be used
  AliTRDparameter *parameter = new AliTRDparameter("TRDparameter"
						  ,"TRD parameter class");
  digitizer->SetParameter(parameter);

  // Create the digits
  digitizer->MakeDigits();

  // Write the digits into the input file
  digitizer->WriteDigits();

  // Save the parameter object in the AliROOT file
  AliRunLoader* rl = AliRunLoader::GetRunLoader(AliConfig::GetDefaultEventFolderName());
  rl->CdGAFile();
  parameter->Write();

}
