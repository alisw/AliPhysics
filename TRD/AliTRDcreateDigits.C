Int_t AliTRDcreateDigits()
{
  //
  // Creates the digits from the hits of the slow simulator
  //

  Int_t rc = 0;

  if (!gAlice) {
    cout << "<AliTRDcreateDigits> No AliRun object found" << endl;
    rc = 1;
    return rc;
  }
  gAlice->GetEvent(0);

  // Create the TRD digitzer 
  AliTRDdigitizer *digitizer = new AliTRDdigitizer("digitizer","Digitizer class");
  digitizer->InitDetector();

  // Set the parameter (for TRF ~200ns)
  digitizer->SetGasGain(1600.);
  digitizer->SetChipGain(8.0);
  digitizer->SetNoise(1000.);
  digitizer->SetADCinRange(1000.);
  digitizer->SetADCoutRange(1023.);
  digitizer->SetADCthreshold(0);
  digitizer->SetVerbose(1);

  // Create the digits
  if (!(digitizer->MakeDigits())) {
    rc = 2;
    return rc;
  }

  // Write the digits into the input file
  if (!(digitizer->MakeBranch())) {
    rc = 3;
    return rc;
  }

  // Write the digits into the input file
  if (!(digitizer->WriteDigits())) {
    rc = 4;
    return rc;
  }

  // Save the digitizer class in the AliROOT file
  if (!(digitizer->Write())) {
    rc = 4;
    return rc;
  }

  return rc;

}
