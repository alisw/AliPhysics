Int_t AliTRDcreateDigits()
{
  //
  // Creates the digits from the hits of the slow simulator
  //

  Int_t rc = 0;

  return rc;

  if (!gAlice) {
    cout << "<AliTRDcreateDigits> No AliRun object found" << endl;
    rc = 1;
    return rc;
  }
  gAlice->GetEvent(0);

  // Create the TRD digitzer 
  AliTRDdigitizer *Digitizer = new AliTRDdigitizer("digitizer","Digitizer class");

  // Initialize the TRD and the geometry
  if (!Digitizer->InitDetector()) {
    cout << "<AliTRDcreateDigits> No TRD geometry found" << endl;
    rc = 2;
    return rc;
  }

  // Set the parameter
  Digitizer->SetDiffusion();
  Digitizer->SetVerbose(1);

  // Create the digits
  if (!Digitizer->MakeDigits()) {
    rc = 3;
    return rc;
  }

  return rc;

}
