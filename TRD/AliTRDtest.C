Int_t AliTRDtest() 
{
  //
  // Test macro for the TRD code
  //

  Int_t rc = 0;

  // Initialize the test setup 
  gAlice->Init("$(ALICE_ROOT)/TRD/AliTRDconfig.C");

  // Run one event and create the hits
  gAlice->Run(1);

  // Analyze the TRD hits
  gROOT->LoadMacro("$(ALICE_ROOT)/TRD/AliTRDanalyzeHits.C");
  if (rc = AliTRDanalyzeHits()) return rc;

  // Run the digitization
  gROOT->LoadMacro("$(ALICE_ROOT)/TRD/AliTRDcreateDigits.C");
  if (rc = AliTRDcreateDigits()) return rc;

  // Analyze the digits
  gROOT->LoadMacro("$(ALICE_ROOT)/TRD/AliTRDanalyzeDigits.C");
  if (rc = AliTRDanalyzeDigits()) return rc;

  return rc;

}
