Int_t AliTOFtest(Int_t nevents=1) 
{
  //
  // Test macro for the TOF code
  // report bug to Fabrizio.Pierella@cern.ch
  // Use case:
  // start aliroot
  // root [0] .L AliTOFtest.C
  // root [1] AliTOFtest()
  //
  // Updated to the new I/O: A. De Caro, C. Zampolli
  //

  Int_t rc = 0;

  // Initialize the test setup 

  gAlice->Init("$ALICE_ROOT/TOF/AliTOFconfig.C");

  // Run one central Hijing event and create the hits
  // (time required: some minuts)

  gAlice->SetDebug(2);
  gAlice->Run(nevents);
  
  if (gAlice)
    {
      delete AliRunLoader::GetRunLoader();
      delete gAlice;
      gAlice = 0x0;
    }

  gROOT->LoadMacro("$(ALICE_ROOT)/TOF/AliTOFanalyzeHits.C");
  if (rc=AliTOFanalyzeHits()) return rc;
  
  if (gAlice)
    {
      delete AliRunLoader::GetRunLoader();
      delete gAlice;
      gAlice = 0x0;
    }

  gROOT->LoadMacro("$(ALICE_ROOT)/TOF/AliTOFhits2sdigits.C");
  if (rc=AliTOFhits2sdigits()) return rc;

    if (gAlice)
    {
      delete AliRunLoader::GetRunLoader();
      delete gAlice;
      gAlice = 0x0;
    }

  gROOT->LoadMacro("$(ALICE_ROOT)/TOF/AliTOFanalyzeSDigitsV2.C");
  if (rc=AliTOFanalyzeSDigitsV2()) return rc;
  
  if (gAlice)
    {
      delete AliRunLoader::GetRunLoader();
      delete gAlice;
      gAlice = 0x0;
    }  

  gROOT->LoadMacro("$(ALICE_ROOT)/TOF/AliTOFSDigits2Digits.C");
  if (rc=AliTOFSDigits2Digits()) return rc;
  //gROOT->LoadMacro("$(ALICE_ROOT)/TOF/AliTOFtestDigitizer.C");
  //if (rc=AliTOFtestDigitizer()) return rc;

    if (gAlice)
    {
      delete AliRunLoader::GetRunLoader();
      delete gAlice;
      gAlice = 0x0;
    }

  gROOT->LoadMacro("$(ALICE_ROOT)/TOF/AliTOFanalyzeDigits.C");
  if (rc=AliTOFanalyzeDigits()) return rc;
  if (gAlice)
    {
      delete AliRunLoader::GetRunLoader();
      delete gAlice;
      gAlice = 0x0;
    }

}
