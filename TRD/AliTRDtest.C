Int_t AliTRDtest() 
{
  //
  // Test macro for the TRD code
  //

  Int_t rc = 0;

  // Initialize the test setup 
  gAlice->Init("$(ALICE_ROOT)/TRD/AliTRDconfig.C");

  // Run one event and create the hits
  gAlice->SetDebug(2);
  gAlice->Run(1);

  if (gAlice) delete gAlice;
  TFile *file = (TFile *) gROOT->GetListOfFiles()->FindObject("TRD_test.root");
  gAlice = (AliRun *) file->Get("gAlice");

  // Analyze the TRD hits
  gROOT->LoadMacro("$(ALICE_ROOT)/TRD/AliTRDanalyzeHits.C");
  if (rc = AliTRDanalyzeHits()) return rc;

  // Run the digitization
  gROOT->LoadMacro("$(ALICE_ROOT)/TRD/AliTRDcreateDigits.C");
  if (rc = AliTRDcreateDigits()) return rc;

  if (gAlice) delete gAlice;
  TFile *file = (TFile *) gROOT->GetListOfFiles()->FindObject("TRD_test.root");
  gAlice = (AliRun *) file->Get("gAlice");

  // Analyze the digits
  gROOT->LoadMacro("$(ALICE_ROOT)/TRD/AliTRDanalyzeDigits.C");
  if (rc = AliTRDanalyzeDigits()) return rc;

  // Create the cluster
  gROOT->LoadMacro("$(ALICE_ROOT)/TRD/AliTRDcreateCluster.C");
  if (rc = AliTRDcreateCluster()) return rc;

  if (gAlice) delete gAlice;
  TFile *file = (TFile *) gROOT->GetListOfFiles()->FindObject("TRD_test.root");
  gAlice = (AliRun *) file->Get("gAlice");

  // Analyze the digits
  gROOT->LoadMacro("$(ALICE_ROOT)/TRD/AliTRDanalyzeCluster.C");
  if (rc = AliTRDanalyzeCluster()) return rc;

  return rc;

}
