void run(Int_t runNumber)
{
  gEnv->SetValue("XSec.GSI.DelegProxy","2");
  // Select ROOT version
  TProof::Mgr("XXX@alicecaf")->SetROOTVersion("v5-26-00");
  // Login to CAF
  TProof::Open("XXX@alicecaf");

  // Enable AliRoot
  gProof->UploadPackage("/afs/cern.ch/alice/caf/sw/ALICE/PARs/v4-18-Release.rec/AF-v4-18-rec.par");
  gProof->EnablePackage("AF-v4-18-rec.par");

  // Run reconstruction
  gROOT->LoadMacro("rec.C");
  gROOT->ProcessLine(Form("rec(%d);",runNumber));

  TProof::Mgr("XXX@alicecaf")->GetSessionLogs()->Save("*",Form("run%d.log",runNumber));

  // Check the produced dataset
  TFileCollection *coll = gProof->GetDataSet(Form("run%d",runNumber));
  if (coll) {
    Int_t nEvents = coll->GetTotalEntries("/esdTree");
    if (nEvents > 0) {
      cout << "===========================================================================" << endl;
      cout << nEvents << " events reconstructed and stored in the dataset run" << runNumber << endl;
      cout << "===========================================================================" << endl;
      cout << "The dataset is:" << endl;
      coll->Print();
      cout << "===========================================================================" << endl;
    }
  }
}
