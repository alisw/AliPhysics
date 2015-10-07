// Use PROOF Lite to verify TDataSetManagerAliEn functionalities.
// AliEn credentials must be available on the system.
// Blame: Dario Berzano <dario.berzano@cern.ch>
void TestDataSets() {
  TString val;
  val.Form("alien cache:%s/datasetcache "
           "urltemplate:file://%s/remotestore<path> "
           "cacheexpiresecs:86400",
           gSystem->pwd(), gSystem->pwd());
  gEnv->SetValue("Proof.DataSetManager", val.Data());
  TProof::Open("workers=1");
  if (!gProof) gSystem->Exit(1);

  // Datasets are loaded from an external file
  TString ds = gSystem->GetFromPipe("cat input_datasets.txt");
  TFileCollection *fc = gProof->GetDataSet(ds.Data());
  if (fc == 0x0) {
    cout << "No file collection returned!" << endl;
    gSystem->Exit(1);
  }
  fc->Print("filter:SsCc");

  // Check default tree
  if (strcmp(fc->GetDefaultTreeName(), "/aodTree") != 0) {
    cout << "Default tree is not /aodTree!" << endl;
    gSystem->Exit(1);
  }

  // Check number of entries
  if (fc->GetNFiles() <= 0) {
    cout << "Dataset is empty!" << endl;
    gSystem->Exit(1);
  }

  TIter nxt(fc->GetList());
  TFileInfo *inf;
  TString m;
  m.Form("^%s/remotestore/alice/data/2013/LHC13e/000(196203|195949)"
         "/ESDs/muon_pass2/AOD134/[0-9]{4}/root_archive.zip$",
         gSystem->pwd());
  TPMERegexp r(m.Data());
  while ((inf = (TFileInfo *)nxt())) {
    if (strcmp(inf->GetCurrentUrl()->GetAnchor(), "AliAOD.root") != 0) {
      cout << "This entry's anchor is not AliAOD.root!" << endl;
      inf->Print();
      gSystem->Exit(1);
    }
    if (r.Match(inf->GetCurrentUrl()->GetFile()) == 0) {
      cout << "This entry's file does not match " << m.Data() << "!" << endl;
      inf->Print();
      gSystem->Exit(1);
    }
  }

  cout << "OK: dataset seem to match expectations" << endl;
  gSystem->Exit(0);
}
