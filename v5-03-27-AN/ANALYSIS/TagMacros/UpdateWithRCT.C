void UpdateWithRCT(const char *inrun, const char *inrct)
{
  TFile *inrunfile = new TFile(inrun);
  TTree *runtree = (TTree *) inrunfile->Get("T");
  AliRunTag *runt = new AliRunTag();
  runtree->SetBranchAddress("AliTAG", &runt);

  TFile *inrctfile = new TFile(inrct);
  TTree *rcttree = (TTree *) inrctfile->Get("T");
  AliRunTag *rctt = new AliRunTag();
  rcttree->SetBranchAddress("AliTAG", &rctt);

  runtree->GetEntry(0);

  cout << "Looking for RCT match for run " << runt->GetRunId() << endl;

  for (int iter=0; iter<rcttree->GetEntries(); iter++) {
    rcttree->GetEntry(iter);
    if (rctt->GetRunId() == runt->GetRunId()) {
      cout << "Found match in RCT for run " << rctt->GetRunId() << endl;
      cout << "Updating " << endl;

      runt->UpdateFromRunTable(rctt);
    }

  }

  TFile* ftag = TFile::Open("Updated.root", "recreate");

  TTree * ttag = new TTree("T","A Tree with event tags");
  TBranch * btag = ttag->Branch("AliTAG", &runt);
  ttag->Fill();
  
  btag->SetCompressionLevel(9);

  ftag->cd();
  ttag->Write();
  ftag->Close();


}
