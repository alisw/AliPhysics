void cloneAODTreeAndRemoveObject(const char *newFileName = "AliAOD_new.root", const char *orgFileName = "AliAOD.root") {
  // This little macro takes an already created AOD file and clones it.
  // After removing an old brach, the new TTree is written to a new file.

  // open input file and get the TTree
  TFile orgFile(orgFileName, "READ");
  // get original TTree
  TTree *orgAodTree = (TTree*)orgFile.Get("aodTree");
  // do your gymnastics with the old TTree
  AliAODEvent *evOrg = new AliAODEvent();
  evOrg->ReadFromTree(orgAodTree);

  // switch off one branch (and its subbranches!)
  orgAodTree->SetBranchStatus("tracks*", 0);
  // remove TObject from the list
  evOrg->RemoveObject(evOrg->GetTracks());

  // open new output file
  TFile *newFile = new TFile(newFileName, "RECREATE");
  // clone old TTree (only clones branches that are switched on)
  TTree *newAodTree = orgAodTree->CloneTree();

  // write new TTree to file
  newAodTree->Write();

  // close files
  newFile->Close();
  delete newFile;

  orgFile.Close();
}
