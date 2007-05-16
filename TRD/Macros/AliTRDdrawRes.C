void AliTRDdrawRes(const char *path) {
  
  TGaxis::SetMaxDigits(3);
  gStyle->SetPadGridX(1);

  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  AliCDBManager::Instance()->SetRun(0);
  
  AliTRDtrackingAnalysis *analysis = new AliTRDtrackingAnalysis();
  analysis->SetPath(path);
  
  //analysis->DrawTrackletResolution(0, 1);
  analysis->DrawResolutionPt(0, 100);
  //analysis->DrawRecPointResolution(0, 100);
}
