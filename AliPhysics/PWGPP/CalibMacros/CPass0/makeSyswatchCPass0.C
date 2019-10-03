//
// Make syswatch resource usage tree
//
void makeSyswatchCPass0(const char * fname){
  //
  //
  //
  gSystem->Exec("sed -i 's|\t| |g;s|[[:space:]]*$||' syswatch*.log");
  TFile * fout= new TFile(fname,"update");
  TTree * treeRec = AliSysInfo::MakeTree("syswatch_rec.log");
  TTree * treeCalib = AliSysInfo::MakeTree("syswatch_calib.log");
  if (treeRec) {
    treeRec->SetName("syswatchRec");
    treeRec->Write("syswatchRec");
  }
  if (treeCalib){
    treeCalib->SetName("syswatchCalib");
    treeCalib->Write("syswatchCalib");
  }
  fout->Close();
}
