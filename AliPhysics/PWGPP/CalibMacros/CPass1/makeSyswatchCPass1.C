//
// Make syswatch resource usage tree
//
void makeSyswatchCPass1(const char * fname){
  //
  //
  //
  TFile * fout= new TFile(fname,"update");
  TTree * treeRec = AliSysInfo::MakeTree("syswatch_rec_Barrel.log");
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
