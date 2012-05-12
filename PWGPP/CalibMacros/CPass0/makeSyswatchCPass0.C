//
// Make syswatch resource usage tree
//
void makeSyswatchCPass0(const char * fname){
  //
  //
  //
  TFile * fout= new TFile(fname,"update");
  TTree * treeRec = AliSysInfo::MakeTree("syswatch_rec.log");
  TTree * treeCalib = AliSysInfo::MakeTree("syswatch_calib.log");
  if (treeRec) treeRec->Write("syswatchRec");
  if (treeCalib) treeCalib->Write("syswatchCalib");
  fout->Close();
}
