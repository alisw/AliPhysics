void Create_xml(TString aliendir, int stage,TString pattern="*root_archive.zip"){

  if (gGrid && gGrid->IsConnected()) return kTRUE;
  if (!gGrid) {
    Info("Connect", "Trying to connect to AliEn ...");
    TGrid::Connect("alien://");
  }
  if (!gGrid || !gGrid->IsConnected()) {
    Error("Connect", "Did not managed to connect to AliEn. Make sure you have a valid token.");
    return kFALSE;
  } 
  //int stage=1;
  //  = "*root_archive.zip";
  //  TString pattern = "*DxHFE_eD0Corr.root";

  //  if (stage>1) pattern = Form("Stage_%d/*root_archive.zip", stage-1);
  //  const char *aliendir = "MinBias_020113/output";
  //  TString aliendir="MinBias_020513_3/output";
  //TString aliendir="LHC2010_020513_2/output";
  if(stage>1) aliendir+=Form("/Stage_%d",stage-1);
  cout << aliendir.Data() << endl;
  TGridResult *res = gGrid->Command(Form("find -x Stage_%d %s %s", stage, aliendir.Data(), pattern.Data()));
  if (res) delete res;
  // Write standard output to file
  gROOT->ProcessLine(Form("gGrid->Stdout(); > %s", Form("Stage_%d.xml",stage)));
}
