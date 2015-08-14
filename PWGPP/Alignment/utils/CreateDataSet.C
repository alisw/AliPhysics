void CreateDataSet(TString path="/alice/data/2012/LHC12g/000188503/cpass1_pass2", 
		   int nTest=-1, 
		   TString file="AliESDs_Barrel.root",
		   TString out="coll.xml")
{
  // Create dataset for the grid data directory + run number.
  const Int_t maxEntries = 9999999;
  if (nTest<0||nTest>maxEntries) nTest = maxEntries;
  if (!gGrid) {
    TGrid::Connect("alien://");
    if (!gGrid) {exit(1);}
  }   
  // Compose the 'find' command arguments
  TString command = Form("find -x collection -l%d %s *%s", nTest , path.Data(), file.Data());
  printf("Call %s\n",command.Data());
  TGridResult *res = gGrid->Command(command.Data());
  if (res) delete res;
  gROOT->ProcessLine(Form("gGrid->Stdout(); > %s",out.Data()));
}


